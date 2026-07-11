cmake_minimum_required(VERSION 3.14.0)
set(CMAKE_POLICY_VERSION_MINIMUM "3.5" CACHE STRING "Minimum policy version to suppress errors in older dependencies")
set(ABSL_PROPAGATE_CXX_STD ON)

project(normal_coordinates)

# Enforce C++17 across the project to match vendored Eigen expectations
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wno-invalid-specialization)
endif()

# Add cmake recipe directory
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

# Build the interactive Polyscope viewer? Turn OFF for a lightweight, headless
# CLI-only build (no OpenGL/GLFW dependency). The `subgrid` and `dualSubgrid`
# executables still write their output meshes; only the on-screen viewer and the
# GUI-only demos (singleTetSubgrid / ringTetSubgrid) require this.
option(SUBGRID_POLYSCOPE_VIEWER "Build with the interactive Polyscope viewer" ON)

# Dependencies
include(geometry-central)
if(SUBGRID_POLYSCOPE_VIEWER)
    include(polyscope)
endif()
include(args)
# include(igl)
include(fcpw)

# Add sdf-dataset as a subdirectory. In a headless build it must not pull in
# Polyscope either (it otherwise builds its own viewer and sources glm from it).
if(NOT SUBGRID_POLYSCOPE_VIEWER)
    set(SDF_WITH_VIEWER OFF CACHE BOOL "sdf-dataset viewer follows SUBGRID_POLYSCOPE_VIEWER" FORCE)
endif()
add_subdirectory(deps/sdf-dataset)

# Library source files
set(SRC_FILES
    src/common/utils.cpp
    src/common/visual_utils.cpp
    src/common/file_IO.cpp
    src/common/npz_reader.cpp
    src/grid/mesh_factory.cpp
    src/assembly/mesh_processing.cpp
    src/assembly/dual_mesh_assembly.cpp
    src/query/sdf_queries.cpp
    src/nc/NC_solver.cpp
    src/query/fcpw_robust_queries.cpp
    src/query/intersection_query.cpp
    src/query/mod2_reduction.cpp
    src/subgrid_MT/dual/qef_solver.cpp
    src/subgrid_MT/dual/dual_construction.cpp
    src/query/input_query_handler.cpp
    src/query/precomputed_query_handler.cpp
    src/subgrid_MT/boundary_curve.cpp
    src/common/triangle_soup.cpp
    src/subgrid_MT/cv_interpolation.cpp
    src/subgrid_MT/primal_reconstruction.cpp
    src/subgrid_MT/normal_reconstruction.cpp
    src/subgrid_MT/evensum_construction.cpp
    src/subgrid_MT/simplicial_embedding.cpp
    src/subgrid_pipeline.cpp
    src/dual_subgrid_pipeline.cpp
  )

# Add normal_coordinates library, which will be shared by all executables.
add_library(normal_coordinates ${SRC_FILES})
# target_include_directories(normal_coordinates PUBLIC include "${GUROBI_INCLUDE_DIRS}") #
target_include_directories(normal_coordinates PUBLIC include headerDeps) #
target_link_libraries(normal_coordinates
  PUBLIC
  geometry-central::geometry-central
  # igl_copyleft::cgal
  sdf_lib
)
if(SUBGRID_POLYSCOPE_VIEWER)
  target_link_libraries(normal_coordinates PUBLIC polyscope::polyscope)
  target_compile_definitions(normal_coordinates PUBLIC HAVE_POLYSCOPE=1)
endif()
# Link fcpw into the core library so its include dirs propagate to any compile unit that includes your headers.
if (TARGET fcpw::fcpw)
  target_link_libraries(normal_coordinates PUBLIC fcpw::fcpw)
  target_compile_definitions(normal_coordinates PUBLIC HAVE_FCPW=1)
  target_compile_definitions(normal_coordinates PUBLIC GLM_ENABLE_EXPERIMENTAL GLM_FORCE_SWIZZLE)
endif()

# ChatGPT'd the hell outta this part
# Ensure Enoki's include dirs are also propagated (CPM may create 'enoki' or 'enoki::enoki').
if (TARGET enoki::enoki OR TARGET enoki)
  if (TARGET enoki::enoki)
    target_link_libraries(normal_coordinates PUBLIC enoki::enoki)
  else()
    target_link_libraries(normal_coordinates PUBLIC enoki)
    add_library(enoki::enoki ALIAS enoki)
  endif()
else()
  # Fallback: CPM/CPM-fetched enoki lives under build/_deps; export its include dir so <enoki/array.h> is found.
  if (EXISTS "${CMAKE_BINARY_DIR}/_deps/enoki-src/include")
    target_include_directories(normal_coordinates PUBLIC "${CMAKE_BINARY_DIR}/_deps/enoki-src/include")
  elseif (EXISTS "${CMAKE_BINARY_DIR}/_deps/fcpw-src/deps/enoki/include")
    target_include_directories(normal_coordinates PUBLIC "${CMAKE_BINARY_DIR}/_deps/fcpw-src/deps/enoki/include")
  endif()
endif()

add_library(normal_coordinates::normal_coordinates ALIAS normal_coordinates)

enable_testing()
add_subdirectory(tests)

# NOTE: The dependency-free standalone port and its verification harnesses live
# entirely under standalone/ and are NOT wired into this build. To run the
# verifications, configure standalone/ as its own project (it pulls this repo in
# as the core library on its own); see standalone/README.md.


# GUI-only debugging demos: only built with the Polyscope viewer.
if(SUBGRID_POLYSCOPE_VIEWER)
  # single/double tet tests
  add_executable(singleTetSubgrid apps/main_single_tet.cpp)
  target_link_libraries(singleTetSubgrid
    PRIVATE
    normal_coordinates::normal_coordinates
    taywee::args
  )
  set_property(TARGET singleTetSubgrid PROPERTY CXX_STANDARD 17)

  # 4-tet ring for comb-merge debugging
  add_executable(ringTetSubgrid apps/main_ring_tet.cpp)
  target_link_libraries(ringTetSubgrid
    PRIVATE
    normal_coordinates::normal_coordinates
    taywee::args
  )
  set_property(TARGET ringTetSubgrid PROPERTY CXX_STANDARD 17)
endif()


# Implicit tet mesh + mesh intersection query
add_executable(subgrid apps/main_subgrid.cpp)
target_link_libraries(subgrid
  PRIVATE
  normal_coordinates::normal_coordinates
  taywee::args
)
set_property(TARGET subgrid PROPERTY CXX_STANDARD 17)


# Implicit tet mesh + mesh intersection query
add_executable(dualSubgrid apps/main_dual_subgrid.cpp)
target_link_libraries(dualSubgrid
  PRIVATE
  normal_coordinates::normal_coordinates
  taywee::args
)
set_property(TARGET dualSubgrid PROPERTY CXX_STANDARD 17)
