include(catch2)
include(CTest)

list(APPEND CMAKE_MODULE_PATH "${Catch2_SOURCE_DIR}/extras")
include(Catch)

add_executable(tests
    test_nc_solver.cpp
    test_subgrid_pipeline.cpp
    test_boundary_curves.cpp
    test_comb_merge.cpp
    test_dual_subgrid_pipeline.cpp
    test_self_intersection.cpp
)

target_link_libraries(tests
    PRIVATE
    normal_coordinates::normal_coordinates
    Catch2::Catch2WithMain
)

set_property(TARGET tests PROPERTY CXX_STANDARD 17)

set(TEST_DATA_DIR "${CMAKE_SOURCE_DIR}/tests/data/meshes"
    CACHE PATH "Directory containing mesh files for pipeline tests")

set(TEST_NPZ_DIR "${CMAKE_SOURCE_DIR}/tests/data/npz"
    CACHE PATH "Directory containing .npz files for explicit-input pipeline tests")

target_compile_definitions(tests PRIVATE
    TEST_DATA_DIR="${TEST_DATA_DIR}"
    TEST_NPZ_DIR="${TEST_NPZ_DIR}"
)

catch_discover_tests(tests)
