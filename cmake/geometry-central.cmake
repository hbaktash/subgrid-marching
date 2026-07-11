if (TARGET geometry-central::geometry-central)
    return()
endif()

# Use a git submodule checkout of geometry-central instead of CPM.
# Expected location: deps/geometry-central at the project root.

get_filename_component(NC_PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(GEOMETRY_CENTRAL_SOURCE_DIR "${NC_PROJECT_ROOT}/deps/geometry-central")

if (NOT EXISTS "${GEOMETRY_CENTRAL_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "geometry-central submodule not found at ${GEOMETRY_CENTRAL_SOURCE_DIR}. Did you run 'git submodule update --init --recursive'?")
endif()

add_subdirectory("${GEOMETRY_CENTRAL_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/geometry-central" EXCLUDE_FROM_ALL)

if (TARGET geometry-central)
    set_target_properties(geometry-central PROPERTIES FOLDER third_party)
    add_library(geometry-central::geometry-central ALIAS geometry-central)
elseif (NOT TARGET geometry-central::geometry-central)
    message(FATAL_ERROR "geometry-central target was not created by the submodule CMakeLists.txt")
endif()


# Original CPM-based version for reference:
# include(CPM)
# CPMAddPackage(
#     NAME geometry-central
#     GITHUB_REPOSITORY MarkGillespie/geometry-central
#     GIT_TAG NonmanifoldTwinConstructor
# )
# set_target_properties(geometry-central PROPERTIES FOLDER third_party)
# add_library(geometry-central::geometry-central ALIAS geometry-central)

