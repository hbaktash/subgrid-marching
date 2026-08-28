# Python bindings (pybind11). Only entered when -DSUBGRID_BUILD_PYTHON=ON, so the
# plain C++ build never looks for Python or pybind11.

# A previous `pip install .` caches pybind11_DIR inside pip's temporary build
# environment, which is deleted right after. Drop the stale entry so the search
# below resolves against the environment being built in now, instead of quietly
# falling through to the FetchContent copy.
if(pybind11_DIR AND NOT pybind11_DIR MATCHES "NOTFOUND$" AND NOT EXISTS "${pybind11_DIR}")
    message(STATUS "subgrid: cached pybind11_DIR no longer exists; searching again.")
    unset(pybind11_DIR CACHE)
endif()

# Prefer an installed pybind11 (pip / conda / system, and the build isolation
# environment used by `pip install .`), fall back to fetching a pinned copy.
find_package(pybind11 CONFIG QUIET)

if(NOT pybind11_FOUND)
    # Ask the active interpreter where its pybind11 keeps its CMake package.
    find_package(Python 3.8 COMPONENTS Interpreter QUIET)
    if(Python_EXECUTABLE)
        execute_process(
            COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
            OUTPUT_VARIABLE _pybind11_cmakedir
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(_pybind11_cmakedir)
            find_package(pybind11 CONFIG QUIET HINTS "${_pybind11_cmakedir}")
        endif()
    endif()
endif()

if(NOT pybind11_FOUND)
    message(STATUS "subgrid: pybind11 not found; fetching it.")
    include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.13.6
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(pybind11)
else()
    message(STATUS "subgrid: using pybind11 ${pybind11_VERSION} from ${pybind11_INCLUDE_DIRS}")
endif()

pybind11_add_module(_subgrid_marching_core bindings.cpp)

set_target_properties(_subgrid_marching_core PROPERTIES OUTPUT_NAME _core)

target_link_libraries(_subgrid_marching_core PRIVATE normal_coordinates::normal_coordinates)
set_property(TARGET _subgrid_marching_core PROPERTY CXX_STANDARD 17)

# `pip install .` (scikit-build-core) picks the module up from here; a plain CMake
# build just leaves it in the build tree next to the executables.
if(SKBUILD OR DEFINED SKBUILD_PROJECT_NAME)
    install(TARGETS _subgrid_marching_core LIBRARY DESTINATION subgrid_marching)
else()
    set_target_properties(_subgrid_marching_core PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/subgrid_marching")
    # Mirror the pure-Python package next to the extension so the build tree is
    # importable directly: PYTHONPATH=build/python python -c "import subgrid_marching"
    add_custom_command(TARGET _subgrid_marching_core POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                "${CMAKE_CURRENT_SOURCE_DIR}/subgrid_marching"
                "${CMAKE_CURRENT_BINARY_DIR}/subgrid_marching"
        COMMENT "Staging the subgrid_marching Python package in the build tree")
endif()
