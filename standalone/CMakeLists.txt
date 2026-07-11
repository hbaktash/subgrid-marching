# Self-contained build for the standalone port's verification harnesses.
#
# This is a top-level CMake project of its own: it does NOT need to be wired into
# the repo-root CMakeLists.txt. It pulls the repo root in as a subdirectory
# (EXCLUDE_FROM_ALL, so none of the root's apps or tests build) only to obtain
# the core library target `normal_coordinates::normal_coordinates`, which is the
# ground truth the port is checked against.
#
# Usage (from the repo root):
#   cmake -S standalone -B standalone/build -DSUBGRID_POLYSCOPE_VIEWER=OFF
#   cmake --build standalone/build --target verify_port verify_port_dual
#   ./standalone/build/verify_port
#   ./standalone/build/verify_port_dual
#
# The harnesses compare the copy-pastable ports in cpp/ (and, via dumped
# fixtures, python/) against that core, in-process. They are verification-only
# and are never part of the shipped, dependency-free port.

cmake_minimum_required(VERSION 3.14.0)
project(subgrid_standalone_verify CXX)

# Pull in the repo root to build the core library. Guarded so this file also
# works if it is ever add_subdirectory()'d from a build that already defines it.
if(NOT TARGET normal_coordinates::normal_coordinates)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/..
                   ${CMAKE_CURRENT_BINARY_DIR}/normal_coordinates_core
                   EXCLUDE_FROM_ALL)
endif()

add_executable(verify_port verify/verify_primal.cpp)
target_link_libraries(verify_port PRIVATE normal_coordinates::normal_coordinates)
target_include_directories(verify_port PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../tests
  ${CMAKE_CURRENT_SOURCE_DIR}/cpp)
set_property(TARGET verify_port PROPERTY CXX_STANDARD 17)

add_executable(verify_port_dual verify/verify_dual.cpp)
target_link_libraries(verify_port_dual PRIVATE normal_coordinates::normal_coordinates)
target_include_directories(verify_port_dual PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../tests
  ${CMAKE_CURRENT_SOURCE_DIR}/cpp)
set_property(TARGET verify_port_dual PROPERTY CXX_STANDARD 17)
