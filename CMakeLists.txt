cmake_minimum_required(VERSION 3.4)

project(RajaProxies LANGUAGES CXX)

option(ENABLE_LULESH_ONE "Build LULESH v1.0" On)
option(ENABLE_LULESH_TWO "Build LULESH v2.0" Off)

################################
# BLT
################################
if (NOT EXISTS ${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)
  message(FATAL_ERROR "\
  The BLT submodule is not present. \
  If in git repository run the following two commands:\n \
  git submodule init\n \
  git submodule update")
endif()
include(${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)

add_subdirectory(tpl/RAJA)

message (STATUS "ENABLE_LULESH_ONE = ${ENABLE_LULESH_ONE}")
if (ENABLE_LULESH_ONE)
  message(STATUS "Building LULESH v1.0")
  add_subdirectory(lulesh-v1.0)
endif ()

if (ENABLE_LULESH_2)
  add_subdirectory()
endif ()
