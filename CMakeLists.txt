cmake_minimum_required(VERSION 3.4)

project(RajaProxies LANGUAGES CXX)

option(ENABLE_LULESH_1 On "Build LULESH v1.0")
option(ENABLE_LULESH_2 On "Build LULESH v2.0")

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

if (ENABLE_LULESH_1)
  add_subdirectory()
endif ()

if (ENABLE_LULESH_2)
  add_subdirectory()
endif ()
