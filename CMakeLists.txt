cmake_minimum_required(VERSION 3.4)

project(RajaProxies LANGUAGES CXX)

option(ENABLE_LULESH_ONE "Build LULESH v1.0" On)
option(ENABLE_LULESH_TWO "Build LULESH v2.0" On)
option(ENABLE_KRIPKE "Build Kripke" On)

set (ENABLE_EXAMPLES Off CACHE Bool "Don't build examples")

if (PROXIES_ENABLE_WARNINGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")
endif()

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

if (ENABLE_CHAI)
  set(ENABLE_COPY_HEADERS On CACHE BOOL "")
endif ()

include(${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)

set(CUDA_NVCC_FLAGS ${RAJA_NVCC_FLAGS})

if (ENABLE_CHAI)
  add_subdirectory(chai)
endif ()

add_subdirectory(tpl/RAJA)

if (ENABLE_DOCUMENTATION)
  add_subdirectory(docs)
endif ()

if (ENABLE_LULESH_ONE)
  message(STATUS "Building LULESH v1.0")
  add_subdirectory(lulesh-v1.0)
endif ()

if (ENABLE_LULESH_TWO)
  message(STATUS "Building LULESH v2.0")
  add_subdirectory(lulesh-v2.0)
endif ()

if (ENABLE_KRIPKE)
  if (ENABLE_CUDA)
    message(FATAL_ERROR "Cannot build Kripke with CUDA")
  endif()
  if (ENABLE_OPENMP)
    set(KRIPKE_ARCH "OpenMP")
  else()
    set(KRIPKE_ARCH "Sequential")
  endif()
  message(STATUS "Building Kripke")
  add_subdirectory(kripke)
endif ()
