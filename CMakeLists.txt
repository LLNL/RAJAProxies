###############################################################################
#
# Copyright (c) 2017-22, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# All rights reserved.
#
# This file is part of the RAJA Proxy App Suite
#
###############################################################################

cmake_minimum_required(VERSION 3.4)

project(RajaProxies LANGUAGES CXX)

option(ENABLE_LULESH_ONE "Build LULESH v1.0" On)
option(ENABLE_LULESH_TWO "Build LULESH v2.0" On)
option(ENABLE_KRIPKE "Build Kripke" On)
option(ENABLE_COMD "Build CoMD" On)

set (ENABLE_TESTS Off CACHE BOOL "Don't build tests")
set (ENABLE_EXAMPLES Off CACHE BOOL "Don't build examples")
set (ENABLE_EXERCISES Off CACHE BOOL "Don't build exercises")

if (PROXIES_ENABLE_WARNINGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")
endif()

################################
# BLT
################################
if (ENABLE_KRIPKE)
  set(BLT_CXX_STD "c++14" CACHE STRING "Version of C++ Standard")
endif()
if (NOT EXISTS ${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)
  message(FATAL_ERROR "\
  The BLT submodule is not present. \
  If in git repository run the following two commands:\n \
  git submodule init\n \
  git submodule update")
endif()
include(${PROJECT_SOURCE_DIR}/blt/SetupBLT.cmake)

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
    set(KRIPKE_ARCH "Cuda")
    set(kripke_compile_options --expt-extended-lambda -arch=${CUDA_ARCH})
  elseif (ENABLE_HIP)
    set(ENABLE_CHAI On CACHE BOOL "")
    set(KRIPKE_ARCH "HIP")
  elseif (ENABLE_OPENMP)
    set(KRIPKE_ARCH "OpenMP")
  else()
    set(KRIPKE_ARCH "Sequential")
  endif()
  message(STATUS "Building Kripke")
  add_subdirectory(kripke)

  target_compile_options(kripke
    PUBLIC
    ${kripke_compile_options})
endif ()


if (ENABLE_COMD)
  message(STATUS "Building CoMD v1.0")
  add_subdirectory(CoMD)
endif ()
