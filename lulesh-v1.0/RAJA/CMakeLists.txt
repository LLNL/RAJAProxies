###############################################################################
#
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# All rights reserved.
#
# This file is part of the RAJA Proxy App Suite
#
###############################################################################

blt_add_executable(
  NAME lulesh-v1.0-RAJA-seq.exe
  SOURCES lulesh.cpp
  DEPENDS_ON RAJA)

target_compile_definitions(
  lulesh-v1.0-RAJA-seq.exe
  PRIVATE
  -DUSE_CASE=1
)

target_compile_options(
  lulesh-v1.0-RAJA-seq.exe
  PRIVATE
  -Wno-unknown-pragmas)

target_include_directories(
  lulesh-v1.0-RAJA-seq.exe
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR})

if (ENABLE_OPENMP)
  blt_add_executable(
    NAME lulesh-v1.0-RAJA-omp.exe
    SOURCES lulesh.cpp
    DEPENDS_ON RAJA openmp)

  target_include_directories(
    lulesh-v1.0-RAJA-omp.exe
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})
endif ()

if (ENABLE_CUDA)
  blt_add_executable(
    NAME lulesh-v1.0-RAJA-cuda.exe
    SOURCES lulesh.cpp
    DEPENDS_ON RAJA openmp cuda)

  target_compile_definitions(
    lulesh-v1.0-RAJA-cuda.exe
    PRIVATE
    -DUSE_CASE=9
  )

  target_include_directories(
    lulesh-v1.0-RAJA-cuda.exe
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})
endif ()
