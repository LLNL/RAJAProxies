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

set (lulesh_sources
  lulesh.cpp
  lulesh-comm.cpp
  lulesh-init.cpp
  lulesh-util.cpp
  lulesh-viz.cpp)

blt_add_executable(
  NAME lulesh-v2.0-RAJA-seq.exe
  SOURCES ${lulesh_sources}
  DEPENDS_ON RAJA
)

target_compile_options(
  lulesh-v2.0-RAJA-seq.exe
  PRIVATE
  -Wno-unknown-pragmas)

target_compile_definitions(
  lulesh-v2.0-RAJA-seq.exe
  PRIVATE 
  -DUSE_MPI=0)

target_compile_definitions(
  lulesh-v2.0-RAJA-seq.exe
  PRIVATE
  -DUSE_CASE=1)

if (ENABLE_OPENMP)
  blt_add_executable(
    NAME lulesh-v2.0-RAJA-omp.exe
    SOURCES ${lulesh_sources}
    DEPENDS_ON RAJA openmp)

  target_compile_definitions(
    lulesh-v2.0-RAJA-omp.exe
    PRIVATE 
    -DUSE_MPI=0)
endif ()

if (ENABLE_CUDA)
  blt_add_executable(
    NAME lulesh-v2.0-RAJA-cuda.exe
    SOURCES ${lulesh_sources}
    DEPENDS_ON RAJA openmp cuda)

  target_compile_definitions(
    lulesh-v2.0-RAJA-cuda.exe
    PRIVATE
    -DUSE_CASE=9)

  target_compile_definitions(
    lulesh-v2.0-RAJA-cuda.exe
    PRIVATE 
    -DUSE_MPI=0)
endif ()
