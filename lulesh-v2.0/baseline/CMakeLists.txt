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

set (lulesh_depends
  RAJA)

# We need to include RAJA with the same features it was built with, even if
# they are not used.
if (ENABLE_OPENMP)
  set (lulesh_depends
    ${lulesh_depends}
    openmp)
endif ()

set (lulesh_sources
  lulesh.cpp
  lulesh-comm.cpp
  lulesh-init.cpp
  lulesh-util.cpp
  lulesh-viz.cpp)

blt_add_executable(
  NAME lulesh-v2.0-baseline-seq.exe
  SOURCES ${lulesh_sources}
  DEPENDS_ON ${lulesh_depends})

target_compile_definitions(
  lulesh-v2.0-baseline-seq.exe
  PRIVATE 
  -DUSE_OMP=0)

target_compile_options(
  lulesh-v2.0-baseline-seq.exe
  PRIVATE
  -Wno-unknown-pragmas)

target_compile_definitions(
  lulesh-v2.0-baseline-seq.exe
  PRIVATE 
  -DUSE_MPI=0)

if(ENABLE_OPENMP)
  blt_add_executable(
    NAME lulesh-v2.0-baseline-omp.exe
    SOURCES ${lulesh_sources}
    DEPENDS_ON RAJA openmp)

  target_compile_definitions(
    lulesh-v2.0-baseline-omp.exe
    PRIVATE 
    -DUSE_OMP=1)

  target_compile_definitions(
    lulesh-v2.0-baseline-omp.exe
    PRIVATE 
    -DUSE_MPI=0)
endif()
