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
  NAME lulesh-v2.0-baseline-seq.exe
  SOURCES lulesh.cc
          lulesh-comm.cc
          lulesh-init.cc
          lulesh-util.cc
          lulesh-viz.cc
  DEPENDS_ON RAJA)

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
    NAME lulesh-v2.0-baseline-OMP.exe
    SOURCES lulesh.cc
            lulesh-comm.cc
            lulesh-init.cc
            lulesh-util.cc
            lulesh-viz.cc
    DEPENDS_ON RAJA openmp)

  target_compile_definitions(
    lulesh-v2.0-baseline-OMP.exe
    PRIVATE 
    -DUSE_OMP=1)

  target_compile_definitions(
    lulesh-v2.0-baseline-OMP.exe
    PRIVATE 
    -DUSE_MPI=0)
endif()
