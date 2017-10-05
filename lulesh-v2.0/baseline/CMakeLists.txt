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

if(ENABLE_OPENMP)
  add_definitions(-DUSE_OMP=1)
  add_definitions(-DUSE_MPI=0)
  blt_add_executable(
    NAME lulesh-v2.0-baseline.exe
    SOURCES lulesh.cpp
            lulesh-comm.cc
            lulesh-init.cc
            lulesh-util.cc
            lulesh-viz.cc
    DEPENDS_ON RAJA openmp)
endif()
