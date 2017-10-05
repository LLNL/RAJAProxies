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

if (ENABLE_CUDA)
  add_definitions(-DUSE_CASE=9)
endif()

if (ENABLE_OPENMP)
  blt_add_executable(
    NAME lulesh-v1.0-RAJA.exe
    SOURCES lulesh.cpp
    DEPENDS_ON RAJA openmp)

  target_include_directories(
    lulesh-v1.0-RAJA.exe
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})
endif ()
