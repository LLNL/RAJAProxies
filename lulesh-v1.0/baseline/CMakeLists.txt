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

blt_add_executable(
  NAME lulesh-v1.0-baseline-seq.exe
  SOURCES lulesh.cpp
  DEPENDS_ON ${lulesh_depends})

target_compile_options(
  lulesh-v1.0-baseline-seq.exe
  PRIVATE
  -Wno-unknown-pragmas)

if(ENABLE_OPENMP)
  blt_add_executable(
    NAME lulesh-v1.0-baseline-omp.exe
    SOURCES lulesh.cpp
    DEPENDS_ON RAJA openmp)
endif()

