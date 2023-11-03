###############################################################################
# Copyright (c) 2016-21, Lawrence Livermore National Security, LLC
# and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
###############################################################################

set(RAJA_COMPILER "RAJA_COMPILER_XLC" CACHE STRING "")

set(CMAKE_C_COMPILER "xlc_r" CACHE PATH "")
set(CMAKE_CXX_COMPILER "xlc++_r" CACHE PATH "")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++1y --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1/" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g9 -std=c++1y --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1/" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -qsmp=omp:noopt -std=c++1y --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1/" CACHE STRING "")
#set(CMAKE_EXE_LINKER_FLAGS "-qpic,-Wl,-z,muldefs" CACHE STRING "")

set(HOST_OPT_FLAGS "-Xcompiler -O3 -Xcompiler -qxlcompatmacros -Xcompiler -qalias=noansi -Xcompiler -qsmp=omp -Xcompiler -qhot -Xcompiler -qnoeh -Xcompiler -qsuppress=1500-029 -Xcompiler -std=c++1y -Xcompiler --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1/ -Xcompiler -qsuppress=1500-036")

set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -std=c++14 --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1/ ${HOST_OPT_FLAGS}" CACHE STRING "")
set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -std=c++14 --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1/" CACHE STRING "")
set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo -O3 -std=c++14 --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1/ ${HOST_OPT_FLAGS}" CACHE STRING "")

# Suppressed XLC warnings:
# - 1500-029 cannot inline
# - 1500-036 nostrict optimizations may alter code semantics
#   (can be countered with -qstrict, with less optimization)

set(RAJA_DATA_ALIGN 64 CACHE STRING "")

set(RAJA_HOST_CONFIG_LOADED On CACHE BOOL "")

