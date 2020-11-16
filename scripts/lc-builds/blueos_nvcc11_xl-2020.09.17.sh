#!/usr/bin/env bash

###############################################################################
# Copyright (c) 2016-20, Lawrence Livermore National Security, LLC
# and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
###############################################################################

#
## NOTE: After building code, you need to load the cuda 11 module to run
##       your code or RAJA tests
#

BUILD_SUFFIX=lc_blueos-nvcc11-xl_2020.09.17
RAJA_HOSTCONFIG=../tpl/RAJA/host-configs/lc-builds/blueos/nvcc_xl_X.cmake

rm -rf build_${BUILD_SUFFIX} 2>/dev/null
mkdir build_${BUILD_SUFFIX} && cd build_${BUILD_SUFFIX}

module load cmake/3.14.5

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/usr/tce/packages/xl/xl-2020.09.17/bin/xlc++_r \
  -C ${RAJA_HOSTCONFIG} \
  -DENABLE_OPENMP=On \
  -DENABLE_CUDA=On \
  -DCUDA_TOOLKIT_ROOT_DIR=/usr/tce/packages/cuda/cuda-11.1.0 \
  -DCMAKE_CUDA_COMPILER=/usr/tce/packages/cuda/cuda-11.1.0/bin/nvcc \
  -DCUDA_ARCH=sm_70 \
  -DCMAKE_CUDA_STANDARD="14" \
  -DCMAKE_INSTALL_PREFIX=../install_${BUILD_SUFFIX} \
  "$@" \
  ..
