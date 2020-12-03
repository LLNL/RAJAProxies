#!/usr/bin/env bash

###############################################################################
# Copyright (c) 2016-20, Lawrence Livermore National Security, LLC
# and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
###############################################################################

BUILD_SUFFIX=lc_toss3-hipcc-3.8.0
RAJA_HOSTCONFIG=../tpl/RAJA/host-configs/lc-builds/toss3/hip.cmake

rm -rf build_${BUILD_SUFFIX} >/dev/null
mkdir build_${BUILD_SUFFIX} && cd build_${BUILD_SUFFIX}


module load cmake/3.14.5

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DRAJA_HIPCC_FLAGS="--amdgpu-target=gfx906" \
  -DHIP_ROOT_DIR=/opt/rocm-3.8.0/hip \
  -C ${RAJA_HOSTCONFIG} \
  -DENABLE_HIP=ON \
  -DENABLE_OPENMP=OFF \
  -DENABLE_CUDA=OFF \
  -DCMAKE_INSTALL_PREFIX=../install_${BUILD_SUFFIX} \
  "$@" \
  ..
