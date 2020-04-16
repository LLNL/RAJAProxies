#!/bin/bash

##
## Copyright (c) 2017-20, Lawrence Livermore National Security, LLC.
##
## Produced at the Lawrence Livermore National Laboratory.
##
## All rights reserved.
## 
## This file is part of the RAJA Proxy App Suite.
##

rm -rf build_blueos_nvcc8.0_clang-coral-2017.09.06 >/dev/null
mkdir build_blueos_nvcc8.0_clang-coral-2017.09.06 && cd build_blueos_nvcc8.0_clang-coral-2017.09.06

module load cmake/3.7.2

RAJA_PROXIES_DIR=$(git rev-parse --show-toplevel)

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -C ${RAJA_PROXIES_DIR}/host-configs/blueos/nvcc_clang_coral_2017_09_06.cmake \
  -DENABLE_OPENMP=On \
  -DENABLE_CUDA=On \
  -DCUDA_TOOLKIT_ROOT_DIR=/usr/tcetmp/packages/cuda-8.0 \
  -DPROXIES_ENABLE_WARNINGS=Off \
  -DENABLE_ALL_WARNINGS=Off \
  -DCMAKE_INSTALL_PREFIX=../install_blueos_nvcc8.0_clang-coral-2017.09.06 \
  "$@" \
  ${RAJA_PROXIES_DIR}
