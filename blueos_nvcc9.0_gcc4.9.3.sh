#!/bin/bash

##
## Copyright (c) 2017, Lawrence Livermore National Security, LLC.
##
## Produced at the Lawrence Livermore National Laboratory.
##
## LLNL-CODE-738930
##
## All rights reserved.
## 
## This file is part of the RAJA Performance Suite.
##
## For details about use and distribution, please read raja-perfsuite/LICENSE.
##

rm -rf build_blueos_nvcc9.0_gcc4.9.3 >/dev/null
mkdir build_blueos_nvcc9.0_gcc4.9.3 && cd build_blueos_nvcc9.0_gcc4.9.3

module load cmake/3.7.2
module load gcc/4.9.3

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -C ../host-configs/blueos/nvcc_gcc_4_9_3.cmake \
  -DENABLE_OPENMP=Off \
  -DENABLE_CUDA=On \
  -DENABLE_KRIPKE=Off \
  -DENABLE_LULESH_ONE=Off \
  -DENABLE_LULESH_TWO=Off \
  -DCUDA_TOOLKIT_ROOT_DIR=/usr/tce/packages/cuda/cuda-9.2.148 \
  -DENABLE_ALL_WARNINGS=Off \
  -DCMAKE_INSTALL_PREFIX=../install_blueos_nvcc9.0_gcc4.9.3 \
  -DENABLE_MPI=On \
  -DCMAKE_EXE_LINKER_FLAGS="-L/g/g92/meng3/software/caliper_install_ray/lib64 -lcaliper" \
  -DCUDA_ARCH=sm_60 \
  -DCUDA_NVCC_FLAGS_RELEASE="-I/g/g92/meng3/software/caliper_install_ray/include" \
  "$@" \
  ..
