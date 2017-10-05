#!/bin/bash

##
## Copyright (c) 2017, Lawrence Livermore National Security, LLC.
##
## Produced at the Lawrence Livermore National Laboratory.
##
## All rights reserved.
## 
## This file is part of the RAJA Proxy App Suite.
##

rm -rf build_blueos-gcc-4.9.3 2>/dev/null
mkdir build_blueos-gcc-4.9.3 && cd build_blueos-gcc-4.9.3

module load cmake/3.7.2

RAJA_PROXIES_DIR=$(git rev-parse --show-toplevel)

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -C ${RAJA_PROXIES_DIR}/host-configs/blueos/gcc_4_9_3.cmake \
  -DENABLE_OPENMP=On \
  -DPROXIES_ENABLE_WARNING=Off \
  -DENABLE_ALL_WARNINGS=Off \
  -DCMAKE_INSTALL_PREFIX=../install_blueos-gcc-4.9.3 \
  "$@" \
  ${RAJA_PROXIES_DIR}
