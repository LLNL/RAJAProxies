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

rm -rf build_blueos-clang-coral-2017.09.06 2>/dev/null
mkdir build_blueos-clang-coral-2017.09.06 && cd build_blueos-clang-coral-2017.09.06

module load cmake/3.7.2

RAJA_PROXIES_DIR=$(git rev-parse --show-toplevel)

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -C ${RAJA_PROXIES_DIR}/host-configs/blueos/clang_coral_2017_09_06.cmake \
  -DENABLE_OPENMP=On \
  -DPROXIES_ENABLE_WARNINGS=Off \
  -DENABLE_ALL_WARNINGS=Off \
  -DCMAKE_INSTALL_PREFIX=../install_blueos-clang-coral-2017.09.06 \
  "$@" \
  ${RAJA_PROXIES_DIR}
