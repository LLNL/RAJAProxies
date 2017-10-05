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

rm -rf build_blueos-clang-coral-2017.08.31 2>/dev/null
mkdir build_blueos-clang-coral-2017.08.31 && cd build_blueos-clang-coral-2017.08.31

module load cmake/3.7.2

RAJA_PROXIES_DIR=$(git rev-parse --show-toplevel)

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -C ${RAJA_PROXIES_DIR}/host-configs/blueos/clang_coral_2017_08_31.cmake \
  -DENABLE_OPENMP=On \
  -DPROXIES_ENABLE_WARNING=Off \
  -DENABLE_ALL_WARNINGS=Off \
  -DCMAKE_INSTALL_PREFIX=../install_blueos-clang-coral-2017.08.31 \
  "$@" \
  ${RAJA_PROXIES_DIR}
