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

rm -rf build_bgqos-clang-4.0.0 2>/dev/null
mkdir build_bgqos-clang-4.0.0 && cd build_bgqos-clang-4.0.0
. /usr/local/tools/dotkit/init.sh && use cmake-3.4.3

RAJA_PROXIES_DIR=$(git rev-parse --show-toplevel)

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -C ${RAJA_PROXIES_DIR}/host-configs/bgqos/clang_4_0_0.cmake \
  -DENABLE_OPENMP=On \
  -DPROXIES_ENABLE_WARNING=Off \
  -DENABLE_ALL_WARNINGS=Off \
  -DCMAKE_INSTALL_PREFIX=../install_bgqos_clang-4.0.0 \
  "$@" \
  ${RAJA_PROXIES_DIR}
