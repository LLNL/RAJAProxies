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

rm -rf build_chaos-clang-4.0.0_debug 2>/dev/null
mkdir build_chaos-clang-4.0.0_debug && cd build_chaos-clang-4.0.0_debug
. /usr/local/tools/dotkit/init.sh && use cmake-3.4.1

RAJA_PROXIES_DIR=$(git rev-parse --show-toplevel)

cmake \
  -DCMAKE_BUILD_TYPE=Debug \
  -C ${RAJA_PROXIES_DIR}/host-configs/chaos/clang_4_0_0.cmake \
  -DENABLE_OPENMP=On \
  -DPROXIES_ENABLE_WARNING=On \
  -DENABLE_ALL_WARNINGS=Off \
  -DCMAKE_INSTALL_PREFIX=../install_chaos-clang-4.0.0_debug \
  "$@" \
  ${RAJA_PROXIES_DIR}
