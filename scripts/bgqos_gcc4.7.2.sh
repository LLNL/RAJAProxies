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

rm -rf build-bgqos_gcc-4.7.2 2>/dev/null
mkdir build-bgqos_gcc-4.7.2 && cd build-bgqos_gcc-4.7.2
. /usr/local/tools/dotkit/init.sh && use cmake-3.4.3

RAJA_PROXIES_DIR=$(git rev-parse --show-toplevel)

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -C ${RAJA_PROXIES_DIR}/host-configs/bgqos/gcc_4_7_2.cmake \
  -DENABLE_OPENMP=On \
  -DPROXIES_ENABLE_WARNINGS=Off \
  -DENABLE_ALL_WARNINGS=Off \
  -DCMAKE_INSTALL_PREFIX=../install-bgqos_gcc-4.7.2 \
  "$@" \
  ${RAJA_PROXIES_DIR}
