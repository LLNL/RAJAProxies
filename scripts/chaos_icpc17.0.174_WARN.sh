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

rm -rf build_chaos-icpc-17.0.174 2>/dev/null
mkdir build_chaos-icpc-17.0.174 && cd build_chaos-icpc-17.0.174
. /usr/local/tools/dotkit/init.sh && use cmake-3.4.1 && use gcc-4.9.3p

RAJA_PROXIES_DIR=$(git rev-parse --show-toplevel)

cmake \
  -DCMAKE_BUILD_TYPE=Debug \
  -C ${RAJA_PROXIES_DIR}/host-configs/chaos/icpc_17_0_174.cmake \
  -DENABLE_OPENMP=On \
  -DPROXIES_ENABLE_WARNINGS=On \
  -DENABLE_ALL_WARNINGS=Off \
  -DCMAKE_INSTALL_PREFIX=../install_chaos-icpc-17.0.174 \
  "$@" \
  ${RAJA_PROXIES_DIR}
