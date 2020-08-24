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

rm -rf build_toss3-icpc-18.0-beta 2>/dev/null
mkdir build_toss3-icpc-18.0-beta && cd build_toss3-icpc-18.0-beta

module load cmake/3.5.2
module load gcc/4.9.3

RAJA_PROXIES_DIR=$(git rev-parse --show-toplevel)

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -C ${RAJA_PROXIES_DIR}/host-configs/toss3/icpc_18_0_beta.cmake \
  -DENABLE_OPENMP=On \
  -DPROXIES_ENABLE_WARNINGS=Off \
  -DENABLE_ALL_WARNINGS=Off \
  -DCMAKE_INSTALL_PREFIX=../install_toss3-icpc-18.0-beta \
  "$@" \
  ${RAJA_PROXIES_DIR}
