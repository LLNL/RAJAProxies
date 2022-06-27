###############################################################################
#
# Copyright (c) 2017-22, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# All rights reserved.
#
# This file is part of the RAJA Proxy App Suite
#
###############################################################################

set (lulesh_depends
  RAJA)

if (ENABLE_OPENMP)
  set (lulesh_depends
    ${lulesh_depends}
    openmp)
endif ()

if (ENABLE_CUDA)
  set (lulesh_depends
    ${lulesh_depends}
    cuda)
  set (lulesh_compile_options
    --expt-extended-lambda -arch=${CUDA_ARCH})
elseif (ENABLE_HIP)
  set (lulesh_depends
    ${lulesh_depends}
    blt::hip)
  set (lulesh_depends
    ${lulesh_depends}
    blt::hip_runtime)
else()
  set (lulesh_compile_options
    -Wno-unknown-pragmas)
endif ()

blt_add_executable(
  NAME lulesh-v1.0-RAJA-seq.exe
  SOURCES lulesh.cpp
  DEPENDS_ON ${lulesh_depends})

target_compile_definitions(
  lulesh-v1.0-RAJA-seq.exe
  PUBLIC
  -DUSE_CASE=1
)

target_compile_definitions(
  lulesh-v1.0-RAJA-seq.exe
  PUBLIC
  -DLULESH_DEVICE= 
)

target_compile_options(
  lulesh-v1.0-RAJA-seq.exe
  PUBLIC
  ${lulesh_compile_options})

target_include_directories(
  lulesh-v1.0-RAJA-seq.exe
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR})

if (ENABLE_OPENMP)
  blt_add_executable(
    NAME lulesh-v1.0-RAJA-omp.exe
    SOURCES lulesh.cpp
    DEPENDS_ON ${lulesh_depends})

  target_include_directories(
    lulesh-v1.0-RAJA-omp.exe
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

  target_compile_definitions(
    lulesh-v1.0-RAJA-omp.exe
    PUBLIC
    -DUSE_CASE=2)

  target_compile_definitions(
    lulesh-v1.0-RAJA-omp.exe
    PUBLIC
    -DLULESH_DEVICE= 
  )
endif ()

if (ENABLE_CUDA)
  configure_file(
    lulesh.cpp
    lulesh_cuda.cpp COPYONLY)

  blt_add_executable(
    NAME lulesh-v1.0-RAJA-cuda.exe
    SOURCES ${CMAKE_CURRENT_BINARY_DIR}/lulesh_cuda.cpp
    DEPENDS_ON ${lulesh_depends})

  target_compile_options(
    lulesh-v1.0-RAJA-cuda.exe
    PRIVATE
    ${lulesh_compile_options})

  target_compile_definitions(
    lulesh-v1.0-RAJA-cuda.exe
    PUBLIC
    -DUSE_CASE=8
  )

  target_compile_definitions(
    lulesh-v1.0-RAJA-cuda.exe
    PUBLIC
    -DLULESH_DEVICE=__device__
  )

  target_include_directories(
    lulesh-v1.0-RAJA-cuda.exe
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})
endif ()

if (ENABLE_HIP)
  configure_file(
    lulesh.cpp
    lulesh_hip.cpp COPYONLY)

  blt_add_executable(
    NAME lulesh-v1.0-RAJA-hip.exe
    SOURCES ${CMAKE_CURRENT_BINARY_DIR}/lulesh_hip.cpp
    DEPENDS_ON ${lulesh_depends})

  target_compile_options(
    lulesh-v1.0-RAJA-hip.exe
    PRIVATE
    ${lulesh_compile_options})

  target_compile_definitions(
    lulesh-v1.0-RAJA-hip.exe
    PUBLIC
    -DUSE_CASE=10
  )

  target_compile_definitions(
    lulesh-v1.0-RAJA-hip.exe
    PUBLIC
    -DLULESH_DEVICE=__device__
  )

  target_include_directories(
    lulesh-v1.0-RAJA-hip.exe
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})
endif ()
