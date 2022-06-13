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

set (lulesh_sources
  lulesh.cpp
  lulesh-comm.cpp
  lulesh-init.cpp
  lulesh-util.cpp
  lulesh-viz.cpp)

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
  NAME lulesh-v2.0-RAJA-seq.exe
  SOURCES ${lulesh_sources}
  DEPENDS_ON ${lulesh_depends})

target_compile_options(
  lulesh-v2.0-RAJA-seq.exe
  PRIVATE
  ${lulesh_compile_options})

target_compile_definitions(
  lulesh-v2.0-RAJA-seq.exe
  PRIVATE 
  -DUSE_MPI=0)

target_compile_definitions(
  lulesh-v2.0-RAJA-seq.exe
  PRIVATE
  -DUSE_CASE=1)

target_compile_definitions(
  lulesh-v2.0-RAJA-seq.exe
  PUBLIC
  -DLULESH_DEVICE= 
)

if (ENABLE_OPENMP)
  blt_add_executable(
    NAME lulesh-v2.0-RAJA-omp.exe
    SOURCES ${lulesh_sources}
    DEPENDS_ON ${lulesh_depends})

  target_compile_definitions(
    lulesh-v2.0-RAJA-omp.exe
    PRIVATE 
    -DUSE_MPI=0)

  target_compile_definitions(
    lulesh-v2.0-RAJA-omp.exe
    PUBLIC
    -DLULESH_DEVICE= 
  )
endif ()

if (ENABLE_CUDA)
  configure_file(
    lulesh.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-cuda.cpp COPYONLY)

  configure_file(
    lulesh-comm.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-comm-cuda.cpp COPYONLY)
  configure_file(
    lulesh-init.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-init-cuda.cpp COPYONLY)
  configure_file(
    lulesh-util.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-util-cuda.cpp COPYONLY)
  configure_file(
    lulesh-viz.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-viz-cuda.cpp COPYONLY)

  blt_add_executable(
    NAME lulesh-v2.0-RAJA-cuda.exe
    SOURCES 
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-cuda.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-comm-cuda.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-init-cuda.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-util-cuda.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-viz-cuda.cpp
    DEPENDS_ON ${lulesh_depends})

  target_compile_options(
    lulesh-v2.0-RAJA-cuda.exe
    PRIVATE
    ${lulesh_compile_options})

  target_compile_definitions(
    lulesh-v2.0-RAJA-cuda.exe
    PRIVATE
    -DUSE_CASE=9)

  target_compile_definitions(
    lulesh-v2.0-RAJA-cuda.exe
    PRIVATE 
    -DUSE_MPI=0)

  target_compile_definitions(
    lulesh-v2.0-RAJA-cuda.exe
    PUBLIC
    -DLULESH_DEVICE=__device__
  )

  target_include_directories(
    lulesh-v2.0-RAJA-cuda.exe
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})
endif ()

if (ENABLE_HIP)
  configure_file(
    lulesh.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-hip.cpp COPYONLY)

  configure_file(
    lulesh-comm.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-comm-hip.cpp COPYONLY)
  configure_file(
    lulesh-init.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-init-hip.cpp COPYONLY)
  configure_file(
    lulesh-util.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-util-hip.cpp COPYONLY)
  configure_file(
    lulesh-viz.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-viz-hip.cpp COPYONLY)

  blt_add_executable(
    NAME lulesh-v2.0-RAJA-hip.exe
    SOURCES 
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-hip.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-comm-hip.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-init-hip.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-util-hip.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/lulesh-viz-hip.cpp
    DEPENDS_ON ${lulesh_depends})

  target_compile_options(
    lulesh-v2.0-RAJA-hip.exe
    PRIVATE
    ${lulesh_compile_options})

  target_compile_definitions(
    lulesh-v2.0-RAJA-hip.exe
    PRIVATE
    -DUSE_CASE=10)

  target_compile_definitions(
    lulesh-v2.0-RAJA-hip.exe
    PRIVATE 
    -DUSE_MPI=0)

  target_compile_definitions(
    lulesh-v2.0-RAJA-hip.exe
    PUBLIC
    -DLULESH_DEVICE=__device__
  )

  target_include_directories(
    lulesh-v2.0-RAJA-hip.exe
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})
endif ()
