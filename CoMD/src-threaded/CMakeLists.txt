add_definitions(-DDOUBLE)

set (CoMD_depends RAJA)

if (ENABLE_MPI)
  add_definitions(-DENABLE_MPI)
  set (CoMD_depends ${CoMD_depends} MPI)
endif()

if (ENABLE_OPENMP)
  add_definitions(-DENABLE_OPENMP)
  set (CoMD_depends ${CoMD_depends} openmp)
endif ()

if (ENABLE_CUDA)
  add_definitions(-DENABLE_CUDA)
  set (CoMD_depends ${CoMD_depends} cuda)
elseif (ENABLE_HIP)
  add_definitions(-DENABLE_HIP)
  set (CoMD_depends ${CoMD_depends} blt::hip)
  set (CoMD_depends ${CoMD_depends} blt::hip_runtime)
endif ()

set_source_files_properties( CoMD_info.h PROPERTIES GENERATED true)
add_custom_command ( OUTPUT CoMD_info.h
  COMMAND ./generate_info_header CoMD-local ${CMAKE_CXX_COMPILER}" "
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

set_source_files_properties(cmdLineParser.cpp
  CoMD.cpp
  decomposition.cpp
  eam.cpp
  haloExchange.cpp
  initAtoms.cpp
  linkCells.cpp
  ljForce.cpp
  mycommand.cpp
  parallel.cpp
  performanceTimers.cpp
  random.cpp
  timestep.cpp
  yamlOutput.cpp
  PROPERTIES
  LANGUAGE CXX)

blt_add_executable(
  NAME coMD_threaded
  SOURCES
  cmdLineParser.cpp
  CoMD.cpp
  decomposition.cpp
  eam.cpp
  haloExchange.cpp
  initAtoms.cpp
  linkCells.cpp
  ljForce.cpp
  mycommand.cpp
  parallel.cpp
  performanceTimers.cpp
  random.cpp
  timestep.cpp
  yamlOutput.cpp
  DEPENDS_ON ${CoMD_depends})

if(ENABLE_CUDA)
  set (comd_compile_options
    --expt-extended-lambda -arch=${CUDA_ARCH})
target_compile_options(
  coMD_threaded
  PUBLIC
  ${comd_compile_options})
endif()

