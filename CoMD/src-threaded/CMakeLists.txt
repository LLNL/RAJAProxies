add_definitions(-DDOUBLE)

set (CoMD_depends RAJA)

if (ENABLE_MPI)
  add_definitions(-DDO_MPI)
  set (CoMD_depends ${CoMD_depends} MPI)
endif()

if (ENABLE_OPENMP)
  add_definitions(-DENABLE_OPENMP)
  set (CoMD_depends ${CoMD_depends} openmp)
endif ()

if (ENABLE_CUDA)
  add_definitions(-DDO_CUDA)
  set (CoMD_depends ${CoMD_depends} cuda)
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

raja_add_executable(
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
  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_BINARY_DIR}/src
    ${CUDA_INCLUDE_PATH}
    ${CUDA_CXXFLAGS})
endif()

if(ENABLE_MPI)
  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_BINARY_DIR}/src
    ${MPI_INCLUDE_PATH}
    ${MPI_CXXFLAGS})
endif()

