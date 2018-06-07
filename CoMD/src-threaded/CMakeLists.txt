add_definitions(-DDOUBLE)

set (CoMD_depends
  RAJA)

if (ENABLE_OPENMP)
  set (CoMD_depends
    ${CoMD_depends}
    openmp)
endif ()

if (ENABLE_CUDA)
  set (CoMD_depends
    ${CoMD_depends}
    cuda)
endif ()



set_source_files_properties( CoMD_info.h PROPERTIES GENERATED true)
add_custom_command ( OUTPUT CoMD_info.h
  COMMAND ./generate_info_header CoMD-local ${CMAKE_CXX_COMPILER}" "
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

set_source_files_properties(cmdLineParser.c
  CoMD.c
  decomposition.c
  eam.c
  haloExchange.c
  initAtoms.c
  linkCells.c
  ljForce.c
  mycommand.c
  parallel.c
  performanceTimers.c
  random.c
  timestep.c
  yamlOutput.c
  PROPERTIES
  LANGUAGE CXX)

if (RAJA_USE_CUDA)
  cuda_add_executable(coMD_threaded.exe
    cmdLineParser.c
    CoMD.c
    decomposition.c
    eam.c
    haloExchange.c
    initAtoms.c
    linkCells.c
    ljForce.c
    mycommand.c
    parallel.c
    performanceTimers.c
    random.c
    timestep.c
    yamlOutput.c)
else()
  add_executable(coMD_threaded.exe
    cmdLineParser.c
    CoMD.c
    decomposition.c
    eam.c
    haloExchange.c
    initAtoms.c
    linkCells.c
    ljForce.c
    mycommand.c
    parallel.c
    performanceTimers.c
    random.c
    timestep.c
    yamlOutput.c)
endif()

target_link_libraries(coMD_threaded.exe
  RAJA)