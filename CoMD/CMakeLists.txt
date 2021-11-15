cmake_minimum_required(VERSION 3.1)

project(COMD LANGUAGES CXX)

#
# Initialize the BLT build system
#

set(ENABLE_GTEST Off CACHE BOOL "")
set(ENABLE_TESTS Off CACHE BOOL "")
set(ENABLE_TBB Off CACHE BOOL "")
set(ENABLE_EXAMPLES Off CACHE BOOL "")
set(ENABLE_DOCUMENTATION Off CACHE BOOL "")

#
# Configure RAJA
#

#set(RAJA_ENABLE_TESTS Off CACHE Bool "")
#add_subdirectory(tpl/raja)# EXCLUDE_FROM_ALL)

get_property(RAJA_INCLUDE_DIRS TARGET RAJA PROPERTY INCLUDE_DIRECTORIES)
include_directories(${RAJA_INCLUDE_DIRS})

set(COMD_VERSION "1.0.0")

#configure_file(${PROJECT_SOURCE_DIR}/src/CoMDConfig.h.in
#               ${PROJECT_BINARY_DIR}/include/CoMDConfig.h)


#
# Configure MPI
#

function(notify_package name)
   if (${name}_FOUND)
     message(STATUS "Found ${name} in ${${name}_DIR}.")
   else()
     message(STATUS "${name} package not found.")
   endif()
endfunction()
# Find the MPI library and set some definitions
# This line ensures that we skip C++ headers altogether, avoiding unnecessary symbols in the .o files.
if(ENABLE_MPI)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
   set(MPI_CXX_SKIP_FLAGS
     -DMPI_NO_CPPBIND -DOMPI_SKIP_MPICXX -D_MPICC_H -DMPICH_SKIP_MPICXX -DMPIBULL_SKIP_MPICXX
     CACHE STRING "Default skip-flags used (for SGI, OpenMPI, HP-MPI, MPICH, BULL-MPI)")
   find_package(MPI REQUIRED)
   if (MPI_FOUND)
      set(COMD2_HAVE_MPI TRUE)
      notify_package(MPI)
   endif()
endif()

#
# Define build targets for CoMD
#


include_directories(./src ${PROJECT_BINARY_DIR}/include)

#add_subdirectory(local-update)
#add_subdirectory(symmetric-update)
add_subdirectory(src-threaded)


#add_dependencies(CoMD RAJA)
