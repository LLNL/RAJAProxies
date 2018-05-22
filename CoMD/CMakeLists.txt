cmake_minimum_required(VERSION 3.1)

project(COMD LANGUAGES CXX)

#
# Initialize the BLT build system
#

set(ENABLE_GTEST Off CACHE Bool "")
set(ENABLE_TESTS Off CACHE Bool "")
set(ENABLE_TBB Off CACHE Bool "")
set(ENABLE_EXAMPLES Off CACHE Bool "")
set(ENABLE_DOCUMENTATION Off CACHE Bool "")


#
# Configure RAJA
#

#set(RAJA_ENABLE_TESTS Off CACHE Bool "")
#add_subdirectory(tpl/raja)# EXCLUDE_FROM_ALL)

get_property(RAJA_INCLUDE_DIRS TARGET RAJA PROPERTY INCLUDE_DIRECTORIES)
include_directories(${RAJA_INCLUDE_DIRS})

set(COMD_VERSION "1.0.0")

configure_file(${PROJECT_SOURCE_DIR}/src/CoMDConfig.h.in
               ${PROJECT_BINARY_DIR}/include/CoMDConfig.h)




#
# Define build targets for CoMD
#


include_directories(./src ${PROJECT_BINARY_DIR}/include)

#add_subdirectory(local-update)
#add_subdirectory(symmetric-update)
add_subdirectory(src-threaded)


#add_dependencies(CoMD RAJA)
