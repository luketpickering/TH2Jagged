cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

project(TH2Jagged)

include(${CMAKE_SOURCE_DIR}/cmake/StringAndListUtils.cmake)

set(BUILD_DEBUG_MSGS TRUE)
include(${CMAKE_SOURCE_DIR}/cmake/cmessage.cmake)

set(VERBOSE TRUE)

#Changes default install path to be a subdirectory of the build dir.
#Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE DEBUG)
elseif(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE DEBUG)
endif()

################################################################################
#                            Check Dependencies
################################################################################

SET(EXTRA_CXX_FLAGS)
SET(EXTRA_LINK_DIRS)
SET(EXTRA_LIBS)

##################################  ROOT  ######################################
include(${CMAKE_SOURCE_DIR}/cmake/ROOTSetup.cmake)

################################## COMPILER ####################################
include(${CMAKE_SOURCE_DIR}/cmake/c++CompilerSetup.cmake)

add_subdirectory(src)

if(DEFINED DO_TEST AND NOT "${DO_TEST} " STREQUAL " ")
  add_subdirectory(test)
endif()

configure_file(cmake/setup.sh.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" DESTINATION
  ${CMAKE_INSTALL_PREFIX})
