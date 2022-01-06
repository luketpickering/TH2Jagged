cmake_minimum_required (VERSION 3.14 FATAL_ERROR)

project(TH2Jagged)
set(TH2Jagged_VERSION 0.9)

if(NOT DEFINED CMAKE_INSTALL_PREFIX OR 
    "${CMAKE_INSTALL_PREFIX}x" STREQUAL "x" OR 
    "${CMAKE_INSTALL_PREFIX}x" STREQUAL "/usr/localx")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
endif()

if(NOT DEFINED CMAKE_BUILD_TYPE OR 
    "${CMAKE_BUILD_TYPE}x" STREQUAL "x")
  set(CMAKE_BUILD_TYPE Release)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/TH2JaggedConfig.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/c++CompilerSetup.cmake)

add_subdirectory(src)

if(DEFINED DO_TEST AND NOT "${DO_TEST} " STREQUAL " ")
  add_compile_definitions(NDEBUG)
  add_compile_options(-Wno-error=unused-but-set-variable -Wno-error=unused-variable)
  add_subdirectory(test)
endif()

configure_file(cmake/setup.sh.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" DESTINATION
  ${CMAKE_INSTALL_PREFIX})

install(FILES
  ${CMAKE_SOURCE_DIR}/cmake/TH2JaggedConfig.cmake DESTINATION
  ${CMAKE_INSTALL_PREFIX}/cmake)

install(EXPORT th2jagged-targets
  FILE TH2JaggedTargets.cmake
  NAMESPACE TH2Jagged::
  DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake
)

FILE(WRITE ${CMAKE_INSTALL_PREFIX}/cmake/TH2JaggedVersion.cmake 
  "SET(TH2Jagged_VERSION ${TH2Jagged_VERSION})\n" )
