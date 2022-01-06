if(NOT "${CMAKE_PROJECT_NAME} " STREQUAL "TH2Jagged ")
  get_filename_component(TH2Jagged_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
endif()

include(CMakeFindDependencyMacro)

if(NOT COMMAND find_package_or_dependency)
  macro(find_package_or_dependency ARGS)
    if("${CMAKE_PROJECT_NAME} " STREQUAL "TH2Jagged ")
      find_package(${ARGS})
    else()
      find_dependency(${ARGS})
    endif()
  endmacro()
endif()

if(NOT COMMAND cmessage)
  if(NOT WIN32)
    string(ASCII 27 Esc)
    set(CM_ColourReset "${Esc}[m")
    set(CM_ColourBold "${Esc}[1m")
    set(CM_Red "${Esc}[31m")
    set(CM_Green "${Esc}[32m")
    set(CM_Yellow "${Esc}[33m")
    set(CM_Blue "${Esc}[34m")
    set(CM_Magenta "${Esc}[35m")
    set(CM_Cyan "${Esc}[36m")
    set(CM_White "${Esc}[37m")
    set(CM_BoldRed "${Esc}[1;31m")
    set(CM_BoldGreen "${Esc}[1;32m")
    set(CM_BoldYellow "${Esc}[1;33m")
    set(CM_BoldBlue "${Esc}[1;34m")
    set(CM_BoldMagenta "${Esc}[1;35m")
    set(CM_BoldCyan "${Esc}[1;36m")
    set(CM_BoldWhite "${Esc}[1;37m")
  endif()

  message(STATUS "Setting up colored messages...")

  function(cmessage)
    list(GET ARGV 0 MessageType)
    if(MessageType STREQUAL FATAL_ERROR OR MessageType STREQUAL SEND_ERROR)
      list(REMOVE_AT ARGV 0)
      message(${MessageType} "${CM_BoldRed}${ARGV}${CM_ColourReset}")
    elseif(MessageType STREQUAL WARNING)
      list(REMOVE_AT ARGV 0)
      message(${MessageType} "${CM_BoldYellow}${ARGV}${CM_ColourReset}")
    elseif(MessageType STREQUAL AUTHOR_WARNING)
      list(REMOVE_AT ARGV 0)
      message(${MessageType} "${CM_BoldCyan}${ARGV}${CM_ColourReset}")
    elseif(MessageType STREQUAL STATUS)
      list(REMOVE_AT ARGV 0)
      message(${MessageType} "${CM_Green}[INFO]:${CM_ColourReset} ${ARGV}")
    elseif(MessageType STREQUAL CACHE)				
      list(REMOVE_AT ARGV 0)
      message(-- "${CM_Blue}[CACHE]:${CM_ColourReset} ${ARGV}")
    elseif(MessageType STREQUAL DEBUG)
      list(REMOVE_AT ARGV 0)
      if(BUILD_DEBUG_MSGS)
        message("${CM_Magenta}[DEBUG]:${CM_ColourReset} ${ARGV}")
      endif()
    else()
      message(${MessageType} "${CM_Green}[INFO]:${CM_ColourReset} ${ARGV}")
    endif()
  endfunction()
endif()

##################################  ROOT  ######################################
if(NOT TARGET ROOT::ROOT) # Only ROOT if we really need to
  find_package_or_dependency(ROOT REQUIRED)
  include(${ROOT_USE_FILE})

  STRING(STRIP ROOT_CXX_FLAGS ${ROOT_CXX_FLAGS})
  STRING(REPLACE " " ";" ROOT_CXX_FLAGS ${ROOT_CXX_FLAGS})

  if("${CMAKE_PROJECT_NAME} " STREQUAL "TH2Jagged ")
    list (FIND ROOT_CXX_FLAGS "-std=c++11" CPP11_INDEX)
    list (FIND ROOT_CXX_FLAGS "-std=c++14" CPP14_INDEX)
    list (FIND ROOT_CXX_FLAGS "-std=c++17" CPP17_INDEX)
    list (FIND ROOT_CXX_FLAGS "-std=c++20" CPP20_INDEX)
    if (${CPP17_INDEX} GREATER -1)
      SET(CMAKE_CXX_STANDARD 17)
    elseif (${CPP20_INDEX} GREATER -1)
      SET(CMAKE_CXX_STANDARD 20)
    else()
      SET(CMAKE_CXX_STANDARD 14)
    endif()
  endif()
  list(FILTER ROOT_CXX_FLAGS EXCLUDE REGEX "-std=.*")

  execute_process (COMMAND root-config
    --version OUTPUT_VARIABLE ROOT_CONFIG_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

  #We should let CMake set this
  list(REMOVE_ITEM ROOT_CXX_FLAGS "-fPIC")

  add_library(ROOT::ROOT INTERFACE IMPORTED)

  LIST(APPEND ROOT_LIB_NAMES 
    Core
    RIO
    Hist
    MathCore)

  set(ROOT_LIBRARIES)
  foreach(LN ${ROOT_LIB_NAMES})
    if(NOT DEFINED ROOT_${LN}_LIBRARY)
      cmessage(FATAL_ERROR "ROOTConfig.cmake did not define ROOT_${LN}_LIBRARY, which is required.")
    endif()
    LIST(APPEND ROOT_LIBRARIES ${ROOT_${LN}_LIBRARY})
  endforeach()

  set_target_properties(ROOT::ROOT PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIRS}"
      INTERFACE_COMPILE_OPTIONS "${ROOT_CXX_FLAGS}"
      INTERFACE_LINK_LIBRARIES "${ROOT_LIBRARIES}"
  )
  cmessage(STATUS "Built ROOT::ROOT Imported target")
  cmessage(STATUS "        ROOT_INCLUDE_DIRS: ${ROOT_INCLUDE_DIRS}: ")
  cmessage(STATUS "        ROOT_CXX_FLAGS: ${ROOT_CXX_FLAGS}")
  cmessage(STATUS "        ROOT_LIBRARIES: ${ROOT_LIBRARIES}")
endif()

##########################  For External Users Only  ###############################

if(NOT "${CMAKE_PROJECT_NAME} " STREQUAL "TH2Jagged ")
  include(${TH2Jagged_CMAKE_DIR}/TH2JaggedTargets.cmake)
  include(${TH2Jagged_CMAKE_DIR}/TH2JaggedVersion.cmake)

  find_path(TH2Jagged_INCLUDE_DIR
    NAMES TH2Jagged.h
    PATHS ${TH2Jagged_CMAKE_DIR}/../include
  )

  find_path(TH2Jagged_LIB_DIR
    NAMES libTH2Jagged.so
    PATHS ${TH2Jagged_CMAKE_DIR}/../lib
  )

  find_path(TH2Jagged_PREFIX
    NAMES setup.sh
    PATHS ${TH2Jagged_CMAKE_DIR}/../
  )

  cmessage(STATUS "TH2Jagged_INCLUDE_DIR: ${TH2Jagged_INCLUDE_DIR}")
  cmessage(STATUS "TH2Jagged_LIB_DIR: ${TH2Jagged_LIB_DIR}")
  cmessage(STATUS "TH2Jagged_PREFIX: ${TH2Jagged_PREFIX}")
  cmessage(STATUS "TH2Jagged_VERSION: ${TH2Jagged_VERSION}")

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(TH2Jagged
      REQUIRED_VARS 
        TH2Jagged_INCLUDE_DIR 
        TH2Jagged_LIB_DIR
        TH2Jagged_PREFIX
      VERSION_VAR 
        TH2Jagged_VERSION
  )
  if(NOT TARGET TH2Jagged::All)
      add_library(TH2Jagged::All INTERFACE IMPORTED)
      set_target_properties(TH2Jagged::All PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES ${TH2Jagged_INCLUDE_DIR}
          INTERFACE_LINK_LIBRARIES TH2Jagged::TH2Jagged
      )
  endif()
endif()