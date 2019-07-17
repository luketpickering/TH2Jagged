include_directories(${INCDIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#rootcling refuses to generate most code if given a -I so we have to copy them here.
configure_file(TH2Jagged.h ${CMAKE_CURRENT_BINARY_DIR}/TH2Jagged.h COPYONLY)
configure_file(TH2JaggedLinkdef.h ${CMAKE_CURRENT_BINARY_DIR}/TH2JaggedLinkdef.h COPYONLY)

add_custom_command(
  OUTPUT "TH2JaggedDict.cxx" "TH2JaggedDict_rdict.pcm"
  COMMAND rootcling
  ARGS -f TH2JaggedDict.cxx TH2Jagged.h TH2JaggedLinkdef.h
  DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/TH2Jagged.h
    ${CMAKE_CURRENT_BINARY_DIR}/TH2JaggedLinkdef.h)

add_library(TH2Jagged SHARED
  TH2Jagged.cxx
  TH2JaggedDict.cxx
)

install(TARGETS TH2Jagged DESTINATION lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/TH2JaggedDict_rdict.pcm DESTINATION lib)

if(NOT "${CMAKE_LINK_FLAGS} " STREQUAL " ")
  set_target_properties(TH2Jagged PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
endif()

if(NOT "${TH2J_LINK_DIRS} " STREQUAL " ")
  target_link_libraries(TH2Jagged ${TH2J_LINK_DIRS})
endif()

if(NOT "${TH2J_DEPEND_LIBS} " STREQUAL " ")
  target_link_libraries(TH2Jagged ${TH2J_DEPEND_LIBS})
endif()
