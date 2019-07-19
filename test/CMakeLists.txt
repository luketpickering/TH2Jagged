SET(TESTS
  BinningTest
  ErrorTest
  WriteTest
  ReadTest
  Test
  TH2DTest
)

include_directories(${CMAKE_SOURCE_DIR}/src)

foreach(appname ${TESTS})

  add_executable(${appname} ${appname}.cxx)

  target_link_libraries(${appname} TH2Jagged)

  if(NOT "${CMAKE_LINK_FLAGS} " STREQUAL " ")
    set_target_properties(${appname} PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
  endif()

  if(NOT "${TH2J_LINK_DIRS} " STREQUAL " ")
    target_link_libraries(${appname} ${TH2J_LINK_DIRS})
  endif()

  if(NOT "${TH2J_DEPEND_LIBS} " STREQUAL " ")
    target_link_libraries(${appname} ${TH2J_DEPEND_LIBS})
  endif()

  install(TARGETS ${appname} DESTINATION test)

endforeach()
