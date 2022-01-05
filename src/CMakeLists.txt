ROOT_GENERATE_DICTIONARY(TH2JaggedDict
            ${CMAKE_CURRENT_SOURCE_DIR}/TH2Jagged.h
            LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/TH2JaggedLinkDef.h
            OPTIONS -p)

add_library(TH2Jagged SHARED
  TH2Jagged.cxx
  TH2JaggedDict.cxx
)

set_target_properties(TH2Jagged PROPERTIES 
  PUBLIC_HEADER TH2Jagged.h
  EXPORT_NAME TH2Jagged)

target_link_libraries(TH2Jagged ROOT::ROOT)

install(TARGETS TH2Jagged 
    EXPORT th2jagged-targets
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libTH2JaggedDict_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/libTH2JaggedDict.rootmap
    DESTINATION lib/)
