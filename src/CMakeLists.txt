add_library(TH2Jagged SHARED
  TH2Jagged.cxx
)

set_target_properties(TH2Jagged PROPERTIES 
  PUBLIC_HEADER TH2Jagged.h
  EXPORT_NAME TH2Jagged)

target_link_libraries(TH2Jagged ROOT::ROOT)

configure_file(TH2Jagged.h ${CMAKE_CURRENT_BINARY_DIR}/TH2Jagged.h COPYONLY)
configure_file(TH2JaggedLinkdef.h ${CMAKE_CURRENT_BINARY_DIR}/TH2JaggedLinkdef.h COPYONLY)

ROOT_GENERATE_DICTIONARY(TH2Jagged_dict
            TH2Jagged.h
            LINKDEF TH2JaggedLinkDef.h
            MODULE TH2Jagged)

install(TARGETS TH2Jagged 
    EXPORT th2jagged-targets
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libTH2Jagged_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/libTH2Jagged.rootmap
    DESTINATION lib/)