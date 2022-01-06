add_library(TH2Jagged SHARED
  TH2Jagged.cxx
)

set_target_properties(TH2Jagged PROPERTIES 
  PUBLIC_HEADER TH2Jagged.h
  EXPORT_NAME TH2Jagged)

target_link_libraries(TH2Jagged ROOT::ROOT)

ROOT_GENERATE_DICTIONARY(TH2Jagged_dict
            ${CMAKE_CURRENT_SOURCE_DIR}/TH2Jagged.h
            LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/TH2JaggedLinkDef.h
            MODULE TH2Jagged)

install(TARGETS TH2Jagged 
    EXPORT th2jagged-targets
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libTH2Jagged_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/libTH2Jagged.rootmap
    DESTINATION lib/)
