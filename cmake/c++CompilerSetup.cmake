if (NOT DEFINED CMAKE_CXX_STANDARD OR "${CMAKE_CXX_STANDARD} " STREQUAL " ")
  SET(CMAKE_CXX_STANDARD 14)
endif()

cmessage(STATUS "CMAKE CXX Standard: ${CMAKE_CXX_STANDARD}")

add_compile_options(-Werror -Wall -fvisibility=default)