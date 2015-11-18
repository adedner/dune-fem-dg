# the main files for all the tests are in the same dir
set(FEMDG_MAIN_DIR "${CMAKE_SOURCE_DIR}/dune/fem-dg/main")

set(GENERATE_FEMDG_PARAMETER_FILE TRUE)
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  set(GENERATE_FEMDG_PARAMETER_FILE FALSE)
endif()

include(Codegen)
include(TargetDistclean)
