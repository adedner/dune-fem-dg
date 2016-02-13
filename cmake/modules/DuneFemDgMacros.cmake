# the main files for all the tests are in the same dir
set(FEMDG_MAIN_DIR "${CMAKE_SOURCE_DIR}/dune/fem-dg/main")

set(OUT_OF_SOURCE_FEMDG_PARAMETER_FILE TRUE)
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  set(OUT_OF_SOURCE_FEMDG_PARAMETER_FILE FALSE)
endif()

# macro for configuring the parameter files from parameter.in
macro(configure_parameter_file)
  set(CURRENT_PARAMETER_PATH .)
  if(OUT_OF_SOURCE_FEMDG_PARAMETER_FILE)
    set(CURRENT_PARAMETER_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
endmacro(configure_parameter_file)

# do a fast test build by default,
# i.e. only build the most important tests
# when calling 'make test' and 'make build_tests', respectively
set(FEMDG_FAST_TESTBUILD ON CACHE BOOL
    "only build the most important tests when calling 'make test'" )

include(Codegen)
include(TargetDistclean)
