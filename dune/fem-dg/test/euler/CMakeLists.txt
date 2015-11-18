dune_install( models.hh  problemcreator.hh  problems.hh )

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/euler/")
if(GENERATE_FEMDG_PARAMETER_FILE)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
endif()

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )

dune_add_test( NAME eulerall
               SOURCES main.cc problems/chorjo.cc
               COMPILE_DEFINITIONS LIMITER )

dune_add_test( NAME euler
               SOURCES main.cc problems/chorjo.cc
               COMPILE_DEFINITIONS POLORDER=2 )

add_code_generate_targets( euler eulerall )
