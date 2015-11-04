dune_install( models.hh problemcreator.hh steppertraits.hh)
     
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/advdiff/")
configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )

dune_add_test( NAME advdiff 
               SOURCES main.cc
               COMPILE_DEFINITIONS POLORDER=2 )
dune_add_test( NAME advdiffall 
               SOURCES main.cc )

add_code_generate_targets( advdiff advdiffall )
