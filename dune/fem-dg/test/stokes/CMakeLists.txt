dune_install(models.hh  problemcreator.hh  problems.hh  stokesalgorithm.hh  stokesassembler.hh  tensorhelper.hh)
     
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/stokes/")
configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)

add_definitions( "-DALUGRID_CUBE" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DDGSCHEME" )
add_definitions( "-DONB" )

dune_add_test( NAME stokes
               SOURCES main.cc
               COMPILE_DEFINITIONS POLORDER=2 )
dune_add_test( NAME stokesall 
               SOURCES main.cc )

add_code_generate_targets( stokes stokesall )
