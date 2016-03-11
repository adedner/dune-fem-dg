include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/stokes/")

configure_parameter_file()

# specify type of grid implemention, dimension and polynmomial order
set( GRIDTYPE ALUGRID_CUBE )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )

configure_parameter_file( stokes )
dune_add_test( NAME stokes
               SOURCES ../main.cc
               COMPILE_DEFINITIONS "POLORDER=${POLORDER}" 
               CMD_ARGS parameters/stokes)
add_code_generate_targets( stokes )

if( NOT FEMDG_FAST_TESTBUILD )
  configure_parameter_file( stokesall )
  dune_add_test( NAME stokesall 
                 SOURCES ../main.cc 
                 CMD_ARGS parameters/stokesall )
  add_code_generate_targets( stokesall )
endif()


