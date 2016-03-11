include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/navierstokes/")

configure_parameter_file()

# specify type of grid implemention, dimension and polynmomial order
set( GRIDTYPE ALUGRID_CUBE )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )

configure_parameter_file( navierstokes )
dune_add_test( NAME navierstokes
               SOURCES ../main.cc
               COMPILE_DEFINITIONS "POLORDER=${POLORDER}"
               CMD_ARGS parameters/navierstokes )

add_code_generate_targets( navierstokes )

if( NOT FEMDG_FAST_TESTBUILD )
  configure_parameter_file( navierstokesall )
  dune_add_test( NAME navierstokesall 
                 SOURCES ../main.cc 
                 CMD_ARGS parameters/navierstokesall )

  add_code_generate_targets( navierstokesall )
endif()

