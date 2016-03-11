include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/poisson/")

configure_parameter_file()

# specify type of grid implemention, dimension and polynmomial order
set( GRIDTYPE ALUGRID_CUBE )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DPADAPTSPACE" )

configure_parameter_file( poisson )
dune_add_test( NAME poisson 
               SOURCES ../main.cc 
               COMPILE_DEFINITIONS "POLORDER=${POLORDER}" 
               CMD_ARGS parameters/poisson )
add_code_generate_targets( poisson )

if( NOT FEMDG_FAST_TESTBUILD )
  configure_parameter_file( poissonall )
  dune_add_test( NAME poissonall 
                 SOURCES ../main.cc
                 CMD_ARGS parameters/poissonall )
  
  add_code_generate_targets( poissonall )
endif()

