include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/poisson/")

# specify type of grid implemention, dimension and polynmomial order
set( GRIDTYPE ALUGRID_CUBE )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DPADAPTSPACE" )

dune_add_test_case( NAME poisson
                    SOURCES ../main.cc
                    COMPILE_DEFINITIONS "POLORDER=${POLORDER}" )
add_code_generate_targets( poisson )

if( NOT FEMDG_FAST_TESTBUILD )
  dune_add_test_case( NAME poissonall
                      SOURCES ../main.cc )

  add_code_generate_targets( poissonall )
endif()

