include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/stokes/")

# specify type of grid implemention, dimension and polynmomial order
set( GRIDTYPE ALUGRID_CUBE )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )

if( NOT FEMDG_FAST_TESTBUILD )
  dune_add_test_case( NAME stokes
                      SOURCES ../main.cc
                      COMPILE_DEFINITIONS "POLORDER=${POLORDER}" )
  add_code_generate_targets( stokes )

  dune_add_test_case( NAME stokesall
                      SOURCES ../main.cc )
  add_code_generate_targets( stokesall )
endif()


