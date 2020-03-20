include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/navierstokes/")

# specify type of grid implemention, dimension and polynmomial order
set( GRIDTYPE ALUGRID_CUBE )
#set( GRIDTYPE SPGRID_COUNT_FLOPS )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
#add_definitions( "-DCOUNT_FLOPS" )

if( NOT FEMDG_FAST_TESTBUILD )
  dune_add_test_case( NAME navierstokes
                      SOURCES ../main.cc
                      COMPILE_DEFINITIONS "POLORDER=${POLORDER}" )
  add_code_generate_targets( navierstokes )

  dune_add_test_case( NAME navierstokesall
                      SOURCES ../main.cc )

  add_code_generate_targets( navierstokesall )
endif()
