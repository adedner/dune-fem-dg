include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/euler/")

# specify type of grid implemention, dimension and polynomial order
set( GRIDTYPE ALUGRID_CUBE )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )

set( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc )
set( CHORJO ${CMAKE_CURRENT_SOURCE_DIR}/../problems/chorjo.cc )

dune_add_test_case( NAME euler
                    SOURCES ${MAIN} ${CHORJO}
                    COMPILE_DEFINITIONS "POLORDER=${POLORDER}" )

add_code_generate_targets( euler )

if( NOT FEMDG_FAST_TESTBUILD )
  dune_add_test_case( NAME eulerall
                      SOURCES ${MAIN} ${CHORJO}
                      COMPILE_DEFINITIONS LIMITER )

  add_code_generate_targets( eulerall )
endif()



