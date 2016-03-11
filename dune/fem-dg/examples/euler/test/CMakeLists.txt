include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/euler/")

configure_parameter_file()

# specify type of grid implemention, dimension and polynomial order
set( GRIDTYPE ALUGRID_CUBE )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )

set( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc )
set( CHORJO ${CMAKE_CURRENT_SOURCE_DIR}/../problems/chorjo.cc )

configure_parameter_file( euler )
dune_add_test( NAME euler
               SOURCES ${MAIN} ${CHORJO}
               COMPILE_DEFINITIONS "POLORDER=${POLORDER}"
               CMD_ARGS parameters/euler )

add_code_generate_targets( euler )

if( NOT FEMDG_FAST_TESTBUILD )
  configure_parameter_file( eulerall )
  dune_add_test( NAME eulerall
                 SOURCES ${MAIN} ${CHORJO}
                 COMPILE_DEFINITIONS LIMITER
                 CMD_ARGS parameters/eulerall )

  add_code_generate_targets( eulerall )
endif()



