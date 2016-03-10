include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/advdiff/")

# specify type of grid implemention, dimension and polynomial order
set( GRIDTYPE ALUGRID_CUBE )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DDIMRANGE=1" )

configure_parameter_file( advdiff )
dune_add_test( NAME advdiff 
               SOURCES ../main.cc
               COMPILE_DEFINITIONS "POLORDER=${POLORDER}"
               TIMEOUT 3000
               CMD_ARGS parameters/advdiff )

add_code_generate_targets( advdiff )


if( NOT FEMDG_FAST_TESTBUILD )
  configure_parameter_file( advdiffall )
  dune_add_test( NAME advdiffall 
                 SOURCES ../main.cc
                 TIMEOUT 3000 
                 CMD_ARGS parameters/advdiffall )

  dune_add_test_case( advdiffall polorder1 )
  dune_add_test_case( advdiffall polorder2 )
  dune_add_test_case( advdiffall polorder3 )
  dune_add_test_case( advdiffall polorder4 )

  dune_add_test_case( advdiff pulse )
  dune_add_test_case( advdiff sin )
  dune_add_test_case( advdiff heat )
  dune_add_test_case( advdiff quasi )
  
  add_code_generate_targets( advdiffall )
endif()


