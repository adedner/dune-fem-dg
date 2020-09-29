include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/advdiff/")

# specify type of grid implementation, dimension and polynomial order
set( GRIDTYPE ALUGRID_CUBE )
#set( GRIDTYPE SPGRID )
set( GRIDDIM 2 )
set( POLORDER 4 )

add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DDIMRANGE=3" )

dune_add_test_case( NAME advdiff
                    SOURCES ../main.cc
                    COMPILE_DEFINITIONS "${GRIDTYPE};POLORDER=${POLORDER}"
                    TIMEOUT 3000 )
add_code_generate_targets( advdiff )

## flops counting program
#if( dune-spgrid_FOUND )
#  dune_add_test_case( NAME advdiffflops
#                      SOURCES ../main.cc
#                      COMPILE_DEFINITIONS "SPGRID_COUNT_FLOPS;COUNT_FLOPS;POLORDER=${POLORDER}" )
#  add_code_generate_targets( advdiffflops )
#endif()


if( NOT FEMDG_FAST_TESTBUILD )
  dune_add_test_case( advdiff pulse )
  dune_add_test_case( advdiff sin )
  dune_add_test_case( advdiff heat )
  dune_add_test_case( advdiff quasi )


  dune_add_test_case( NAME advdiffall
                      SOURCES ../main.cc
                      COMPILE_DEFINITIONS "${GRIDTYPE}"
                      TIMEOUT 3000 )
  add_code_generate_targets( advdiffall )
  dune_add_test_case( advdiffall polorder1 )
  dune_add_test_case( advdiffall polorder2 )
  dune_add_test_case( advdiffall polorder3 )
  dune_add_test_case( advdiffall polorder4 )


  dune_add_test_case( NAME advdiff_legendre
                      SOURCES ../main.cc
                      COMPILE_DEFINITIONS "${GRIDTYPE};POLORDER=${POLORDER};DISCRETEFUNCTIONSPACESENUM=2"
                      TIMEOUT 3000 )

  dune_add_test_case( NAME advdiff_hierarchic_legendre
                      SOURCES ../main.cc
                      COMPILE_DEFINITIONS "${GRIDTYPE};POLORDER=${POLORDER};DISCRETEFUNCTIONSPACESENUM=3"
                      TIMEOUT 3000 )

  dune_add_test_case( NAME advdiff_orthonormal
                      SOURCES ../main.cc
                      COMPILE_DEFINITIONS "${GRIDTYPE};POLORDER=${POLORDER};DISCRETEFUNCTIONSPACESENUM=4"
                      TIMEOUT 3000 )

endif()
