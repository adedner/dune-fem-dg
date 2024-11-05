include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/euler/")

# specify type of grid implemention, dimension and polynomial order

#set( GRIDTYPE YASPGRID )
set( GRIDTYPE ALUGRID_CUBE )
#set( GRIDTYPE SPGRID_COUNT_FLOPS )
#set( GRIDTYPE SPGRID )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
#add_definitions( "-DCOUNT_FLOPS" )
#add_definitions( "-DCODEGEN_SIMD_WIDTH=4" )

set( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc )

dune_add_test_case( NAME euler
                    SOURCES ${MAIN}
                    COMPILE_DEFINITIONS "POLORDER=${POLORDER}" )

add_code_generate_targets( euler )

if( NOT FEMDG_FAST_TESTBUILD )
set(EFFCOMP_POLORD 4)
# efficiency comparison with Python
dune_add_test_case( NAME euler_effcomp
                    SOURCES ${MAIN}
                    COMPILE_DEFINITIONS DGRIDDIM=${GRIDDIM} POLORDER=${EFFCOMP_POLORD} )

# efficiency comparison with Python
dune_add_test_case( NAME euler_effcomp_ufl
                    SOURCES ${MAIN}
                    COMPILE_DEFINITIONS GRIDDIM=${GRIDDIM} POLORDER=${EFFCOMP_POLORD} EULER_WRAPPER_TEST )

  dune_add_test_case( NAME eulerall
                      SOURCES ${MAIN}
                      COMPILE_DEFINITIONS LIMITER )

  add_code_generate_targets( eulerall )
endif()
