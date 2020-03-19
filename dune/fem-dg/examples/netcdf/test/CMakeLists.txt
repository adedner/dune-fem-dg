include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/netcdf/")

# specify type of grid implemention, dimension and polynomial order
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

dune_add_test_case( NAME netcdf
                    SOURCES ${MAIN}
                    COMPILE_DEFINITIONS "POLORDER=${POLORDER}"
                    CMAKE_GUARD netcdf_FOUND)
