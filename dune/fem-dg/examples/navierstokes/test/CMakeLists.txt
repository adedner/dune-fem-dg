include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/navierstokes/")
# configure parameter file from parameter.in
configure_parameter_file()

# specify type of grid implemention, dimension and polynmomial order
set( GRIDTYPE ALUGRID_CUBE )
#set( GRIDTYPE SPGRID_COUNT_FLOPS )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
#add_definitions( "-DCOUNT_FLOPS" )

if( NOT FEMDG_FAST_TESTBUILD )
  if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
    add_executable( navierstokesall ../main.cc )
    dune_target_enable_all_packages( navierstokesall )
    dune_add_test( navierstokesall )
  else()
    dune_add_test( NAME navierstokesall 
                   SOURCES ../main.cc )
  endif()
  add_code_generate_targets( navierstokesall )
endif()
