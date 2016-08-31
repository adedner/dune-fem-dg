include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/advdiff/")

# configure parameter file from parameter.in
configure_parameter_file()

# specify type of grid implementation, dimension and polynomial order
set( GRIDTYPE ALUGRID_CUBE )
#set( GRIDTYPE SPGRID_COUNT_FLOPS )
#set( GRIDTYPE SPGRID )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DDIMRANGE=1" )
#add_definitions( "-DCOUNT_FLOPS" )

if( NOT FEMDG_FAST_TESTBUILD )
  if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
    add_executable( advdiffall ../main.cc )
    dune_target_enable_all_packages( advdiffall )
    dune_add_test( advdiffall )
    set_tests_properties( advdiffall PROPERTIES TIMEOUT 3000 )
  else()
    dune_add_test( NAME advdiffall 
                   SOURCES ../main.cc
                   TIMEOUT 3000 ) 
  endif()
  add_code_generate_targets( advdiffall )
endif()

if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  add_executable( advdiff ../main.cc )
  set_property(TARGET advdiff APPEND PROPERTY 
    COMPILE_DEFINITIONS "POLORDER=${POLORDER}" )
  dune_target_enable_all_packages( advdiff )
  dune_add_test( advdiff )
  set_tests_properties( advdiff PROPERTIES TIMEOUT 3000 )
else()
  dune_add_test( NAME advdiff 
                 SOURCES ../main.cc
                 COMPILE_DEFINITIONS "POLORDER=${POLORDER}"
                 TIMEOUT 3000 )
endif()
add_code_generate_targets( advdiff )
