include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/advdiff/")

# configure parameter file from parameter.in
configure_parameter_file()

# specify type of grid implementation, dimension and polynomial order
set( GRIDTYPE ALUGRID_CUBE )
#set( GRIDTYPE SPGRID_COUNT_FLOPS )
#set( GRIDTYPE SPGRID )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DDIMRANGE=5" )

if( NOT FEMDG_FAST_TESTBUILD )
  if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
    add_executable( advdiffall ../main.cc )
    set_property(TARGET advdiffall APPEND PROPERTY
      COMPILE_DEFINITIONS "${GRIDTYPE}" )
    dune_target_enable_all_packages( advdiffall )
    dune_add_test( advdiffall )
    set_tests_properties( advdiffall PROPERTIES TIMEOUT 3000 )
  else()
    dune_add_test( NAME advdiffall
                   SOURCES ../main.cc
                   COMPILE_DEFINITIONS "${GRIDTYPE}"
                   TIMEOUT 3000 )
  endif()
  add_code_generate_targets( advdiffall )
endif()

if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  add_executable( advdiff ../main.cc )
  set_property(TARGET advdiff APPEND PROPERTY
    COMPILE_DEFINITIONS "${GRIDTYPE};POLORDER=${POLORDER}" )
  dune_target_enable_all_packages( advdiff )
  dune_add_test( advdiff )
  set_tests_properties( advdiff PROPERTIES TIMEOUT 3000 )

  ## flops counting program
  if( dune-spgrid_FOUND )
    add_executable( advdiffflops ../main.cc )
    set_property(TARGET advdiffflops APPEND PROPERTY
      COMPILE_DEFINITIONS "SPGRID_COUNT_FLOPS;COUNT_FLOPS;POLORDER=${POLORDER}" )
    dune_target_enable_all_packages( advdiffflops )
  endif()
else()
  dune_add_test( NAME advdiff
                 SOURCES ../main.cc
                 COMPILE_DEFINITIONS "${GRIDTYPE};POLORDER=${POLORDER}"
                 TIMEOUT 3000 )
endif()

add_code_generate_targets( advdiff )
if( dune-spgrid_FOUND )
  add_code_generate_targets( advdiffflops )
endif()
