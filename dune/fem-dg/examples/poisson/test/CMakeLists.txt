include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/poisson/")
# configure parameter file from parameter.in
configure_parameter_file()


# specify type of grid implemention, dimension and polynmomial order
set( GRIDTYPE ALUGRID_CUBE )
set( GRIDDIM 2 )
set( POLORDER 2 )

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DPADAPTSPACE" )


if( NOT FEMDG_FAST_TESTBUILD )
  if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
    add_executable( poissonall ../main.cc )
    dune_target_enable_all_packages( poissonall )
    dune_add_test( poissonall )
  else()
    dune_add_test( NAME poissonall 
                   SOURCES ../main.cc )
  endif()
  add_code_generate_targets( poissonall )
endif()

if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  add_executable( poisson ../main.cc )
  set_property(TARGET poisson APPEND PROPERTY 
               COMPILE_DEFINITIONS "POLORDER=${POLORDER}" )
  dune_target_enable_all_packages( poisson )
  dune_add_test( poisson )
else()
  dune_add_test( NAME poisson 
                 SOURCES ../main.cc 
                 COMPILE_DEFINITIONS "POLORDER=${POLORDER}" )
endif()
add_code_generate_targets( poisson )
