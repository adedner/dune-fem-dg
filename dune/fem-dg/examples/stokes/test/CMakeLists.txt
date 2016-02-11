include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/stokes/")
# configure parameter file from parameter.in
configure_parameter_file()

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )


if( NOT FEMDG_FAST_TESTBUILD )
  if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
    add_executable( stokesall ../main.cc )
    dune_target_enable_all_packages( stokesall )
    dune_add_test( stokesall )
  else() 
    dune_add_test( NAME stokesall 
                   SOURCES ../main.cc )
  endif()
  add_code_generate_targets( stokesall )
endif()

if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  add_executable( stokes ../main.cc )
  set_property(TARGET stokes APPEND PROPERTY 
               COMPILE_DEFINITIONS "POLORDER=2" )
  dune_target_enable_all_packages( stokes )
  dune_add_test( stokes )
else() 
  dune_add_test( NAME stokes
                 SOURCES ../main.cc
                 COMPILE_DEFINITIONS POLORDER=2 )
endif()

add_code_generate_targets( stokes )

