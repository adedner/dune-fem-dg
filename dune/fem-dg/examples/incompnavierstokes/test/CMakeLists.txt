include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/incompnavierstokes/")
configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter )

# default settings for the test case
add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )

if( NOT FEMDG_FAST_TESTBUILD )
  if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
    #message( STATUS "Using old test build system (2.4 compability mode) " )
    #add_executable( incompnavierstokesall ../main.cc )
    #dune_target_enable_all_packages( incompnavierstokesall )
    #dune_add_test( incompnavierstokesall )
  else()
    dune_add_test( NAME incompnavierstokesall
                   SOURCES ../main.cc )
    add_code_generate_targets( incompnavierstokesall )
  endif()
endif()

if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  #message( STATUS "Using old test build system (2.4 compability mode) " )
  #add_executable( incompnavierstokes ../main.cc )
  #set_property(TARGET incompnavierstokes APPEND PROPERTY 
  #             COMPILE_DEFINITIONS "POLORDER=2" )
  #dune_target_enable_all_packages( incompnavierstokes )
  #dune_add_test( incompnavierstokes )
else()
  dune_add_test( NAME incompnavierstokes
                 SOURCES ../main.cc 
                 COMPILE_DEFINITIONS POLORDER=2 )
  add_code_generate_targets( incompnavierstokes )
endif()


