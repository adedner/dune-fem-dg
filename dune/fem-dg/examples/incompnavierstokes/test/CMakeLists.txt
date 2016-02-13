include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/incompnavierstokes/")
# configure parameter file from parameter.in
configure_parameter_file()

# default settings for the test case
add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )

if( NOT FEMDG_FAST_TESTBUILD )
  if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
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
  # fake empty test
  configure_file( BuildTests.cmake BuildTests.cmake @ONLY)
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


