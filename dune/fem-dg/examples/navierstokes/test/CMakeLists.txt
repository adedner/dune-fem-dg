include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/navierstokes/")
# configure parameter file from parameter.in
configure_parameter_file()

# default settings for the test case
add_definitions( "-DALUGRID_CUBE" )
add_definitions( "-DGRIDDIM=2" )
#add_definitions( "-DUSE_PTHREADS" )

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

if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  add_executable( navierstokes ../main.cc )
  set_property(TARGET navierstokes APPEND PROPERTY 
               COMPILE_DEFINITIONS "POLORDER=2" )
  dune_target_enable_all_packages( navierstokes )
  dune_add_test( navierstokes )
else()
  dune_add_test( NAME navierstokes
                 SOURCES ../main.cc
                 COMPILE_DEFINITIONS POLORDER=2 )
endif()
add_code_generate_targets( navierstokes )
