dune_install( ns_model.hh problemcreator.hh nswaves.hh ns_model_spec.hh nssmooth.hh thermodynamics.hh passtraits.hh)

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/incompnavierstokes/")
configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter )

# default settings for the test case
add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )

if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  message( STATUS "Using old test build system (2.4 compability mode) " )
  add_executable( incompnavierstokes main.cc )
  add_executable( incompnavierstokesall main.cc )
  set_property(TARGET incompnavierstokesall APPEND PROPERTY 
               COMPILE_DEFINITIONS "POLORDER=2" )
  dune_target_enable_all_packages( incompnavierstokes incompnavierstokesall )
  dune_add_test( incompnavierstokes incompnavierstokesall )
else()
  dune_add_test( NAME incompnavierstokesall
                 SOURCES main.cc )

  dune_add_test( NAME incompnavierstokes
                 SOURCES main.cc 
                 COMPILE_DEFINITIONS POLORDER=2 )
endif()

             #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")

add_code_generate_targets( incompnavierstokes incompnavierstokesall )
