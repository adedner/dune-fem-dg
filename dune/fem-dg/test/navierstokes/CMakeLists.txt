dune_install( ns_model.hh problemcreator.hh ns_model_spec.hh thermodynamics.hh passtraits.hh)

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/navierstokes/")
if(GENERATE_FEMDG_PARAMETER_FILE)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter )
endif()

# default settings for the test case
add_definitions( "-DALUGRID_SIMPLEX" )
add_definitions( "-DGRIDDIM=2" )

dune_add_test( NAME navierstokes
               SOURCES main.cc
               COMPILE_DEFINITIONS POLORDER=2 )
dune_add_test( NAME navierstokesall 
               SOURCES main.cc )

             #set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0" )

add_code_generate_targets( navierstokes navierstokesall )
