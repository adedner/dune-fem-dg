include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/incompnavierstokes/")

configure_parameter_file()

# default settings for the test case
add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )

configure_file( BuildTests.cmake BuildTests.cmake @ONLY)
#configure_parameter_file( incompnavierstokes )
#dune_add_test( NAME incompnavierstokes
#               SOURCES ../main.cc 
#               COMPILE_DEFINITIONS POLORDER=2
#               CMD_ARGS parameters/incompnavierstokes )
#add_code_generate_targets( incompnavierstokes )
#
#if( NOT FEMDG_FAST_TESTBUILD )
#  configure_parameter_file( incompnavierstokesall )
#  dune_add_test( NAME incompnavierstokesall
#                 SOURCES ../main.cc 
#                 CMD_ARGS parameters/incompnavierstokesall )
#  add_code_generate_targets( incompnavierstokesall )
#endif()



