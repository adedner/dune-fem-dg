include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/incompnavierstokes/")

# default settings for the test case
add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )

#set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0" )

#configure_file( BuildTests.cmake BuildTests.cmake @ONLY)
#dune_add_test_case( NAME incompnavierstokes
#                    SOURCES ../main.cc
#                    COMPILE_DEFINITIONS POLORDER=2 )
#add_code_generate_targets( incompnavierstokes )



#if( NOT FEMDG_FAST_TESTBUILD )
#  dune_add_test_case( NAME incompnavierstokesall
#                      SOURCES ../main.cc )
#  add_code_generate_targets( incompnavierstokesall )
#endif()



