include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/euler/")
if(GENERATE_FEMDG_PARAMETER_FILE)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
endif()

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )

if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  message( STATUS "Using old test build system (2.4 compability mode) " )
  add_executable( euler ../main.cc ../problems/chorjo.cc )
  add_executable( eulerall ../main.cc ../problems/chorjo.cc )
  set_property(TARGET euler APPEND PROPERTY 
               COMPILE_DEFINITIONS "POLORDER=2" )
  dune_target_enable_all_packages( euler eulerall )
  dune_add_test( euler eulerall )
else()
  dune_add_test( NAME eulerall
                 SOURCES ../main.cc ../problems/chorjo.cc
                 COMPILE_DEFINITIONS LIMITER )

  dune_add_test( NAME euler
                 SOURCES ../main.cc ../problems/chorjo.cc
                 COMPILE_DEFINITIONS POLORDER=2 )
endif()

add_code_generate_targets( euler eulerall )
