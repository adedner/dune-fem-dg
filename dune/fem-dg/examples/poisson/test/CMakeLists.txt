include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/poisson/")
if(GENERATE_FEMDG_PARAMETER_FILE)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
endif()

add_definitions( "-DALUGRID_CUBE" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DDIMRANGE=1" )

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
               COMPILE_DEFINITIONS "POLORDER=2" )
  dune_target_enable_all_packages( poisson )
  dune_add_test( poisson )
else()
  dune_add_test( NAME poisson 
                 SOURCES ../main.cc 
                 COMPILE_DEFINITIONS POLORDER=2 )
endif()
add_code_generate_targets( poisson )
