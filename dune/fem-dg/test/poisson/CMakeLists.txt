dune_install( models.hh  problemcreator.hh  problems.hh )

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/poisson/")
if(GENERATE_FEMDG_PARAMETER_FILE)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
endif()

add_definitions( "-DALUGRID_CUBE" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DDIMRANGE=1" )

if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  message( STATUS "Using old test build system (2.4 compability mode) " )
  add_executable( poisson main.cc )
  add_executable( poissonall main.cc )
  set_property(TARGET poissonall APPEND PROPERTY 
               COMPILE_DEFINITIONS "POLORDER=2" )
  dune_target_enable_all_packages( poisson poissonall )
  dune_add_test( poisson poissonall )
else()
  dune_add_test( NAME poissonall 
                 SOURCES main.cc )
  dune_add_test( NAME poisson 
                 SOURCES main.cc 
                 COMPILE_DEFINITIONS POLORDER=2)
endif()

             #set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0" )             

add_code_generate_targets( poisson poissonall )
