dune_install(models.hh  problemcreator.hh  problems.hh  stokesalgorithm.hh  stokesassembler.hh  tensorhelper.hh)
     
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/stokes/")
if(GENERATE_FEMDG_PARAMETER_FILE)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
endif()

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )

if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  message( STATUS "Using old test build system (2.4 compability mode) " )
  add_executable( stokes main.cc )
  add_executable( stokesall main.cc )
  set_property(TARGET stokesall APPEND PROPERTY 
               COMPILE_DEFINITIONS "POLORDER=2" )
  dune_target_enable_all_packages( stokes stokesall )
  dune_add_test( stokes stokesall )
else() 
  dune_add_test( NAME stokes
                 SOURCES main.cc
                 COMPILE_DEFINITIONS POLORDER=2 )
  dune_add_test( NAME stokesall 
                 SOURCES main.cc )
endif()

             #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")

add_code_generate_targets( stokes stokesall )
