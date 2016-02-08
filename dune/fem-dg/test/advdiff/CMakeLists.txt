dune_install( models.hh problemcreator.hh steppertraits.hh)
     
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/advdiff/")
if(GENERATE_FEMDG_PARAMETER_FILE)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
endif()

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )

if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  message( STATUS "Using old test build system (2.4 compability mode) " )
  add_executable( advdiff main.cc )
  add_executable( advdiffall main.cc )
  set_property(TARGET advdiffall APPEND PROPERTY 
               COMPILE_DEFINITIONS "POLORDER=2" )
  dune_add_test( advdiff advdiffall )
else()
  dune_add_test( NAME advdiff 
                 SOURCES main.cc
                 COMPILE_DEFINITIONS POLORDER=2 )
  dune_add_test( NAME advdiffall 
                 SOURCES main.cc )
endif()

add_code_generate_targets( advdiff advdiffall )
