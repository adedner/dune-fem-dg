dune_install( dataconvert.cc )

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/dataio/")

if(GENERATE_FEMDG_PARAMETER_FILE)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
endif()

add_definitions( "-DALUGRID_SIMPLEX" )
add_definitions( "-DGRIDDIM=3" )
add_definitions( "-DUSE_SIONLIB=0" )
add_definitions( "-DPOLORDER=2" )
add_definitions( "-DDIMRANGE=5" )
add_definitions( "-DENABLE_ADAPTIVELEAFINDEXSET_FOR_YASPGRID" ) 


if( DUNE_FEM_DG_VERSION VERSION_LESS "2.4.1" )
  message( STATUS "Using old test build system (2.4 compability mode) " )
  add_executable( dataconvert dataconvert.cc )
  add_executable( checkpointing main.cc )
  dune_target_enable_all_packages( dataconvert checkpointing )
  dune_add_test( dataconvert checkpointing )
else()
  dune_add_test( NAME dataconvert
                 SOURCES dataconvert.cc )

  dune_add_test( NAME checkpointing
                 SOURCES main.cc )
endif()
