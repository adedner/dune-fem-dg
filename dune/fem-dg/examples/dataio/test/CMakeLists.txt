include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/examples/dataio/")

add_definitions( "-DALUGRID_SIMPLEX" )
add_definitions( "-DGRIDDIM=3" )
add_definitions( "-DUSE_SIONLIB=0" )
add_definitions( "-DPOLORDER=2" )
add_definitions( "-DDIMRANGE=5" )
add_definitions( "-DENABLE_ADAPTIVELEAFINDEXSET_FOR_YASPGRID" ) 

  #dune_add_test_case( NAME dataconvert
  #                    SOURCES ../dataconvert.cc )

dune_add_test_case( NAME checkpointing
                    SOURCES ../main.cc )
