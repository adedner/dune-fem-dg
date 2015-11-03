dune_install( dataconvert.cc )

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/dataio/")

set(GRIDTYPE YASPGRID)
set(GRIDDIM  3)
set(POLORDER  2)
set(DIMRANGE 5)
set(USE_SIONLIB 0)

configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=3" )
add_definitions( "-DUSE_SIONLIB=0" )
add_definitions( "-DPOLORDER=2" )
add_definitions( "-DDIMRANGE=5" )
add_definitions( "-DENABLE_ADAPTIVELEAFINDEXSET_FOR_YASPGRID" ) 

add_executable( dataconvert dataconvert.cc )
add_executable( checkpointing main.cc )

dune_target_link_libraries( checkpointing "${DUNE_LIBS}")
dune_target_link_libraries( dataconvert "${DUNE_LIBS}")

add_dune_superlu_flags( dataconvert checkpointing )
add_dune_umfpack_flags( dataconvert checkpointing )

#set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g" )

dune_add_test( checkpointing dataconvert )

