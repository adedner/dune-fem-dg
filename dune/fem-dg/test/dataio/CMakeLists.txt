set(MAININCLUDE_PATH  "${dune-fem-dg_SOURCE_DIR}/dune/fem-dg/main/")

include_directories(${MAININCLUDE_PATH})
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/dataio/")


set(GRIDTYPE YASPGRID)
set(GRIDDIM  3)
set(POLORDER  2)
set(DIMRANGE 5)
set(USE_SIONLIB 0)


add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DUSE_SIONLIB=${USE_SIONLIB}" )
add_definitions( "-DPOLORDER=${POLORDER}" )
add_definitions( "-DONLY_ONE_P" )
add_definitions( "-DDIMRANGE=${DIMRANGE}" )
add_definitions( "-DENABLE_ADAPTIVELEAFINDEXSET_FOR_YASPGRID" ) 

add_executable( dataconvert dataconvert.cc )

add_executable( checkpointing "${MAININCLUDE_PATH}/main.cc" 
                             "${MAININCLUDE_PATH}/main_pol.cc" )


dune_target_link_libraries( checkpointing "${DUNE_LIBS}")
dune_target_link_libraries( dataconvert "${DUNE_LIBS}")

dune_add_test( advdiff advdiffall )

