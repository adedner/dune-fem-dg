dune_install( benchmarkproblems.hh models.hh passtraits.hh problemcreator.hh poissonproblem.hh )
                                      
set(GRIDTYPE  YASPGRID)
set(POLORDER  2)
set(GRIDDIM  2)
set(DIMRANGE  1)
set(FLUX 1) #2

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DDIMRANGE=${DIMRANGE}" )
add_definitions( "-DFLUX=${FLUX}" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DUSE_ASSERT_THROW" ) 

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/poisson")

add_executable( poissonall ${SOURCEALL})
dune_target_link_libraries(poissonall ${DUNE_LIBS})

add_executable( poisson ${SOURCEONE})
dune_target_link_libraries(poisson ${DUNE_LIBS})
set_property(TARGET poisson APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=${POLORDER}")

dune_add_test( poisson poissonall )
