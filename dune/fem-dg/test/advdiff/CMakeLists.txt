dune_install( models.hh problemcreator.hh steppertraits.hh)
                                      
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

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/advdiff")

add_executable( advdiffall ${SOURCEALL})
dune_target_link_libraries(advdiffall ${DUNE_LIBS})

add_executable( advdiff ${SOURCEONE})
dune_target_link_libraries(advdiff ${DUNE_LIBS})
set_property(TARGET advdiff APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=${POLORDER}")
           
dune_add_test( advdiff advdiffall )
