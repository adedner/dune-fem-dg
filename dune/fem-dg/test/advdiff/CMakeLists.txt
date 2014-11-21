# install these headers
set(advdiffdir ${CMAKE_INSTALL_INCLUDEDIR}/test/advdiff)
set(advdiff_HEADERS  models.hh  
                     problemcreator.hh 
                     steppertraits.hh)
                                      
set(GRIDTYPE  YASPGRID)
set(POLORDER  2)
set(GRIDDIM  2)
set(DIMRANGE  1)
set(FLUX 1) #2

set(TESTS advdiff
          advdiffall) 

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

# We do not want want to build the tests during make all
# but just build them on demand
install(FILES ${advdiff_HEADERS} DESTINATION ${advdiffdir})
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
