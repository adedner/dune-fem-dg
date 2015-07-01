dune_install( models.hh problemcreator.hh steppertraits.hh)
     
set(MAININCLUDE_PATH  "${dune-fem-dg_SOURCE_DIR}/dune/fem-dg/main/")

include_directories(${MAININCLUDE_PATH})
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/advdiff/")

configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DUSE_ASSERT_THROW" ) 

add_executable( advdiffall "${MAININCLUDE_PATH}/main.cc"  
                            "${MAININCLUDE_PATH}/main_0.cc" 
                            "${MAININCLUDE_PATH}/main_1.cc"
                            "${MAININCLUDE_PATH}/main_2.cc" 
                            "${MAININCLUDE_PATH}/main_3.cc")

dune_target_link_libraries(advdiffall "${DUNE_LIBS}")

add_executable( advdiff "${MAININCLUDE_PATH}/main.cc" "${MAININCLUDE_PATH}/main_pol.cc")
dune_target_link_libraries(advdiff "${DUNE_LIBS}")
set_property(TARGET advdiff APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=2")

dune_add_test( advdiff advdiffall )
