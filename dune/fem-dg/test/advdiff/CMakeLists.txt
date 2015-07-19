dune_install( models.hh problemcreator.hh steppertraits.hh)
     
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/advdiff/")
configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DUSE_ASSERT_THROW" ) 

add_executable( advdiffall main.cc )
dune_target_link_libraries(advdiffall "${DUNE_LIBS}")

add_executable( advdiff main.cc )
dune_target_link_libraries(advdiff "${DUNE_LIBS}")
set_property(TARGET advdiff APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=2")

add_code_generate_targets( advdiff main.cc )
add_code_generate_targets( advdiffall main.cc )

dune_add_test( advdiff advdiffall )
