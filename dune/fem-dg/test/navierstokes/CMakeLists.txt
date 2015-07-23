dune_install( ns_model.hh problemcreator.hh nswaves.hh ns_model_spec.hh nssmooth.hh thermodynamics.hh passtraits.hh)

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/navierstokes/")
configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter )

# default settings for the test case
set(POLORDER 2)
add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DPROBLEM=2" )

add_executable( navierstokesall main.cc )
dune_target_link_libraries(navierstokesall ${DUNE_LIBS})

add_executable( navierstokes main.cc )
dune_target_link_libraries(navierstokes ${DUNE_LIBS})
set_property(TARGET navierstokes APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=${POLORDER}")

add_code_generate_targets( navierstokes navierstokesall )

dune_add_test( navierstokes navierstokesall )
