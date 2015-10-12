dune_install( ns_model.hh problemcreator.hh nswaves.hh ns_model_spec.hh nssmooth.hh thermodynamics.hh passtraits.hh)

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/incompnavierstokes/")
configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter )

# default settings for the test case
set(POLORDER 2)
add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DPROBLEM=2" )

add_executable( incompnavierstokesall main.cc )
dune_target_link_libraries(incompnavierstokesall ${DUNE_LIBS})

add_executable( incompnavierstokes main.cc )
dune_target_link_libraries(incompnavierstokes ${DUNE_LIBS})
set_property(TARGET incompnavierstokes APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=${POLORDER}")

add_dune_superlu_flags( incompnavierstokes incompnavierstokesall )
add_dune_umfpack_flags( incompnavierstokes incompnavierstokesall )

add_code_generate_targets( incompnavierstokes incompnavierstokesall )

dune_add_test( incompnavierstokes incompnavierstokesall )
