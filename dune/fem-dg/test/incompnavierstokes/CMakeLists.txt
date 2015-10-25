dune_install( ns_model.hh problemcreator.hh nswaves.hh ns_model_spec.hh nssmooth.hh thermodynamics.hh passtraits.hh)

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/incompnavierstokes/")
if(EXISTS parameter.in)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter )
endif()

# default settings for the test case
set(POLORDER 2)
add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DPROBLEM=2" )

add_executable( incompnavierstokesall main.cc )
dune_target_enable_all_packages( incompnavierstokesall )

add_executable( incompnavierstokes main.cc )
dune_target_enable_all_packages( incompnavierstokes )
set_property(TARGET incompnavierstokes APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=${POLORDER}")

add_code_generate_targets( incompnavierstokes incompnavierstokesall )

dune_fem_add_test( incompnavierstokes incompnavierstokesall )
