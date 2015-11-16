dune_install( models.hh  problemcreator.hh  problems.hh )

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/poisson/")
if(GENERATE_FEMDG_PARAMETER_FILE)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
endif()

set(POLORDER 2)
add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DUSE_ASSERT_THROW" ) 

add_executable( poissonall main.cc ) 
dune_target_enable_all_packages( poissonall )
set_property(TARGET poissonall APPEND PROPERTY 
             COMPILE_DEFINITIONS "DGSCHEME;ONB")

add_executable( poisson main.cc )
dune_target_enable_all_packages( poisson )
set_property(TARGET poisson APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=${POLORDER}")

add_code_generate_targets( poisson poissonall )

dune_fem_add_test( poisson poissonall )
