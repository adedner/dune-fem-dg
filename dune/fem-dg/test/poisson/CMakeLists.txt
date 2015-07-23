dune_install( benchmarkproblems.hh models.hh passtraits.hh problemcreator.hh poissonproblem.hh )

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/poisson/")
configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
                                
set(POLORDER 2)
add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DUSE_ASSERT_THROW" ) 

add_executable( poissonall main.cc ) 
dune_target_link_libraries(poissonall ${DUNE_LIBS})
set_property(TARGET poissonall APPEND PROPERTY 
             COMPILE_DEFINITIONS "DGSCHEME;ONB")

add_executable( poisson main.cc )
dune_target_link_libraries(poisson ${DUNE_LIBS})
set_property(TARGET poisson APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=${POLORDER}")

add_code_generate_targets( poisson poissonall )

dune_add_test( poisson poissonall )
