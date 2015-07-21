dune_install(problem.hh corner.hh problemcreator.hh stokesassembler.hh)
     
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/stokes/")
configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)

add_definitions( "-DALUGRID_CUBE" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DUSE_ASSERT_THROW" ) 
add_definitions( "-DDGSCHEME" )
add_definitions( "-DONB" )

add_executable( stokesall main.cc )
dune_target_link_libraries(stokesall "${DUNE_LIBS}")

add_executable( stokes main.cc )
dune_target_link_libraries(stokes "${DUNE_LIBS}")
set_property(TARGET stokes APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=2")

add_code_generate_targets( stokes stokesall )

dune_add_test( stokes stokesall )
