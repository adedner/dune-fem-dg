dune_install(models.hh  problemcreator.hh  problems.hh  stokesalgorithm.hh  stokesassembler.hh  tensorhelper.hh)
     
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
dune_target_enable_all_packages( stokesall )

add_executable( stokes main.cc )
dune_target_enable_all_packages( stokes )
set_property(TARGET stokes APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=2")

add_code_generate_targets( stokes stokesall )

dune_add_test( stokes stokesall )
