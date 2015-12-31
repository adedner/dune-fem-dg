dune_install(models.hh  problemcreator.hh  problems.hh  stokesalgorithm.hh  stokesassembler.hh  tensorhelper.hh)
     
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/stokes/")
if(GENERATE_FEMDG_PARAMETER_FILE)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
endif()

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )

dune_add_test( NAME stokes
               SOURCES main.cc
               COMPILE_DEFINITIONS POLORDER=2 )
dune_add_test( NAME stokesall 
               SOURCES main.cc )

             #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")

add_code_generate_targets( stokes stokesall )
