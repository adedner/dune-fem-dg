dune_install( models.hh problemcreator.hh steppertraits.hh)
     
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/advdiff/")
if(GENERATE_FEMDG_PARAMETER_FILE)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
endif()

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )

add_executable( advdiffall main.cc )
dune_target_enable_all_packages( advdiffall )

add_executable( advdiff main.cc )
dune_target_enable_all_packages( advdiff )
set_property(TARGET advdiff APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=2")

add_code_generate_targets( advdiff advdiffall )

dune_fem_add_test( advdiff advdiffall )
