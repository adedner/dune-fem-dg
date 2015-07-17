dune_install(problem.hh corner.hh problemcreator.hh stokesassembler.hh)
     
set(MAININCLUDE_PATH  "${dune-fem-dg_SOURCE_DIR}/dune/fem-dg/main/")

include_directories(${MAININCLUDE_PATH})
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/stokes/")

configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)

add_definitions( "-DALUGRID_CUBE" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DUSE_ASSERT_THROW" ) 
add_definitions( "-DWANT_ISTL" ) 
add_definitions( "-DDGSCHEME" )
add_definitions( "-DONB" )
add_definitions( "-DSTOKES" )

add_executable( stokesall "${MAININCLUDE_PATH}/main.cc"  
                          "${MAININCLUDE_PATH}/main_0.cc" 
                          "${MAININCLUDE_PATH}/main_1.cc"
                          "${MAININCLUDE_PATH}/main_2.cc" 
                          "${MAININCLUDE_PATH}/main_3.cc")

dune_target_link_libraries(stokesall "${DUNE_LIBS}")

add_executable( stokes "${MAININCLUDE_PATH}/main.cc" "${MAININCLUDE_PATH}/main_pol.cc")
dune_target_link_libraries(stokes "${DUNE_LIBS}")
set_property(TARGET stokes APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=2")

dune_add_test( stokes stokesall )
