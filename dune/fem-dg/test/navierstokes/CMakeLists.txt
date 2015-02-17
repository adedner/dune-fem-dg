dune_install( ns_model.hh problemcreator.hh nswaves.hh ns_model_spec.hh nssmooth.hh thermodynamics.hh passtraits.hh)

set(MAININCLUDE_PATH  "${dune-fem-dg_SOURCE_DIR}/dune/fem-dg/main/")

include_directories(${MAININCLUDE_PATH})
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/navierstokes/")


set(GRIDTYPE YASPGRID)
set(GRIDDIM  3)
set(POLORDER 2)
set(PROBLEM 2)
set(FLUX 1)
set(DIFFFLUX PRIMALDG)

configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DFLUX=${FLUX}" )
add_definitions( "-D${DIFFFLUX}" )
add_definitions( "-DPROBLEM=${PROBLEM}" )

add_executable( navierstokesall "${MAININCLUDE_PATH}/main.cc" 
                                "${MAININCLUDE_PATH}/main_0.cc" 
                                "${MAININCLUDE_PATH}/main_1.cc"
                                "${MAININCLUDE_PATH}/main_2.cc" 
                                "${MAININCLUDE_PATH}/main_3.cc")
dune_target_link_libraries(navierstokesall ${DUNE_LIBS})

add_executable( navierstokes "${MAININCLUDE_PATH}/main.cc" "${MAININCLUDE_PATH}/main_pol.cc" )
dune_target_link_libraries(navierstokes ${DUNE_LIBS})
set_property(TARGET navierstokes APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=${POLORDER}")


dune_add_test( navierstokes navierstokesall )
