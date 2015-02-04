dune_install( ns_model.hh problemcreator.hh nswaves.hh ns_model_spec.hh nssmooth.hh thermodynamics.hh passtraits.hh)

set(GRIDTYPE YASPGRID)
set(GRIDDIM  3)
set(POLORDER 2)
set(PROBLEM 2)
set(FLUX 1)
set(DIFFFLUX PRIMALDG)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DFLUX=${FLUX}" )
add_definitions( "-D${DIFFFLUX}" )
add_definitions( "-DPROBLEM=${PROBLEM}" )

include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/navierstokes")

add_executable( navierstokes ${SOURCEONE})
dune_target_link_libraries(navierstokes ${DUNE_LIBS})
set_property(TARGET navierstokes APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=${POLORDER}")

add_executable( navierstokesall ${SOURCEALL})
dune_target_link_libraries(navierstokes ${DUNE_LIBS})
set_property(TARGET navierstokesall APPEND PROPERTY 
  COMPILE_DEFINITIONS "LIMITER")

dune_add_test( navierstokes navierstokesall )
