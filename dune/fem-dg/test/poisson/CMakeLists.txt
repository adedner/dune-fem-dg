dune_install( benchmarkproblems.hh models.hh passtraits.hh problemcreator.hh poissonproblem.hh )

set(MAININCLUDE_PATH  "${dune-fem-dg_SOURCE_DIR}/dune/fem-dg/main/")

include_directories(${MAININCLUDE_PATH})
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/poisson/")

                                
set(GRIDTYPE  YASPGRID)
set(POLORDER  2)
set(GRIDDIM  2)
set(DIMRANGE  1)
set(FLUX 1) #2

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" )
add_definitions( "-DDIMRANGE=${DIMRANGE}" )
add_definitions( "-DFLUX=${FLUX}" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DUSE_ASSERT_THROW" ) 

add_executable( poissonall "${MAININCLUDE_PATH}/main.cc" 
                           "${MAININCLUDE_PATH}/main_0.cc" 
                           "${MAININCLUDE_PATH}/main_1.cc"
                           "${MAININCLUDE_PATH}/main_2.cc" 
                           "${MAININCLUDE_PATH}/main_3.cc")
dune_target_link_libraries(poissonall ${DUNE_LIBS})
set_property(TARGET poissonall APPEND PROPERTY 
  COMPILE_DEFINITIONS "LAG")

add_executable( poisson "${MAININCLUDE_PATH}/main.cc" "${MAININCLUDE_PATH}/main_pol.cc")
dune_target_link_libraries(poisson ${DUNE_LIBS})
set_property(TARGET poisson APPEND PROPERTY 
             COMPILE_DEFINITIONS "ONLY_ONE_P;POLORDER=${POLORDER}")

dune_add_test( poisson poissonall )
