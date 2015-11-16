dune_install( problemcreator.hh combinedsteppertraits.hh)

add_definitions( "-DYASPGRID" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DDIMRANGE=1" )
add_definitions( "-DFLUX=1" )
add_definitions( "-DPRIMALDG" )
add_definitions( "-DUSE_ASSERT_THROW" ) 
add_definitions( "-DCOMBINEDSTEPPER" ) 

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O0 -g" )

set( ${DUNE_FEM_DG_ROOT}/dune/fem-dg/main)

set(MAININCLUDE_PATH  "${dune-fem-dg_SOURCE_DIR}/dune/fem-dg/main/")

include_directories(${MAININCLUDE_PATH})
include_directories("${CMAKE_SOURCE_DIR}/dune/fem-dg/test/combined/")

if(GENERATE_FEMDG_PARAMETER_FILE)
  configure_file(parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter)
endif()

add_executable(combinedall "${MAININCLUDE_PATH}/main.cc" 
                            "${MAININCLUDE_PATH}/main_0.cc" 
                            "${MAININCLUDE_PATH}/main_1.cc"
                            "${MAININCLUDE_PATH}/main_2.cc" 
                            "${MAININCLUDE_PATH}/main_3.cc")
dune_target_link_libraries(combinedall ${DUNE_LIBS})
add_dune_papi_flags(combinedall)
add_executable(combined "${MAININCLUDE_PATH}/main.cc" 
                            "${MAININCLUDE_PATH}/main_pol.cc")
dune_target_link_libraries(combined ${DUNE_LIBS})
add_dune_papi_flags(combined)

set_property(TARGET combined APPEND PROPERTY 
  COMPILE_DEFINITIONS "POLORDER=${POLORDER};" )

add_dune_alugrid_flags( combined combinedall )


