# install these headers
set(nseqdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem-dg/test/navierstokes)
set(nseq_HEADERS  ns_model.hh 
                  problemcreator.hh 
                  nswaves.hh
                  ns_model_spec.hh 
                  nssmooth.hh 
                  thermodynamics.hh 
                  passtraits.hh)

set(GRIDTYPE YASPGRID)
set(GRIDDIM  3)
set(POLORDER 2)
set(PROBLEM 2)
set(FLUX 1)
set(DIFFFLUX PRIMALDG)

set(TESTS navierstokes
          navierstokesall)

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

# We do not want want to build the tests during make all
# but just build them on demand
install(FILES ${nseq_HEADERS} DESTINATION ${nseqdir})
foreach(i ${TESTS})
  add_dune_alugrid_flags(${i})
  add_test(${i} ${i})
endforeach(i ${TESTS})
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
