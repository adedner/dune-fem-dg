# the main files for all the tests are in the same dir
set(FEMDG_MAIN_DIR "${CMAKE_SOURCE_DIR}/dune/fem-dg/main")



# macro for configuring the parameter files from parameter.in
function(configure_parameter_file)
  set(CURRENT_PARAMETER_PATH ${CMAKE_CURRENT_SOURCE_DIR})

  #default value
  set( TESTCASE_OUTPUT "data" )
  set( TESTCASE_INPUT "parameter" )
  
  if( ARGC EQUAL 1 )
    set( TESTCASE_OUTPUT "${ARGV0}" )
    set( TESTCASE_INPUT "${ARGV0}" )
  endif()
 
  if( ARGC EQUAL 2 )
    set( TESTCASE_OUTPUT "${ARGV0}" )
    set( TESTCASE_INPUT "${ARGV1}" )
  endif()

  # write start parameter file
  # This file points to the real parameter file contained in the directory "parameters"
  configure_file(${CMAKE_SOURCE_DIR}/cmake/scripts/parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter )

endfunction(configure_parameter_file)

function(dune_add_test_case target paramfile )
  set( abbr "${CMAKE_CURRENT_SOURCE_DIR}/" )
  set( default_params "fem.verboserank:0" "fem.prefix:${abbr}data/${paramfile}" "fem.prefix.input:${abbr}" "fem.eoc.outputpath:${abbr}data/${paramfile}" )

  if( "${target}" STREQUAL NAME )
    ##copy configure file for target calls without testing tools
    #configure_parameter_file( ${paramfile} )
    configure_parameter_file()
    ##we are creating a real new target
    dune_add_test( ${target} ${paramfile} ${ARGN} CMD_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/parameters/${paramfile} ${default_params} )
  else()
    if(NOT TARGET ${target})
      message( ERROR "You have tried to create a test case depending on a non existing target. The missing target name is '${target}'" )
    endif()
    foreach( testcase ${ARGN})
      add_test( NAME ${base_name}_${testcase} 
                COMMAND ./${target} ${CMAKE_CURRENT_SOURCE_DIR}/parameters/${target}_${testcase} ${default_params} )
    endforeach()
  endif()
endfunction()


# do a fast test build by default,
# i.e. only build the most important tests
# when calling 'make test' and 'make build_tests', respectively
set(FEMDG_FAST_TESTBUILD ON CACHE BOOL
    "only build the most important tests when calling 'make test'" )

include(Codegen)
include(TargetDistclean)
