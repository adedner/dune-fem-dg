# Module that provides tools for testing the Dune-Fem-DG way.
#
# What is the difference between a 'test' and a 'test case' (in Dune-Fem-Dg)?
# Or: What is the difference between :code:`dune_add_test()` and :code:`dune_add_test_case()`?
#
# Creating a test with the :code:`dune_add_test()` function creates an executable
# and adds this executable to the testing framework. Running the target can be either done via
#
# 1. Using the _testing framework_: call :code:`make build_tests` and calling :code:`make test` or
#
# 2. _Directly_: call :code:`make <target>` and running :code:`./<target>`.
#
# Once specified the test, it is often unclear which parameters (given by a parameter file)
# has to be read to run the test _properly_ and where to write the data.
#
# But what is a 'test case'? A 'test case' is simply said a 'test' which knows
# where to write data and knows the parameter file to run the test properly.
#
# Of course, it is also possible to run :code:`dune_add_test()` and use the
# :code:`CMD_ARG` argument to bind the parameter file 'by hand' to the test.
#
# One drawback is that this parameter file is only added to the testing framework
# and not to a direct call of the target.
#
# Giving up some responsibility for generating tests does not come for free:
# In order to use the :code:`dune_add_test_case()` framework the user has to
# stick to some basic simple rules:
#
# * Write a CMakeList.txt and use the :code:`dune_add_test_case(NAME <target>)` version
#   In this version every parameter which can be added to :ref:`dune_add_test()`
#   can be used. Nevertheless, using :code:`CMD_ARGS` to bind the parameter file
#   to the test is not necessary anymore (and should be avoided...)
#
# * create a folder 'parameters' where the CMakeList.txt is located
#
# * Inside the parameters directory: create a parameter file 'parameter'.
#   This parameter file is called when you call the target directly.
#
# * Inside the parameters directory: create a parameter file :code:`<target>`.
#   This parameter file is called when you call the target via the testing framework.
#
# Optionally, you can add test cases depending on this existing target :code:`<target>`.
# This is done in the following way:
#
# * Use the :code:`dune_add_test_case(<target> <paramfile>)` version where
#   :code:`<target>` is the already existing target.
#
# * Inside the parameters directory: Create a parameter file called <target>_<paramfile>.
#
# All data is written to the directory 'data/<target>/' (testing framework) and
# 'data/data' (direct call)
#
# WARNING: Do not edit or create parameter files called 'parameter' in the directory
# where the executable is located. These file will be overwritten automatically.
# The location parameters/parameter is the proper way to manipulate parameters
# in a parameter file!
#
# .. cmake_variable:: DUNE_FEMDG_FAST_TESTBUILD
#
#    You may set this variable through your opts file or on a per module level (in the toplevel
#    :code:`CMakeLists.txt` to have the Dune build system to build all test builds.
#
#

function(dune_add_test_case target paramfile )
  set( abbr "${CMAKE_CURRENT_SOURCE_DIR}/" )

  if( "${target}" STREQUAL NAME )
    #First version of this function: we are creating a real new target
    set( newTarget "${paramfile}" )
    #set( default_params "fem.verboserank:1" "fem.prefix:${abbr}data/${newTarget}" "fem.prefix.input:${abbr}" "fem.eoc.outputpath:${abbr}data/${newTarget}" )
    set( default_params "testcase.path:${abbr}" )

    # default directory name for direct call (i.e. withouch testing tools)
    set( TESTCASE_OUTPUT "data" )
    # default parameter name for direct call (i.e. withouch testing tools)
    set( TESTCASE_INPUT "parameter" )
    # copy default parameter file to location of executable
    configure_file(${CMAKE_SOURCE_DIR}/cmake/scripts/parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter )
    #dune_add_test( NAME ${newTarget} ${ARGN} CMD_ARGS ${default_params} paramfile:parameters/${newTarget} )
    dune_add_test( NAME ${newTarget} ${ARGN} CMD_ARGS ${default_params} )
  else()
    set( newTarget "${target}_${paramfile}"  )
    #set( default_params "fem.verboserank:1" "fem.prefix:${abbr}data/${newTarget}" "fem.prefix.input:${abbr}" "fem.eoc.outputpath:${abbr}data/${newTarget}" )
    set( default_params "testcase.path:${abbr}" "testcase:${paramfile}" )

    #Section version of this function: We just append another parameter file to an existing target
    if(NOT TARGET ${target})
      message( ERROR "You have tried to create a test case depending on a non existing target '${target}'." )
    endif()
    #add_test( NAME ${newTarget} COMMAND ./${target} ${default_params} paramfile:parameters/${newTarget} )
    add_test( NAME ${newTarget} COMMAND ./${target} ${default_params} )
  endif()
endfunction()


# do a fast test build by default,
# i.e. only build the most important tests
# when calling 'make test' and 'make build_tests', respectively
set(FEMDG_FAST_TESTBUILD ON CACHE BOOL
    "only build the most important tests when calling 'make test'" )


macro(add_header_listing)
    # header
    file( GLOB_RECURSE header "${CMAKE_CURRENT_SOURCE_DIR}/*.hh" "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.ui")
    set( COMMON_HEADER ${header} ${DUNE_HEADERS} )

    # add header of dependent modules for header listing
    foreach(_mod ${ALL_DEPENDENCIES})
        file(GLOB_RECURSE HEADER_LIST "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/*.hh" "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/*.ui")
        list(APPEND COMMON_HEADER ${HEADER_LIST})
    endforeach(_mod DEPENDENCIES)
    #set_source_files_properties(${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1)
    #add_custom_target(common_header SOURCES ${COMMON_HEADER})

    #needed at least one official target for finding header files in source code
    file(WRITE ${CMAKE_BINARY_DIR}/.qtcreator/main.cc "//this is a dummy file for including files into the project file structure of qtcreator\nint main(){}" )
    add_executable( qtcreator ${CMAKE_BINARY_DIR}/.qtcreator/main.cc ${COMMON_HEADER} )
    set_target_properties( qtcreator PROPERTIES LINKER_LANGUAGE CXX)
endmacro(add_header_listing)

macro(make_dependent_modules_sys_included)
    #disable most warnings from dependent modules
    foreach(_mod ${ALL_DEPENDENCIES})
        if(${_mod}_INCLUDE_DIRS)
            foreach( _idir ${${_mod}_INCLUDE_DIRS} )
                add_definitions("-isystem ${_idir}")
            endforeach( _idir )
        endif(${_mod}_INCLUDE_DIRS)
    endforeach(_mod DEPENDENCIES)
endmacro(make_dependent_modules_sys_included)


include(Codegen)
include(TargetDistclean)

find_package(NETCDF QUIET)
