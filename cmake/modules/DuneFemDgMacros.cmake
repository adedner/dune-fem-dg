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
  set( default_params "fem.verboserank:0" "fem.prefix:${abbr}data/${paramfile}" "fem.prefix.input:${abbr}" "fem.eoc.outputpath:${abbr}data/${paramfile}" )

  if( "${target}" STREQUAL NAME )
    #First version of this function: we are creating a real new target

    # default directory name for direct call (i.e. withouch testing tools)
    set( TESTCASE_OUTPUT "data" )
    # default parameter name for direct call (i.e. withouch testing tools)
    set( TESTCASE_INPUT "parameter" )
    # copy default parameter file to location of executable
    configure_file(${CMAKE_SOURCE_DIR}/cmake/scripts/parameter.in ${CMAKE_CURRENT_BINARY_DIR}/parameter )
    dune_add_test( ${target} ${paramfile} ${ARGN} CMD_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/parameters/${paramfile} ${default_params} )
  else()
    #Section version of this function: We just append another parameter file to an existing target
    if(NOT TARGET ${target})
      message( ERROR "You have tried to create a test case depending on a non existing target '${target}'." )
    endif()
    add_test( NAME ${base_name}_${paramfile}
              COMMAND ./${target} ${CMAKE_CURRENT_SOURCE_DIR}/parameters/${target}_${paramfile} ${default_params} )
  endif()
endfunction()


# do a fast test build by default,
# i.e. only build the most important tests
# when calling 'make test' and 'make build_tests', respectively
set(FEMDG_FAST_TESTBUILD ON CACHE BOOL
    "only build the most important tests when calling 'make test'" )

include(Codegen)
include(TargetDistclean)
