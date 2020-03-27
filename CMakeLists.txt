# set up project
project("dune-fem-dg" C CXX)

set(DUNE_REENABLE_ADD_TEST TRUE)

# this option is mandatory for dune-fem-dg
option(DUNE_GRID_GRIDTYPE_SELECTOR "Grid selector definition added to config.h" ON)

# general stuff
cmake_minimum_required(VERSION 2.8.12)

#find dune-common and set the module path
find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
     "${CMAKE_SOURCE_DIR}/cmake/modules")# make sure our own modules are found

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

include(DuneMPI)

if(dune-python_FOUND)
  set(PYDEMO_DIR pydemo)
endif()

#add sub directories
dune_add_subdirs( dune lib doc cmake/modules ${PYDEMO_DIR})

#add fem-dg test, i.e. tests which are not included in a directory called "test"
#set( _dirnew dune/fem-dg/examples/advdiff_test )
#set_property(DIRECTORY ${_dirnew} PROPERTY TEST_INCLUDE_FILE ${CMAKE_BINARY_DIR}/${_dirnew}/BuildTests.cmake)

make_dependent_modules_sys_included()
add_header_listing()

add_subdirectory("python")
dune_python_install_package(PATH "python")

# finalize the dune project, e.g., generate config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
