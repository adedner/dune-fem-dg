# set up project
project("dune-fem-dg" C CXX)

#circumvent not building docs
set(BUILD_DOCS 1)

# general stuff
cmake_minimum_required(VERSION 2.8.6)

#find dune-common and set the module path
find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
     "${CMAKE_SOURCE_DIR}/cmake/modules")# make sure our own modules are found

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# $Id: duneproject 5842 2010-01-20 18:48:34Z joe $
# we need the module file to be able to build via dunecontrol
set(EXTRA_DIST dune.module)
set(SUBDIRS  dune m4)
if( BUILD_DOCS)
# TODO: set up documentation tree automatically
#SUBDIRS += doc
endif( BUILD_DOCS )
# don't follow the full GNU-standard
# we need automake 1.5
set(AUTOMAKE_OPTIONS  foreign 1.5)
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-fem=${DUNE_FEM_ROOT}  CXX="${CXX}" CC="${CC}")
# include not needed for CMake
# include $(top_srcdir)/am/top-rules
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})

add_subdirectory(cmake/modules)
# finalize the dune project, e.g., generate config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
