#!/bin/bash

WORKDIR=`pwd`

echo "This script will download and build all DUNE modules"
echo "necessary to run the examples in dune-fem-dg."
echo
echo "The installation directory is: $WORKDIR"
echo "Some third party libraries have to be downloaded manually."
echo "Please take a look at this script for parameters and options."
echo
read -p "Install DUNE modules to $WORKDIR? (Y/N) " YN
if [ "$YN" != "Y" ] ;then
  exit 1
fi

# this script downloads the necessary set of DUNE modules
# to build and run the examples in dune-fem-dg
# NOTE: Zoltan has to be downloaded separately from

#change appropriately, i.e. 2.4 or empty (which refers to master)
DUNEVERSION=2.4

# your favorite compiler optimization flags
FLAGS="-O3 -DNDEBUG"
MAKE_FLAGS="-j4"

# PETSc from http://www.mcs.anl.gov/petsc/
# set the environment variable PETSC_DIR for PETSc to be found
# export PETSC_DIR=$WORKDIR/petsc

# Zoltan has to be downloaded from
# http://www.cs.sandia.gov/zoltan/
#WITH_ZOLTAN="-DZOLTAN_ROOT=$WORKDIR/zoltan"
WITH_ZOLTAN=

# SIONlib 1.5p1 (or later) has to be downloaded from
# http://www.fz-juelich.de/ias/jsc/EN/Expertise/Support/Software/SIONlib/sionlib-download_node
#WITH_SIONLIB="-DSIONLIB_ROOT=$WORKDIR/sionlib"
WITH_SIONLIB=

# dune modules needed to build dune-fem-dg
DUNEMODULES="dune-common dune-geometry dune-grid dune-istl dune-alugrid dune-fem dune-fem-dg"

# build dir for cmake (. for in-source build or empty then this will default to build-cmake)
BUILDDIR=.

# build flags for all DUNE modules
# change according to your needs
if ! test -f config.opts ; then
  echo "MAKE_FLAGS=\"$MAKE_FLAGS\" \\
  USE_CMAKE=yes \\
  BUILDDIR=$BUILDDIR \\
  CMAKE_FLAGS=\"-DCMAKE_CXX_FLAGS=\\\"$FLAGS\\\" \\
   -DDUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS=TRUE \\
   -DCMAKE_DISABLE_FIND_PACKAGE_LATEX=TRUE \\
   $WITH_ZOLTAN\" " > config.opts
fi

DUNEBRANCH=
ALUGRIDBRANCH=
FEMBRANCH=
FEMDGBRANCH=
if [ "$DUNEVERSION" != "" ] ; then
  DUNEBRANCH="-b releases/$DUNEVERSION"
  ALUGRIDBRANCH="-b releases/2.4"
  FEMBRANCH="-b releases/2.4-dune-fem-dg"
  FEMDGBRANCH="-b releases/2.4"
fi

# get all dune modules necessary
for MOD in $DUNEMODULES ; do
  if [ "$MOD" == "dune-alugrid" ] ; then
    # use the special branch for dune-alugrid
    git clone $DUNEALUGRIDBRANCH https://gitlab.dune-project.org/extensions/dune-alugrid.git
  elif [ "$MOD" == "dune-fem" ] ; then
    # use the special branch for dune-fem
    git clone $FEMBRANCH https://gitlab.dune-project.org/dune-fem/dune-fem.git
  elif [ "$MOD" == "dune-fem-dg" ] ; then
    # use the special branch for dune-fem-dg
    git clone $FEMDGBRANCH https://gitlab.dune-project.org/dune-fem/dune-fem-dg.git
  else
    git clone $DUNEBRANCH https://gitlab.dune-project.org/core/$MOD.git
  fi
done

# build all DUNE modules in the correct order
./dune-common/bin/dunecontrol --opts=config.opts all

if [ "$BUILDDIR" == "" ]; then
  cd dune-fem-dg/build-cmake
else
  cd dune-fem-dg/$BUILDDIR/
fi
TARGET=test

# compile and run tests
make $MAKE_FLAGS $TARGET
