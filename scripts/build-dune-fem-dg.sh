#!/bin/bash

WORKDIR=${PWD}

if [ "$DUNE_CONTROL_PATH" != "" ]; then
  if [ "$WORKDIR" != "$DUNE_CONTROL_PATH" ]; then
    echo "DUNE_CONTROL_PATH is already set to $DUNE_CONTROL_PATH"
    exit 1
  fi
fi

CMAKE_VERSION=`cmake --version | head -1 | cut -d " " -f 3 | cut -d " " -f 1`
REQUIRED_VERSION="3.13.3"
# check if cmake version is ok
if awk 'BEGIN {exit !('$CMAKE_VERSION' < '$REQUIRED_VERSION')}'; then
  CMAKEPIP=cmake
fi
echo $CMAKEPIP

# create necessary python virtual environment
if ! test -d $WORKDIR/dune-venv ; then
  python3 -m venv $WORKDIR/dune-venv
  source $WORKDIR/dune-venv/bin/activate
  pip install --upgrade pip
  pip install $CMAKEPIP ufl numpy matplotlib mpi4py
  #pip install numpy matplotlib mpi4py
else
  source $WORKDIR/dune-venv/bin/activate
fi

#change appropriately, i.e. 2.6 or empty which refers to master
DUNEVERSION=

FLAGS="-O3 -DNDEBUG -funroll-loops -finline-functions -Wall -ftree-vectorize -fno-stack-protector -mtune=native"

DUNECOREMODULES="dune-common dune-istl dune-geometry dune-grid dune-localfunctions"
DUNEEXTMODULES="dune-alugrid dune-spgrid"
DUNEFEMMODULES="dune-fem dune-fem-dg"
EXTRAFEMMODULES="dune-fempy"

# build flags for all DUNE modules
# change according to your needs
if ! test -f config.opts ; then
echo "\
DUNEPATH=`pwd`
BUILDDIR=build-cmake
USE_CMAKE=yes
MAKE_FLAGS=-j4
CMAKE_FLAGS=\"-DCMAKE_CXX_FLAGS=\\\"$FLAGS\\\"  \\
 -DDUNE_ENABLE_PYTHONBINDINGS=ON \\
 -DADDITIONAL_PIP_PARAMS="-upgrade" \\
 -DCMAKE_LD_FLAGS=\\\"$PY_LDFLAGS\\\" \\
 -DCMAKE_POSITION_INDEPENDENT_CODE=TRUE \\
 -DDISABLE_DOCUMENTATION=TRUE \\
 -DCMAKE_DISABLE_FIND_PACKAGE_LATEX=TRUE\" " > config.opts
fi

if ! test -f activate.sh ; then
echo "#!/bin/bash

# set current main working directory
export DUNE_CONTROL_PATH=\${PWD}

source \$DUNE_CONTROL_PATH/dune-venv/bin/activate

# defines CMAKE_FLAGS
source \${DUNE_CONTROL_PATH}/config.opts

# CRITICAL, DEBUG or NOTSET
export DUNE_LOG_LEVEL=CRITICAL
export DUNE_LOG_FORMAT='%(asctime)s - %(name)s - %(levelname)s - %(message)s'

export DUNE_PY_DIR=\${DUNE_CONTROL_PATH}/cache/

export DUNE_CMAKE_FLAGS="\${CMAKE_FLAGS}"

MODULES=\"$DUNECOREMODULES $DUNEEXTMODULES $DUNEFEMMODULES\"
for MOD in \$MODULES; do
  MODPATH=\"\${PWD}/\${MOD}/build-cmake/python\"
  MODFOUND=\`echo \$PYTHONPATH | grep \$MODPATH\`
  if [ \"\$MODFOUND\" == \"\" ]; then
    export PYTHONPATH=\$PYTHONPATH:\$MODPATH
  fi
done
" > activate.sh
fi

DUNEBRANCH=
if [ "$DUNEVERSION" != "" ] ; then
  DUNEBRANCH="-b releases/$DUNEVERSION"
fi

# get all dune modules necessary
for MOD in $DUNECOREMODULES ; do
  git clone $DUNEBRANCH https://gitlab.dune-project.org/core/$MOD.git
done

# get all dune extension modules necessary
for MOD in $DUNEEXTMODULES ; do
  if [ "$MOD" == "dune-alugrid" ]; then
    git clone $DUNEBRANCH https://gitlab.dune-project.org/extensions/$MOD.git
  elif [ "$MOD" == "dune-spgrid" ]; then
    git clone $DUNEBRANCH https://gitlab.dune-project.org/extensions/$MOD.git
  else
    git clone $DUNEBRANCH https://gitlab.dune-project.org/staging/$MOD.git
  fi
done

# get all dune extension modules necessary
for MOD in $DUNEFEMMODULES ; do
  git clone $DUNEBRANCH https://gitlab.dune-project.org/dune-fem/$MOD.git
done

# get all dune extension modules necessary
for MOD in $EXTRAFEMMODULES ; do
  git clone $DUNEBRANCH https://gitlab.dune-project.org/dune-fem/$MOD.git
done

# load environment variables
source activate.sh

# build all DUNE modules using dune-control
./dune-common/bin/dunecontrol --opts=config.opts all

#cd dune-fem-dg/build-cmake
#make build_tests
#make test
