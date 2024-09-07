#!/bin/bash

#change appropriately, i.e. 2.8 or leave empty which refers to master
# use latest-stable to get a recent stable version
DUNEVERSION=latest-stable
UFLVERSION=2022.2.0

# use ON or OFF
USEVENV=ON

WORKDIR=${PWD}
echo "Placing DUNE modules in $WORKDIR"

if [ "$DUNE_CONTROL_PATH" != "" ]; then
  if [ "$WORKDIR" != "$DUNE_CONTROL_PATH" ]; then
    echo "DUNE_CONTROL_PATH is already set to $DUNE_CONTROL_PATH"
    exit 1
  fi
fi

# check for system installation
if test -f /usr/share/dune/cmake/modules/DuneMacros.cmake || test -f /usr/bin/dunecontrol ; then
  echo "DUNE system installation seems to exists, remove it first by running:"
  echo ""
  echo "sudo apt remove libdune*"
  echo "sudo apt autoremove"
  echo ""
  exit 1
fi

# create necessary python virtual environment
# this script assumes the name venv.
# Otherwise copy the instructions from the script
# to build you own

CMAKE_NOT_FOUND=`command -v cmake`

if [ "$CMAKE_NOT_FOUND" == "" ]; then
  CMAKEPIP=cmake
  echo "Installing cmake since no cmake was found!"
else
  CMAKE_VERSION=`cmake --version | head -1 | cut -d " " -f 3 | cut -d " " -f 1`
  REQUIRED_VERSION="3.13.3"
  # check if cmake version is ok
  if awk 'BEGIN {exit !('$CMAKE_VERSION' < '$REQUIRED_VERSION')}'; then
    CMAKEPIP=cmake
    echo "Installing cmake since current version is not new enough!"
  fi
fi

if [ "$USEVENV" == "ON" ]; then
  # create necessary python virtual environment
  VENVDIR=$WORKDIR/venv
  if ! test -d $VENVDIR ; then
    python3 -m venv --system-site-packages $VENVDIR
    source $VENVDIR/bin/activate
    pip install --upgrade pip
    pip install $CMAKEPIP fenics-ufl==$UFLVERSION numpy matplotlib mpi4py
  else
    source $VENVDIR/bin/activate
  fi
fi

FLAGS="-O3 -DNDEBUG -funroll-loops -finline-functions -Wall -ftree-vectorize -fno-stack-protector -mtune=native"

DUNECOREMODULES="dune-common dune-istl dune-geometry dune-grid dune-localfunctions"
DUNEEXTMODULES="dune-alugrid dune-spgrid dune-polygongrid"
DUNEFEMMODULES="dune-fem dune-fempy dune-fem-dg"

# build flags for all DUNE modules
# change according to your needs
if test -f $WORKDIR/config.opts ; then
  read -p "Found config.opts. Overwrite with default? (y,n) " YN
  if [ "$YN" == "y" ] ;then
    echo "Overwriting config.opts!"
    rm -f $WORKDIR/config.opts
  fi
fi

if ! test -f $WORKDIR/config.opts ; then
echo "\
DUNEPATH=`pwd`
BUILDDIR=build-cmake
USE_CMAKE=yes
MAKE_FLAGS=-j4
CMAKE_FLAGS=\"-DCMAKE_CXX_FLAGS=\\\"$FLAGS\\\"  \\
 -DDUNE_ENABLE_PYTHONBINDINGS=ON \\
 -DDUNE_PYTHON_USE_VENV=$USEVENV \\
 -DADDITIONAL_PIP_PARAMS="-upgrade" \\
 -DCMAKE_LD_FLAGS=\\\"$PY_LDFLAGS\\\" \\
 -DCMAKE_POSITION_INDEPENDENT_CODE=TRUE \\
 -DDISABLE_DOCUMENTATION=TRUE \\
 -DCMAKE_DISABLE_FIND_PACKAGE_Vc=TRUE \\
 -DCMAKE_DISABLE_FIND_PACKAGE_LATEX=TRUE\" " > $WORKDIR/config.opts
fi

ACTIVATE=$WORKDIR/activate.sh
if [ "$USEVENV" == "ON" ]; then
  ACTIVATE=$VENVDIR/bin/activate
else
  DEACTIVATEFUNCTION="deactivate() {
  echo \"\"
}
"
fi

FOUND_DUNE_ACTIVATE=`grep "DUNE_VENV_SPECIFIC_SETUP" $ACTIVATE`

if [ "$FOUND_DUNE_ACTIVATE" == "" ]; then
echo "
## DUNE_VENV_SPECIFIC_SETUP

# setVariable( varname newvalue )
setVariable() {
  if [ -n \"\${!1}\" ]; then
    export _OLD_VIRTUAL_\${1}=\"\${!1}\"
  fi
  export \${1}=\"\${2}\"
}

# set current main working directory
setVariable DUNE_CONTROL_PATH \"$WORKDIR\"
setVariable DUNE_LOG_LEVEL \"info\"

# defines CMAKE_FLAGS
source \${DUNE_CONTROL_PATH}/config.opts

setVariable DUNE_PY_DIR \${DUNE_CONTROL_PATH}/cache/
setVariable DUNE_CMAKE_FLAGS \"\$CMAKE_FLAGS\"

DUNEPYTHONPATH=
MODULES=\`\$DUNE_CONTROL_PATH/dune-common/bin/dunecontrol --print 2> /dev/null\`
for MOD in \$MODULES; do
  MODPATH=\"\$DUNE_CONTROL_PATH/\${MOD}/build-cmake/python\"
  MODFOUND=\`echo \$DUNEPYTHONPATH | grep \$MODPATH\`
  if [ \"\$MODFOUND\" == \"\" ]; then
    DUNEPYTHONPATH=\$DUNEPYTHONPATH:\$MODPATH
  fi
done

setVariable PYTHONPATH \$DUNEPYTHONPATH

echo \"DUNE_LOG_LEVEL=\$DUNE_LOG_LEVEL\"
echo \"Change with 'export DUNE_LOG_LEVEL={none,critical,info,debug}' if necessary!\"

save_function() {
    local ORIG_FUNC=\$(declare -f \$1)
    local NEWNAME_FUNC=\"\$2\${ORIG_FUNC#\$1}\"
    eval \"\$NEWNAME_FUNC\"
}

restoreVariable() {
  VARNAME=_OLD_VIRTUAL_\$1
  if [ -n \"\${!VARNAME}\" ]; then
    export \${1}=\${!VARNAME}
    unset \${VARNAME}
  else
    unset \${1}
  fi
}

$DEACTIVATEFUNCTION

save_function deactivate venv_deactivate
deactivate() {
  restoreVariable DUNE_CONTROL_PATH
  restoreVariable DUNE_PY_DIR
  restoreVariable DUNE_CMAKE_FLAGS
  restoreVariable DUNE_LOG_LEVEL
  restoreVariable PYTHONPATH

  # call original deactivate
  venv_deactivate
  # self destroy
  unset venv_deactivate
  unset deactivate
}
" >> $ACTIVATE
fi


# load environment variables
source $ACTIVATE

#################################################################
##
## Obtain DUNE modules from servers
##
#################################################################

DUNEBRANCH=
URL=https://gitlab.dune-project.org
EXT=core

if [ "$DUNEVERSION" == "latest-stable" ] ; then
  # remove extension, this is not used on Lund server
  EXT=
  URL=https://gitlab.maths.lu.se/dune
elif [ "$DUNEVERSION" != "" ] ; then
  # assume that some version was selected (make sure that version is correct)
  DUNEBRANCH="-b releases/$DUNEVERSION"
fi

if [ "$DUNEVERSION" == "latest-stable" ]; then
  REPO=$URL/dune-common.git
  # get latest stable tag starting with s20, e.g. s2024.08
  TAG=`git ls-remote --refs --tags $REPO | grep "s20" | cut --delimiter='/' --fields=3 | tr '-' '~'  | sort --version-sort  | tail --lines=1`
  DUNEBRANCH="-b $TAG"
  echo "Selecting stable tag $TAG"
fi

# get all dune modules necessary
for MOD in $DUNECOREMODULES ; do
  git clone --depth 1 $DUNEBRANCH $URL/$EXT/$MOD.git
done

if [ "$EXT" != "" ]; then
  EXT=extensions
fi

# get all dune extension modules necessary
for MOD in $DUNEEXTMODULES ; do
  git clone --depth 1 $DUNEBRANCH $URL/$EXT/$MOD.git
done

if [ "$EXT" != "" ]; then
  EXT=dune-fem
fi
# get all dune extension modules necessary
for MOD in $DUNEFEMMODULES ; do
  git clone --depth 1 $DUNEBRANCH $URL/$EXT/$MOD.git
done

# build all DUNE modules using dune-control
./dune-common/bin/dunecontrol --opts=config.opts all

echo "####################################################

Build finished (hopefully successful). Use

source $ACTIVATE

to activate the virtual environment!"
