#!/bin/bash

WORKDIR=${PWD}

if [ "$DUNE_CONTROL_PATH" != "" ]; then
  if [ "$WORKDIR" != "$DUNE_CONTROL_PATH" ]; then
    echo "DUNE_CONTROL_PATH is already set to $DUNE_CONTROL_PATH"
    exit 1
  fi
fi

# update all DUNE modules using dune-control
./dune-common/bin/dunecontrol git pull

# NOTE: if conflicts arise due to diverging repositories then fix these
# conflicts and re-run the script
# Optional: clean the build directories of all modules to invoke a build from
# scratch, i.e.
# MODULES=`./dune-common/bin/dunecontrol print 2>/dev/null`
MODULES=
for MOD in $MODULES; do
  rm -rf $MOD/build-cmake
done

# build all DUNE modules using dune-control
./dune-common/bin/dunecontrol --opts=config.opts all
