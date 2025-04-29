#!/bin/bash

WDIR=`pwd`
DATE=`date`
echo "DUNE Modules used at $DATE"
if [ "$DUNE_CONTROL_PATH" != "" ]; then
  MODS=`$DUNE_CONTROL_PATH/dune-common/bin/dunecontrol print`

  for MOD in $MODS; do
    COMMIT=`cd $DUNE_CONTROL_PATH/$MOD && git rev-parse HEAD`
    printf "%-20s%-4s"  $MOD
    printf "%-30s%-4s\n"  $COMMIT
  done
fi

echo "*****************************************************************"
echo ""

#DUNECOREMODULES="dune-common dune-istl dune-geometry dune-grid dune-localfunctions"
#DUNEEXTMODULES="dune-alugrid dune-spgrid dune-polygongrid"
#DUNEFEMMODULES="dune-fem dune-fempy dune-fem-dg"
