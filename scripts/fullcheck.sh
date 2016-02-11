#!/bin/bash

# check for parameter pointing to DUNE base directory
# ---------------------------------------------------

DUNECONTROL="dune-common/bin/dunecontrol"

if test \( $# -lt 1 \) -o ! -e $1/$DUNECONTROL ; then
  echo "Usage: $0 <dune-base-dir>"
  exit 1
fi

echo "Full Check of dune-fem-dg"
echo "-------------------------"

echo
echo "Host Name: $HOSTNAME"
echo "Host Type: $HOSTTYPE"

# set up some variables
# ---------------------

WORKINGDIR=`pwd`
cd $1
DUNEDIR=`pwd`
FEMDIR="$DUNEDIR/dune-fem-dg"
SCRIPTSDIR="$FEMDIR/scripts"
OPTSDIR="$SCRIPTSDIR/opts"

errors=0

NONOPTSFILE=""
if test -e $OPTSDIR/config_none.opts; then
  NONOPTSFILE="config_none.opts"
else
  echo "config_none.opts not found."
  exit 1
fi

echo ""
echo "Performing checks for config_none.opts..."

cd $WORKINGDIR
# perform headercheck only ones
if ! $SCRIPTSDIR/check-opts.sh $DUNEDIR $NONOPTSFILE true; then
  errors=$((errors+1))
fi

if test $errors -gt 0 ; then
  exit 1
else
  exit 0
fi
