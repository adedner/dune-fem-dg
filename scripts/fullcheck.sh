#!/bin/bash

# check for parameter pointing to DUNE base directory
# ---------------------------------------------------

DUNECONTROL="dune-common/bin/dunecontrol"

if test \( $# -lt 1 \) -o ! -e $1/$DUNECONTROL ; then
  echo "Usage: $0 <dune-base-dir>"
  exit 1
fi

MODNAME="dune-fem-dg"

echo "Full Check of $MODNAME"
echo "-------------------------"

echo
echo "Host Name: $HOSTNAME"
echo "Host Type: $HOSTTYPE"

# set up some variables
# ---------------------

WORKINGDIR=`pwd`
cd $1
DUNEDIR=`pwd`
FEMDIR="$DUNEDIR/$MODNAME"
BUILDDIR="$DUNEDIR"
SCRIPTSDIR="$FEMDIR/scripts"
OPTSDIR="$SCRIPTSDIR/opts"

# search for all existing dune modules
# ------------------------------------
MODULES=""
for modctrl in $(find -name dune.module -print); do
  MODULES+="$(dirname $modctrl | sed -e 's@^[.]/@@g') " 
done

errors=0

# configure with minimal options
# ------------------------------

MINIMALOPTS="$OPTSDIR/minimal.opts"
BUILDDIR="$BUILDDIR/$(source $MINIMALOPTS; echo $BUILD_DIR)"

MAKE_FLAGS=""
MAKE_FLAGS="$(source $MINIMALOPTS; echo $MAKE_FLAGS)"

if test ! -e $MINIMALOPTS ; then
  echo "Error: $MINIMALOPTS not found."
  exit 1
fi


minimal_configure()
{
  local check=`mktemp -p $WORKINGDIR check.XXXXXX`
  {
    $DUNECONTROL --opts=$MINIMALOPTS --builddir=$BUILDDIR all
    echo $? > $check
  } 2>&1 | dd conv=notrunc > $WORKINGDIR/minimal-conf.out 2>/dev/null
  local return_value=`cat $check`
  rm $check
  return $return_value
}


echo
echo "Configuring with minimal options..."
cd $DUNEDIR
if ! minimal_configure ; then
  echo "Fatal: Cannot configure with minimal options (see $WORKINGDIR/minimal-conf.out)."
  exit 1
fi


# check headers
# -------------

for module in $MODULES;
do
  CHECKLOG="$WORKINGDIR/minimal-hc-$module.out"
  echo
  echo "Checking headers in $module ..."
  cd $BUILDDIR/$module
  make headercheck $MAKE_FLAGS -i &> $CHECKLOG
  hc_errors=$(grep error: $CHECKLOG)
  if test -z "$hc_errors" ; then
    rm $CHECKLOG
  else
    if test "x$module" == "x$MODNAME"; then
      echo "Error: headercheck for module $module failed (see $CHECKLOG )"
      errors=$((errors+1))
    else
      echo "Warning: headercheck for module $module failed (see $CHECKLOG )"
    fi
  fi
done


# perform make test
# ------------------

cd $WORKINGDIR 
echo
echo "Checking for minimal options ..."

CHECKLOG="$WORKINGDIR/minimal-check.out"

if ! $SCRIPTSDIR/check-tests.sh $BUILDDIR/$MODNAME "$MAKE_FLAGS"; then
  echo "Error: make test failed with minimal options (see $CHECKLOG)"
  errors=$((errors+1))
fi
mv $WORKINGDIR/check-tests.out $CHECKLOG

if test $errors -gt 0 ; then
  exit 1
else
  exit 0
fi
