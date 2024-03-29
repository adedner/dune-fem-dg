#!/bin/bash

# check for parameter pointing to DUNE base directory
# ---------------------------------------------------

DUNECONTROL="dune-common/bin/dunecontrol"

if test \( $# -lt 1 \) -o ! -e $1/$DUNECONTROL ; then
  echo "Usage: $0 <dune-base-dir>"
  exit 1
fi

# set up some variables
# ---------------------

WORKINGDIR=`pwd`
cd $1
DUNEDIR=`pwd`
FEMDIR="$DUNEDIR/dune-fem"
SCRIPTSDIR="$FEMDIR/scripts"
OPTSDIR="$SCRIPTSDIR/opts"

MODULES="dune-common dune-geometry dune-grid dune-istl dune-localfunctions dune-spgrid dune-fem dune-alugrid dune-fem-dg"
EXISTINGMODULES=""

# fetch missing tarballs from website
# -----------------------------------

for MODULE in $MODULES ; do
  # ignore missing modules since istl may be missing
  if ! test -d $DUNEDIR/$MODULE ; then
    continue;
  fi

  cd $DUNEDIR/$MODULE
  if test x`find -maxdepth 1 -name "*.tar.gz"` != x ; then
    continue
  fi

  if test "$MODULE" != "dune-fem" ; then
    echo
    cached="/usr/people/dune/src/tarballs/$MODULE-1.2svn.tar.gz"
    if test -e $cached ; then
      echo "Copying tarball from backup ($cached)..."
      cp $cached .
    fi
  fi

  if test x`find -maxdepth 1 -name "*.tar.gz"` == x ; then
    echo "Fatal: No tarball available for $MODULE"
    exit 1
  fi

  EXISTINGMODULES+=" $MODULE"
done

# perform tests
# -------------

TESTDIR=`mktemp -d -p $WORKINGDIR dune-tmp-XXXXXX`

do_configure()
{
  local check=`mktemp -p $WORKINGDIR check.XXXXXX`
  {
    $DUNECONTROL --opts=$OPTSDIR/$OPTS all
    echo $? > $check
  } 2>&1 | dd conv=notrunc > $CONFIGLOG 2>/dev/null
  local return_value=`cat $check`
  rm $check
  return $return_value
}

errors=0
for OPTS in `cd $OPTSDIR ; ls *.opts` ; do
  echo
  echo "Checking $OPTS..."

  cd $TESTDIR
  rm -rf *
 
  find $DUNEDIR -maxdepth 2 -name "*.tar.gz" -exec tar -xzf \{\} \;

  for MODULE in $MODULES ; do
    if ls -d ${MODULE}* &> /dev/null ; then
      mv `ls -d ${MODULE}*` $MODULE
    fi
  done

  CONFIGLOG="$WORKINGDIR/${OPTS%.opts}-conf.out"
  if ! do_configure ; then
    echo "Error: Cannot configure with $OPTS (see $CONFIGLOG)"
    errors=$((errors+1))
    continue
  fi

  cd $WORKINGDIR
  CHECKLOG="$WORKINGDIR/${OPTS%.opts}-check.out"
  MAKE_CHECK_FLAGS=""
  MAKE_CHECK_FLAGS="$(source $OPTSDIR/$OPTS; echo $MAKE_CHECK_FLAGS)"

  # check for dependencies
  MAKE_CHECK_DEPS="$(source $OPTSDIR/$OPTS; echo $MAKE_CHECK_DEPS)"

  MISSINGDEPS=""
  for dep in $MAKE_CHECK_DEPS ; do
    if ! echo $EXISTINGMODULES | grep $dep; then
      MISSINGDEPS+=" $dep"
    fi
  done
 
  if test -z "$MISSINGDEPS" ; then 
    if ! $SCRIPTSDIR/check-tests.sh $TESTDIR/dune-fem "$MAKE_CHECK_FLAGS"; then
      echo "Error: Check failed with $OPTS (see $CHECKLOG)"
      errors=$((errors+1))
    fi
    mv $WORKINGDIR/check-tests.out $CHECKLOG
  else
    echo "Skipping $OPTS due to missing dependencies:$MISSINGDEPS"
  fi
done

# clean up
# --------

rm -rf $TESTDIR

if test $errors -gt 0 ; then
  exit 1
else
  exit 0
fi
