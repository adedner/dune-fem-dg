# searches for bgq_l1prefetch-headers and libs

AC_DEFUN([DUNE_PATH_BGQ_L1PREFETCH],[
  AC_REQUIRE([AC_PROG_CC])

  AC_ARG_WITH(bgq_l1prefetch,
    AC_HELP_STRING([--with-bgq_l1prefetch=PATH],[directory with BGQ_L1PREFETCH inside]))

# store old values
ac_save_LDFLAGS="$LDFLAGS"
ac_save_CPPFLAGS="$CPPFLAGS"
ac_save_LIBS="$LIBS"

BGQ_L1PREFETCHYES=0
## do nothing if no --with-bgq_l1prefetch was supplied
if test x$with_bgq_l1prefetch != x && test x$with_bgq_l1prefetch != xno ; then
  BGQ_L1PREFETCHYES=1
fi

if test x$BGQ_L1PREFETCHYES = x1 ; then

  # is --with-bgq_l1prefetch=bla used?
  if test "x$with_bgq_l1prefetch" != x ; then
    BGQ_L1PREFETCHROOT=`cd $with_bgq_l1prefetch && pwd`
	  if ! test -d $BGQ_L1PREFETCHROOT;  then
      AC_MSG_WARN([BGQ_L1PREFETCH directory $with_bgq_l1prefetch does not exist])
  	fi

    if test "x$BGQ_L1PREFETCHROOT" = x; then
      # use some default value...
      BGQ_L1PREFETCHROOT="/usr/local/bgq_l1prefetch"
    fi

    BGQ_L1PREFETCH_LIB_PATH="$BGQ_L1PREFETCHROOT/lib"
    BGQ_L1PREFETCH_INCLUDE_PATH="$BGQ_L1PREFETCHROOT/include"
  fi

  # set variables so that tests can use them
  REM_CPPFLAGS=$CPPFLAGS

  LDFLAGS="$LDFLAGS -L$BGQ_L1PREFETCH_LIB_PATH"
  BGQ_L1PREFETCH_INC_FLAG="-I$BGQ_L1PREFETCH_INCLUDE_PATH -DENABLE_BGQ_L1PREFETCH=1"
  CPPFLAGS="$CPPFLAGS $BGQ_L1PREFETCH_INC_FLAG"

  # check for header
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([l1p/pprefetch.h], 
     [BGQ_L1PREFETCH_CPPFLAGS="$BGQ_L1PREFETCH_INC_FLAG"
    HAVE_BGQ_L1PREFETCH="1"],
    AC_MSG_WARN([pprefetch.h not found in $BGQ_L1PREFETCH_INCLUDE_PATH]))
   
  CPPFLAGS="$REM_CPPFLAGS"
  REM_CPPFLAGS=

  REM_LDFLAGS=$LDFLAGS

  # if header is found...
  if test x$HAVE_BGQ_L1PREFETCH = x1 ; then
    LIBNAME="SPI_l1p"
    AC_CHECK_LIB($LIBNAME,[BGQ_L1PREFETCH_flops],
    [BGQ_L1PREFETCH_LIBS="-l$LIBNAME"
           BGQ_L1PREFETCH_LDFLAGS="-L$BGQ_L1PREFETCH_LIB_PATH"],
	  [HAVE_BGQ_L1PREFETCH="0"
	  AC_MSG_WARN(lib$LIBNAME not found!)])
  fi

  LDFLAGS=$REM_LDFLAGS
  AC_LANG_POP
fi

# survived all tests?
if test x$HAVE_BGQ_L1PREFETCH = x1 ; then
  AC_SUBST(BGQ_L1PREFETCH_LIBS, $BGQ_L1PREFETCH_LIBS)
  AC_SUBST(BGQ_L1PREFETCH_LDFLAGS, $BGQ_L1PREFETCH_LDFLAGS)
  AC_SUBST(BGQ_L1PREFETCH_CPPFLAGS, $BGQ_L1PREFETCH_CPPFLAGS)
  AC_DEFINE(HAVE_BGQ_L1PREFETCH, ENABLE_BGQ_L1PREFETCH,
    [This is only true if bgq_l1prefetch-library was found by configure 
     _and_ if the application uses the BGQ_L1PREFETCH_CPPFLAGS])

  # add to global list
  DUNE_ADD_ALL_PKG([BGQ_L1PREFETCH], [\${BGQ_L1PREFETCH_CPPFLAGS}],
                   [\${BGQ_L1PREFETCH_LDFLAGS}], [\${BGQ_L1PREFETCH_LIBS}])

  # set variable for summary
  with_bgq_l1prefetch="yes"

else
  AC_SUBST(BGQ_L1PREFETCH_LIBS, "")
  AC_SUBST(BGQ_L1PREFETCH_LDFLAGS, "")
  AC_SUBST(BGQ_L1PREFETCH_CPPFLAGS, "")

  # set variable for summary
  with_bgq_l1prefetch="no"
fi
  
# also tell automake
AM_CONDITIONAL(BGQ_L1PREFETCH, test x$HAVE_BGQ_L1PREFETCH = x1)

# reset old values
LIBS="$ac_save_LIBS"
CPPFLAGS="$ac_save_CPPFLAGS"
LDFLAGS="$ac_save_LDFLAGS"

DUNE_ADD_SUMMARY_ENTRY([BGQ_L1PREFETCH],[$with_bgq_l1prefetch])
])
