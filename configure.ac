dnl
dnl Process this file with autoconf to produce a configure script.
dnl 

dnl
dnl Setup autoconf
dnl
AC_PREREQ([2.69])
AC_INIT([omf5],[1.7.3],[https://github.com/assaferan/omf5],[omf5],[https://github.com/assaferan/omf5])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([cnf])

# AC_CANONICAL_HOST is needed to access the 'host_os' variable    
AC_CANONICAL_HOST

build_linux=no
build_windows=no
build_mac=no

# Detect the target system
case "${host_os}" in
    linux*)
        build_linux=yes
        ;;
    cygwin*|mingw*)
        build_windows=yes
        ;;
    darwin*)
        build_mac=yes
        ;;
    *)
        AC_MSG_ERROR(["OS $host_os is not supported"])
        ;;
esac

AM_CONDITIONAL([LINUX], [test "$build_linux" = "yes"])
AM_CONDITIONAL([WINDOWS], [test "$build_windows" = "yes"])
AM_CONDITIONAL([OSX], [test "$build_mac" = "yes"])

# x86-specific optimization flags (-msse3 -mavx) are not valid elsewhere (e.g. arm64)
case "${host_cpu}" in
    x86_64|i?86)
        build_x86=yes
        ;;
    *)
        build_x86=no
        ;;
esac

AM_CONDITIONAL([X86], [test "$build_x86" = "yes"])

dnl
dnl C is the language
dnl
: ${CFLAGS=""}
AC_PROG_CC
AC_LANG([C])

dnl
dnl Setup automake
dnl
AM_INIT_AUTOMAKE([1.11 -Wall -Werror foreign subdir-objects])
AM_SILENT_RULES([yes])
AM_PROG_AR
LT_INIT
LT_INIT

dnl
dnl Setting the locations of flint, antic and gmp
dnl
omf5_flint_path_set=yes
omf5_antic_path_set=yes
omf5_gmp_path_set=yes
omf5_carat_path_set=yes
if test "${FLINT_PATH+set}" != set; then
     export FLINT_PATH="/usr/local/"
     omf5_flint_path_set=no
fi
if test "${ANTIC_PATH+set}" != set; then
     export ANTIC_PATH="/usr/local/"
     omf5_antic_path_set=no
fi
if test "${GMP_PATH+set}" != set; then
     export GMP_PATH="/usr/local/"
     omf5_gmp_path_set=no
fi
if test "${CARAT_PATH+set}" != set; then
     export CARAT_PATH="/usr/local/"
     omf5_carat_path_set=no
fi

echo FLINT_PATH=${FLINT_PATH}
echo ANTIC_PATH=${ANTIC_PATH}
echo GMP_PATH=${GMP_PATH}
echo CARAT_PATH=${CARAT_PATH}

AC_ARG_VAR(FLINT_PATH, [location of flint installation])
AC_ARG_VAR(ANTIC_PATH, [location of antic installation])
AC_ARG_VAR(GMP_PATH, [location of gmp installation])
AC_ARG_VAR(CARAT_PATH, [location of carat installation])

dnl
dnl Finally, generate the Makefiles and output everything
dnl
AC_CONFIG_FILES([Makefile src/Makefile])
AM_PROG_CC_C_O

dnl this doesn't seem to work
dnl AC_SEARCH_LIBS([], [flint], [
dnl AC_MSG_ERROR([unable to find the flint library])
dnl ])

dnl This doesn't seem to work
dnl AC_CHECK_HEADERS([flint.h],
dnl        [found_flint_headers=yes; break;])

AC_SEARCH_LIBS([gmp])

dnl
dnl Include and library paths of the dependencies, flint first, so that
dnl the build compiles against the flint headers checked here (see
dnl src/Makefile.am). The flint library is checked below.
dnl
DEPS_CPPFLAGS="-I${FLINT_PATH}/include -I${GMP_PATH}/include"
DEPS_LDFLAGS="-L${FLINT_PATH}/lib -L${GMP_PATH}/lib"

save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $DEPS_CPPFLAGS"

AC_CHECK_HEADER([flint/flint.h], [],
  [AC_MSG_ERROR([flint/flint.h not found under ${FLINT_PATH}/include, set FLINT_PATH to the flint installation])])

dnl omf5 uses functions added in flint 2.9, such as fmpq_set_fmpz
AC_MSG_CHECKING([whether flint is version 2.9 or later])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <flint/flint.h>]],
  [[#if __FLINT_RELEASE < 20900
    #error flint older than 2.9
    #endif]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([omf5 requires flint 2.9 or later, but the flint under ${FLINT_PATH}/include is older])])

dnl
dnl flint >= 3.0 includes the number field modules (nf, nf_elem)
dnl that used to live in antic, so antic is only needed for flint 2.x
dnl
AC_MSG_CHECKING([whether flint is version 3.0 or later])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <flint/flint.h>]],
  [[#if __FLINT_RELEASE < 30000
    #error flint 2.x
    #endif]])],
  [AC_MSG_RESULT([yes]); need_antic=no],
  [AC_MSG_RESULT([no, linking with antic]); need_antic=yes])

ANTIC_LIBS=""
if test "$need_antic" = "yes"; then
   DEPS_CPPFLAGS="$DEPS_CPPFLAGS -I${ANTIC_PATH}/include"
   DEPS_LDFLAGS="$DEPS_LDFLAGS -L${ANTIC_PATH}/lib"
   ANTIC_LIBS="-lantic"
   CPPFLAGS="$save_CPPFLAGS $DEPS_CPPFLAGS"
   AC_CHECK_HEADER([antic/nf.h], [],
     [AC_MSG_ERROR([flint 2.x requires antic, but antic/nf.h was not found under ${ANTIC_PATH}/include, set ANTIC_PATH to the antic installation])])
fi

dnl
dnl The flint library found by the linker has to be the one of these
dnl headers: a -L of the paths above, which the build puts before the
dnl user's LDFLAGS (src/Makefile.am), can put another libflint first,
dnl and omf5 then crashes. The test links as the build does, with the
dnl -L of CARAT too. If the versions differ, the -L of the paths that
dnl were not set are left out, so that the linker searches its own
dnl directories; if they still differ, configure stops.
dnl
CARAT_LDFLAGS="-L${CARAT_PATH}/lib -L${CARAT_PATH}"
m4_define([OMF5_FLINT_VERSION_PROGRAM],
  [AC_LANG_PROGRAM([[#include <string.h>
#include <flint/flint.h>]],
    [[return strcmp(flint_version, FLINT_VERSION) != 0;]])])
save_LDFLAGS="$LDFLAGS"
save_LIBS="$LIBS"
LIBS="-lflint -lgmp -lm $save_LIBS"
LDFLAGS="$DEPS_LDFLAGS $CARAT_LDFLAGS $save_LDFLAGS"
AC_MSG_CHECKING([whether a program links with -lflint -lgmp -lm])
AC_LINK_IFELSE([OMF5_FLINT_VERSION_PROGRAM],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([cannot link a program with -lflint -lgmp -lm and $LDFLAGS; set FLINT_PATH and GMP_PATH to the installations of flint and gmp])])
AC_MSG_CHECKING([whether the flint library linked matches the flint headers])
AC_RUN_IFELSE([OMF5_FLINT_VERSION_PROGRAM],
  [AC_MSG_RESULT([yes]); flint_link=yes],
  [AC_MSG_RESULT([no]); flint_link=no],
  [AC_MSG_RESULT([not checked, cross compiling]); flint_link=yes])
if test "$flint_link" = "no"; then
   DEPS_LDFLAGS=""
   test "$omf5_flint_path_set" = "yes" && DEPS_LDFLAGS="$DEPS_LDFLAGS -L${FLINT_PATH}/lib"
   test "$omf5_gmp_path_set" = "yes" && DEPS_LDFLAGS="$DEPS_LDFLAGS -L${GMP_PATH}/lib"
   test "$need_antic" = "yes" && test "$omf5_antic_path_set" = "yes" && DEPS_LDFLAGS="$DEPS_LDFLAGS -L${ANTIC_PATH}/lib"
   DEPS_LDFLAGS=`echo $DEPS_LDFLAGS`
   test "$omf5_carat_path_set" = "yes" || CARAT_LDFLAGS=""
   omf5_retry_ldflags=`echo $DEPS_LDFLAGS $CARAT_LDFLAGS`
   AC_MSG_CHECKING([whether it does without the -L of the default paths (${omf5_retry_ldflags:-none})])
   LDFLAGS="$DEPS_LDFLAGS $CARAT_LDFLAGS $save_LDFLAGS"
   AC_RUN_IFELSE([OMF5_FLINT_VERSION_PROGRAM],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_ERROR([the flint library found by the linker is not the flint of ${FLINT_PATH}/include, or it cannot be loaded: another libflint may come first on the library path; set FLINT_PATH, GMP_PATH and CARAT_PATH to installations without another libflint])],
     [AC_MSG_RESULT([not checked, cross compiling])])
fi
LDFLAGS="$save_LDFLAGS"
LIBS="$save_LIBS"

CPPFLAGS="$save_CPPFLAGS"
AC_SUBST([DEPS_CPPFLAGS])
AC_SUBST([DEPS_LDFLAGS])
AC_SUBST([CARAT_LDFLAGS])
AC_SUBST([ANTIC_LIBS])
AC_OUTPUT()
