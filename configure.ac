dnl
dnl Process this file with autoconf to produce a configure script.
dnl 

dnl
dnl Setup autoconf
dnl
AC_PREREQ([2.69])
AC_INIT([omf5],[1.0.0-beta.1],[https://github.com/assaferan/omf5],[omf5],[https://github.com/assaferan/omf5])
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
if test "${FLINT_PATH+set}" != set; then
     export FLINT_PATH="/usr/local/"
     echo FLINT_PATH=${FLINT_PATH}
fi
if test "${ANTIC_PATH+set}" != set; then
     export ANTIC_PATH="/usr/local/"
     echo ANTIC_PATH=${ANTIC_PATH}
fi
if test "${GMP_PATH+set}" != set; then
     export GMP_PATH="/usr/local/"
     echo GMP_PATH=${GMP_PATH}
fi

AC_ARG_VAR(FLINT_PATH, [location of flint installation])
AC_ARG_VAR(ANTIC_PATH, [location of antic installation])
AC_ARG_VAR(GMP_PATH, [location of gmp installation])

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
AC_OUTPUT()
