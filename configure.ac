dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(create_pox, 1.0)


AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE()

AC_PROG_CXX

AC_CHECK_HEADER_STDBOOL
AC_CHECK_HEADERS([stdlib.h unistd.h sys/param.h])
AC_CHECK_FUNCS([getcwd strstr memset pow sqrt])

AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT

