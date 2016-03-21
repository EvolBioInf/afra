AC_INIT([afra], [m4_esyscmd_s([git describe --always])])
AM_INIT_AUTOMAKE([-Wall foreign subdir-objects])

AC_PROG_CC
AC_PROG_CPP

AC_LANG(C)

AC_CHECK_HEADERS([stdlib.h string.h])
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
