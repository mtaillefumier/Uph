dnl #                                               -*- Autoconf -*-
dnl # Process this file with autoconf to produce a configure script.



AC_INIT([SPH], [0.1.0], [kalale.chola@oist.jp])
AC_PREREQ([2.52])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADERS([config.h])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_ISC_POSIX
AM_INIT_AUTOMAKE([1.6 dist-bzip2])
LT_INIT([disable_shared])
AC_DISABLE_SHARED

dnl # Checks for programs.
AC_PROG_LIBTOOL
AC_PROG_CXX(g++ icc)
AC_PROG_CC(gcc icc)
AC_OPENMP

dnl #AX_MPI([],[AC_MSG_ERROR([could not find mpi library])])
dnl #AC_CHECK_PROG(MPIRUN, mpirun, mpirun)
dnl #AC_SUBST(MPIRUN)
PKG_PROG_PKG_CONFIG
VMAJ=`echo $PACKAGE_VERSION | awk -F. '{printf("%s", $1);}'`
VMIN=`echo $PACKAGE_VERSION | awk -F. '{printf("%s", $2);}'`
VMIC=`echo $PACKAGE_VERSION | awk -F. '{printf("%s", $3);}'`
version_info=`expr $VMAJ + $VMIN`":$VMIC"
AC_SUBST(version_info)
release_info="-release ver-pre-svn"
AC_SUBST(release_info)
dnl # Checks for libraries.
if $PKG_CONFIG --atleast-pkgconfig-version 0.22; then
   pkgconfig_requires_private="Requires.private"
else
   pkgconfig_requires_private="Requires"
fi
AC_SUBST(pkgconfig_requires_private)
dnl # Checks for header files.
AC_CHECK_HEADERS([malloc.h stdlib.h string.h math.h])

dnl # Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_PROG_CC_C99

dnl # Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset pow sqrt strcasecmp strdup])

have_gsl="no"
PKG_CHECK_MODULES(GSL, [gsl >= 1.12])

requirement="gsl"

dnl # mpi support

dnl # cuda support

have_cuda="no"

AC_ARG_WITH([cuda],
[AS_HELP_STRING([--with-cuda=DIR],[Use GPU to do some of the calculations. You can specify the directory where cuda is installed. No need if you have a cuda.pc])],
[],
[with_cuda="no"])

AM_CONDITIONAL(HAVE_CUDA, false)
AS_IF([test "x$with_cuda" != "xno"],
	    [if test "x$with_cuda" = "xyes"; then
	    PKG_CHECK_EXISTS([cuda >= 4.0], [cudaexists="yes"],[cudaexists="no"])
    	       if test "x$cudaexists" = "xyes" ; then
       	           PKG_CHECK_MODULES(CUDA, [cuda >= 4.0], [have_cuda="yes"])
	           AC_DEFINE(HAVE_CUDA, 1, [cuda is found])
		   AM_CONDITIONAL(HAVE_CUDA, true)
	       fi
	    fi])

AC_SUBST(CUDA_CFLAGS)
AC_SUBST(CUDA_LIBS)
AC_SUBST(OPENMP_CFLAGS)


dnl # mpi transpose support (included in fftw3)
dnl #tmp_libs="$LIBS"
dnl #AC_CHECK_LIB([fftw3_mpi],[fftw_mpi_local_size_many_transposed],[FFTW_MPI_LIBS=-lfftw3_mpi], [AC_MSG_ERROR([For now the mpi transpose from fftw3 should be used. Please check if the mpi option for fftw3.3 was activated during compilation. If not recompile the lib with --enable-mpi])], [-lfftw3 $MPI_LIBS -lstdc++])
dnl #AC_SUBST(FFTW_MPI_LIBS)
dnl #LIBS="$tmp_libs"
dnl #blas support

mkl_found="no"
check_for_other_blas="yes" # by default if mkl is not found
PKG_CHECK_EXISTS([mkl],mkl_found="yes")

if test "x${mkl_found}" = "xyes"; then
   PKG_CHECK_MODULES(BLAS, [mkl >= 10.0.0.0])
   check_for_other_blas="no"
   mkl_found="yes"
   AC_DEFINE(HAVE_MKL, 1, [mkl library found])
   AC_DEFINE(HAVE_BLAS, 1, [blas library found])
   AC_DEFINE(HAVE_VML, 1, [vml library found])
fi
ipp_found="no"
PKG_CHECK_EXISTS([ipp],ipp_found="yes")
if test "x${ipp_found}" = "xyes"; then
make   PKG_CHECK_MODULES(IPP, [ipp >= 6.0.0.0])
   AC_DEFINE(HAVE_IPP, 1, [ipp library found])
fi

if test "x$check_for_other_blas" = "xyes" ; then
   ACX_BLAS()
fi
   
AC_OUTPUT([Makefile
	src/Makefile
	src/bin/Makefile
	src/lib/Makefile],[])
