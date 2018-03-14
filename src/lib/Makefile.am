MAINTAINERCLEANFILES = Makefile.in
lib_LTLIBRARIES = libSPH.la


include_HEADERS = derivatives.h \
	inoutflow.h \
	inputbp.h \
	input.h \
	output.h \
	parameterfile.h \
	param.h \
	System.h


SPH_SOURCES =  derivatives.c \
	inoutflow.c \
	inputbp.c \
	input.c \
	output.c \
	parameterfile.c \
	System.c \
	Search.cc


libSPH_la_SOURCES = $(SPH_SOURCES)

libSPH_la_CFLAGS = -I$(top_srcdir) -Ideplib/CPU -Ideplib/GPU \
@GSL_CFLAGS@ \
@BLAS_CFLAGS@ \
@IPP_CFLAGS@ \
@OPENMP_CFLAGS@ -I/usr/include/eigen3

libSPH_la_CXXFLAGS = -I$(top_srcdir) -Ideplib/CPU -Ideplib/GPU \
@GSL_CFLAGS@ \
@BLAS_CFLAGS@ \
@IPP_CFLAGS@ \
@OPENMP_CFLAGS@ -I/usr/include/eigen3

libSPH_la_LIBADD = @GSL_LIBS@ \
@BLAS_LIBS@ \
@IPP_LIBS@ \
@OPENMP_CFLAGS@ \
-lnabo -lpthread -lm

libSPH_la_LDFLAGS = -no-undefined -version-info @version_info@
