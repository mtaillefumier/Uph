MAINTAINERCLEANFILES = Makefile.in
lib_LTLIBRARIES = libMC.la 

if HAVE_CUDA
.cu.lo:
	$(LIBTOOL) --tag=CC --mode=compile  nvcc  -v -O2 -g  -arch sm_21 -lineinfo --use_fast_math --ptxas-options=-v -ccbin $(CXX) --compiler-options `echo -fPIC -DPIC $(libMC_la_CXXFLAGS)  -DMPICH_SKIP_MPICXX -DMPICH_IGNORE_CXX_SEEK | sed -e "s/ /,/g"` -c $^ -static  
endif

include_HEADERS = MonteCarlo.hpp \
SpinConfiguration.hpp \
Transpose.hpp \
StructFactor.hpp \
Graph.hpp \
AutoCorrelation.hpp \
DynamicMeasurement.hpp \
StaticMeasurement.hpp \
TimeIntegration.hpp \
Hamiltonian.hpp \
fft.hpp \
utils.hpp \
functions.h \
Parameters.hpp \
Coplanarity.hpp \
Spin.hpp \
MonteCarloSteps.hpp \
StructFactorPowder.hpp \
Cell.hpp \
Neighbors.hpp \
Quadripole.hpp \
SpinSpinCorrelations.hpp \
Coloring.hpp \
Lattice.hpp \
OrderParameters.hpp

MC_SOURCES = MonteCarlo.cc \
MonteCarlo2.cc \
SpinConfiguration.cc \
StructFactor.cc \
Graph.cc \
AutoCorrelation.cc \
DynamicMeasurement.cc \
StaticMeasurement.cc \
TimeIntegration.cc \
Hamiltonian.cc \
fft.cc \
utils.cc \
Parameters.cc \
Coplanarity.cc \
Spin.cc \
MonteCarloSteps.cc \
StructFactorPowder.cc \
Cell.cc \
Neighbors.cc \
Transpose.cc \
Quadripole.cc \
SpinSpinCorrelations.cc \
Coloring.cc \
RungeKutta.cc \
Math.cc \
functions.c \
dSFMT.c \
Lattice.cc \
OrderParameters.cc 

OTHERS_SOURCES = Hamiltonian.cu \
		 Math_cuda.cu \
		 MonteCarlo_cuda.cu \
	         fft.cu \
		 Coloring.cu

if HAVE_CUDA
MC_SOURCES += cuda.cu
endif

libMC_la_SOURCES = $(MC_SOURCES)

if HAVE_CUDA
libMC_la_CFLAGS = -I$(top_srcdir) \
@GSL_CFLAGS@ \
@BLAS_CFLAGS@ \
@MPI_CFLAGS@ \
@HDF5_CFLAGS@ \
@FFTW_CFLAGS@ \
@IPP_CFLAGS@ \
@XML2_CFLAGS@ \
@CUDA_CFLAGS@
libMC_la_CXXFLAGS = -I$(top_srcdir) \
@GSL_CFLAGS@ \
@BLAS_CFLAGS@ \
@MPI_CFLAGS@ \
@HDF5_CFLAGS@ \
@FFTW_CFLAGS@ \
@IPP_CFLAGS@ \
@XML2_CFLAGS@ \
@CUDA_CFLAGS@ 
libMC_la_LIBADD = @GSL_LIBS@ \
@FFTW_LIBS@ \
@FFTW_MPI_LIBS@ \
@HDF5_LIBS@  \
@MPI_LIBS@ \
@BLAS_LIBS@ \
@IPP_LIBS@ \
@XML2_LIBS@ \
@CUDA_LIBS@ \
-lm 
else
libMC_la_CFLAGS = -I$(top_srcdir) \
@GSL_CFLAGS@ \
@BLAS_CFLAGS@ \
@MPI_CFLAGS@ \
@HDF5_CFLAGS@ \
@FFTW_CFLAGS@ \
@IPP_CFLAGS@ \
@XML2_CFLAGS@
libMC_la_CXXFLAGS =  -I$(top_srcdir) \
@GSL_CFLAGS@ \
@BLAS_CFLAGS@ \
@MPI_CFLAGS@ \
@HDF5_CFLAGS@ \
@FFTW_CFLAGS@ \
@IPP_CFLAGS@ \
@XML2_CFLAGS@
libMC_la_LIBADD = @GSL_LIBS@ \
@FFTW_LIBS@ \
@FFTW_MPI_LIBS@ \
@MPI_LIBS@ \
@HDF5_LIBS@  \
@BLAS_LIBS@ \
@IPP_LIBS@ \
@XML2_LIBS@ \
-lm 
endif

libMC_la_LDFLAGS = -no-undefined -version-info @version_info@ 

cuda.cu : $(OTHERS_SOURCES)
	rm -f cuda.cu
	echo "#include <mpi.h>" > cuda.cu
	echo "#undef _GLIBCXX_ATOMIC_BUILTINS" >> cuda.cu
	echo "#undef _GLIBCXX_USE_INT128" >> cuda.cu
	echo "#include \"Math_cuda.cu\" " >> cuda.cu
	echo "#include \"Hamiltonian.cu\" " >> cuda.cu
	echo "#include \"MonteCarlo_cuda.cu\" " >> cuda.cu
	echo "#include \"fft.cu\" " >> cuda.cu
	echo "#include \"Coloring.cu\" " >> cuda.cu
	echo "#include \"header.hu\" " >> cuda.cu
