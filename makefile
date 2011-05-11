#################################################################################
# Copyright (C) 2009-2011 Technische Universitaet Muenchen                      #
# This file is part of the SG++ project. For conditions of distribution and     #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp            #
#                                                                               #
# author Alexander Heinecke (Alexander.Heinecke@mytum.de)                       #
#################################################################################

###################################################################
# Needed Pathes
###################################################################
SRCDIR=./../../../src/sgpp
#only for extensions:
#####################
# Intel Array Building Blocks
ARBBINCLUDE = /opt/intel/arbb/1.0.0.015/include
ARBBLIB = /opt/intel/arbb/1.0.0.015/lib/intel64
# NVidia OpenCL
OCLINCLUDE = /opt/cuda/include
OCLLIB = /usr/lib64

###################################################################
# Default Variables, overwirtten by CLI
###################################################################	
# use OpenMP Version 3
OMP=0
# use the TR1 Implementations for Hashmaps
TR1=0
# default compiler: g++; possible values: g++, icpc (Intel Compiler)
CC=g++
#CC=icpc
# vectorization option for intel compiler
VEC=sse3
# extensions, manages extensions to be included, possible values (only when using Intel Compiler):
#	ArBB - Intel Array Building Blocks support
#	OCL - NVIDIA OpenCL support
#	MPI - Intel MPI support
#	NO - no extensions, default
EXT=NO

###################################################################
# Compiler Flags
###################################################################	
CFLAGS_GCC:=-Wall -pedantic -ansi -c -Wno-long-long -fno-strict-aliasing -O3 -funroll-loops -ffloat-store -I$(SRCDIR) -DSG_BASE -DSG_PDE -DSG_DATADRIVEN -DSG_SOLVER -DSG_FINANCE -DSG_PARALLEL -DSG_COMBIGRID
LFLAGS_GCC:=-Wall -pedantic -ansi -O3

CFLAGS_ICC:=-Wall -ansi -c -fno-strict-aliasing -ipo -ip -ansi-alias -O3 -funroll-loops -I$(SRCDIR) -DSG_BASE -DSG_PDE -DSG_DATADRIVEN -DSG_SOLVER -DSG_FINANCE -DSG_PARALLEL -DSG_COMBIGRID
LFLAGS_ICC:=-Wall -ansi -O3 -static-intel -ipo -ip

ifeq ($(CC),g++)
CFLAGS:=$(CFLAGS_GCC)
LFLAGS:=$(LFLAGS_GCC)
EXT=NO
ifeq ($(OMP),1)
CFLAGS:=$(CFLAGS) -fopenmp
LFLAGS:=$(LFLAGS) -fopenmp
endif
ifeq ($(TR1),1)
CFLAGS:=$(CFLAGS) -DUSETRONE -std=c++0x
endif
ifeq ($(EXT), OCL)
CFLAGS:=$(CFLAGS) -I$(OCLINCLUDE) -DUSEOCL -fopenmp
LFLAGS:=$(LFLAGS) -L$(OCLLIB) -lOpenCL -fopenmp
endif
endif

ifeq ($(CC),icpc)
CFLAGS:=$(CFLAGS_ICC)
LFLAGS:=$(LFLAGS_ICC)
ifeq ($(VEC),sse3)
CFLAGS:=$(CFLAGS) -xSSE3
endif
ifeq ($(VEC),sse4)
CFLAGS:=$(CFLAGS) -xSSE4.2
endif
ifeq ($(VEC),avx)
CFLAGS:=$(CFLAGS) -xAVX -DUSEAVX
endif
ifeq ($(OMP),1)
CFLAGS:=$(CFLAGS) -openmp
LFLAGS:=$(LFLAGS) -openmp
endif
ifeq ($(TR1),1)
CFLAGS:=$(CFLAGS) -DUSETRONE -std=c++0x
endif
ifeq ($(EXT), ArBB)
CFLAGS:=$(CFLAGS) -I$(ARBBINCLUDE) -DUSEARBB
LFLAGS:=$(LFLAGS) -L$(ARBBLIB) -larbb -ltbb
endif
ifeq ($(EXT), OCL)
CFLAGS:=$(CFLAGS) -I$(OCLINCLUDE) -DUSEOCL -openmp
LFLAGS:=$(LFLAGS) -L$(OCLLIB) -lOpenCL -openmp
endif
endif

ifeq ($(CC),mpiicpc)
CFLAGS:=$(CFLAGS_ICC)
LFLAGS:=$(LFLAGS_ICC)
CFLAGS:=$(CFLAGS) -DUSE_MPI
EXT=MPI
ifeq ($(VEC),sse3)
CFLAGS:=$(CFLAGS) -xSSE3
endif
ifeq ($(VEC),sse4)
CFLAGS:=$(CFLAGS) -xSSE4.2
endif
ifeq ($(VEC),avx)
CFLAGS:=$(CFLAGS) -xAVX -DUSEAVX
endif
ifeq ($(OMP),1)
CFLAGS:=$(CFLAGS) -openmp
LFLAGS:=$(LFLAGS) -openmp
endif
endif

###################################################################
# Builds a lib containing all SG Algorithms
###################################################################	
default:
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/sgpplib_gcc
	make -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/sgpplib_icc
	make -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "EXT=$(EXT)"
endif
ifeq ($(CC),mpiicpc)
	mkdir -p tmp/build_native/sgpplib_mpiicc
	make -f ./../../../src/makefileSGppLIB --directory=./tmp/build_native/sgpplib_mpiicc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpiicc.a" "EXT=$(EXT)"
endif


###################################################################
# Builds a Balck Scholes Solver
###################################################################	
BSSolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/BSSolver_gcc
	make -f ./../../../src/makefileNativeBlackScholesSolver --directory=./tmp/build_native/BSSolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=BSSolver_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/BSSolver_icc
	make -f ./../../../src/makefileNativeBlackScholesSolver --directory=./tmp/build_native/BSSolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=BSSolver_ICC" "EXT=$(EXT)"
endif


###################################################################
# Builds a Black Scholes Solver with Stretching
###################################################################	
BSSolverWithStretching: default

ifeq ($(CC),g++)
	mkdir -p tmp/build_native/BSSolverWithStretching_gcc
	make -f ./../../../src/makefileNativeBlackScholesSolverWithStretching --directory=./tmp/build_native/BSSolverWithStretching_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=BSSolverWithStretching_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/BSSolverWithStretching_icc
	make -f ./../../../src/makefileNativeBlackScholesSolverWithStretching --directory=./tmp/build_native/BSSolverWithStretching_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=BSSolverWithStretching_ICC" "EXT=$(EXT)"
endif

###################################################################
# Builds a Hull White Solver
###################################################################	
HWSolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/HWSolver_gcc
	make -f ./../../../src/makefileNativeHullWhiteSolver --directory=./tmp/build_native/HWSolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=HWSolver_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/HWSolver_icc
	make -f ./../../../src/makefileNativeHullWhiteSolver --directory=./tmp/build_native/HWSolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=HWSolver_ICC" "EXT=$(EXT)"
endif

###################################################################
# Builds a Hull White combine Black Scholes Solver
###################################################################	
BSHWSolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/BSHWSolver_gcc
	make -f ./../../../src/makefileNativeBSHWSolver --directory=./tmp/build_native/BSHWSolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=BSHWSolver_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/BSHWSolver_icc
	make -f ./../../../src/makefileNativeBSHWSolver --directory=./tmp/build_native/BSHWSolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=BSHWSolver_ICC" "EXT=$(EXT)"
endif

###################################################################
# Builds a simple Heat Equation Solver
###################################################################	
HESolver: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/HESolver_gcc
	make -f ./../../../src/makefileNativeHeatEquationSolver --directory=./tmp/build_native/HESolver_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=HESolver_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/HESolver_icc
	make -f ./../../../src/makefileNativeHeatEquationSolver --directory=./tmp/build_native/HESolver_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=HESolver_ICC" "EXT=$(EXT)"
endif
ifeq ($(CC),mpiicpc)
	mkdir -p tmp/build_native/HESolver_mpiicc
	make -f ./../../../src/makefileNativeHeatEquationSolverMPI --directory=./tmp/build_native/HESolver_mpiicc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_mpiicc.a" "BINNAME=HESolver_ICC_MPI" "EXT=$(EXT)"
endif

###################################################################
# Builds a simple Heat Equation Solver (rotating Laser test case)
###################################################################	
LaserHESolver2D: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/LaserHESolver2D_gcc
	make -f ./../../../src/makefileNativeLaserHeatEquationSolver --directory=./tmp/build_native/LaserHESolver2D_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=LaserHESolver2D_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/LaserHESolver2D_icc
	make -f ./../../../src/makefileNativeLaserHeatEquationSolver --directory=./tmp/build_native/LaserHESolver2D_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=LaserHESolver2D_ICC" "EXT=$(EXT)"
endif

###################################################################
# Builds a simple Heat Equation Solver with Stretching
###################################################################	
HESolverWithStretching: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/HESolverWithStretching_gcc
	make -f ./../../../src/makefileNativeHeatEquationSolverWithStretching --directory=./tmp/build_native/HESolverWithStretching_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=HESolverWithStretching_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/HESolverWithStretching_icc
	make -f ./../../../src/makefileNativeHeatEquationSolverWithStretching --directory=./tmp/build_native/HESolverWithStretching_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=HESolverWithStretching_ICC" "EXT=$(EXT)"
endif

###################################################################
# Builds a ClassifyBenchmark Application
###################################################################	
ClassifyBenchmark: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/ClassifyBenchmark_gcc
	make -f ./../../../src/makefileNativeClassifyBenchmark --directory=./tmp/build_native/ClassifyBenchmark_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=ClassifyBenchmark_GCC" "EXT=$(EXT)"
endif
ifeq ($(CC),icpc)
	mkdir -p tmp/build_native/ClassifyBenchmark_icc
	make -f ./../../../src/makefileNativeClassifyBenchmark --directory=./tmp/build_native/ClassifyBenchmark_icc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_icc.a" "BINNAME=ClassifyBenchmark_ICC" "EXT=$(EXT)"
endif

###################################################################
# Builds a Up/Down Test Application
###################################################################	
UpDownTest: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/UpDownTest_gcc
	make -f ./../../../src/makefileUpDownTest --directory=./tmp/build_native/UpDownTest_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=UpDownTest_GCC" "EXT=$(EXT)"
endif

###################################################################
# Builds a Up/Down Test Application
###################################################################	
UpDownTestForStretching: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/UpDownTestForStretching_gcc
	make -f ./../../../src/makefileUpDownTestForStretching --directory=./tmp/build_native/UpDownTestForStretching_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=UpDownTestForStretching_GCC" "EXT=$(EXT)"
endif

###################################################################
# Builds a Refine/Coarsen Test Application
###################################################################	
RefineCoarsenTest: default
ifeq ($(CC),g++)
	mkdir -p tmp/build_native/RefineCoarsen_gcc
	make -f ./../../../src/makefileRefineCoarsenTest --directory=./tmp/build_native/RefineCoarsen_gcc "CC=$(CC)" "CFLAGS=$(CFLAGS)" "LFLAGS=$(LFLAGS)" "LIBNAME=libsgpp_gcc.a" "BINNAME=RefineCoarsen_gcc" "EXT=$(EXT)"
endif

###################################################################
# test Balck Scholes Solver
###################################################################	
		
test_BS_1d:
	cd bin; \
	./copyBSSolverToTest.sh; \
	cd ./../tests/CPP_Apps/BSSolver/1d; \
	./test_BSSolver_1d.sh;
	
test_BS_2d:
	cd bin; \
	./copyBSSolverToTest.sh; \
	cd ./../tests/CPP_Apps/BSSolver/2d; \
	./test_BSSolver_2d.sh;
		
test_BS_3d:
	cd bin; \
	./copyBSSolverToTest.sh; \
	cd ./../tests/CPP_Apps/BSSolver/3d; \
	./test_BSSolver_3d.sh;
	
test_BS_all: test_BS_1d test_BS_2d test_BS_3d
	echo "executed all BS tests!"

###################################################################
# test Black Scholes Solver with Stretching
###################################################################	
		
test_BSS_1d:
	cd bin; \
	./copyBSSolverWithStretchingToTest.sh; \
	cd ./../tests/CPP_Apps/BSSolverWithStretching/1d; \
	./test_BSSolverWithStretching_1d.sh;
	
test_BSS_2d:
	cd bin; \
	./copyBSSolverWithStretchingToTest.sh; \
	cd ./../tests/CPP_Apps/BSSolverWithStretching/2d; \
	./test_BSSolverWithStretching_2d.sh;
		
test_BSS_3d:
	cd bin; \
	./copyBSSolverWithStretchingToTest.sh; \
	cd ./../tests/CPP_Apps/BSSolverWithStretching/3d; \
	./test_BSSolverWithStretching_3d.sh;
	
test_BSS_all: test_BSS_1d test_BSS_2d test_BSS_3d
	echo "executed all BS tests!"
		
###################################################################
# test Combined Hull Wihte Solver Solver
###################################################################			

test_BSHW:
	cd bin; \
	./copyBSHWSolverToTest.sh; \
	cd ./../tests/CPP_Apps/BSHWSolver; \
	./test_BSHWSolver_2d_cart.sh;
		
###################################################################
# test ClassifyBenchmark
###################################################################	

test_ClassifyBenchmark:
	cd bin; \
	./copyClassifyBenchmarkToTest.sh; \
	cd ./../tests/CPP_Apps/ClassifyBenchmark; \
	./test_ClassifyBenchmark.sh;

clean:
	rm -rdfv tmp/build_native