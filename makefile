#Change upon moving folder
PETSC_DIR=../PSC-FMM/petsc-3.2-p2
PETSC_ARCH=arch-linux2-cxx-debug

# Hack for linker... (ARMADILLO)
#export LD_LIBRARY_PATH=/home/pletou1/Dropbox/Stanford/ResearchPSC/C++/PSC-FMM-ELASTIC/Armadillo/usr/lib64/
export LD_LIBRARY_PATH=/home/pletou1/Dropbox/Stanford/ResearchPSC/C++/PSC-FMM-ELASTIC/Armadillo/usr/lib64/

#PETSC_DIR=/usr/local/petsc-3.2-p2/

CFLAGS	         = -O3 -Wall -pedantic -I/home/pletou1/Dropbox/Stanford/ResearchPSC/C++/gsl-1.16/include
FFLAGS	         = 
CPPFLAGS         = $(CFLAGS) -pg
FPPFLAGS         =
LOCDIR           = 

#CLINKER          = g++
FC=gfortran

CLEANFILES       =
NP               = 1

OBJFILES         = Coordinates.o Indexing.o Scatterer.o Distributions.o PSCEnv.o Transfer.o Translation_Function.o Quadrature.o SphericalHarmonicsTransform.o TransferUtils.o Transfer_Function.o Reterpolator.o Level.o NTree.o MLFMM_Env.o IncWave.o Solver.o L_choice.o Finalize.o Imaging.o

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

##################################################
# Paths and Directories
LIBDIR	= -L../AmosBessel/ -L../gsl-1.16/lib/ -L./Armadillo/usr/lib64/
INCLDIR	= -I/home/pletou1/Dropbox/Stanford/ResearchPSC/C++/gsl-1.16/include -I./libAmosBessel/ -I./Armadillo/usr/include/ -I./LAPACK++/include/

##################################################
# Libraries
LIBS 	= -lgfortran -g77libs -lfftw3 -larmadillo -llapack -lblas -lm ../gsl-1.16/lib/libgsl.a ../gsl-1.16/lib/libgslcblas.a 

#../AmosBessel/libAmosBessel.a

##################################################
# Executable 

Main: Main.o  chkopts
	-${CLINKER} -o Main $(OBJFILES) Main.o   $(LIBDIR) $(INCLDIR) ${LIBS}
	${RM} Main.o

HK: HelmholtzKirchoff.o chkopts
	-${CLINKER} -o HelmholtzKirchoff $(OBJFILES) HelmholtzKirchoff.o   $(LIBDIR) $(INCLDIR) ${LIBS}
	${RM} HelmholtzKirchoff.o

DI: DisplacementImaging.o chkopts
	-${CLINKER} -o DisplacementImaging $(OBJFILES) DisplacementImaging.o   $(LIBDIR) $(INCLDIR) ${LIBS}
	${RM} DisplacementImaging.o

SR: SuperRes.o chkopts
	-${CLINKER} -o SuperRes $(OBJFILES) SuperRes.o   $(LIBDIR) $(INCLDIR) ${LIBS}
	${RM} SuperRes.o


# Core files
Coordinates: Coordinates.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Coordinates.cpp

Indexing: Indexing.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Indexing.cpp

Scatterer: Scatterer.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Scatterer.cpp

Distributions: Distributions.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Distributions.cpp

Transfer: Transfer.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Transfer.cpp

IncWave: IncWave.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) IncWave.cpp

PSCEnv: PSCEnv.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) PSCEnv.cpp

Solver: Solver.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Solver.cpp

L_choice: L_choice.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) L_choice.cpp

Finalize: Finalize.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Finalize.cpp

Imaging: Imaging.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Imaging.cpp



# FMMPS files
Translation_Function: ./FMMPS/Translation_Function.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/Translation_Function.cpp

Quadrature: ./FMMPS/Quadrature.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/Quadrature.cpp

SHT: ./FMMPS/SphericalHarmonicsTransform.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/SphericalHarmonicsTransform.cpp

TransferUtils: ./FMMPS/TransferUtils.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/TransferUtils.cpp

Transfer_Function: ./FMMPS/Transfer_Function.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/Transfer_Function.cpp

Reterpolator: ./FMMPS/Reterpolator.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/Reterpolator.cpp

Level: ./FMMPS/Level.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/Level.cpp

NTree: ./FMMPS/NTree.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/NTree.cpp

MLFMM_Env: ./FMMPS/MLFMM_Env.o  chkopts
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/MLFMM_Env.cpp



# Groups
Core: Coordinates Indexing Distributions Scatterer IncWave Transfer PSCEnv Solver L_choice Finalize Imaging

FMM: Translation_Function Quadrature SHT TransferUtils Transfer_Function Reterpolator Level NTree MLFMM_Env

all: Core FMM Main