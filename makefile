export LD_LIBRARY_PATH=/home/pletou1/Dropbox/Stanford/ResearchPSC/C++/PSC-FMM-ELASTIC/Armadillo/usr/lib64/

CFLAGS	         = -O3 -Wall -pedantic -I/home/pletou1/Dropbox/Stanford/ResearchPSC/C++/gsl-1.16/include -I/usr/include/mpi
FFLAGS	         = 
CPPFLAGS         = $(CFLAGS) -pg
FPPFLAGS         =
LOCDIR           = 

CLINKER          = mpicxx #g++
FC               = gfortran

CLEANFILES       =
NP               = 1

OBJFILES         = Coordinates.o Indexing.o Scatterer.o Distributions.o PSCEnv.o Transfer.o Translation_Function.o Quadrature.o SphericalHarmonicsTransform.o TransferUtils.o Transfer_Function.o Reterpolator.o Level.o NTree.o MLFMM_Env.o IncWave.o Solver.o L_choice.o Finalize.o Imaging.o


##################################################
# Paths and Directories
LIBDIR	= -L../AmosBessel/ -L../gsl-1.16/lib/ -L./Armadillo/usr/lib64/
INCLDIR	= -I/home/pletou1/Dropbox/Stanford/ResearchPSC/C++/gsl-1.16/include -I./libAmosBessel/ -I./Armadillo/usr/include/ -I./LAPACK++/include/ -I/usr/include/mpi

##################################################
# Libraries
LIBS 	= -lgfortran -g77libs -lfftw3 -larmadillo -llapack -lblas -lm ./gsl-1.16/lib/libgsl.a ./gsl-1.16/lib/libgslcblas.a 
#-lmpichcxx -lmpich -lopa -lmpl -lrt -lcr -lpthread

##################################################
# Executable 

Main: $(OBJFILES) Main.o
	${CLINKER} ${CFLAGS} $(OBJFILES) Main.o $(LIBDIR) $(INCLDIR) ${LIBS} -o Main
	${RM} Main.o

HK: $(OBJFILES) HelmholtzKirchoff.o 
	${CLINKER} $(OBJFILES) HelmholtzKirchoff.o   $(LIBDIR) $(INCLDIR) ${LIBS} -o HelmholtzKirchoff
	${RM} HelmholtzKirchoff.o

DI: $(OBJFILES) DisplacementImaging.o 
	${CLINKER} ${CFLAGS} $(OBJFILES) DisplacementImaging.o   $(LIBDIR) $(INCLDIR) ${LIBS} -o DisplacementImaging
	${RM} DisplacementImaging.o

MF: $(OBJFILES) MultiFreq.o
	${CLINKER} ${CFLAGS} $(OBJFILES) MultiFreq.o $(LIBDIR) $(INCLDIR) ${LIBS} -o MultiFreq
	${RM} MultiFreq.o

SR: SuperRes.o -${CLINKER} -o SuperRes $(OBJFILES) SuperRes.o   $(LIBDIR) $(INCLDIR) ${LIBS}
	${RM} SuperRes.o

EC: $(OBJFILES) ErrorCheck.o
	${CLINKER} ${CFLAGS} $(OBJFILES) ErrorCheck.o $(LIBDIR) $(INCLDIR) ${LIBS} -o ErrorCheck
	${RM} ErrorCheck.o

FFH: $(OBJFILES) FarFieldHomogenization.o
	${CLINKER} ${CFLAGS} $(OBJFILES) FarFieldHomogenization.o $(LIBDIR) $(INCLDIR) ${LIBS} -o FarFieldHomogenization
	${RM} FarFieldHomogenization.o

LS: $(OBJFILES) LS_check.o
	${CLINKER} ${CFLAGS} $(OBJFILES) LS_check.o $(LIBDIR) $(INCLDIR) ${LIBS} -o LS_check
	${RM} LS_check.o

FLC: $(OBJFILES) FoldyLaxCheck.o
	${CLINKER} ${CFLAGS} $(OBJFILES) FoldyLaxCheck.o $(LIBDIR) $(INCLDIR) ${LIBS} -o FoldyLaxCheck
	${RM} FoldyLaxCheck.o



# Core files
Coordinates: Coordinates.o
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Coordinates.cpp

Indexing: Indexing.o
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Indexing.cpp

Scatterer: Scatterer.o
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Scatterer.cpp

Distributions: Distributions.o
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Distributions.cpp

Transfer: Transfer.o
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Transfer.cpp

IncWave: IncWave.o
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) IncWave.cpp

PSCEnv: PSCEnv.o
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) PSCEnv.cpp

Solver: Solver.o
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Solver.cpp

L_choice: L_choice.o
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) L_choice.cpp

Finalize: Finalize.o
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Finalize.cpp

Imaging: Imaging.o
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) Imaging.cpp



# FMMPS files
Translation_Function: ./FMMPS/Translation_Function.cpp
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/Translation_Function.cpp

Quadrature: ./FMMPS/Quadrature.cpp
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/Quadrature.cpp

SHT: ./FMMPS/SphericalHarmonicsTransform.cpp
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/SphericalHarmonicsTransform.cpp

TransferUtils: ./FMMPS/TransferUtils.cpp
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/TransferUtils.cpp

Transfer_Function: ./FMMPS/Transfer_Function.cpp
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/Transfer_Function.cpp

Reterpolator: ./FMMPS/Reterpolator.cpp
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/Reterpolator.cpp

Level: ./FMMPS/Level.cpp
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/Level.cpp

NTree: ./FMMPS/NTree.cpp
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/NTree.cpp

MLFMM_Env: ./FMMPS/MLFMM_Env.cpp
	-${CLINKER} -c $(LIBDIR) $(INCLDIR) ./FMMPS/MLFMM_Env.cpp



# Groups
Core: Coordinates Indexing Distributions Scatterer IncWave Transfer PSCEnv Solver L_choice Finalize Imaging

FMM: Translation_Function Quadrature SHT TransferUtils Transfer_Function Reterpolator Level NTree MLFMM_Env

all: Core FMM Main
