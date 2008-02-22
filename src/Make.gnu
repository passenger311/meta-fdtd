#-*- Makefile -*-======================================================
#
# Make.default
#
# for gnu >= 4.1 compiler set

# ---------------------------------------------------------------------

M4_ARCH=GNU

# --- compiler set

ifeq ($(MPI),)

F90=gfortran
CXX=g++
CC=gcc

else
 
F90=mpif90 
CXX=mpicxx
CC=mpicc

endif

# --- release flags

RELFLAGS=-static

# --- openmp flags

ifneq ($(OMP),)

OMPFLAGS=-fopenmp

endif 

# --- debug flags

ifneq ($(DBG),)

DBGFLAGS=-g -O0 -ffpe-trap="overflow,invalid,zero"
F90DBGFLAGS= -fbounds-check 

endif 

# --- optimization flags

ifeq ($(DBG),)

OPTFLAGS=-O3 -ffast-math -funroll-loops -mtune=i686
NOOPTFLAGS=-O1
LDOPTFLAGS=

endif

# --- flags

F90RELFLAGS= 

FLAGS=$(DBGFLAGS) $(RELFLAGS) $(OMPFLAGS)

F90FLAGS+= $(FLAGS) $(F90RELFLAGS) $(F90DBGFLAGS)
CXXFLAGS+= $(FLAGS)
CFLAGS+= $(FLAGS)
LDFLAGS+= $(FLAGS) $(LDOPTFLAGS)

F90FLAGS_NOOPT=$(NOOPTFLAGS) $(F90FLAGS)
F90FLAGS_OPT=$(OPTFLAGS) $(F90FLAGS)
CXXFLAGS_NOOPT=$(NOOPTFLAGS) $(CXXFLAGS)
CXXFLAGS_OPT=$(OPTFLAGS) $(CXXFLAGS)
CFLAGS_NOOPT=$(NOOPTFLAGS) $(CFLAGS)
CFLAGS_OPT=$(OPTFLAGS) $(CFLAGS)


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================

