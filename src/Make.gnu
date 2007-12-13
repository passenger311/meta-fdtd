#-*- Makefile -*-======================================================
#
# Make.default
#
# for gnu >= 4.1 compiler set

# ---------------------------------------------------------------------

M4_ARCH=GCC

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

# --- openmp flags

ifneq ($(OMP),)

OMPFLAGS=-fopenmp

endif 

# --- debug flags

ifneq ($(DBG),)

DBGFLAGS=-g -O0 -ffpe-trap="overflow,invalid,zero"
F90DBGFLAGS= -fbounds-check 

endif 

# --- release flags

RELFLAGS= 

# --- optimization flags

ifeq ($(DBG),)

OPTFLAGS=-O3
NOOPTFLAGS=-O1

endif

# --- flags

FLAGS=$(DBGFLAGS) $(RELFLAGS) $(OMPFLAGS)

F90FLAGS+=$(NOOPTFLAGS) $(FLAGS) $(F90DBGFLAGS)
F90FLAGS_OPT=$(OPTFLAGS) $(FLAGS) $(F90DBGFLAGS) 
CXXFLAGS+=$(NOOPTFLAGS) $(FLAGS)
CXXFLAGS_OPT=$(OPTFLAGS) $(FLAGS)
CCFLAGS+=$(NOOPTFLAGS) $(FLAGS)
CCFLAGS_OPT=$(OPTFLAGS) $(FLAGS)
LDFLAGS+=$(FLAGS)


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================

