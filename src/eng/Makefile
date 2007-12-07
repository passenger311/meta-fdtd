
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------


.SUFFIXES: .f90 .F90


# Compilation
# -----------

ifeq ($(ARCH),)
include ../Make.default
else 
include ../Make.$(ARCH)
endif 

ifneq ($(MPI),)
MPI=-mpi
M4DEFINE+= -DM4_MPI=1
endif

ifneq ($(OMP),)
OMP=-omp
M4DEFINE+= -DM4_OMP=1
endif

ifneq ($(CT),)
CT=ct
M4DEFINE+= -DM4_CT=1
endif

ifneq ($(NG),)
NG=ng
M4DEFINE+=-DM4_NG=1
endif

ifneq ($(TE),)
TE=te
M4DEFINE+=-DM4_TE=1
endif

# ENGINE=$(CT)$(NG)$(TE)$(MPI)$(OMP)

M4= m4 $(M4DEFINE) meta3.m4

MATLIST=$(shell cat mat.def | tr "," " " 2>/dev/null)
DIAGLIST=$(shell cat diag.def | tr "," " " 2>/dev/null)

MATOBJ= $(patsubst %,%.o,$(MATLIST))
DIAGOBJ= $(patsubst %,%.o,$(DIAGLIST))


FOBJ= strings.o\
constant.o\
mpiworld.o\
grid.o\
reglist.o\
buflist.o\
outlist.o\
list.o\
fdtd.o\
pec.o\
pml.o\
fdtd_outgpl.o\
$(MATOBJ)\
mat.o\
$(DIAGOBJ)\
diag.o\
outgpl.o\
out.o\
meta3.o


EXE = meta3

default: $(EXE)

%.o : %.f90
	$(M4) $< >$*_m4.f90 
	$(F90) $(F90FLAGS) $(INCLUDE) -c $*_m4.f90 -o $*.o

%.o : %.F90
	$(M4) $< >$*_m4.F90 
	$(F90) $(F90FLAGS) $(INCLUDE) -c $*_m4.F90 -o $*.o

%.o : %.c
	 $(CC) $(CCFLAGS) $(INCLUDE) -c $<


$(EXE) : $(FOBJ)
	$(F90) $(LDFLAGS) $(FOBJ) -o $(EXE)

.PHONY: clean proper

clean:
	 rm -f  core* a.out *.L *.par *.log *.list *.diag *.txt *.mod *_m4.f90 *_m4.F90

proper: clean
	 rm -f $(EXE) $(PREEXE) *.o *.in *.out 


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
