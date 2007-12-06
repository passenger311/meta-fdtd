
#======================================================================
#
# Makefile
#

.SUFFIXES: .f90 .F90


# Compilation
# -----------

ifeq ($(ARCH),)
include Make.default
else 
include Make.$(ARCH)
endif 

ifneq ($(MPI),)
MPI=-mpi
M4DEFINE+= -DM4_MPI=1
endif

ifneq ($(OMP),)
OMP=-omp
M4DEFINE+= -DM4_OMP=1
endif

ifneq ($(CMPLX),)
CT=ct
M4DEFINE+= -DM4_CMPLX=1
endif

ifneq ($(NGRID),)
NG=ng
M4DEFINE+=-DM4_NGRID=1
endif

ifneq ($(TEPS),)
TE=te
M4DEFINE+=-DM4_TEPS=1
endif

ENGINE=$(CT)$(NG)$(TE)$(MPI)$(OMP)

M4= m4 $(M4DEFINE) meta3.m4

SRC= strings.f90\
constant.f90\
reglist.f90\
buflist.f90\
outlist.f90\
list.f90\
mpiworld.F90\
grid.f90\
fdtd.f90\
pec.f90\
pml.f90\
fdtd_outgpl.f90\
tfsf.f90\
matsource.f90\
pre-obj-1.f90\
mat.f90\
outgpl.f90\
out.f90\
meta3.F90

FOBJ= strings.o\
constant.o\
reglist.o\
buflist.o\
outlist.o\
list.o\
mpiworld.o\
grid.o\
fdtd.o\
pec.o\
pml.o\
fdtd_outgpl.o\
matsource.o\
tfsf.o\
mat.o\
outgpl.o\
out.o\
meta3.o

PREOBJ= strings.o\
constant.o\
mpiworld.o\
grid.o\
fdtd.o\
pre-obj-1.o

PREEXE = pre-obj-1

EXE = meta3-$(ENGINE)

#
default: $(PREEXE) $(EXE)
#


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

$(PREEXE) : $(PREOBJ)
	$(F90) $(LDFLAGS) $(PREOBJ) -o $(PREEXE)


.PHONY: clean proper

clean:
	 rm -f  core* a.out *.L *.par *.log *.list *.diag *.txt *.mod *_m4.f90 *_m4.F90

proper: clean
	 rm -f $(EXE) $(PREEXE) *.o *.in *.out 


#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
