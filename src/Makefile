
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

ifeq ($(CMPLX),)
M4FT="real(kind=8)"
FT=r
else
M4FT="complex(kind=8)"
FT=c
F90FLAGS+= -DCMPLX
endif


M4= m4 max3d.m4 -DM4_FT=$(M4FT) 

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
outgpl-fdtd.f90\
tfsf.f90\
matsource.f90\
pre-obj-1.f90\
mat.f90\
outgpl.f90\
out.f90\
max3d.F90

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
outgpl-fdtd.o\
matsource.o\
tfsf.o\
mat.o\
outgpl.o\
out.o\
max3d.o

PREOBJ= strings.o\
constant.o\
mpiworld.o\
grid.o\
fdtd.o\
pre-obj-1.o

PREEXE = pre-obj-1

EXE = max3d

#
default: $(PREEXE) $(EXE)
#


%.o : %.f90
	@echo $(F90) $(F90FLAGS) $(INCLUDE) -c $<
	@$(M4) $< >$*_m4.f90 
	@$(F90) $(F90FLAGS) $(INCLUDE) -c $*_m4.f90 -o $*.o

%.o : %.F90
	@echo $(F90) $(F90FLAGS) $(INCLUDE) -c $<
	@$(M4) $< >$*_m4.F90 
	@$(F90) $(F90FLAGS) $(INCLUDE) -c $*_m4.F90 -o $*.o

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
