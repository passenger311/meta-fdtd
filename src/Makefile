
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

ifeq ($(ENGINE),)
ENGINE=r
endif



M4= m4 max3d.m4

SRC= strings.f90\
constant.f90\
reglist.f90\
buflist.f90\
outlist.f90\
mpiworld.F90\
grid-$(ENGINE).f90\
fdtd-$(ENGINE).f90\
pec-$(ENGINE).f90\
upml-$(ENGINE).f90\
outgpl-fdtd-$(ENGINE).f90\
tfsf-$(ENGINE).f90\
matsource-$(ENGINE).f90\
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
mpiworld.o\
grid-$(ENGINE).o\
fdtd-$(ENGINE).o\
pec-$(ENGINE).o\
upml-$(ENGINE).o\
outgpl-fdtd-$(ENGINE).o\
matsource-$(ENGINE).o\
tfsf-$(ENGINE).o\
mat.o\
outgpl.o\
out.o\
max3d.o

PREOBJ= strings.o\
constant.o\
mpiworld.o\
grid-$(ENGINE).o\
fdtd-$(ENGINE).o\
pre-obj-1.o

PREEXE = pre-obj-1

EXE = max3d-$(ENGINE)

#
default: $(PREEXE) $(EXE)
#


%.o : %.f90
	@echo $(F90) $(F90FLAGS) $(INCLUDE) -c $<
	@$(M4) $< >$*_m4.f90 
	@$(F90) $(F90FLAGS) $(INCLUDE) -c $*_m4.f90 -o $*.o
	@rm $*_m4.f90

%.o : %.F90
	@echo $(F90) $(F90FLAGS) $(INCLUDE) -c $<
	@$(M4) $< >$*_m4.F90 
	@$(F90) $(F90FLAGS) $(INCLUDE) -c $*_m4.F90 -o $*.o
	@rm $*_m4.F90

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
