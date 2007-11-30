
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


SRC= strings.f90\
constant.f90\
region.f90\
mpiworld.F90\
outgpl.f90\
grid-$(ENGINE).f90\
fdtd-$(ENGINE).f90\
pec-$(ENGINE).f90\
upml-$(ENGINE).f90\
fdtd-outgpl-$(ENGINE).f90\
tfsf-$(ENGINE).f90\
src-point-$(ENGINE).f90\
pre-obj-1.f90\
max3d.F90

FOBJ= strings.o\
constant.o\
region.o\
mpiworld.o\
outgpl.o\
grid-$(ENGINE).o\
fdtd-$(ENGINE).o\
pec-$(ENGINE).o\
upml-$(ENGINE).o\
fdtd-outgpl-$(ENGINE).o\
src-point-$(ENGINE).o\
tfsf-$(ENGINE).o\
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
	 $(F90) $(F90FLAGS) $(INCLUDE) -c $<

%.o : %.F90
	 $(F90) $(F90FLAGS) $(INCLUDE) -c $<

%.o : %.c
	 $(CC) $(CCFLAGS) $(INCLUDE) -c $<


$(EXE) : $(FOBJ)
	$(F90) $(LDFLAGS) $(FOBJ) -o $(EXE)

$(PREEXE) : $(PREOBJ)
	$(F90) $(LDFLAGS) $(PREOBJ) -o $(PREEXE)


.PHONY: clean proper

clean:
	 rm -f  core* a.out *.L *.par *.log *.list *.diag *.txt *.mod

proper: clean
	 rm -f $(EXE) $(PREEXE) *.o *.in *.out 

