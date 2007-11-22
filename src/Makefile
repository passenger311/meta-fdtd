
.SUFFIXES: .f90 .F90


# Compilation
# -----------

ifeq ($(ARCH),)
include Make.default
else 
include Make.$(ARCH)
endif 

SRC= strings.f90\
constant.f90\
grid-r.f90\
bound-r.f90\
fdtd3d-r.f90\
upml-r.f90\
config.f90\
output-r.f90\
tfsf-r.f90\
pre-obj-1.f90\
pointsource.f90\
mpistart.F90\
max3d.F90

FOBJ= strings.o\
constant.o\
mpistart.o\
grid-r.o\
bound-r.o\
fdtd3d-r.o\
upml-r.o\
config.o\
output-r.o\
pointsource.o\
tfsf-r.o\
max3d.o

PREOBJ= strings.o\
constant.o\
mpistart.o\
grid-r.o\
upml-r.o\
config.o\
pre-obj-1.o

PREEXE = pre-obj-1

EXE = max3d

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
	 rm -f $(EXE) *.o 

