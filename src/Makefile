.SUFFIXES: .f90
#
# Compilation
# -----------
#
F90 = mpif90 
CC= mpicc

F90FLAGS= -g -fpp -O3 -tpp7 -xW -r8
LDFLAGS= -g -O3 -xW -r8
CCFLAGS=


SRC= strings.f90\
constant.f90\
mpistart.f90\
grid-r.f90\
bound-r.f90\
fdtd3d-r.f90\
upml-r.f90\
config.f90\
output-r.f90\
tfsf-r.f90\
pre-obj-1.f90\
second_wall.c\
max3d.f90

FOBJ= strings.o\
constant.o\
mpistart.o\
grid-r.o\
bound-r.o\
fdtd3d-r.o\
upml-r.o\
config.o\
output-r.o\
tfsf-r.o\
second_wall.o\
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

%.o : %.c
	 $(CC) $(CCFLAGS) $(INCLUDE) -c $<


$(EXE) : $(FOBJ)
	$(F90) $(LDFLAGS) $(FOBJ) -o $(EXE) second_wall.o

$(PREEXE) : $(PREOBJ)
	$(F90) $(LDFLAGS) $(PREOBJ) -o $(PREEXE)

#
clean:
	 rm -f *.o???? core* a.out *.L *.par *.log *.list *.diag *.txt

proper:
	 rm -f $(EXE) *.mod *.o core* a.out *.o???? *.L *.par *.log *.list *.diag *.txt *~
#
cleanobj: clean
	 rm -f $(FOBJ) $(EXE) $(PREEXE)
