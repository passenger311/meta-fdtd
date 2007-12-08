
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------

ifeq ($(ARCH),)
include ../Make.default
else 
include ../Make.$(ARCH)
endif 
include ../define.mk

# ---------------------------------------------------------------------

.PHONY: clean proper default install

M4= m4 $(M4DEFINE) meta3.m4

MATLIST=$(shell mlist.pl -s mat.def 2>/dev/null)
DIAGLIST=$(shell mlist.pl -s diag.def 2>/dev/null)
OUTGPLLIST=$(shell mlist.pl -s outgpl.def 2>/dev/null)
MATOBJ= $(patsubst %,%.o,$(MATLIST))
DIAGOBJ= $(patsubst %,%.o,$(DIAGLIST))
OUTGPLOBJ= $(patsubst %,%_outgpl.o,$(OUTGPLLIST))

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

INSTALLFILES=\
$(EXE)

default: $(EXE)

install: 
	@echo [INSTALL] $(INSTALLFILES)
	@ [ -d ../$(BUILDDIR) ] || mkdir ../$(BUILDDIR)
	@ cp $(INSTALLFILES) ../$(BUILDDIR)

%.o : %.f90
	@echo [F90] $*
	@$(M4) $< >$*_m4.f90 
	@$(F90) $(F90FLAGS) $(INCLUDE) -c $*_m4.f90 -o $*.o

%.o : %.F90
	@echo [F90] $*
	@$(M4) $< >$*_m4.F90 
	@$(F90) $(F90FLAGS) $(INCLUDE) -c $*_m4.F90 -o $*.o

%.o : %.c
	@echo [CC ] $*
	@$(CC) $(CCFLAGS) $(INCLUDE) -c $<


$(EXE) : $(FOBJ)
	$(F90) $(LDFLAGS) $(FOBJ) -o $(EXE)

clean:
	 rm -f  core* a.out *.L *.par *.log *.list *.diag *.txt *.mod *_m4.f90 *_m4.F90

proper: clean
	 rm -f $(EXE) $(PREEXE) *.o *.in *.out *~


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
