
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

.PHONY: clean proper default install all info

M4= m4 $(M4DEFINE) meta.m4

MATLIST=$(shell ./mlist.pl -s mat.def 2>/dev/null)
DIAGLIST=$(shell ./mlist.pl -s diag.def 2>/dev/null)
OUTGPLLIST=$(shell ./mlist.pl -s outgpl.def 2>/dev/null)
MATOBJ= $(patsubst %,%.o,$(MATLIST))
DIAGOBJ= $(patsubst %,%.o,$(DIAGLIST))
OUTGPLOBJ= $(patsubst %,%_outgpl.o,$(OUTGPLLIST))

BOBJ=\
strings.o\
constant.o\
mpiworld.o\
reglist.o\
buflist.o\
outlist.o\
list.o\
grid.o

FOBJ=\
fdtd.o\
pec.o\
pml.o\
bound.o\
fdtd_outgpl.o\
$(MATOBJ)\
mat.o\
$(DIAGOBJ)\
diag.o\
$(OUTGPLOBJ)\
outgpl.o\
out.o\
config.o\
mpicomms.o\
meta.o

EXE = meta

ENGINEFILES=\
$(EXE)

default: all

all: info $(EXE)

info:
	@echo "+---- BUILDING ENGINE: $(ENGINE)"

install: build-dir
	@echo [INSTALL] $(ENGINEFILES)
	@cp $(ENGINEFILES) ../$(ENGINEDIR)

$(BOBJ) : %.o : %.f90
	@echo [F90] $@
	@$(M4) $< >$*_m4.f90 
	@$(strip $(F90) $(F90FLAGS) -c $*_m4.f90 -o $*.o)

$(FOBJ) : %.o : %.f90 $(BOBJ)
	@echo [F90] $@
	@$(M4) $< >$*_m4.f90 
	@$(strip $(F90) $(F90FLAGS_OPT) -c $*_m4.f90 -o $*.o)


$(EXE) : $(FOBJ) $(BOBJ)
	@echo [LD] $@
	@$(strip @$(F90) $(LDFLAGS) $(FOBJ) $(BOBJ) -o $(EXE))

clean:
	@rm -f  core* a.out *.L *.par *.log *.list *.diag *.txt *.mod *_m4.f90 *_m4.F90

proper: clean
	@rm -f $(EXE) $(PREEXE) *.o *~ *.gpl


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================