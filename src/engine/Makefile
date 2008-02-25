
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

M4= m4 $(M4DEFINE)
M4= m4 $(M4DEFINE)

MATLIST=$(shell ./mlist.pl -s mat.def 2>/dev/null)
DIAGLIST=$(shell ./mlist.pl -s diag.def 2>/dev/null)
OUTGPLLIST=$(shell ./mlist.pl -s outgpl.def 2>/dev/null)
MATOBJ= $(patsubst %,%.o,$(MATLIST))
DIAGOBJ= $(patsubst %,%.o,$(DIAGLIST))
OUTGPLOBJ= $(patsubst %,%_outgpl.o,$(OUTGPLLIST))
OUTVTKOBJ= $(patsubst %,%_outgpl.o,$(OUTVTKLIST))
ifneq ($(HD5),)
OUTHD5OBJ= $(patsubst %,%_outgpl.o,$(OUTHD5LIST))
OUTHD5= outhd5.o
FDTD_OUTHD5= fdtd_outhd5.o
endif

# basic objects
BOBJ=\
constant.o\
strings.o\
parse.o\
utils.o\
mpiworld.o\
reglist.o\
vallist.o\
buflist.o\
outlist.o\
list.o\
grid.o

# core objects
FOBJ=\
numerics.o\
fdtd.o\
fdtd_calc.o\
pec.o\
sbc.o\
pbc.o\
pml.o\
out_calc.o\
fdtd_outset.o\
fdtd_outgpl.o\
fdtd_outvtk.o\
$(FDTD_OUTHD5)\
$(MATOBJ)\
mat.o\
$(DIAGOBJ)\
diag.o\
$(OUTGPLOBJ)\
$(OUTVTKOBJ)\
$(OUTHD5OBJ)\
outset.o\
outgpl.o\
outvtk.o\
$(OUTHD5)\
mpicomms.o

# top objects
TOBJ=\
out.o\
bound.o\
config.o\
manifest.o\
meta.o


OBJ=$(BOBJ) $(FOBJ)
M4OBJ=$(patsubst %.o,%_m4.f90,$(OBJ))

EXE = meta

ENGINEFILES=\
$(EXE)\
MANIFEST

default: all

all: info check $(EXE)

info:
	@echo ". ENGINE [$(BUILDDIR)/$(ENGINE)]"

check:
	@if [ "`cat BUILD 2>/dev/null`" != "$(BUILDDIR)" ]; then \
	echo "$(BUILDDIR)" >BUILD; \
	fi
	@if [ "`cat FLAVOUR 2>/dev/null`" != "$(ENGINE)" ]; then \
	echo "$(ENGINE)" >FLAVOUR; \
	fi

install: MANIFEST build-dir
	@echo [INSTALL] $(ENGINEFILES)
	@cp $(ENGINEFILES) ../$(ENGINEDIR)

MANIFEST: manifest.inf 
	@echo [SIGN] $@
	@$(M4) $< | tail -n 4 >$@

$(BOBJ): %.o : %.f90 BUILD
	@echo [F90] $@
	@if [ -f $*.m4 ]; then \
	$(M4) $*.m4 $< >$*_m4.f90; \
	else \
	cp $*.f90 $*_m4.f90; \
	fi
	@$(strip $(F90) $(F90FLAGS_NOOPT) -c $*_m4.f90 -o $@)

$(FOBJ) : %.o : %.f90 $(BOBJ) BUILD FLAVOUR
	@echo [F90] $@
	@if [ -f $*.m4 ]; then \
	$(M4) $*.m4 $< >$*_m4.f90; \
	else \
	cp $*.f90 $*_m4.f90; \
	fi
	@$(strip $(F90) $(F90FLAGS_OPT) -c $*_m4.f90 -o $@)

$(TOBJ): %.o : %.f90 $(BOBJ) $(FOBJ) BUILD FLAVOUR
	@echo [F90] $@
	@if [ -f $*.m4 ]; then \
	$(M4) $*.m4 $< >$*_m4.f90; \
	else \
	cp $*.f90 $*_m4.f90; \
	fi
	@$(strip $(F90) $(F90FLAGS_NOOPT) -c $*_m4.f90 -o $@)

$(EXE) : $(TOBJ) $(FOBJ) $(BOBJ) second_wall.o
	@echo [LD] $@
	@$(strip @$(F90) $(LDFLAGS) $^ -o $(EXE))

%.o : %.c
	@echo [CC] $@
	@$(strip @$(CC) $(CCFLAGS) -c $<)

clean:
	@rm -f  core* a.out *.L *.par *.log *.list *.diag *.txt *.mod *_m4.f90 *_m4.F90 *_m4.tmp BUILD FLAVOUR

proper: clean
	@rm -f $(EXE) $(PREEXE) *.o *~ *.gpl *.vtk *.set


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
