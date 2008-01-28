
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
OUTVTKOBJ= $(patsubst %,%_outgpl.o,$(OUTVTKLIST))
ifneq ($(HD5),)
OUTHD5OBJ= $(patsubst %,%_outgpl.o,$(OUTHD5LIST))
OUTHD5= outhd5.o
FDTD_OUTHD5= fdtd_outhd5.o
endif

BOBJ=\
constant.o\
strings.o\
utils.o\
mpiworld.o\
reglist.o\
vallist.o\
buflist.o\
outlist.o\
list.o\
grid.o


FOBJ=\
fdtd.o\
fdtd_calc.o\
pec.o\
sbc.o\
pbc.o\
pml.o\
bound.o\
out_calc.o\
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
outgpl.o\
outvtk.o\
$(OUTHD5)\
out.o\
config.o\
mpicomms.o\
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
	@echo ". BUILD ENGINE [$(BUILDDIR)/$(ENGINE)]"

check:
	@if [ "`cat BUILD 2>/dev/null`" != "$(BUILDDIR)/$(ENGINE)" ]; then \
	echo "$(BUILDDIR)/$(ENGINE)" >BUILD; \
	else \
	touch BUILD; \
	fi

install: MANIFEST build-dir
	@echo [INSTALL] $(ENGINEFILES)
	@cp $(ENGINEFILES) ../$(ENGINEDIR)

MANIFEST: manifest.inf 
	@echo [SIGN] $@
	@$(M4) $< | tail -n 4 >$@


$(M4OBJ): %_m4.f90: %.f90 BUILD
	@$(M4) $< >$*_m4.tmp; \
	if [ -f $@ ]; then \
	if [ "`diff $*_m4.tmp $@`x" != "x" ]; then \
	mv $*_m4.tmp $@; \
	fi \
	else \
	mv $*_m4.tmp $@; \
	fi

$(BOBJ) : %.o : %_m4.f90
	@echo [F90] $@
	@$(strip $(F90) $(F90FLAGS) -c $< -o $@)

$(FOBJ) : %.o : %_m4.f90 $(BOBJ)
	@echo [F90] $@
	@$(strip $(F90) $(F90FLAGS_OPT) -c $< -o $@)

$(EXE) : $(FOBJ) $(BOBJ) second_wall.o
	@echo [LD] $@
	@$(strip @$(F90) $(LDFLAGS) $^ -o $(EXE))

%.o : %.c
	@echo [CC] $@
	@$(strip @$(CC) $(CCFLAGS) -c $<)

clean:
	@rm -f  core* a.out *.L *.par *.log *.list *.diag *.txt *.mod *_m4.f90 *_m4.F90 *_m4.tmp

proper: clean
	@rm -f $(EXE) $(PREEXE) *.o *~ *.gpl *.vtk BUILD


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
