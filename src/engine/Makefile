
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

.PHONY: clean proper default install all info sign check

M4= m4 $(M4DEFINE)
M4= m4 $(M4DEFINE)


MATLIST=$(shell ./scripts/mlist.pl -s mat.def 2>/dev/null)
DIAGLIST=$(shell ./scripts/mlist.pl -s diag.def 2>/dev/null)
OUTGPLLIST=$(shell ./scripts/mlist.pl -s outgpl.def 2>/dev/null)
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
signal.o\
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
timer.o\
meta.o


OBJ=$(BOBJ) $(FOBJ)
M4OBJ=$(patsubst %.o,%_m4.f90,$(OBJ))

EXE = meta

default: all

all: info check $(ENGINE)/$(EXE) sign

info:
	@echo ". ENGINE [$(BUILDDIR)/$(ENGINE)]"

check: 
	@if [ "`cat BUILD 2>/dev/null`" != "$(BUILDDIR)" ]; then \
	echo "$(BUILDDIR)" >BUILD; \
	fi
	@if [ "`cat FLAVOUR 2>/dev/null`" != "$(ENGINE)" ]; then \
	echo "$(ENGINE)" >FLAVOUR; \
	fi
	@echo "`date +%y%m%d`" >VERSION

install: build-dir
	@echo [INSTALL] $(ENGINE)
	@cp -a $(ENGINE) ../$(BUILDDIR)


install-tests: build-dir
	@echo [INSTALL] testsuite
	@rm -rf ../$(BUILDDIR)/engine-tests
	@cp -a tests ../$(BUILDDIR)/engine-tests


$(ENGINE):
	@mkdir $(ENGINE)

sign: $(ENGINE)
	@echo [SIGN]  $(ENGINE)/MANIFEST
	@$(M4) manifest.m4 manifest.inf | sed '/^ *$$/d' >$(ENGINE)/MANIFEST

$(BOBJ): %.o : %.f90
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

$(ENGINE)/$(EXE) : $(TOBJ) $(FOBJ) $(BOBJ) second_wall.o $(ENGINE)
	@echo [LD] $@
	@$(strip $(F90) $(LDFLAGS) $(TOBJ) $(FOBJ) $(BOBJ) second_wall.o -o $(ENGINE)/$(EXE))

%.o : %.c
	@echo [CC] $@
	@$(strip @$(CC) $(CCFLAGS) -c $<)

clean:
	@rm -f  core* a.out *.L *.par *.log *.list *.diag *.txt *.mod *_m4.f90 *_m4.F90 *_m4.tmp BUILD VERSION FLAVOUR MANIFEST

proper: clean
	@cd testsuite; $(MAKE) -s proper
	@rm -f $(PREEXE) *.o *~ *.gpl *.vtk *.set
	@rm -rf meta1* meta2* meta3* 

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
