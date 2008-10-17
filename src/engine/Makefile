
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

.PHONY: clean default install all info sign check 

M4= m4 $(M4DEFINE)
M4= m4 $(M4DEFINE)


SRCLIST=$(shell ./scripts/mlist.lua -s src.def 2>/dev/null)
MATLIST=$(shell ./scripts/mlist.lua -s mat.def 2>/dev/null)
DIAGLIST=$(shell ./scripts/mlist.lua -s diag.def 2>/dev/null)
OUTGPLLIST=$(shell ./scripts/mlist.lua -s outgpl.def 2>/dev/null)
SRCOBJ= $(patsubst %,%.o,$(SRCLIST))
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
BASIC_OBJ=\
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

# third party objects

THIRD_OBJ=\
fftpack5.o

# core objects
CORE_OBJ=\
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
$(SRCOBJ)\
src.o\
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
TOP_OBJ=\
out.o\
bound.o\
config.o\
manifest.o\
timer.o\
meta.o

EXE = meta

default: all

all: check $(ENGINE)/$(EXE) sign

info:
	@echo "----- ENGINE [$(BUILDDIR)/$(ENGINE)]"

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
	@echo [INSTALL] engine-tests
	@rm -rf ../$(BUILDDIR)/engine-tests
	@cp -a engine-tests ../$(BUILDDIR)/engine-tests

$(ENGINE):
	@mkdir $(ENGINE)

sign: $(ENGINE)
	@echo [SIGN]  $(ENGINE)/MANIFEST
	@$(M4) manifest.m4 manifest.inf | sed '/^ *$$/d' >$(ENGINE)/MANIFEST

fftpack5/libfftpack.a:
	@ cd fftpack5; $(MAKE) -s all

$(BASIC_OBJ): %.o : %.f90
	@echo [F90] $@
	@if [ -f $*.m4 ]; then \
	$(M4) $*.m4 $< >$*_m4.f90; \
	else \
	cp $*.f90 $*_m4.f90; \
	fi
	@$(strip $(F90) $(F90FLAGS) -c $*_m4.f90 -o $@)

$(THIRD_OBJ) : %.o : %.f90
	@echo [F90] $@
	@$(strip $(F90) $(F90FLAGS_HOPT) -c $*.f90 -o $@)

$(CORE_OBJ) : %.o : %.f90 $(BASIC_OBJ) BUILD FLAVOUR
	@echo [F90] $@
	@if [ -f $*.m4 ]; then \
	$(M4) $*.m4 $< >$*_m4.f90; \
	else \
	cp $*.f90 $*_m4.f90; \
	fi
	@$(strip $(F90) $(F90FLAGS_HOPT) -c $*_m4.f90 -o $@)

$(TOP_OBJ): %.o : %.f90 $(BASIC_OBJ) $(CORE_OBJ) BUILD FLAVOUR
	@echo [F90] $@
	@if [ -f $*.m4 ]; then \
	$(M4) $*.m4 $< >$*_m4.f90; \
	else \
	cp $*.f90 $*_m4.f90; \
	fi
	@$(strip $(F90) $(F90FLAGS) -c $*_m4.f90 -o $@)

$(ENGINE)/$(EXE) : $(TOP_OBJ) $(CORE_OBJ) $(THIRD_OBJ) $(BASIC_OBJ) second_wall.o $(ENGINE)
	@echo [LD] $@
	@$(strip $(F90) $(LDFLAGS) $(TOP_OBJ) $(CORE_OBJ) $(THIRD_OBJ) $(BASIC_OBJ) second_wall.o -o $(ENGINE)/$(EXE))

%.o : %.c
	@echo [CC] $@
	@$(strip @$(CC) $(CCFLAGS) -c $<)

clean:
	@rm -f  core* a.out *.L *.par *.log *.list *.diag *.txt *_m4.f90 *_m4.F90 *_m4.tmp BUILD VERSION FLAVOUR MANIFEST
	@cd engine-tests; $(MAKE) -s clean
	@rm -f $(PREEXE) *.o *~ *.gpl *.vtk *.set *.mod
	@rm -f fftpack5/*.o fftpack5/*.a
	@rm -rf meta1* meta2* meta3* 

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
