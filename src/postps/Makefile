
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

BUILDFILES=\
stitch-outgpl

default: all

info:
	@echo "----- POSTPS"

all: info install
 
install: build-dir
	@echo [INSTALL] $(BUILDFILES)
	#@ cp $(BUILDFILES) ../$(BUILDDIR)

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

clean:
	@rm -f  core* a.out *.L *.par *.log *.list *.diag *.txt *.mod *_m4.f90 *_m4.F90

proper: clean
	@rm -f $(EXE) $(PREEXE) *.o *.in *.out *~


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
