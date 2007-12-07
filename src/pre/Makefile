
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------


.SUFFIXES: .f90 .F90


# Compilation
# -----------

ifeq ($(ARCH),)
include ../Make.default
else 
include ../Make.$(ARCH)
endif 

ifneq ($(MPI),)
MPI=-mpi
M4DEFINE+= -DM4_MPI=1
endif

ifneq ($(OMP),)
OMP=-omp
M4DEFINE+= -DM4_OMP=1
endif

ifneq ($(CT),)
CT=ct
M4DEFINE+= -DM4_CT=1
endif

ifneq ($(NG),)
NG=ng
M4DEFINE+=-DM4_NG=1
endif

ifneq ($(TE),)
TE=te
M4DEFINE+=-DM4_TE=1
endif


ENGINE=$(CT)$(NG)$(TE)$(MPI)$(OMP)

INSTALLFILES=\
material.pm


.PHONY: clean proper default install

default: 

install: 
	@echo [INSTALL] $(INSTALLFILES)
	@ [ -d ../build$(ENGINE) ] || mkdir ../build$(ENGINE)
	@ cp $(INSTALLFILES) ../build$(ENGINE)

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

.PHONY: clean proper

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
