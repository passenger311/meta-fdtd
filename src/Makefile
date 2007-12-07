
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------


.SUFFIXES: .f90 .F90


ifneq ($(MPI),)
MPI=-mpi
M4DEFINE+= -DM4_MPI=1
endif

ifneq ($(OMP),)
OMP=-omp
M4DEFINE+= -DM4_OMP=1
endif

ifneq ($(CT),)
CT=-ct
M4DEFINE+= -DM4_CT=1
endif

ifneq ($(NG),)
NG=-ng
M4DEFINE+=-DM4_NG=1
endif

ifneq ($(TE),)
TE=-te
M4DEFINE+=-DM4_TE=1
endif

ENGINE=$(CT)$(NG)$(TE)$(MPI)$(OMP)


.PHONY: default clean proper eng pre post install

default: eng pre post run install

eng pre post:
	@cd $@ && $(MAKE) $(M4DEFINE) $(ENGINE)


clean proper install: 
	@cd eng;  $(MAKE) $@; cd ..
	@cd pre;  $(MAKE) $@; cd ..
	@cd post; $(MAKE) $@; cd ..
	@cd run; $(MAKE) $@; cd ..


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
