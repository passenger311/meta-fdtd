
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------


.SUFFIXES: .f90 .F90


ifneq ($(MPI),)
MPI=-mpi
endif

ifneq ($(OMP),)
OMP=-omp
endif

ifneq ($(CT),)
CT=-ct
endif

ifneq ($(NG),)
NG=-ng
endif

ifneq ($(TE),)
TE=-te
endif

ENGINE=$(CT)$(NG)$(TE)$(MPI)$(OMP)


.PHONY: default clean proper eng pre post install

default: eng pre post run install

eng pre post:
	@cd $@ && $(MAKE) 


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
