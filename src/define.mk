#======================================================================
#
# define.mk
#

# ---------------------------------------------------------------------


ifneq ($(DBG),)
ENG_DBG=-dbg
M4DEFINE+= -DDBG=1
endif

ifneq ($(MPI),)
ENG_MPI=-mpi
M4DEFINE+= -DMPI=1
endif

ifneq ($(OMP),)
ENG_OMP=-omp
M4DEFINE+= -DOMP=1
endif

ifneq ($(CF),)
ENG_CF=-cf
M4DEFINE+= -DCF=1
endif

ifneq ($(NG),)
ENG_NG=-ng
M4DEFINE+=-DNG=1
endif

ifneq ($(TE),)
ENG_TE=-te
M4DEFINE+=-DTE=1
endif

ENGINE=$(ENG_CF)$(ENG_NG)$(ENG_TE)$(ENG_MPI)$(ENG_OMP)
BUILDPFX=meta3-build
BUILDDIR=$(BUILDPFX)$(ENGINE)

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
