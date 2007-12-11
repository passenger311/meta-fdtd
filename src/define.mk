#======================================================================
#
# define.mk
#

# ---------------------------------------------------------------------


ifneq ($(DBG),)
ENG_DBG=-dbg
M4DEFINE+= -DM4_DBG=1
endif

ifneq ($(MPI),)
ENG_MPI=-mpi
M4DEFINE+= -DM4_MPI=1
endif

ifneq ($(OMP),)
ENG_OMP=-omp
M4DEFINE+= -DM4_OMP=1
endif

ifneq ($(CF),)
ENG_CF=-cf
M4DEFINE+= -DM4_CF=1
endif

ifneq ($(NG),)
ENG_NG=-ng
M4DEFINE+=-DM4_NG=1
endif

ifneq ($(TE),)
ENG_TE=-te
M4DEFINE+=-DM4_TE=1
endif

ENGINE=$(ENG_CF)$(ENG_NG)$(ENG_TE)$(ENG_MPI)$(ENG_OMP)$(ENG_DBG)
BUILDPFX=meta3-build
BUILDDIR=$(BUILDPFX)$(ENGINE)

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
