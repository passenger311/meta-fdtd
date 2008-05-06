#======================================================================
#
# define.mk
#
# ---------------------------------------------------------------------

ENG_ARCH=-$(ARCH)
ifneq ($(M4_ARCH),)
M4DEFINE+= -DM4_ARCH=$(M4_ARCH)
endif

ifeq ($(DIM),2)
ENG_DIM=2
M4DEFINE+= -DM4_DIM=2
else
ifeq ($(DIM),1)
ENG_DIM=1
M4DEFINE+= -DM4_DIM=1
else
ENG_DIM=3
M4DEFINE+= -DM4_DIM=3
endif
endif

ifneq ($(DBG),)
ENG_DBG=-dbg
M4DEFINE+= -DM4_DBG=1
endif

ifneq ($(HD5),)
ENG_DBG=-hd5
M4DEFINE+= -DM4_HD5=1
endif

ifneq ($(MPI),)
ENG_MPI=-mpi
M4DEFINE+= -DM4_MPI=1
endif

ifneq ($(MPELOG),)
ifneq ($(MPI),)
ENG_MPELOG=-mpelog
M4DEFINE+= -DM4_MPELOG=1
endif
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

ifneq ($(WMU),)
ENG_WMU=-wmu
M4DEFINE+=-DM4_WMU=1
endif

ifneq ($(TM),)
ifeq ($(DIM),2)
ENG_TE=-tm
M4DEFINE+=-DM4_TMLBL=1
M4DEFINE+=-DM4_TM=1
endif
endif

ifneq ($(TE),)
ifeq  ($(DIM),2)
ifeq ($(TM),)
ENG_TE=-te
M4DEFINE+=-DM4_TELBL=1
M4DEFINE+=-DM4_TE=1
endif
endif
endif

ifeq ($(TE)$(TM),)
M4DEFINE+=-DM4_TE=1
M4DEFINE+=-DM4_TM=1
endif

ifneq ($(TEPS),)
ENG_TE=-teps
M4DEFINE+=-DM4_TEPS=1
endif


ENGINE=meta$(ENG_DIM)$(ENG_CF)$(ENG_TM)$(ENG_TE)$(ENG_WMU)$(ENG_TEPS)$(ENG_NG)
BUILDDIR=build-$(ARCH)$(ENG_MPI)$(ENG_MPELOG)$(ENG_OMP)$(ENG_DBG)
ENGINEDIR=$(BUILDDIR)/$(ENGINE)

.PHONY: build-dir

default: all

build-dir:
	@if [ ! -d ../$(BUILDDIR) ]; then \
	echo "[CREATE] $(BUILDDIR)"; \
	mkdir ../$(BUILDDIR); \
	fi



# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
