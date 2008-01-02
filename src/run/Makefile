
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

.PHONY: clean proper default install

INSTALLFILES=\
qscript.mpi\
mpiwrap.sh\
config.0.in

.PHONY: clean proper default install

default: all

all: install

install: 
	@echo [INSTALL] $(INSTALLFILES)
	@ [ -d ../$(BUILDDIR) ] || mkdir ../$(BUILDDIR)
	@ cp $(INSTALLFILES) ../$(BUILDDIR)

.PHONY: clean proper

clean:

proper: clean
	@rm -f *~


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
