
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

.PHONY: clean proper default install info

BUILDFILES=\
qscript.mpi\
mpiwrap.sh


.PHONY: clean proper default install

default: all

info: 
	@echo ". RUN SCRIPTS"

all: info install

install: build-dir
	@echo [INSTALL] $(BUILDFILES)
	@ cp $(BUILDFILES) ../$(BUILDDIR)

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
