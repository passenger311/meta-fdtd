
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

.PHONY: clean proper all tests default install

BUILDFILES=\
qscript.mpi\
mpiwrap.sh


.PHONY: clean proper default install

default: all

all: tests proper

tests: 
	@echo "+---- RUNNING TESTS FOR: $(BUILDDIR)"
	@./run-all

install:

.PHONY: clean proper

clean:
	@echo "[CLEANUP]"
	@rm -f */*.gpl

proper: clean


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
