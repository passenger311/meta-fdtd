
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
TEST.LOG


.PHONY: clean proper default install

default: all

all: tests

tests: 
	@echo ". RUNNING TESTS FOR [$(BUILDDIR)]"
	@echo "."
	@./run-all $(BUILDDIR)

install:
	@echo [INSTALL] $(BUILDFILES)
	@cp $(BUILDFILES) ../$(BUILDDIR)

.PHONY: clean proper

clean:
	@rm -f */*.gpl */log TEST.LOG

proper: clean


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
