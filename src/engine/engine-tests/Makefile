
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------


#ifeq ($(ARCH),)
#include ../../Make.default
#else 
#include ../../Make.$(ARCH)
#endif 
include ../../define.mk

# ---------------------------------------------------------------------

.PHONY: clean proper all tests info default

.PHONY: clean proper default

default: all

all: tests

info:
	@echo "."
	@echo ". RUNNING ENGINE TESTS"
	@echo "."

tests: 
	@./run-all ../$(ENGINE)

.PHONY: clean proper

clean:
	@rm -f */*.gpl */log */log.err */*_sh

proper: clean
	@rm -f TESTED

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 6/5/2008
#
#======================================================================
