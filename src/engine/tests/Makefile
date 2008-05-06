
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

.PHONY: clean proper all tests default

.PHONY: clean proper default

default: all

all: tests clean

tests: 
	@echo "."
	@echo ". RUNNING TESTS"
	@echo "."
	@./run-all ../$(ENGINE)
	@echo "."

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
