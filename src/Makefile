
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------

ifeq ($(ARCH),)
ARCH=gnu
endif
include ./define.mk

.PHONY: default clean proper distclean engine config post run all build tests

build: 
	@echo "+==== CREATE FULL BUILD: $(BUILDDIR)"
	-@$(MAKE) config post run
	-@$(MAKE) engine DIM=1 
	-@$(MAKE) engine DIM=2 
	-@$(MAKE) engine DIM=3
	-@$(MAKE) engine DIM=1 CF=1
	-@$(MAKE) engine DIM=2 CF=1
	-@$(MAKE) engine DIM=3 CF=1
	-@$(MAKE) tests
	@echo "+==== FULL BUILD COMPLETED: $(BUILDDIR)"

all: config post run engine 


engine:
	@cd $@ && $(MAKE) proper all install

config post run tests:
	@cd $@ && $(MAKE) all install

clean proper: 
	@cd engine; $(MAKE) $@; cd ..
	@cd config; $(MAKE) $@; cd ..
	@cd post; $(MAKE) $@; cd ..
	@cd run; $(MAKE) $@; cd ..
	@cd tests; $(MAKE) $@; cd ..

distclean: proper
	@rm -rf build-* build

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
