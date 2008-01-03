
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
	@echo ". CREATE FULL BUILD [$(BUILDDIR)]"
	-@$(MAKE) -s config post run
	-@$(MAKE) -s engine DIM=1 
	-@$(MAKE) -s engine DIM=2 
	-@$(MAKE) -s engine DIM=3
	-@$(MAKE) -s engine DIM=1 CF=1
	-@$(MAKE) -s engine DIM=2 CF=1
	-@$(MAKE) -s engine DIM=3 CF=1
	-@$(MAKE) -s tests
	-@$(MAKE) -s proper
	@echo ". FULL BUILD COMPLETED [$(BUILDDIR)]"

all: config post run engine 


engine:
	@cd $@ && $(MAKE) -s proper all install

config post run tests:
	@cd $@ && $(MAKE) -s all install

clean proper: 
	@cd engine; $(MAKE) -s $@; cd ..
	@cd config; $(MAKE) -s $@; cd ..
	@cd post; $(MAKE) -s $@; cd ..
	@cd run; $(MAKE) -s $@; cd ..
	@cd tests; $(MAKE) -s $@; cd ..

distclean: proper
	@rm -rf build-* build

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
