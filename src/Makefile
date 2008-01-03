
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------


include ./define.mk

.PHONY: default clean proper distclean engine config post run all build

build: 
	-$(MAKE) config post run
	-$(MAKE) proper engine DIM=1 
	-$(MAKE) proper engine DIM=2 
	-$(MAKE) proper engine DIM=3
	-$(MAKE) proper engine DIM=1 CF=1
	-$(MAKE) proper engine DIM=2 CF=1
	-$(MAKE) proper engine DIM=3 CF=1

all: config post run engine

engine config post run:
	@cd $@ && $(MAKE) all install

clean proper: 
	@cd engine; $(MAKE) $@; cd ..
	@cd config; $(MAKE) $@; cd ..
	@cd post; $(MAKE) $@; cd ..
	@cd run; $(MAKE) $@; cd ..

distclean: proper
	@rm -rf build-* build

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
