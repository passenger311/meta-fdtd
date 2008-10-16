
#======================================================================
#
# Makefile
#
# ---------------------------------------------------------------------

ifeq ($(ARCH),)
ARCH=gcc
endif
include ./define.mk

.PHONY: default clean distclean fftpack5 engine config geopre postps  build tests

build: 
	-@$(MAKE) -s proper
	@echo 
	@echo "+ CREATE BUILD [$(BUILDDIR)]"
	-@$(MAKE) -s config
	-@$(MAKE) -s postps
	-@$(MAKE) -s engine DIM=1 
	-@$(MAKE) -s engine DIM=2 
	-@$(MAKE) -s engine DIM=3
	-@$(MAKE) -s engine DIM=1 CF=1
	-@$(MAKE) -s engine DIM=2 CF=1
	-@$(MAKE) -s engine DIM=3 CF=1
	@cd $(BUILDDIR)/engine-tests; $(MAKE) -s info all; cd .. 
	@echo "."
	@echo "+ BUILD COMPLETED [$(BUILDDIR)]"
	@echo

default: all

all: config postps engine

engine-tests:
	@echo "."
	@echo ". ENGINE TESTS"
	@echo "."
	@cd engine && $(MAKE) -s install-tests

engine:
	@cd $@ && $(MAKE) -s info all install

config postps:
	@cd $@ && $(MAKE) -s MPI= info all install


clean: 
	@cd engine; $(MAKE) -s $@; cd ..
	@cd config; $(MAKE) -s $@; cd ..
	@cd postps; $(MAKE) -s $@; cd ..

distclean: clean
	@rm -rf build-* build

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
