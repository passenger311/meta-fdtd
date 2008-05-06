
#======================================================================
#
# Makefile
#
# ---------------------------------------------------------------------

ifeq ($(ARCH),)
ARCH=gnu
endif
include ./define.mk

.PHONY: default clean proper distclean fftpack5 engine config post run  build tests

build: 
	-@$(MAKE) -s proper
	@echo 
	@echo "+ CREATE BUILD [$(BUILDDIR)]"
	-@$(MAKE) -s config
	-@$(MAKE) -s post
	-@$(MAKE) -s run
	-@$(MAKE) -s engine DIM=1 
	-@$(MAKE) -s engine DIM=2 
	-@$(MAKE) -s engine DIM=3
	-@$(MAKE) -s engine DIM=1 CF=1
	-@$(MAKE) -s engine DIM=2 CF=1
	-@$(MAKE) -s engine DIM=3 CF=1
	-@$(MAKE) -s engine-tests
	@cd $(BUILDDIR)/engine-tests; $(MAKE) -s info all; cd .. 
	@echo "."
	@echo "+ BUILD COMPLETED [$(BUILDDIR)]"
	@echo

default: all

all: fftpack5 config post run engine 

fftpack5: 
ifneq ($(FFT),)
	@cd $@ && $(MAKE) -s all
endif 

engine-tests:
	@echo "."
	@echo ". ENGINE TESTS"
	@echo "."
	@cd engine && $(MAKE) -s install-tests

engine:
	@cd $@ && $(MAKE) -s info all install

config post run:
	@cd $@ && $(MAKE) -s info all install

clean proper: 
	@cd engine; $(MAKE) -s $@; cd ..
	@cd config; $(MAKE) -s $@; cd ..
	@cd post; $(MAKE) -s $@; cd ..
	@cd run; $(MAKE) -s $@; cd ..
ifneq ($(FFT),)
	@cd fftpack5; $(MAKE) -s $@; cd ..
endif 

distclean: proper
	@rm -rf build-* build

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
