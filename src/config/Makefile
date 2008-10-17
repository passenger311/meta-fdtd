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

.PHONY: clean proper default install all info lua cfg

INSTALLFILES=\
src/luacfg/luacfg\
src/luacfg/config.lua

default: all

all: lua cfg

lua: 
	@echo ----- BUILDING VANILLA LUA INTERPRETER
	@cd lua/src; $(MAKE) -s ansi

cfg:
	@echo ----- BUILDING LUACFG PREPROCESSOR
	@cd src; $(MAKE) -s all

install: build-dir
	@echo [INSTALL] $(INSTALLFILES)
	@ cp $(INSTALLFILES) ../$(BUILDDIR)


.PHONY: clean proper

clean:
	@cd lua; $(MAKE) -s $@
	@cd src; $(MAKE) -s $@

proper: clean
	@cd lua; $(MAKE) -s clean
	@cd src; $(MAKE) -s $@


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
