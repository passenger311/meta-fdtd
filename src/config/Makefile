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

.PHONY: clean proper default install all info lua luageo

ifeq ($(LUAPLATFORM),) 
LUAPLATFORM=ansi
endif

INSTALLFILES=\
lua/src/lua\
config.lua


default: all

all: lua luageo

lua: 
	@echo . 
	@echo . BUILDING LUA INTERPRETER
	@echo . 
	@cd lua/src; $(MAKE) $(LUAPLATFORM) 

luageo:
	@echo . 
	@echo . BUILDING GEOMETRY PREPROCESSOR
	@echo . 
	@cd luageo; $(MAKE) all

install: build-dir
	@echo [INSTALL] $(INSTALLFILES)
	@ cp $(INSTALLFILES) ../$(BUILDDIR)


.PHONY: clean proper

clean:
	@cd lua; $(MAKE) -s $@
	@cd luageo; $(MAKE) -s $@

proper: clean
	@cd lua; $(MAKE) -s clean
	@cd luageo; $(MAKE) -s $@


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
