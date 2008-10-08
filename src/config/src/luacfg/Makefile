
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------

ifeq ($(ARCH),)
include ../../../Make.default
else 
include ../../../Make.$(ARCH)
endif 
include ../../../define.mk

# ---------------------------------------------------------------------

.PHONY: clean proper default install all info sign check 

# objects
CCOBJ=\
luageo.o\
cfg_lmod.o\
lua.o

CXXOBJ=\
geo_core.o

default: all

CCFLAGS+= -I../../lua/src
CXXFLAGS+= -I../../lua/src
LDFLAGS+= -L../../lua/src

LUAEXE= luacfg
LIBS+= -lm

all: $(LUAEXE)

cfg_lmod.c: cfg.lua
	@../../lua/src/lua forge.lua

$(CCOBJ): %.o : %.c
	@echo [CC] $@
	@$(strip $(CC) $(CCFLAGS) -c $< -o $@)

$(CXXOBJ): %.o: %.cpp
	@echo [CXX] $@
	@$(strip $(CXX) $(CXXFLAGS) -c $< -o $@)

$(LUAEXE): $(CCOBJ) $(CXXOBJ)
	@echo [LD] $@
	@$(strip $(CXX) $(LDFLAGS) $^ ../../lua/src/liblua.a  $(LIBS) -o $@)

clean:
	@rm -f  *.o cfg_lmod.* luacfg

proper: clean

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/10/2008
#
#======================================================================
