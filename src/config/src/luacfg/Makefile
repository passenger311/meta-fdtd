
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
luacfg_lmod.o\
lua.o

CXXOBJ=\
geo_core.o\
geo_scene.o\
geo_grid.o\
geo_file.o\
geo_utils.o\
geo_objects.o

default: all

CCFLAGS+= -I../../lua/src
CXXFLAGS+= -I../../lua/src
LDFLAGS+= -L../../lua/src

LUAEXE= luacfg
LIBS+= -lm

LIBLUA=../../lua/src/liblua.a
LIBGEO=../libgeo.a

all: $(LUAEXE)

luacfg_lmod.c: cfg.lua
	@../../lua/src/lua forge.lua

$(CCOBJ): %.o : %.c
	@echo [CC] $@
	@$(strip $(CC) $(CCFLAGS) -c $< -o $@)

$(CXXOBJ): %.o: %.cpp
	@echo [CXX] $@
	@$(strip $(CXX) $(CXXFLAGS) -c $< -o $@)

$(LUAEXE): $(CCOBJ) $(CXXOBJ) $(LIBLUA) $(LIBGEO)
	@echo [LD] $@
	@$(strip $(CXX) $(LDFLAGS) $^ $(LIBLUA) $(LIBGEO) $(LIBS) -o $@)

clean:
	@rm -f  *.o *_lmod.* luacfg

proper: clean

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/10/2008
#
#======================================================================
