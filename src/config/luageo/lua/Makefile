
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
OBJ=\
luageo.o\
geo_core.o\
lua.o

default: all

CCFLAGS+= -I../../lua/src
LDFLAGS+= -L../../lua/src

LUAEXE= luacfg
LIBS+= -lm

all: $(LUAEXE)

$(OBJ): %.o : %.c
	@echo [CC] $@
	@$(strip $(CC) $(CCFLAGS) -c $< -o $@)

$(LUAEXE): $(OBJ)
	@echo [LD] $@
	@$(strip $(CC) $(LDFLAGS) $^ ../../lua/src/liblua.a  $(LIBS) -o $@)

clean:
	@rm -f  *.o

proper: clean

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/10/2008
#
#======================================================================
