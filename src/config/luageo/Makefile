
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------

ifeq ($(ARCH),)
include ../../Make.default
else 
include ../../Make.$(ARCH)
endif 
include ../../define.mk

# ---------------------------------------------------------------------

.PHONY: clean proper default install all info sign check objects filehandler lua

# objects
OBJ=\
errorclasses.o\
GridBox.o\
OutputFile.o\
OutputList.o\
Scene.o

default: all

all: objects filehandler $(OBJ) lua 

objects:
	@cd objects; $(MAKE) -s all

filehandler: 
	@cd filehandler; $(MAKE) -s all

$(OBJ): %.o : %.cpp
	@echo [CXX] $@
	@$(strip $(CXX) $(CXXFLAGS) -c $< -o $@)

lua:
	@cd lua; $(MAKE) -s all

clean:
	@cd objects; $(MAKE) -s $@
	@cd filehandler; $(MAKE) -s $@
	@rm -f  *.o

proper: clean
	@cd objects; $(MAKE) -s $@
	@cd filehandler; $(MAKE) -s $@

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/10/2008
#
#======================================================================
