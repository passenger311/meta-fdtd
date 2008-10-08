
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

.PHONY: clean proper default install all info sign check objects filehandler luacfg

# objects
OBJ=\
errorclasses.o\
GridBox.o\
OutputFile.o\
OutputList.o\
Scene.o\
objects/CBinaryObject.o\
objects/CBox.o\
objects/CCylinder.o\
objects/CLogicOrObject.o\
objects/CObject.o\
objects/CSimpleTransform.o\
filehandler/FileHandler.o\
filehandler/FileHandlerFortranIN.o\
filehandler/FileHandlerVTK.o

LIB=libgeo.a

default: all

all: $(LIB) luacfg

$(OBJ): %.o : %.cpp
	@echo [CXX] $@
	@$(strip $(CXX) $(CXXFLAGS) -c $< -o $@)

$(LIB): $(OBJ)
	@echo [AR] $@
	@$(strip $(AR) -cr $@ $^)

luacfg:
	@cd luacfg; $(MAKE) -s all

clean:
	@rm -f  *.o objects/*.o filehandler/*.o $(LIB)
	@cd luacfg; $(MAKE) -s clean


install:
	@cd luacfg; $(MAKE) -s install

proper: clean

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/10/2008
#
#======================================================================
