
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
Grid.o\
objects/CObject.o\
objects/CBinaryObject.o\
objects/CBox.o\
objects/CSphere.o\
objects/CEllipsoid.o\
objects/CCylinder.o\
objects/CConvexPrism.o\
objects/CSimplePrism.o\
objects/CBezierPrism.o\
objects/CLogicOrObject.o\
objects/CLogicAndObject.o\
objects/CLogicAndNotObject.o\
objects/CLogicXOrObject.o\
objects/CLogicOrContainer.o\
objects/CSimpleRotationZ.o\
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
	@rm -f  *.o objects/*.o filehandler/*.o *.ti objects/*.ti filehandler/*.ti  $(LIB)
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
