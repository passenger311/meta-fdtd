
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
FileHandler.o\
FileHandlerFortranIN.o\
FileHandlerVTK.o

default: all

all: $(OBJ)

$(OBJ): %.o : %.cpp
	@echo [CXX] $@
	@$(strip $(CXX) $(CXXFLAGS) -c $< -o $@)

clean:
	@rm -f  *.o

proper: clean

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/10/2008
#
#======================================================================
