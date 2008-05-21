
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

.PHONY: clean proper default install all info

ifeq ($(DBG),1)
GBUILDDIR=src/Release
else
GBUILDDIR= src/Debug
endif


BUILDFILES=\
geopre

default: all

info:
	@echo "."
	@echo ". GEOPRE"
	@echo "."

all: info geopre install


geopre:
	cd $(GBUILDDIR); $(MAKE) -s clean all
	cp $(GBUILDDIR)/geopre .

install: build-dir
	@echo [INSTALL] $(BUILDFILES)
	@cp $(BUILDFILES) ../$(BUILDDIR)


clean:
	@rm -f geopre
	@cd src/Release; make clean
	@cd src/Debug; make clean

proper: clean


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
