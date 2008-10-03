
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

.PHONY: clean proper default install all info libexpat

GBUILDDIR=src/Release
EXPATDIR=expat

BUILDFILES=\
geopre

default: all

info:
	@echo "."
	@echo ". GEOPRE"
	@echo "."

all: info geopre install

libexpat:
	@echo .
	@echo . EXPAT LIBRARY
	@echo .
	@cd $(EXPATDIR); $(MAKE) -s clean all

geopre: libexpat
	@echo "."
	@echo ". GEOPRE EXECUTABLE"
	@echo "."
	@cd $(GBUILDDIR); $(MAKE) MPI= -s clean all
	@cp $(GBUILDDIR)/geopre .

install: build-dir
	@echo [INSTALL] $(BUILDFILES)
	@cp $(BUILDFILES) ../$(BUILDDIR)


clean:
	@rm -f geopre
	@cd expat; $(MAKE) -s clean
	@cd src/Release; $(MAKE) -s clean

proper: clean


# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
