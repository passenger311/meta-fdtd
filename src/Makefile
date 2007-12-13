
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------


include ./define.mk

.PHONY: default clean proper distclean eng pre post install all

default: all

all: eng pre post run install

eng pre post:
	@cd $@ && $(MAKE) 


clean install: 
	@cd eng;  $(MAKE) $@; cd ..
	@cd pre;  $(MAKE) $@; cd ..
	@cd post; $(MAKE) $@; cd ..
	@cd run; $(MAKE) $@; cd ..

proper: 
	@cd eng;  $(MAKE) $@; cd ..
	@cd pre;  $(MAKE) $@; cd ..
	@cd post; $(MAKE) $@; cd ..
	@cd run; $(MAKE) $@; cd ..
	@rm -f *~ \#*

distclean: proper
	@rm -rf meta3*

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
