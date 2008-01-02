
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------


include ./define.mk

.PHONY: default clean proper distclean engine config post install all

default: all

all: engine config post run install

engine config post:
	@cd $@ && $(MAKE) 


clean install: 
	@cd engine;  $(MAKE) $@; cd ..
	@cd config;  $(MAKE) $@; cd ..
	@cd post; $(MAKE) $@; cd ..
	@cd run; $(MAKE) $@; cd ..

proper: 
	@cd engine;  $(MAKE) $@; cd ..
	@cd config;  $(MAKE) $@; cd ..
	@cd post; $(MAKE) $@; cd ..
	@cd run; $(MAKE) $@; cd ..
	@rm -f *~ \#*

distclean: proper
	@rm -rf meta3* meta2*

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
