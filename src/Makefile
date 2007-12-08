
#======================================================================
#
# Makefile
#

# ---------------------------------------------------------------------


include ./define.mk

.PHONY: default clean proper eng pre post install

default: eng pre post run install

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
	@rm *~ \#*
	@rm -rf meta3-build*

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
