
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

.PHONY: clean default install all info sign check 

M4= m4 $(M4DEFINE)
M4= m4 $(M4DEFINE)



# basic objects
BASIC_OBJ=

# core objects
CORE_OBJ=\
fdtd.o\

# top objects
TOP_OBJ=

LIBSO= cfdtd.so
LIBA= cfdtd.a

default: all 

all: $(LIBSO) $(LIBA)


$(BASIC_OBJ): %.o : %.c
	@echo [CC] $@
	@if [ -f $*.m4 ]; then \
	$(M4) $*.m4 $< >$*_m4.c; \
	else \
	cp $*.c $*_m4.c; \
	fi
	@$(strip $(CC) $(CCFLAGS) -c $*_m4.c -o $@)

$(CORE_OBJ) : %.o : %.c $(BASIC_OBJ)
	@echo [F90] $@
	@if [ -f $*.m4 ]; then \
	$(M4) $*.m4 $< >$*_m4.c; \
	else \
	cp $*.c $*_m4.c; \
	fi
	@$(strip $(F90) $(F90FLAGS_HOPT) -c $*_m4.c -o $@)

$(TOP_OBJ): %.o : %.c $(BASIC_OBJ) $(CORE_OBJ)
	@echo [F90] $@
	@if [ -f $*.m4 ]; then \
	$(M4) $*.m4 $< >$*_m4.c; \
	else \
	cp $*.c $*_m4.c; \
	fi
	@$(strip $(F90) $(F90FLAGS) -c $*_m4.c -o $@)

clean:
	@rm -f  core* a.out *.L *.par *.log *.list *.diag *.txt *_m4.c *_m4.c *_m4.tmp
	@rm 

# ---------------------------------------------------------------------

#
# Authors:  J.Hamm 
# Modified: 4/12/2007
#
#======================================================================
