
ifeq ($(ARCH),)
include ../../../Make.default
else 
include ../../../Make.$(ARCH)
endif 


-include ../makefile.init

RM := rm -rf

-include sources.mk
-include subdir.mk
-include xmlparser/subdir.mk
-include geoclasses/subdir.mk
-include filehandler/subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: geopre

geopre: $(OBJS) $(USER_OBJS)
	@echo "[LD] $@"
	@$(CXX) -L../../expat/xmlparse -o "geopre" $(OBJS) $(USER_OBJS) $(LIBS)

# Other Targets
clean:
	@$(RM) $(OBJS)$(C++_DEPS)$(EXECUTABLES)$(CC_DEPS)$(C_DEPS)$(CPP_DEPS)$(CXX_DEPS)$(C_UPPER_DEPS) geopre

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
