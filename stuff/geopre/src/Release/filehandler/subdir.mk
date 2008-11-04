################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../filehandler/FileHandler.cpp \
../filehandler/FileHandlerFortranIN.cpp \
../filehandler/FileHandlerVTK.cpp 

OBJS += \
./filehandler/FileHandler.o \
./filehandler/FileHandlerFortranIN.o \
./filehandler/FileHandlerVTK.o 

# Each subdirectory must supply rules for building sources it contributes
filehandler/%.o: ../filehandler/%.cpp
	@echo "[CXX] $@"
	@$(CXX) $(CXXFLAGS_NOOPT) -I../../expat/xmlparse -c -o"$@" "$<"


