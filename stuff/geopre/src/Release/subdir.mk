################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../ArgumentReader.cpp \
../GridBox.cpp \
../OutputFile.cpp \
../OutputList.cpp \
../Scene.cpp \
../errorclasses.cpp \
../main.cpp 

OBJS += \
./ArgumentReader.o \
./GridBox.o \
./OutputFile.o \
./OutputList.o \
./Scene.o \
./errorclasses.o \
./main.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo "[CXX] $@"
	@$(CXX) $(CXXFLAGS_OPT) -I../../expat/xmlparse -c -o"$@" "$<"


