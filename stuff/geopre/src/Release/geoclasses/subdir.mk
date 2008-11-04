################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../geoclasses/CBinaryObject.cpp \
../geoclasses/CBox.cpp \
../geoclasses/CCylinder.cpp \
../geoclasses/CLogicOrObject.cpp \
../geoclasses/CObject.cpp \
../geoclasses/CSimpleTransform.cpp 

OBJS += \
./geoclasses/CBinaryObject.o \
./geoclasses/CBox.o \
./geoclasses/CCylinder.o \
./geoclasses/CLogicOrObject.o \
./geoclasses/CObject.o \
./geoclasses/CSimpleTransform.o 


# Each subdirectory must supply rules for building sources it contributes
geoclasses/%.o: ../geoclasses/%.cpp
	@echo "[CXX] $@"
	@$(CXX)  $(CXXFLAGS_OPT) -I../../expat/xmlparse -c -o"$@" "$<"


