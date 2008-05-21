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

CPP_DEPS += \
./geoclasses/CBinaryObject.d \
./geoclasses/CBox.d \
./geoclasses/CCylinder.d \
./geoclasses/CLogicOrObject.d \
./geoclasses/CObject.d \
./geoclasses/CSimpleTransform.d 


# Each subdirectory must supply rules for building sources it contributes
geoclasses/%.o: ../geoclasses/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../../expat/xmlparse -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


