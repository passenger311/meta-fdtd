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

CPP_DEPS += \
./ArgumentReader.d \
./GridBox.d \
./OutputFile.d \
./OutputList.d \
./Scene.d \
./errorclasses.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../../expat/xmlparse -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


