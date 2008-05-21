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

CPP_DEPS += \
./filehandler/FileHandler.d \
./filehandler/FileHandlerFortranIN.d \
./filehandler/FileHandlerVTK.d 


# Each subdirectory must supply rules for building sources it contributes
filehandler/%.o: ../filehandler/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../../expat/xmlparse -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


