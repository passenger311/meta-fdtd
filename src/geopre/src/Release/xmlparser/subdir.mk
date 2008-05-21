################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../xmlparser/ClusterParser.cpp \
../xmlparser/GeometryParser.cpp \
../xmlparser/GridBoxParser.cpp \
../xmlparser/OutputParser.cpp \
../xmlparser/ParameterParser.cpp \
../xmlparser/SetupFileParser.cpp \
../xmlparser/UnitParser.cpp \
../xmlparser/expatpp.cpp 

OBJS += \
./xmlparser/ClusterParser.o \
./xmlparser/GeometryParser.o \
./xmlparser/GridBoxParser.o \
./xmlparser/OutputParser.o \
./xmlparser/ParameterParser.o \
./xmlparser/SetupFileParser.o \
./xmlparser/UnitParser.o \
./xmlparser/expatpp.o 

CPP_DEPS += \
./xmlparser/ClusterParser.d \
./xmlparser/GeometryParser.d \
./xmlparser/GridBoxParser.d \
./xmlparser/OutputParser.d \
./xmlparser/ParameterParser.d \
./xmlparser/SetupFileParser.d \
./xmlparser/UnitParser.d \
./xmlparser/expatpp.d 


# Each subdirectory must supply rules for building sources it contributes
xmlparser/%.o: ../xmlparser/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../../expat/xmlparse -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


