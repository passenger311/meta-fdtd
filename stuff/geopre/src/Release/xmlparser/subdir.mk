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

# Each subdirectory must supply rules for building sources it contributes
xmlparser/%.o: ../xmlparser/%.cpp
	@echo "[CXX] $@"
	@$(CXX) $(CXXFLAGS_NOOPT) -I../../expat/xmlparse -c -o"$@" "$<"


