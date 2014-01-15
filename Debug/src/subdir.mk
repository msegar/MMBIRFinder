################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/ConfigDB.cpp \
../src/Fragment.cpp \
../src/alignment.cpp \
../src/birAligner.cpp \
../src/birConsolidate.cpp \
../src/birFinder.cpp \
../src/main.cpp \
../src/printOptions.cpp \
../src/programExecution.cpp \
../src/templateFinder.cpp 

OBJS += \
./src/ConfigDB.o \
./src/Fragment.o \
./src/alignment.o \
./src/birAligner.o \
./src/birConsolidate.o \
./src/birFinder.o \
./src/main.o \
./src/printOptions.o \
./src/programExecution.o \
./src/templateFinder.o 

CPP_DEPS += \
./src/ConfigDB.d \
./src/Fragment.d \
./src/alignment.d \
./src/birAligner.d \
./src/birConsolidate.d \
./src/birFinder.d \
./src/main.d \
./src/printOptions.d \
./src/programExecution.d \
./src/templateFinder.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -Wno-sign-compare -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


