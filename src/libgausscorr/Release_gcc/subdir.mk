################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../libgausscorr.cpp \
../libgaussio.cpp 

OBJS += \
./libgausscorr.o \
./libgaussio.o 

CPP_DEPS += \
./libgausscorr.d \
./libgaussio.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I../include -I../../../src -O3 -U_FORTIFY_SOURCE -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


