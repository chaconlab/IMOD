################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../MarchingCubes.cpp \
../ffourier.cpp \
../fkernel.cpp \
../fmask.cpp \
../fmat.cpp \
../fmod.cpp \
../fpad.cpp \
../frotate.cpp \
../fsurface.cpp \
../ftransform.cpp \
../fwr.cpp \
../vlkernel.cpp \
../vlmrot.cpp \
../vlvoldata.cpp \
../vlvolume.cpp 

OBJS += \
./MarchingCubes.o \
./ffourier.o \
./fkernel.o \
./fmask.o \
./fmat.o \
./fmod.o \
./fpad.o \
./frotate.o \
./fsurface.o \
./ftransform.o \
./fwr.o \
./vlkernel.o \
./vlmrot.o \
./vlvoldata.o \
./vlvolume.o 

CPP_DEPS += \
./MarchingCubes.d \
./ffourier.d \
./fkernel.d \
./fmask.d \
./fmat.d \
./fmod.d \
./fpad.d \
./frotate.d \
./fsurface.d \
./ftransform.d \
./fwr.d \
./vlkernel.d \
./vlmrot.d \
./vlvoldata.d \
./vlvolume.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpx -O3 -inline-level=2 -qmkl=sequential -I../include -I../../../src -I/opt/local/include -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


