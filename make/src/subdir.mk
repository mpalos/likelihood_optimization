################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/functions.c \
../src/hcubature.c \
../src/likelihood.c \
../src/optim.c 

OBJS += \
./src/functions.o \
./src/hcubature.o \
./src/likelihood.o \
./src/optim.o 

C_DEPS += \
./src/functions.d \
./src/hcubature.d \
./src/likelihood.d \
./src/optim.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


