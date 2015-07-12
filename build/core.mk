PWD := $(shell pwd)
BUILD_SYSTEM := $(dir $(lastword $(MAKEFILE_LIST)))
TOP := $(BUILD_SYSTEM)..
CLEAR_VARS:= $(BUILD_SYSTEM)/clear_vars.mk
BUILD_EXECUTABLE:= $(BUILD_SYSTEM)/executable.mk
BUILD_STATIC_LIBRARY:= $(BUILD_SYSTEM)/static_library.mk
BUILD_SHARED_LIBRARY:= $(BUILD_SYSTEM)/shared_library.mk
BUILD_INTERNAL:= $(BUILD_SYSTEM)/build_internal.mk


define my-dir
$(shell pwd)
endef

CC=gcc
CPP=g++
AR=ar
RANLIB=ranlib

#CC=arm-none-eabi-gcc
#CPP=arm-none-eabi-g++
#AR=arm-none-eabi-ar
#RANLIB=arm-none-eabi-ranlib
