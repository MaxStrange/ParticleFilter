# Much stolen from lazyfoo's SDL tutorial

# Set to non-blank for SDL, otherwise we will not have a GUI, just console values.
SDL := yes
ifeq ($(OS),Windows_NT)
	IS_WINDOWS := yes
else
	IS_WINDOWS :=
endif

# Variables
ROOTDIR := ..

# Source
OBJS := pf.cpp
OBJS += cudaparticlefilter.cpp
OBJS += device.cu
OBJS += $(ROOTDIR)/common/graphics.cpp
OBJS += $(ROOTDIR)/common/particlefilter.cpp
OBJS += $(ROOTDIR)/common/robot.cpp

# Compiler
CC := nvcc

# Rando flags
COMPILER_FLAGS := --std=c++14
ifneq ($(and $(IS_WINDOWS),$(SDL)),)
	COMPILER_FLAGS += -Xlinker /subsystem:windows
endif

# Linker Path
ifneq ($(and $(IS_WINDOWS),$(SDL)),)
	LINKER_PATH    := -L$(ROOTDIR)\SDL\lib\win64
else
	LINKER_PATH    :=
endif

# Linker flags
ifneq ($(and $(IS_WINDOWS),$(SDL)),)
	LINKER_FLAGS   := -lSDL2main
	LINKER_FLAGS   += -lSDL2
else ifeq ($(SDL),yes)
	LINKER_FLAGS   := -lSDL2
else
	LINKER_FLAGS   :=
endif

# Includes
INCDIRS  := $(ROOTDIR)/common
INCDIRS  += $(ROOTDIR)/unoptimized
INCDIRS  += $(ROOTDIR)/cudaoptimized
INCLUDES := $(patsubst %, -I%, $(INCDIRS))

ifeq ($(IS_WINDOWS),yes)
	OBJ_NAME := particlefilter.exe
else
	OBJ_NAME := particlefilter.out
endif

# Defines
ifeq ($(SDL),yes)
	DEFINES :=
else
	DEFINES := -DDEBUG
endif

all: $(OBJS)
	$(CC) $(OBJS) $(INCLUDES) $(COMPILER_FLAGS) $(DEFINES) $(LINKER_PATH) $(LINKER_FLAGS) -o $(OBJ_NAME)
