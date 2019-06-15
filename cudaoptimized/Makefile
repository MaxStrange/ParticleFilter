# Much stolen from lazyfoo's SDL tutorial

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
ifeq ($(OS),Windows_NT)
	COMPILER_FLAGS += -Xlinker /subsystem:windows
endif

# Linker Path
ifeq ($(OS),Windows_NT)
	LINKER_PATH    := -L$(ROOTDIR)\SDL\lib\win64
else
	LINKER_PATH    :=
endif

# Linker flags
ifeq ($(OS),Windows_NT)
	LINKER_FLAGS   := -lSDL2main
	LINKER_FLAGS   += -lSDL2
else
	LINKER_FLAGS   := -lSDL2
endif

# Includes
INCDIRS  := $(ROOTDIR)/common
INCDIRS  += $(ROOTDIR)/unoptimized
INCDIRS  += $(ROOTDIR)/cudaoptimized
INCLUDES := $(patsubst %, -I%, $(INCDIRS))

ifeq ($(OS),Windows_NT)
	OBJ_NAME := particlefilter.exe
else
	OBJ_NAME := particlefilter.out
endif

all: $(OBJS)
	$(CC) $(OBJS) $(INCLUDES) $(COMPILER_FLAGS) $(LINKER_PATH) $(LINKER_FLAGS) -o $(OBJ_NAME)
