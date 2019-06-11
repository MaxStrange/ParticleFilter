# Much stolen from lazyfoo's SDL tutorial

# Variables
ROOTDIR := ..

# Source
OBJS := pf.cpp
OBJS += unoptimizedparticlefilter.cpp
OBJS += $(ROOTDIR)/common/graphics.cpp
OBJS += $(ROOTDIR)/common/particlefilter.cpp
OBJS += $(ROOTDIR)/common/robot.cpp

# Compiler
CC := g++

# Warnings
COMPILER_FLAGS := -Wall
COMPILER_FLAGS += -Wextra
COMPILER_FLAGS += -Werror

# Optimization
COMPILER_FLAGS += -Os
COMPILER_FLAGS += -ggdb

# Rando flags
COMPILER_FLAGS += -std=c++14

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
	LINKER_FLAGS   += -lmingw32
else
	LINKER_FLAGS   := -lSDL2
endif

# Includes
INCDIRS  := $(ROOTDIR)/common
INCDIRS  += $(ROOTDIR)/unoptimized
INCLUDES := $(patsubst %, -I%, $(INCDIRS))

ifeq ($(OS),Windows_NT)
	OBJ_NAME := particlefilter.exe
else
	OBJ_NAME := particlefilter.out
endif

all: $(OBJS)
	$(CC) $(OBJS) $(INCLUDES) $(COMPILER_FLAGS) $(LINKER_PATH) $(LINKER_FLAGS) -o $(OBJ_NAME)
