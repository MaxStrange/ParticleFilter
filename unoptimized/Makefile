# Much stolen from lazyfoo's SDL tutorial

# Variables
ROOTDIR := ..

# Source
OBJS := pf.cpp
OBJS += $(ROOTDIR)/common/graphics.cpp

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

# Linker flags
LINKER_FLAGS   += -lSDL2

# Includes
INCDIRS  := $(ROOTDIR)/common
INCDIRS  += $(ROOTDIR)/unoptimized
INCLUDES := $(patsubst %, -I%, $(INCDIRS))

OBJ_NAME := particlefilter.out

all: $(OBJS)
	$(CC) $(OBJS) $(INCLUDES) $(COMPILER_FLAGS) $(LINKER_FLAGS) -o $(OBJ_NAME)
