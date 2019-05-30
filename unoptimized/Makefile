# Much stolen from lazyfoo's SDL tutorial

OBJS = pf.cpp

# Compiler
CC = g++

# Warnings
COMPILER_FLAGS =  -Wall
COMPILER_FLAGS += -Wextra
COMPILER_FLAGS += -Werror

# Optimization
COMPILER_FLAGS += -Os
COMPILER_FLAGS += -ggdb

# Rando flags
COMPILER_FLAGS += -std=c++14

# Linker flags
LINKER_FLAGS   += -lSDL2

OBJ_NAME = particlefilter

all: $(OBJS)
	$(CC) $(OBJS) $(COMPILER_FLAGS) $(LINKER_FLAGS) -o $(OBJ_NAME)
