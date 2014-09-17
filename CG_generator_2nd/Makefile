CC=g++
CFLAGS = -O3 -std=c++0x
LDFLAGS = -lz -lm 
SOURCES = generator.cpp
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = generator

all: generator
 
.c.o:
	$(CC) $(CFLAGS) $< -o $@

generator: function.cpp function.h generator.cpp
	$(CC) $(CFLAGS) function.cpp generator.cpp -o generator

