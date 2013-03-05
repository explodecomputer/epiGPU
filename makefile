CC = gcc
CFLAGS = -c -O2 -Wall
OP = -lOpenCL

EXE1 = exe/linux/epiGPU

UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
OP = -framework opencl
EXE1 = exe/mac/epiGPU
endif

LDFLAGS = -lm -Wall -O2 $(OP)
SOURCES1 = src/main.c src/EpiOpenCL.c src/readplink.c src/timing.c
OBJECTS1 = $(SOURCES1:.c=.o)


all: $(SOURCES1) $(EXE1)


$(EXE1): $(OBJECTS1)
	$(CC) $(LDFLAGS) $(OBJECTS1) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

all:
	rm $(OBJECTS1)

clean:
	rm $(OBJECTS1) $(EXE1)



