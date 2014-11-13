
CC = gcc
CFLAGS = -O3 -g
LIBS = -lm

HEADERS = initfcc.h alpotential.h MD_functions.h
OBJECTS = initfcc.o alpotential.o MD_functions.o MD_main.o
PROGRAM = MD

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

