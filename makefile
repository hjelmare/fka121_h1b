
CC = gcc
CFLAGS = -O3
LIBS = -lm -lgsl -lgslcblas


HEADERS = h1b_func.h 
OBJECTS = h1b_func.o h1b_main.o
PROGRAM = H1

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

