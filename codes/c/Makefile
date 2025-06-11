CC = gcc
CFLAGS = -Wall -Wextra -pedantic -O3 -ggdb -fopenmp
LDFLAGS = -lm -lgmp -lmpfr


all: main

main: main.c
	$(CC) $< $(CFLAGS) $(LDFLAGS) -o $@

clean:
	rm -f main
	
.PHONY: all clean


