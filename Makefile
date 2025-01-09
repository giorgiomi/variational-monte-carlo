CC = gcc
CFLAGS = -I${CURDIR}/include -O2
OBJECTS = src/main.c src/functions.c src/metropolis.c
run: $(OBJECTS)
	$(CC) $(CFLAGS) -o run $(OBJECTS)