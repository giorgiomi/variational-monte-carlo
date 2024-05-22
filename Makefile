CC = gcc
CFLAGS = -I${CURDIR}/include
OBJECTS = src/main.c src/functions.c
run: $(OBJECTS)
	$(CC) $(CFLAGS) -o run $(OBJECTS)