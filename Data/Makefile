CC = gcc
CFLAGS =-Wall -Wextra -pedantic -std=gnu11
LD_FLAGS = -lpthread

all : convert

rebuild : clean all

convert : convert.c
	$(CC) $(CFLAGS) $< -o $@ $(LD_FLAGS)


clean :
	rm -rf *.o convert

../src.zip :
	make clean
	cd .. && zip src.zip -r src
