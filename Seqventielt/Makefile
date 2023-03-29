

CC = gcc
CFLAGS=-Wall -Wextra -pedantic -std=gnu99 -O3

SA = sim_ann
OPT = tsp_opt

$(SA): tsp_data.h
	$(CC) TSP-SA.c -o sa $(CFLAGS) -lm
