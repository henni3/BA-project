

CC = gcc
CFLAGS=-Wall -Wextra -pedantic -std=gnu99 -O3

SA = sim_ann
OPT = tsp_opt
ACO = tsp_ACO

$(SA): tsp_data.h
	$(CC) TSP-SA.c -o sa $(CFLAGS) -lm

$(ACO): tsp_data.h
	$(CC) TSP-ACO.c -o aco $(CFLAGS) -lm

$(OPT): tsp_data.h
	$(CC) TSP-seqRand.c -o 2opt $(CFLAGS) -lm

all: $(SA) $(OPT) $(ACO)

clean:
	rm -f $(SA) $(OPT) $(ACO)
