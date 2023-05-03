#include "tsp-main-helper.cu.h"

int main(int argc, char* argv[]) {
    if (argc != 4) {
        printf("Usage: %s <file-name> <number-of-restarts (multiple of 50)> <cities> \n", argv[0]);
        exit(1);
    }
    // Collect input arguments
    char* file_name = argv[1];
    int restarts = atoi(argv[2]);
    if (restarts % 50 != 0) {
        printf("restarts must be multiple of 50 \n");
        exit(1);
    } 
    if(restarts <= 0){
        printf("Number of restarts has to be a number larger than 0");
        exit(1);
    }
    int cities = atoi(argv[3]);
    int jump = cities;
    initHwd();
    for (int i = 50; i <= restarts; i += jump) {
        runProgram(file_name, i, 4);
        //printf("iteration %d \n", i);
    }

    return 0;
}