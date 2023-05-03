#include "tsp-main-helper.cu.h"

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Usage: %s <file-name> <number-of-restarts (multiple of 50)> \n", argv[0]);
        exit(1);
    }
    // Collect input arguments
    char* file_name = argv[1];
    int restarts = atoi(argv[2]);
    if (restarts % 10 != 0) {
        printf("restarts must be multiple of 10 \n");
        exit(1);
    } 
    if(restarts <= 0){
        printf("Number of restarts has to be a number larger than 0");
        exit(1);
    }
    initHwd();
    for (int i = 50; i <= restarts; i += 50) {
        runProgram(file_name, i, 4);
        //printf("iteration %d \n", i);
    }

    return 0;
}