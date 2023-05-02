#include "tsp-main-helper.cu.h"

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Usage: %s <file-name> <number-of-restarts (multiple of 10)> \n", argv[0]);
        exit(1);
    }
    // Collect input arguments
    char* file_name = argv[1];
    int restarts = atoi(argv[2]);
    if(restarts <= 0){
        printf("Number of restarts has to be a number larger than 0");
        exit(1);
    }
    initHwd();

    runProgram(file_name, restarts, 4);

    return 0;
}