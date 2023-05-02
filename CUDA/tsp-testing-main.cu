#include "tsp-main-helper.cu.h"

int main(int argc, char* argv[]) {
    int restarts, maxRestarts;
    char* file_name;
    if (argc != 4) {
        printf("Usage: %s <file-name> <minimum-number-restarts> <maximum-number-restarts>\n", argv[0]);
        exit(1);
    }
    // Collect input arguments
    file_name = argv[1];
    restarts = atoi(argv[2]);
    if(restarts <= 0){
        printf("Number of minimum restarts has to be a number larger than 0");
        exit(1);
    }
    maxRestarts = atoi(argv[3]);
    if(maxRestarts <= 0){
        printf("Number of maximum restarts has to be a number larger than 0");
        exit(1);
    }else if(maxRestarts < restarts){
        printf("Number of maximum restarts has to be larger than minimum restarts");
        exit(1);
    }
    char* true_name = "../Data/swiss42.tsp";

    initHwd();

    while(restarts <= maxRestarts){
        printf("\n\nResult with %d climbers: \n", restarts);
        runProgram(true_name, restarts, 1);
        restarts = restarts*10;
    }
   

    return 0;
}