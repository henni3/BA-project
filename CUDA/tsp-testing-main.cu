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
    }else if(maxRestarts <= restarts){
        printf("Number of maximum restarts has to be larger than minimum restarts");
        exit(1);
    }

    initHwd();

    while(restarts < maxRestarts){
        printf("\nResults from version 1 w %d climbers: \n", restarts);
        runProgram(file_name, restarts, 1);

        printf("\nResults from version 2 w %d climbers: \n", restarts);
        runProgram(file_name, restarts, 2);

        printf("\nResults from version 3 w %d climbers: \n", restarts);
        runProgram(file_name, restarts, 3);
        
        restarts++;
    }
   

    return 0;
}