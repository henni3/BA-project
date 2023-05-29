#include "tsp-main-helper.cu.h"

/* 
 * To run this program insert example "make test && ./test ../Data/hardcode.tsp 5 50" 
 * to terminal.
 * This main program is created for testing the result with different number 
 * climbers on different datasets. These tests are only performed on version 1.
 */
int main(int argc, char* argv[]) {
    int restarts, maxRestarts;
    char* file_name;
    if (argc != 4) {
        printf("Usage: %s <file-name> <minimum-number-climbers> <maximum-number-climbers>\n", argv[0]);
        exit(1);
    }
    // Collect input arguments
    file_name = argv[1];
    restarts = atoi(argv[2]);
    if(restarts <= 0){
        printf("Number of minimum climbers has to be a number larger than 0");
        exit(1);
    }
    maxRestarts = atoi(argv[3]);
    if(maxRestarts <= 0){
        printf("Number of maximum climbers has to be a number larger than 0");
        exit(1);
    }else if(maxRestarts < restarts){
        printf("Number of maximum climbers has to be larger than minimum restarts");
        exit(1);
    }

    initHwd();

    while(restarts <= maxRestarts){
        printf("\n\nResult with %d climbers: \n", restarts);
        runProgram(file_name, restarts, 1);
        restarts = restarts*10;
    }
   

    return 0;
}