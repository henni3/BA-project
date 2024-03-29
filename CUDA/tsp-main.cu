#include "tsp-main-helper.cu.h"
/* 
 * To run this program insert example "make tsp && ./tsp ../Data/hardcode.tsp 5 1" 
 * to terminal.
 * This main program is created for running the 4 different versions on different number 
 * of climbers on different datasets.
 */
int main(int argc, char* argv[]) {
    if (argc != 4) {
        printf("Usage: %s <file-name> <number-of-climbers> <Which program version? 1 (original), 2 (100Cities), 3 (calculatedIandJ) or 4 (GB/s)>\n", argv[0]);
        exit(1);
    }
    // Collect input arguments
    char* file_name = argv[1];
    int restarts = atoi(argv[2]);
    if(restarts <= 0){
        printf("Number of restarts has to be a number larger than 0");
        exit(1);
    }
    int version = atoi(argv[3]);
    if((version < 1) || (version > 4)){
        printf("Wrong program version. You can choose between 1, 2 or 3, or 4 for byte count");
        exit(1);
    }
    
    initHwd();

    runProgram(file_name, restarts, version);

    return 0;
}