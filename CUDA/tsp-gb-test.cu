#include "tsp-main-helper.cu.h"

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Usage: %s <file-name> <number-of-restarts>\n", argv[0]);
        exit(1);
    }
    // Collect input arguments
    char* file_name = argv[1];
    int restarts = atoi(argv[2]);
    if(restarts <= 0){
        printf("Number of restarts has to be a number larger than 0");
        exit(1);
    }
    //int cities = atoi(argv[3]);
    int jump = restarts / 10;
    initHwd();
    printf("testing for %s, each run repeated 100 times and taken average time \n", file_name);
    for (int i = 10; i <= restarts; i += jump) {
        runProgram(file_name, i, 4);
        //printf("iteration %d \n", i);
    }

    return 0;
}