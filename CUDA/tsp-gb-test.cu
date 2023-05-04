#include "tsp-main-helper.cu.h"

int main(int argc, char* argv[]) {
    if (argc < 3) {
        printf("Usage: %s <number-of-restarts> <list of files>\n", argv[0]);
        exit(1);
    }
    // Collect input arguments
    char* file_names[argc-2];
    for( int i = 2; i < argc; i++) {
        file_names[i-2] = argv[i];
    }
    int restarts = atoi(argv[1]);
    if(restarts <= 0){
        printf("Number of restarts has to be a number larger than 0");
        exit(1);
    }
    //int cities = atoi(argv[3]);
    int jump = restarts / 10;
    initHwd();
    
    for (char* dataset : file_names ) {
         printf("testing for %s, each run repeated 100 times and taken average time \n", dataset);
        for (int i = jump; i <= restarts; i += jump) {
        runProgram(dataset, i, 4);
        //printf("iteration %d \n", i);
        printf("\n\n");
    }
    }

    return 0;
}