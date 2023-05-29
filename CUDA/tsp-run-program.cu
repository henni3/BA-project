#include "tsp-main-helper.cu.h"

/* 
 * Run program by inserting "make pre_pro && ./pre_pro" to terminal. 
 * This will run the program with dataset kroA100.tsp with 75000 
 * climbers on version 1.
 */
int main() {
    int restarts;
    char* file_name;
    restarts = 75000;
    file_name = "../Data/kroA100.tsp";

    initHwd();
    printf("\n\nRun program on dataset %s. Result with %d climbers: \n", file_name, restarts);
    runProgram(file_name, restarts, 1);
   
    return 0;
}