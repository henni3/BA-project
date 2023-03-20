#include "tsp-main-helper.cu.h"

int main() {
    int restarts, maxRestarts;
    const *char file_name[21] = "../Data/hardcode.tsp";
    //const char* file_name[21] = ;
    initHwd();

    restarts = 1;
    maxRestarts = 15;
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