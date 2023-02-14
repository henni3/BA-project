#include <stdio.h>
#include <stdlib.h>

//Two optimization swapping
void twoOptSwap(int* itTour, int i, int j){
    int minI = i+1;
    int minJ = j;
    while ( minI < minJ) {
        int to = itTour[minJ];
        itTour[minJ] = itTour[minI];
        itTour[minI] = to;
        minI++;
        minJ--;
    }
}

//Calculate the change of the swap
int change (int** distM, int* tour, int i, int j){
    int ti, tj, tiplus1, tjplus1;
    ti = tour[i]; tiplus1 = tour[i+1];
    tj = tour[j]; tjplus1 = tour[j+1];
    return (distM[ti][tj] + distM[tiplus1][tjplus1] - (distM[ti][tiplus1] + distM[tj][tjplus1]));
}

//Borrowed from https://stackoverflow.com/questions/1202687/how-do-i-get-a-specific-range-of-numbers-from-rand
//Compute random number in a specific range
int random(int min, int max){
   return min + rand() / (RAND_MAX / (max - min + 1) + 1);
}

// Simulated Annealing
int* simulatedAnnealing(int** distM, int* currTour, int cities, int maxIter){
    //initialising
    int iter, i, j, currChange, optChange;
    double m, temperature;
    int* optTour = malloc((cities + 1 ) * sizeof(int));
    m = (double)(cities * (cities-1))/2; //give another name
    temperature = pow(m,3.0);
    iter = 0;
    //compute simulated annealing
    while(iter < maxIter){
        if(RAND_MAX < cities){ //should we have an alternativ way to handle this problem?
            printf("Usage: number of cities (%d) is larger than RAND_MAX (%d)\n", cities, RAND_MAX);
            exit(1);
        }
        //Tjek efter om denne måde er den bedste måde at lave random på
        i = random(0,(cities - 3)); //select random city
        j = random(i+2,(cities - 1)); //select random city
        if(i =! j) { //Det gør de ikke med ovenstående måde at finde i og j
            memcpy(optTour, currTour, sizeof(int) * (cities + 1));
            twoOptSwap(optTour, i, j);
            currChange = change(distM, currTour, i, j);
            optChange = change(distM, optTour, i, j);
            //double check this is correct
            if(currChange > optChange){
                memcpy(currTour, optTour, sizeof(int) * (cities + 1));
            }else{

            }
            
        }

    }
    free(optPath);

}

int main() {
    int row = 6;
    int column = 6;
    int cities = column - 1;
    int** distM = malloc(row*sizeof(int*));
    for (int i = 0; i < row; i++) {
        distM[i] = malloc(column * sizeof(int));
    }
   
    memcpy(distM[0], (int[6]) {0,1,2,3,4,5}, sizeof(int) * (column));
    memcpy(distM[1], (int[6]) {1,0,4,6,8,3}, sizeof(int) * (column));
    memcpy(distM[2], (int[6]) {2,4,0,4,5,2}, sizeof(int) * (column));
    memcpy(distM[3], (int[6]) {3,6,4,0,2,3}, sizeof(int) * (column));
    memcpy(distM[4], (int[6]) {4,8,5,2,0,4}, sizeof(int) * (column));
    memcpy(distM[5], (int[6]) {5,3,2,3,4,0}, sizeof(int) * (column));

    int* tour = malloc((cities + 1 ) * sizeof(int));
    createTour(tour, cities);
    int * opt_tour = twoOptMove(distM, tour, cities);
    int oldCost = 0;
    int newCost = 0;
    for (int i = 0; i < cities +1 ; i++ ) {
        if(i < cities){
            oldCost += distM[tour[i]][tour[i+1]];
            newCost += distM[opt_tour[i]][opt_tour[i+1]];
            //printf("new_cost is = %d \n",newCost);
            //printf("distM: %d\n", distM[opt_tour[i]][opt_tour[i+1]]);
        }
        printf("%d is City %d \n", i, opt_tour[i]);
    }
    printf("Old cost: %d, new cost: %d\n", oldCost, newCost);
    for (int i = 0; i < row; i++) {
        free(distM[i]);
    }
    free(opt_tour);
    free(distM);
    free(tour);
}