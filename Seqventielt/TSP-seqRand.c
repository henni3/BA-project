#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void createTour(int* tour, int cities, int startingCity){
    for(int t = 1; t < cities; t++){
        tour[t] = t+1;
        //printf ("tour creatinon yeilds %d \n", tour[t]);
    }
    tour[0] = startingCity;
    tour[cities] = tour[0];

    int randNum, to, index;
    for (int i = 1; i < cities; i++) {
        randNum = rand();
        index = randNum % cities;
        if (index) {
            to = tour[i];
            tour[i] = tour[index];
            tour[index] = to;
        }
    }
    /*printf ( "new route: \n") ;
    for(int t = 0; t < cities+1; t++){
        printf(" %d ", tour[t]);
    }
    printf("\n");*/
}


void twoOptSwap(int* itTour, int i , int j){
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

int* twoOptMove(int** distM, int cities, int restarts){
    int bestCost = INT32_MAX;
    int cost = 0;
    int minChange = -1;
    int change, mini, minj, ti, tiplus1, tj, tjplus1;
    int* iterativeTour = malloc((cities + 1 ) * (sizeof(int)));
    int* bestTour = malloc((cities + 1) * (sizeof(int)));
    //memcpy(iterativeTour, tour, sizeof(int) * (cities +1));
    for (int i = 0; i < restarts; i++){
        createTour(iterativeTour, cities, 1);
        while(minChange < 0){
            minChange = 0;
            for(int i = 0; i < cities - 2; i++){
                ti = iterativeTour[i];
                tiplus1 = iterativeTour[i+1] ; 
                for(int j = i+2; j < cities; j++){
                    tj = iterativeTour[j];
                    tjplus1 = iterativeTour[j+1];
                    change = distM[ti][tj] + distM[tiplus1][tjplus1] - (distM[ti][tiplus1] + distM[tj][tjplus1]);
                    //printf("change is %d \n", change);
                    if(change < minChange){
                        minChange = change;
                        //printf("minchange is %d with index i %d and index j %d \n", minChange, i, j);
                        mini = i;
                        minj = j;
                        //printf("after min\n");
                    }   
                }
            }
            if (minChange < 0){
                //printf("before swap\n");
                twoOptSwap(iterativeTour, mini, minj);
                //for (int i = 0; i < cities + 1; i++) {
                //    printf("%d ", iterativeTour[i]);
                //}
                //printf("\n");
                //printf("after swap\n");
            }
        }
        for (int i = 0; i < cities; i++ ) {
            cost += distM[iterativeTour[i]][iterativeTour[i+1]];
        }
        if(cost < bestCost){
            memcpy(bestTour, iterativeTour, sizeof(int) * (cities +1));
            bestCost = cost;
        }
    }
    free(iterativeTour);
    return bestTour;
}


int main(int argc, char *argv[]) {
    if(argc != 2) {
        printf("Expected one argument. Please enter how many climbers");
        return 0;
    }
    srand(time(NULL));
    int row = 6;
    int column = 6;
    int cities = column - 1;
    int restarts = atoi(argv[1]);
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

    /*for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j ++) {
            //printf("enter value for distance between city %d and city %d \n", i, j);
            //scanf("%d", &distM[i][j]);
            printf("%d ",distM[i][j]);
        }
        printf("\n");
    }*/

    //int* tour = malloc((cities + 1 ) * sizeof(int));
    //createTour(tour, cities, 1);
    int * opt_tour = twoOptMove(distM, cities, restarts);
    int newCost = 0;
    for (int i = 0; i < cities +1 ; i++ ) {
        if(i < cities){
            newCost += distM[opt_tour[i]][opt_tour[i+1]];
            //printf("new_cost is = %d \n",newCost);
            //printf("distM: %d\n", distM[opt_tour[i]][opt_tour[i+1]]);
        }
        printf("%d is City %d \n", i, opt_tour[i]);
    }
    printf("Best cost: %d\n", newCost);
    for (int i = 0; i < row; i++) {
        free(distM[i]);
    }
    free(opt_tour);
    free(distM);
    //free(tour);
}