#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "tsp_data.h"
#include <limits.h>

void createTour(int* tour, int cities, int startingCity){
    for(int t = 1; t < cities; t++){
        tour[t] = t;
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

int fitness (int*distM, int* tour, int cities) {
    int cost = 0;
    for (int i = 0; i < cities; i++) {
        cost += distM[cities * tour[i] + tour[i+1]];
    }
    return cost;
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

int* twoOptMove(int* distM, int cities, int restarts){
    int bestCost = INT_MAX;
    int cost = 0;
    int change, mini, minj, ti, tiplus1, tj, tjplus1;
    int* iterativeTour = malloc((cities + 1 ) * (sizeof(int)));
    int* bestTour = malloc((cities + 1) * (sizeof(int)));
    int iteration_count = 0;
    //memcpy(iterativeTour, tour, sizeof(int) * (cities +1));
    for (int i = 0; i < restarts; i++){
        createTour(iterativeTour, cities, 0);
        int minChange = -1;
        int while_count = 0;
        while(minChange < 0){
            iteration_count++;
            while_count ++;
            minChange = 0;
            for(int i = 0; i < cities - 2; i++){
                ti = iterativeTour[i];
                tiplus1 = iterativeTour[i+1] ; 
                for(int j = i+2; j < cities; j++){
                    tj = iterativeTour[j];
                    tjplus1 = iterativeTour[j+1];
                    change = distM[cities * ti + tj] + distM[cities * tiplus1 + tjplus1] - (distM[cities * ti + tiplus1] + distM[cities * tj + tjplus1]);
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
        cost = fitness(distM, iterativeTour,cities);
        if(cost < bestCost){
            memcpy(bestTour, iterativeTour, sizeof(int) * (cities +1));
            bestCost = cost;
        }
        // number of climbers 
        printf("%d\n", iteration_count);
        // tour cost
        //printf("%d\n", bestCost);
    }   
    free(iterativeTour);
    return bestTour;
}


int main(int argc, char *argv[]) {
    if(argc != 3) {
        printf("Expected two arguments. data set,  how many climbers \n");
        return 0;
    }
    srand(123456);
    char* fileName = argv[1];
    int restarts = atoi(argv[2]);
    int* distM = (int*) malloc(sizeof(int) * MAXCITIES * MAXCITIES);
    int cities = fileToDistM(fileName, distM);
    printf("number of cities is %d \n", cities);
    if( cities > MAXCITIES){
        printf("too many cities :( \n");
        exit(1);
    }
    distM = (int*) realloc(distM,sizeof(int) * cities * cities);
    /*for (int i = 0; i < cities; i++ ) {
        for (int j= 0; j < cities; j++) {
            printf("%d,", distM[cities* i + j]);
        }
        printf("\n");
    }*/

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
    int* opt_tour = twoOptMove(distM, cities, restarts);
    int newCost = 0;
    printf("best tour found [");
    for (int i = 0; i < cities +1 ; i++ ) {
        printf(" %d,", opt_tour[i]);
        if(i < cities){
            //printf("edge is edge (%d, %d) \n", opt_tour[i], opt_tour[i + 1]);
            newCost += distM[cities * opt_tour[i] + opt_tour[i+1]];
            //printf("new_cost is = %d \n",newCost);
            //printf("distM: %d\n", distM[cities * opt_tour[i] + opt_tour[i+1]]);
        }
    }
    printf("] \n");
    printf("Best cost: %d\n", newCost);
    free(opt_tour);
    free(distM);
    //free(tour);
}