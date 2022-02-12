#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void twoOptSwap(int* itTour, int i , int j){
    int minI = i-1;
    int minJ = j;
    while ( minI < minJ) {
        int to = itTour[minJ];
        itTour[minJ] = itTour[minI];
        itTour[minI] = to;
        minI++;
        minJ--;
    }
    
}

int* twoOptMove(int** distM, int* tour, int cities){
    int minChange = -1;
    int change, mini, minj, ti, tiplus1, tj, tjplus1;
    int* iterativeTour = malloc((cities + 1 ) * (sizeof(int)));
    memcpy(iterativeTour, tour, sizeof(int) * (cities +1));
    
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
                if(minChange > change){
                    minChange = change;
                    printf("minchange is %d with index i %d and index j %d \n", minChange, i, j);
                    mini = i;
                    minj = j;
                }   
            }
        }
        twoOptSwap(iterativeTour, mini, minj);
        for (int i = 0; i < cities + 1; i++) {
            printf("%d ", iterativeTour[i]);
        }
        printf("\n");
    }
    return iterativeTour;
}



void createTour(int* tour, int cities){
    printf("please input the intial tour \n");
    for(int t = 0; t < cities; t++){
        printf("city %d \n", t+1);
        scanf("%d", &tour[t]);
        //tour[t] = t+1;
    }
    printf("done \n");
    tour[cities] = tour[0];
}



int main(int argc, char *argv[]) {
    if(argc != 2){
        printf("Wrong number of arguments. Write how many restars should be taken: %d\n", argc);
    }
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

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j ++) {
            //printf("enter value for distance between city %d and city %d \n", i, j);
            //scanf("%d", &distM[i][j]);
            printf("%d ",distM[i][j]);
        }
        printf("\n");
    }

    int* tour = malloc((cities + 1 ) * sizeof(int));
    int numRestarts = atoi(argv[1]);
    createTour(tour, cities);
    int * opt_tour = twoOptMove(distM, tour, cities);
    for (int i = 0; i < cities +1 ; i++ ) {
        printf("%d is City %d \n", i, opt_tour[i]);
    }
    for (int i = 0; i < row; i++) {
        free(distM[i]);
    }
    free(opt_tour);
    free(distM);
    free(tour);
}