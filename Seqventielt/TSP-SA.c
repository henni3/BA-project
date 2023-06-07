#include "tsp_data.h"
#include <time.h>

// creating tours
void createTour(int* tour, int cities){
    for(int t = 1; t < cities; t++){
        tour[t] = t;
        //printf ("tour creatinon yeilds %d \n", tour[t]);
    }
    tour[0] = 0;

    int randNumFrom, randNumTo, to;
    for (int i = 0; i < cities*10; i++) {
        randNumFrom = rand() % cities;
        randNumTo = rand() % cities;
        if (randNumFrom && randNumTo){
        to = tour[randNumTo];
        tour[randNumTo] = tour[randNumFrom];
        tour[randNumFrom] = to;
        }
    }
    tour[cities] = tour[0];
}

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

int fitness (int*distM, int* tour, int cities) {
    int cost = 0;
    for (int i = 0; i < cities; i++) {
        cost += distM[cities * tour[i] + tour[i+1]];
    }
    return cost;
}

//Borrowed from https://stackoverflow.com/questions/1202687/how-do-i-get-a-specific-range-of-numbers-from-rand
//Compute random number in a specific range
int randomRange(int min, int max){
   return rand() % (max + 1 - min) + min;
}

// Simulated Annealing
int simulatedAnnealing(int* distM, int* currTour, int cities ){
    //initialising
    int iter, i, j, currChange, optChange;
    double c, temperature;
    int* optTour = malloc((cities + 1 ) * sizeof(int));
    int edges = (cities * (cities-1))/2; //give another name
    int start_temp = pow(edges,3);
    temperature = start_temp;
    //printf("%d \n", start_temp);
    printf("iterations, cost\n");
    iter = 0;
    int equal_counter = 0;
    //compute simulated annealing
    while(temperature > 0.00025){
        //printf("we get here \n");
        if(RAND_MAX < cities){ //should we have an alternativ way to handle this problem?
            printf("Usage: number of cities (%d) is larger than RAND_MAX (%d)\n", cities, RAND_MAX);
            exit(1);
        }
        iter++;
        //Tjek efter om denne måde er den bedste måde at lave random på
        i = randomRange(0,(cities - 1)); //select random city
        j = randomRange(0,(cities - 1)); //select random city
        if(i != j) { //Det gør de ikke med ovenstående måde at finde i og j
            memcpy(optTour, currTour, sizeof(int) * (cities + 1));
            twoOptSwap(optTour, i, j);
            //currChange = change(distM, currTour, i, j);
            //optChange = change(distM, optTour, i, j);
            //double check this is correct
            int currFit = fitness(distM, currTour, cities);
            int newFit = fitness(distM, optTour, cities);
            if( currFit > newFit ){
                memcpy(currTour, optTour, sizeof(int) * (cities + 1));
                equal_counter = 0;
            }
            else if (currFit == newFit)
            {
                equal_counter++;
                if(equal_counter >= cities*10) {
                    //break;
                }
            }
            
            else{
                c = ((double)rand()) / (float)RAND_MAX;
                int costDiff = newFit - currFit;
                double check = exp(- (double)costDiff / temperature);
                if (c < check) {
                    memcpy(currTour, optTour, sizeof(int) * (cities + 1));
                    equal_counter = 0;
                }

            }
            
        }
        
        //double factor = 1.56E5;
        double alpha = 0.89;
        temperature = start_temp * pow(alpha,iter);
        if(iter % 10 == 0) {
            printf("%d, %d \n", iter, fitness(distM,currTour, cities));
            //printf("and temperature %f \n", temperature);
        }

    }
    int finalCost = fitness(distM,currTour, cities);
    printf("cost of tour %d after simulated annelaing with %d iterations \n", finalCost, iter);
    free(optTour);
    return iter;

}

int main(int argc, char** argv) {
    if (argc != 2) {
        exit(1);
    }
    srand(123456);
    char* fileName = argv[1];
    int* distMatrix = (int*) malloc(sizeof(int) * MAXCITIES * MAXCITIES);
    int cities = fileToDistM(fileName, distMatrix);
    printf("number of cities is %d \n", cities);
    if( cities > MAXCITIES){
        printf("too many cities :( \n");
        exit(1);
    }
    distMatrix = (int*) realloc(distMatrix,sizeof(int) * cities * cities);
    /*int row = 6;
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
    */
    int* tour = malloc((cities + 1 ) * sizeof(int));
    createTour(tour, cities);
    printf("tour is:");
    for (int i =0; i < cities + 1 ; i++) {
        printf("%d, ", tour[i]);
    }
    printf("\n done \n");

    //int * opt_tour = twoOptMove(distM, tour, cities);
    int oldCost = 0;
    int newCost = 0;
    simulatedAnnealing(distMatrix, tour, cities);
    printf("tour after annealing is:");
    for (int i =0; i < cities + 1 ; i++) {
        printf("%d, ", tour[i]);
    }
    //for (int i = 0; i < row; i++) {
    //    free(distM[i]);
    //}
    free(distMatrix);
    free(tour);
    return 0;
}