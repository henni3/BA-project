#include "tsp_data.h"
#include <time.h>
#include <float.h>

double beta;

double alpha;

double rho;

double Q;


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
    //tour[cities] = tour[0];
}

int fitness (int*distM, int* tour, int cities) {
    int cost = 0;
    for (int i = 0; i < cities; i++) {
        if (i >= cities - 1) {
            cost += distM[cities * tour[i] + tour[0]];
            break;
        }
        cost += distM[cities * tour[i] + tour[i+1]];
    }
    return cost;
}

typedef struct Ant {
    int* path;
    //int cost;

} Ant;

int findIndex(int* array, int target, int length) {
    for (int i = 0; i < length; i++) {
        if (array[i] == target) {
            return i;
        }
    }
    printf("target could not be found in array \n");
    return -1;
}

double* GenerateMoveProbabilities (int cities, int cityI, int* choseArray, double* pher, int* distM) {

    double* tau = malloc(sizeof(double) * cities);
    double sum = 0.0;

    for (int i = 0; i < cities; i++) {
        if ( i == cityI || choseArray[i]) {
            tau[i] = 0.0;
        }
        else {
            double prob = pow(pher[cities * cityI + i], alpha) * 
            pow((1.0 / (double) distM[cities * cityI + i]), beta);
            if (prob < 0.0001) {

                prob = 0.0001;
            }
            else if (prob > (DBL_MAX / (cities * 100) )) {
                prob = DBL_MAX / (cities * 100);
            }

            tau[i] = prob;
        }
        sum += tau[i];


    }

    double* probs = malloc(sizeof(double) * cities);
    for (int i = 0; i < cities; i++ ) {
        probs[i] = tau[i] / sum;
    }
    free(tau);
    return probs;
}


int chooseCity(int cities, int cityI, int* choseArray, double* pher, int* distM) {
    
    double* probArray = GenerateMoveProbabilities(cities, cityI, choseArray, pher, distM);

    double* culm = malloc(sizeof(double) * (cities + 1));
    culm[0] = 0.0;
    for (int i = 0; i < cities; i++) {
        culm[i+1] = culm[i] + probArray[i];
    }

    double p = ((double) rand()) / ((double ) RAND_MAX) ;

    for (int i = 0; i < cities; i ++) {
        if (p >= culm[i] && p < culm[i+1])
        {
            return i;
        }
    }
    printf("no city could be found \n");
    return -1;


}

int* NewTrail (int start, double* pher, int cities, int* distM) {
    int* trail = malloc(sizeof(int) * cities);
    int* choose = malloc(sizeof(int) * cities);
    for (int i = 0; i < cities; i++) {
        choose[i] = 0;
    }
    choose[start] = 1;
    trail[0] = start;
    for (int i = 0; i < cities-1; i++) {
        int cityX = trail[i];
        int next = chooseCity(cities, cityX, choose, pher, distM);
        trail[i+1]= next;
        choose[next] = 1;
    }
    return trail;

}






Ant* intializeAnts (int numAnts, int numCities, int* distM) {
    Ant* antArray = malloc(sizeof(Ant) * numAnts);
    for (int i = 0; i < numAnts; i++) {
        antArray[i].path = malloc(sizeof(int) * numCities);
        createTour(antArray[i].path,numCities);
        //antArray[i].cost = fitness(distM, antArray[i].path, numCities);
    }
    return antArray;
} 

double* initPhero (int cities) {
    double* pheroArray = malloc(sizeof(double) * cities * cities);
    for (int i = 0; i < cities; i++) {
        for (int j = 0; j < cities; j++) {
            pheroArray[cities * i + j] = 0.01; 
        }
    }
    return pheroArray; 
}

void updateAnts (Ant* antArray, int numAnts, int cities, double* phero, int* distM) {
    for (int i = 0; i < numAnts; i++) {
        int start = rand() % cities;
        //based on all ants start with a path, before creating a new one
        free(antArray[i].path);
        antArray[i].path = NewTrail(start, phero, cities, distM);
    }
}

int isEdge(int cityI, int CityJ, int* trail, int length) {

    int last = length -1;
    int index = findIndex(trail, cityI, length);

    if (index == 0 && trail[1] == CityJ ) {
        return 1;
    }
    else if (index == 0 && trail[last] == CityJ) {
        return 1;
    }
    else if (index == last && trail[last-1] == CityJ) {
        return 1;
    }
    else if (index == last && trail[0] == CityJ ) {
        return 1;
    }
    else if(index == last) {
        return 0;
    }
    else if (trail[index - 1] == CityJ || trail[index + 1] == CityJ){
        return 1;
    }
    else {
        return 0;
    }

}

void updatePheromones(double* phero, Ant* ants, int* distM, int numAnts, int cities) {
    for (int i = 0; i < cities; i++) {
        for (int j = i + 1; j < cities; j++) {
            for (int k = 0; k < numAnts; k++) {
                double tourLength = (double) fitness(distM, ants[k].path, cities);
                double decrease = (1.0 - rho) * phero[cities * i + j];
                double increase = 0.0;
                if (isEdge(i,j,ants[k].path, cities)) {
                    increase = Q / tourLength;
                }

                double change = decrease + increase;
                if ( change < 0.0001) {
                    change = 0.0001;
                }

                if (change > 100000.0) {
                    change = 100000.0;
                }
                phero[cities * i + j] = change;
                phero[cities * j + i] = change;

            }
        }
    }

}

int shortestPath (Ant* ants, int* distM, int numAnts, int cities) {
    int best = fitness(distM, ants[0].path, cities);
    for (int i = 1; i < numAnts; i++) {
        int cost = fitness(distM, ants[i].path, cities);
        if (cost < best) {
            best = cost;
        }
    }
    return best;
}

int main(int argc, char** argv) {

    if (argc != 4) {
        printf("usage: filename, numAnts, Iterations \n");
        exit(0);
    }
    char* fileName = argv[1];
    int numAnts = atoi(argv[2]);
    int* distMatrix= (int*) malloc (sizeof(int) * MAXCITIES * MAXCITIES);
    int cities = fileToDistM(fileName, distMatrix);
    distMatrix = (int*) realloc(distMatrix, sizeof(int) * cities * cities);
    int MaxIter = atoi(argv[3]);
    Q = 2.0;
    alpha = 1.0;
    beta = 1.0;
    rho = 0.05;
    srand(12345);
    printf("we get to intialize the ants \n");
    Ant* ants = intializeAnts(numAnts, cities, distMatrix);
    printf("we get to intialize the ants 2 \n");
    double* pheromones = initPhero(cities);
    int bestCost = shortestPath(ants, distMatrix, numAnts, cities);
    for (int i = 0; i < MaxIter; i++) {
        updateAnts(ants,numAnts,cities,pheromones,distMatrix);
        updatePheromones(pheromones,ants,distMatrix,numAnts,cities);
        //printf("got through updating in iteration %d \n", i);
        int newCost = shortestPath(ants, distMatrix, numAnts, cities);
        if (newCost < bestCost) {
            bestCost = newCost;
        }

    }
    printf("done \n");
    printf("best cost was %d \n", bestCost);
    return 1;    
}