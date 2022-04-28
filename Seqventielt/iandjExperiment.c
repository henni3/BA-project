#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    if(argc != 2){
        printf("Usage: %s <cities>\n", argv[0]);
        exit(1);
    }
    //Get cities
    int cities = atoi(argv[1]);
    //Calculate total number of iterations
    int totIter = ((cities-1) * (cities-2))/2;

    int next, i, j, d;
    float tmp;
    //compute i and j
    for(int gloId = 0; gloId < totIter; gloId++){
        d = 1-(4*(-2*(totIter-gloId)));
        //printf("d %d\n", d);
        tmp = (((-1-(sqrt(d)))/2)*(-1))+0.9999;
        //printf("tmp %f\n", tmp);
        next = (int) tmp;
        //printf("next %d\n", next);
        i = (cities-2) - (next-1);
        j = (i+2) + (gloId-(totIter-((next*(next-1))/2)));
        printf("global id: %d, i: %d, j: %d\n", gloId, i, j);
    }
}
