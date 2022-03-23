#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXLINE 8192

enum DATA_TYPE {EUCLIDIAN,GEO,MATRIX};

int euclidain_float_to_int(float x1, float x2, float y1, float y2) {
    float dx = x1 - x2;
    float dy = y1 - y2;
    float e_dist = sqrt(dx * dx + dy * dy) + 0.5;
    return (int) e_dist;
}

void create_dist_array (int* dist_array, float* xs, float* ys, int type, int cities) {
    if (type == EUCLIDIAN) {
        int euclidian;
        for(int i = 0; i < cities; i++){
            for(int j = 0; j < cities; j++){
                euclidian = euclidain_float_to_int(xs[i], xs[j], ys[i], ys[j]);
                dist_array[i * cities + j] = euclidian;
                dist_array[j * cities + i] = euclidian;
            }
        }      
    }else{
        printf("This type is not supported \n.");
        exit(1);
    }
}

int fileToDistM(char* filename, int* save_array){
    FILE* source = fopen(filename, "r");
    if (source == NULL) {
        printf("Source file not found\n");
        return EXIT_FAILURE;
    }
    char buf[MAXLINE];
    int cities, read3, type, *distM;
    float read1,read2, *X_positions, *Y_positions;
    while(fscanf(source, "%s", buf)) {
        if (strncmp("DIMENSION", buf,9) == 0 ) {
            fscanf(source, "%d", &cities);
            distM = malloc(sizeof(int) * cities * cities);
        } else if(strncmp("EUC_2D", buf,6) == 0){
            type = EUCLIDIAN;
        } else if(strncmp("GEO:", buf,3) == 0){
            type = GEO;
            printf("Does not support GEO currently \n.");
            exit(1);
        } else if(strncmp("EXPLICIT", buf,8) == 0){
            type = MATRIX;
        } else if(strncmp("NODE_COORD_SECTION",buf,18) == 0){
            X_positions = malloc(sizeof(float) * cities);
            Y_positions = malloc(sizeof(float) * cities);
            int i = 0;
            while(fscanf(source,"%d %f %f \n", &read3, &read1, &read2)){
                X_positions[i] = read1;
                Y_positions[i] = read2;
                i++;
            }
        }else if (strncmp("EDGE_WEIGHT_SECTION",buf,19) == 0){
            for (int i = 0; i < cities; i++){
                for (int j = 0; j < cities; j++) {
                    fscanf(source, "%d", &read3);
                    distM[i*cities + j] = read3;
                }
            }
        }else if (strncmp("EOF", buf, 3) == 0) {
            break;
        }
        
    }
    if (type != MATRIX) {
        for (int i = 0; i < cities; i++){
        printf("X , Y  pos value is: %f , %f \n", X_positions[i], Y_positions[i] );
        }
        create_dist_array(distM, X_positions, Y_positions, type, cities);
        free(X_positions);  free(Y_positions);
    }
    memcpy(save_array,distM,sizeof(int) * cities * cities);
    free(distM);
    return cities;
} 