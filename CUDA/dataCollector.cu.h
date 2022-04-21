#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXLINE 8192

enum DATA_TYPE {EUCLIDIAN,GEO,MATRIX};

uint32_t euclidain_float_to_int(float x1, float x2, float y1, float y2) {
    float dx = x1 - x2;
    float dy = y1 - y2;
    float e_dist = sqrt(dx * dx + dy * dy) + 0.5;
    return (uint32_t) e_dist;
}

void create_dist_array (uint32_t* dist_array, float* xs, float* ys, uint32_t type, uint32_t cities) {
    if (type == EUCLIDIAN) {
        uint32_t euclidian;
        for(uint32_t i = 0; i < cities; i++){
            for(uint32_t j = 0; j < cities; j++){
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

uint32_t fileToDistM(char* filename, uint32_t* save_array){
    FILE* source = fopen(filename, "r");
    if (source == NULL) {
        printf("Source file not found\n");
        return EXIT_FAILURE;
    }
    printf("After opening file\n");
    char buf[MAXLINE];
    uint32_t cities, read3, type, *distM;
    float read1,read2, *X_positions, *Y_positions;
    while(fscanf(source, "%s", buf)) {
        if (strncmp("DIMENSION", buf,9) == 0 ) {
            fscanf(source, "%d", &cities);
            distM = (uint32_t*) malloc(sizeof(uint32_t) * cities * cities);
        } else if(strncmp("EUC_2D", buf,6) == 0){
            type = EUCLIDIAN;
        } else if(strncmp("GEO:", buf,3) == 0){
            type = GEO;
            printf("Does not support GEO currently \n.");
            exit(1);
        } else if(strncmp("EXPLICIT", buf,8) == 0){
            type = MATRIX;
        } else if(strncmp("NODE_COORD_SECTION",buf,18) == 0){
            X_positions = (float*) malloc(sizeof(float) * cities);
            Y_positions = (float*) malloc(sizeof(float) * cities);
            uint32_t i = 0;
            while(fscanf(source,"%d %f %f \n", &read3, &read1, &read2)){
                X_positions[i] = read1;
                Y_positions[i] = read2;
                i++;
            }
        }else if (strncmp("EDGE_WEIGHT_SECTION",buf,19) == 0){
            for (uint32_t i = 0; i < cities; i++){
                for (uint32_t j = 0; j < cities; j++) {
                    fscanf(source, "%d", &read3);
                    distM[i*cities + j] = read3;
                }
            }
        }else if (strncmp("EOF", buf, 3) == 0) {
            break;
        }
        
    }
    printf("After While\n");
    if (type != MATRIX) {
        /*for (uint32_t i = 0; i < cities; i++){
            printf("X , Y  pos value is: %f , %f \n", X_positions[i], Y_positions[i] );
        }*/
        create_dist_array(distM, X_positions, Y_positions, type, cities);
        printf("After create_dist_array\n");
        free(X_positions);  free(Y_positions);
    }
    printf("After type!=MATRIX\n");
    fclose(source);
    memcpy(save_array,distM,sizeof(uint32_t) * cities * cities);
    free(distM);
    return cities;
} 