#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXLINE 8192
#define MAXCITIES 10000
#define PI 3.141592
#define RRR 6378.388
enum DATA_TYPE {EUCLIDIAN,GEO,MATRIX};
enum MATRIX_TYPE {UPPER_ROW, UPPER_ROW_DIAG, UPPER_COL, UPPER_COL_DIAG};

double lat_calc(float x){
    int deg = (int) (x + 0.5);
    double min = x - deg;
    double lat = PI * (deg + 5.0 * min / 3.0) / 180.0;
    return lat;
}
int geo_float_to_int(float x1, float x2, float y1, float y2){
    double lat1 = lat_calc(x1);
    double long1 = lat_calc(y1);
    double lat2 = lat_calc(x2);
    double long2 = lat_calc(y2);
    double q1 = cos(long1 - long2);
    double q2 = cos(lat1 - lat2);
    double q3 = cos(lat1 + lat2);
    int dist = (int) (RRR * acos(0.5 * ((1.0+q1)*q2 -(1.0-q1)*q3)) + 1.0);
    return dist;
}

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
    }
    else if (type == GEO){
        int geo;
        for(int i = 0; i < cities; i++){
            for(int j = 0; j < cities; j++){
                geo = geo_float_to_int(xs[i], xs[j], ys[i], ys[j]);
                dist_array[i * cities + j] = geo;
                dist_array[j * cities + i] = geo;
            }
        }      
    }
    else{
        printf("This type is not supported \n.");
        exit(1);
    }
}

int dertimine_matrix_type(char* descript){
    int return_val;
    if (strncmp("UPPER_ROW",descript,9) == 0){
        return_val = UPPER_ROW;
    }
    else if ((strncmp("UPPER_DIAG_ROW",descript,14) == 0)){
        return_val = UPPER_ROW_DIAG;
    }
    else if ((strncmp("UPPER_COL",descript,9) == 0)){
        return_val = UPPER_COL;
    }
    else if ((strncmp("UPPER_DIAG_COL",descript,14) == 0)){
        return_val = UPPER_COL_DIAG;
    }
    else {
        printf("error on match \n");
        exit(1);
    }
    return return_val;
}

int fileToDistM(char* filename, int* save_array){
    FILE* source = fopen(filename, "r");
    if (source == NULL) {
        printf("Source file not found\n");
        return EXIT_FAILURE;
    }
    char buf[MAXLINE];
    int cities, read3, type, *distM, m_type, ch;
    float read1,read2, *X_positions, *Y_positions;
    while(fscanf(source, "%s", buf)) {
        if (strncmp("DIMENSION", buf,9) == 0 ) {
            fscanf(source, "%d", &cities);
            distM = malloc(sizeof(int) * cities * cities);
        } else if(strncmp("EUC_2D", buf,6) == 0){
            type = EUCLIDIAN;
        } else if(strncmp("GEO:", buf,3) == 0){
            type = GEO;
        } else if(strncmp("EXPLICIT", buf,8) == 0){
            type = MATRIX;
        } else if(strncmp("EDGE_WEIGHT_FORMAT", buf,18)== 0){
            ch = fgetc(source);
            while (ch != ':')
            {
                ch = fgetc(source);
            }
            fscanf(source,"%s", buf);
            m_type = dertimine_matrix_type(buf);
        }
         else if(strncmp("NODE_COORD_SECTION",buf,18) == 0){
            X_positions = malloc(sizeof(float) * cities);
            Y_positions = malloc(sizeof(float) * cities);
            int i = 0;
            while(fscanf(source,"%d %f %f \n", &read3, &read1, &read2)){
                X_positions[i] = read1;
                Y_positions[i] = read2;
                i++;
            }
        }else if (strncmp("EDGE_WEIGHT_SECTION",buf,19) == 0){
            if (type = MATRIX){
                for (int i = 0; i < cities; i++){
                    for (int j = 0; j < cities; j++) {
                        fscanf(source, "%d", &read3);
                        distM[i*cities + j] = read3;
                    }
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

void matrixToTxt(int*distM, int cities, char* output_file){
    FILE* source = fopen(output_file, "w+");
    if (source == NULL) {
        printf("Source file not found\n");
        exit(1);
    }
    printf("gets here \n");
    fprintf(source,"%di32\n", cities);
    fprintf(source,"[");
    for (int i = 0; i < cities; i++) {
        for (int j = 0; j < cities; j++){
            if(i == cities-1 && j == cities-1){
                fprintf(source,"%di32]", distM[(cities-1)*cities + (cities-1)]);
            }else{
                fprintf(source, "%di32, ",distM[i*cities+j]);
            }
            
        }
    }
}


int main(int argc, char *argv[]) {
    if (argc != 3){
        printf("Incorrect number of arguments. Expected two argument (address to a datafile and address to return txt file) \n");
        exit(1);
    }
    int *distM = malloc(sizeof(int) * MAXCITIES * MAXCITIES);
    int cities = fileToDistM(argv[1],distM);
    distM = realloc(distM,sizeof(int) * cities * cities);
    for (int i = 0; i < cities; i++) {
        for (int j = 0; j < cities; j++){
            printf("%d ", distM[i* cities + j]);
            }
        printf("\n"); 
    }
    matrixToTxt(distM, cities, argv[2]);
    free(distM);
    return 0;
}

