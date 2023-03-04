#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#define eps 1e-3

typedef struct{ //struct with all data for a particle
    double pos_x;
    double pos_y;
    double mass;
    double vel_x;
    double vel_y;
    double brightness;
}planet_data_t;

void update_particles(planet_data_t *arr[], planet_data_t *arrtemp[], const int N, const int start, const int end, const double dt){
    //variables needed
    double acc_x;
    double acc_y;
    double rx;
    double ry;
    double r;
    double denom;
    const double G = 100/ (double)N; //typecast to get correct result
    double veln_x;
    double veln_y;
    double posn_x;
    double posn_y;

    for(int i = start; i < end; i++){ //loop over each particle, calculating acceleration and new vel & pos
        acc_x = 0;
        acc_y = 0;

        for(int j = 0; j<i; j++){ //second version of inner for-loop, split in two loops
            rx = arr[i]->pos_x - arr[j]->pos_x;
            ry = arr[i]->pos_y - arr[j]->pos_y;
            r = sqrt(rx*rx + ry*ry);

            denom = 1/((r+eps)*(r+eps)*(r+eps)); //computing "denominator" (actually 1/denominator) once instead of twice below

            acc_x +=  denom*arr[j]->mass*rx;
            acc_y +=  denom*arr[j]->mass*ry;
        }

        for(int j=i+1; j<N; j++){
            rx = arr[i]->pos_x - arr[j]->pos_x;
            ry = arr[i]->pos_y - arr[j]->pos_y;
            r = sqrt(rx*rx + ry*ry);

            denom = 1/((r+eps)*(r+eps)*(r+eps)); 

            acc_x +=  denom*arr[j]->mass*rx;
            acc_y +=  denom*arr[j]->mass*ry;
        }
        
        acc_x *= -G;
        acc_y *= -G;

        veln_x = arr[i]->vel_x + dt*acc_x;
        veln_y = arr[i]->vel_y + dt*acc_y;

        posn_x = arr[i]->pos_x + dt*veln_x;
        posn_y = arr[i]->pos_y + dt*veln_y;

        arrtemp[i]->pos_x = posn_x; //store following timesteps values temporarily
        arrtemp[i]->pos_y = posn_y;
        arrtemp[i]->vel_x = veln_x;
        arrtemp[i]->vel_y = veln_y;
    }
}

// MAKING FUNCTION TO WRITE ARRTEMP TO ARR -----------------------
void update_arrrays(planet_data_t *arr[], planet_data_t *arrtemp[], const int N){
    //updating each particle's vel & pos after going through all particles on a timestep
    for(int i=0; i<N; i++){ 
        arr[i]->pos_x = arrtemp[i]->pos_x;
        arr[i]->pos_y = arrtemp[i]->pos_y;
        arr[i]->vel_x = arrtemp[i]->vel_x;
        arr[i]->vel_y = arrtemp[i]->vel_y;
    }
}
// ---------------------------------------------------------------


int main(int argc, char *argv[])  {

    if (argc != 7){
        printf("Wrong number of arguments! \n");
        printf("Expected: './galsim  N  filename  nsteps  delta_t  graphics(1/0)  n_threads' \n");
        return 0;
    }

    //read command line arguments into variables
    int N = atoi(argv[1]);
    char file[100];
    strcpy(file, argv[2]);
    int nsteps = atoi(argv[3]);
    double dt = atof(argv[4]);
    int graphics = atoi(argv[5]);
    int N_THREADS = atoi(argv[6]);
    //printf("%d, %s, %d, %lf, %d\n", N, file, nsteps, dt, graphics, n_threads); //test printing


    //Opening file
    char path[100];
    //strcpy(path, "input_data/");
    //strcat(path, file);
    strcpy(path, file);
    //printf("%s\n", path);

    FILE *stream;

    stream = fopen(path, "rb");

    //reading file into structs in an array

//arrays to store pointers to the data and the following timestep's values
    planet_data_t** arr = (planet_data_t**)malloc(N*sizeof(planet_data_t*));
    planet_data_t** arrtemp = (planet_data_t**)malloc(N*sizeof(planet_data_t*));

    for(int i = 0; i<N; i++){
        arr[i]=(planet_data_t *)malloc(sizeof(planet_data_t));
        arrtemp[i]=(planet_data_t *)malloc(sizeof(planet_data_t));
    }
    
    if(!stream){
        printf("Couldn't read file, exiting!\n");
        return 0;
    }

    for(int i =0; i<N; i++){
        
        fread(&(arr[i]->pos_x), sizeof(double), 1, stream);
        fread(&(arr[i]->pos_y), sizeof(double), 1, stream);
        fread(&(arr[i]->mass), sizeof(double), 1, stream);
        fread(&(arr[i]->vel_x), sizeof(double), 1, stream);
        fread(&(arr[i]->vel_y), sizeof(double), 1, stream);
        fread(&(arr[i]->brightness), sizeof(double), 1, stream);
    
    }

    fclose(stream);

    // THREAD STUFF
    int* THREAD_ID[N_THREADS];
    pthread_t THREADS[N_THREADS];

    int start;
    int end;

    // LOCAL FUNCTION FOR PTHREAD_CREATE
    void* thread_func(void* THREAD_ID){
        int ID = *(int*)(THREAD_ID);    // casting the thread ID into an integer
        start = (ID) * (N / N_THREADS); // calculating start index
        end = start + (N/N_THREADS);    // calculating end index
        update_particles(arr, arrtemp, N, start, end, dt);
       // printf("\nWORK DONE BY THREAD %d\n", ID);
    }


    for (int i = 0; i < nsteps; i++){

        // START WORK
        for (int j = 0; j < N_THREADS; j++){
            THREAD_ID[j] = (int*)malloc(sizeof(int));
            *THREAD_ID[j] = j;
            pthread_create(&THREADS[j], NULL, thread_func, THREAD_ID[j]);
        }

        // WAIT FOR THREADS
        for (int j = 0; j < N_THREADS; j++){
            pthread_join(THREADS[j], NULL);
        }

        // UPDATE ARRAYS
        update_arrrays(arr, arrtemp, N);

        //printf("TIME STEP DONE: %d\n ", i);
    }
  
    //writing to result-file
    char *filename = "result.gal";
    FILE *fp;
    fp = fopen(filename, "wb");

    if(!fp){
        printf("Error writing file!\n");
        return 0;
    }

    

    for (int i = 0; i < N; i++){
        fwrite(&(arr[i]->pos_x), sizeof(double), 1, fp);
        fwrite(&(arr[i]->pos_y), sizeof(double), 1, fp);
        fwrite(&(arr[i]->mass), sizeof(double), 1, fp);
        fwrite(&(arr[i]->vel_x), sizeof(double), 1, fp);
        fwrite(&(arr[i]->vel_y), sizeof(double), 1, fp);
        fwrite(&(arr[i]->brightness), sizeof(double), 1, fp);
    }

    fclose(fp);

    for(int i = 0; i<N; i++){
        free(arr[i]);
        free(arrtemp[i]);
    }

    free(arr);
    free(arrtemp);

    pthread_exit(NULL);
    
    return 0;

}