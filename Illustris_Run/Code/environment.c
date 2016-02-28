/*
 * This program calculates the environment for each halo in Halo.csv
 * The environment is defined as the area projected for nth nearest galaxy
 * conditions for neighbouring galaxies:
 *          -rmag < -19 (already filtered in Galaxy.csv)
 *          -radial velocity relative to the halo < 500km/s
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//------------------------------------------------------------------------------
//   Definitions
//------------------------------------------------------------------------------

// env.csv : the name of the file to write the environments
#define man "./environment.x env.csv"
// Define the name for the halo and galaxy files
#define Hdata "Halo.csv"
#define Gdata "Galaxy.csv"
// Simulation constants
#define Hubb 100.0       // Hubble constant
#define lh 0.7           // Little h (confirm)
#define L 106500.0          // The size of the box in the simulation in kpc
#define L_h 75000.0         // L/lh
// Method parameters
#define nn 3             // Nearest neighbor for environment definition
#define ddd 1            // Defines the range in the grid to count neighbor candidates
#define res 10.0         // The resolution of the 3d spatial grid




//------------------------------------------------------------------------------
//  Structs
//------------------------------------------------------------------------------

/*
 * Chain list
 */
struct List{
  // El hijo proximo
  struct List* next;
  // El ultimo hijo
  struct List* last;
  // El indice que guarda el tetraedro correspondiente
  //(Permite sacar las matrices, los vecref y los volumenes)
  int index;
  // Guarda la longitud de la lista (solo la guarda el primero, o ahi vemos)
  //unsigned int n;
};
typedef struct List List;


//------------------------------------------------------------------------------
//   Methods declaration
//------------------------------------------------------------------------------
void allocate_All();
void read_File();
void get_Env(int j, int ix, int iy, int iz);
void write_Env( char* nameW );
float deltaPBC(float x, float y);
int dPBC(int x, int y);
int rPBC(int x);
List* iniList( int index );
void add(List* lista, unsigned int index);

//------------------------------------------------------------------------------
//   General Variables
//------------------------------------------------------------------------------

// Vectors and matrices for halos and galaxies positions, velocities and masses
float** Hp;
float** Hv;
float** Gp;
float** Gv;

// The environment array
float* env;

// The point of view
float* p;

// The number of halos and galaxies
unsigned int nH, nG = 0;

// The grid of indexes for galaxies and halos
List**** gridH;
List**** gridG;


//------------------------------------------------------------------------------
//   Main
//------------------------------------------------------------------------------

int main(int argc, char **argv){

  // The name of the file to write environments
  char* nameW = argv[1];

  // Allocates memory
  allocate_All();

  // Reads the files Halo and Galaxy
  read_File();

  // Gets the environment for all halos
  time_t start = time(NULL);
  printf("Getting Environments ...\n");

/*
  int i,j,k;
  for( i = 0; i < res; i++){
    for( j = 0; j < res; j++){
      for( k = 0; k < res; k++){

        List* actual = (gridH[i][j][k])->next;

        do{

          //printf("Halo %d\n", i+1);
          get_Env( actual->index, i, j, k );

          actual = actual->next;

        }while(actual != 0);

      }
    }
  }
*/

  int i,j,k;
  for( i = 0; i < 1; i++){
    for( j = 0; j < 1; j++){
      for( k = 1; k < 10; k++){

        List* actual = (gridH[i][j][k])->next;

        int count = 0;
        printf("%d_%d_%d_\n",i,j,k);

        do{

          //printf("Halo %d\n", i+1);
          get_Env( actual->index, i, j, k );
          //printf("HALO_%d__________________________\n",actual->index);
          actual = actual->next;
          count++;

        }while(actual!= 0);

      }
    }
  }

  printf("Time elapsed: %f\n", (float)(time(NULL) - start));

  // Writes the file of environments
  write_Env(nameW);
  return 0;
}



//------------------------------------------------------------------------------
//   Methods
//------------------------------------------------------------------------------


/*
 * Allocate memory for all arrays
 */
void allocate_All(){
  time_t start = time(NULL);
  printf("Allocating...\n");
  // Allocates the grid
  int l,m,n;
  gridH = malloc(res*sizeof(List***));
  gridG = malloc(res*sizeof(List***));
  for( l = 0; l < res; l++){
    gridH[l] = malloc(res*sizeof(List**));
    gridG[l] = malloc(res*sizeof(List**));
    for( m = 0; m < res; m++){
      gridH[l][m] = malloc(res*sizeof(List*));
      gridG[l][m] = malloc(res*sizeof(List*));
      for( n = 0; n < res; n++){
        // Initialises the list with the index 0
        gridH[l][m][n] = iniList(0);
        gridG[l][m][n] = iniList(0);
      }
    }
  }


  // Halos and galaxies files to read
  FILE *Halo;
  FILE *Galaxy;

  // Read the files Halo and Galaxy
  Halo = fopen( Hdata, "r" );
  Galaxy = fopen( Gdata, "r" );

  // Read files to obtain number of halos and subhalos
  int test, test2;
  float tmp;

  do{
    test = fscanf(Halo, "%f,%f,%f,%f,%f,%f\n", &tmp, &tmp, &tmp, &tmp, &tmp, &tmp );
    nH ++;
  }while( test!=EOF );

  printf("Number of Halos: %d\n",nH);

  do{
    test2 = fscanf(Galaxy, "%f,%f,%f,%f,%f,%f\n", &tmp, &tmp, &tmp, &tmp, &tmp, &tmp );
    nG ++;
  }while( test2!=EOF );

  printf("Number of Galaxies: %d\n",nG);

  // Initialises the point of view
  p = malloc(3*sizeof(float));
  p[0] = 0;
  p[1] = 0;
  p[2] = 0;

  // Environment
  env = malloc(nH*sizeof(float));

  // Allocate Halos
  int i;
  Hp = malloc(nH*sizeof(float*));
  Hv = malloc(nH*sizeof(float*));

  for( i = 0; i < nH; i++){
    Hp[i] = malloc(3*sizeof(float));
    Hv[i] = malloc(3*sizeof(float));
  }

  // Allocate galaxies
  Gp = malloc(nG*sizeof(float*));
  Gv = malloc(nG*sizeof(float*));

  for( i = 0; i < nG; i++){
    Gp[i] = malloc(3*sizeof(float));
    Gv[i] = malloc(3*sizeof(float));
  }

  // Close the files
  fclose(Halo);
  fclose(Galaxy);

  printf("Time elapsed: %f\n", (float)(time(NULL) - start));
}

/*
 * Read files with all the info of Halos and Subhalos
 * HaloFormat : x,y,z,vx,vy,vz,mass
 * GalaxyFormat: x,y,z,vx,vy,vz,mass
 * galaxies already filtered by magnitude
 */
void read_File(){

  time_t start = time(NULL);
  printf("Reading data...\n");

  // Halos and galaxies files to read
  FILE *Halo;
  FILE *Galaxy;

  // Read the files Halo and Galaxy
  Halo = fopen( Hdata, "r" );
  Galaxy = fopen( Gdata, "r" );

  // Temporal variables to read
  int test,test2;
  int numH = 0;
  int numG = 0;
  float x,y,z,vx,vy,vz;

  // Reading
  do{

    test = fscanf(Halo, "%f,%f,%f,%f,%f,%f\n", &x, &y, &z, &vx, &vy, &vz );
    Hp[numH][0] = x/lh;
    Hp[numH][1] = y/lh;
    Hp[numH][2] = z/lh;
    Hv[numH][0] = vx;
    Hv[numH][1] = vy;
    Hv[numH][2] = vz;
    add(gridH[(int)floorf(x/(L_h/res))][(int)floorf(y/(L_h/res))][(int)floorf(z/(L_h/res))],numH);
    //printf("%d_%d_%d_\n",(int)floorf(x/(L/res)),(int)floorf(y/(L/res)),(int)floorf(z/(L/res)));
    numH ++;

  }while( test!=EOF );

  printf("Number of Halos read: %d\n",numH);

  int i,j,k;
  for( i = 0; i < res; i++){
    for( j = 0; j < res; j++){
      for( k = 0; k < res; k++){
        if((gridH[i][j][k])->next == 0){
          printf("HVoid__%d_%d_%d____\n",i,j,k);
        }
      }
    }
  }

  do{

    test2 = fscanf(Galaxy, "%f,%f,%f,%f,%f,%f\n", &x, &y, &z, &vx, &vy, &vz );
    Gp[numG][0] = x/lh;
    Gp[numG][1] = y/lh;
    Gp[numG][2] = z/lh;
    Gv[numG][0] = vx;
    Gv[numG][1] = vy;
    Gv[numG][2] = vz;
    add(gridG[(int)floorf(x/(L_h/res))][(int)floorf(y/(L_h/res))][(int)floorf(z/(L_h/res))],numG);
    numG ++;

  }while( test2!=EOF );
  printf("Number of Galaxies read: %d\n",numG);

  // Check how many grid cubes are out of galaxies
  for( i = 0; i < res; i++){
    for( j = 0; j < res; j++){
      for( k = 0; k < res; k++){
        if((gridG[i][j][k])->next == 0){
          printf("GVoid__%d_%d_%d____\n",i,j,k);
        }
      }
    }
  }


  // Close the files
  fclose(Halo);
  fclose(Galaxy);

  printf("Time elapsed: %f\n", (float)(time(NULL) - start));

}

/*
 * Given an index, gets the environment of the corresponding halo and saves it to env
 */
void get_Env(int j, int hx, int hy, int hz){

   // Saves the environment of the candidate galaxies
   int num = 0;
   float* Genv = malloc((nG*0.25)*sizeof(float));

   // Indexes for the loops
   int i,l;
   int ix,iy,iz,cx,cy,cz;
   // Loops in each dimension
   ix = rPBC(hx-ddd);
   for(cx = 0; cx <= 2*ddd; cx++){
     iy = rPBC(hy-ddd);
     for(cy = 0; cy <= 2*ddd; cy++){
       iz = rPBC(hz-ddd);
       for(cz = 0; cz <= 2*ddd; cz++){


         // Actual node of the list
         //printf("OK1.1\n");

         List* actual = (gridG[ix][iy][iz])->next;

         //printf("OK1.2_____%d\n",(gridG[ix][iy][iz])->next != 0);

         // Loop over the nodes
         do{

           // index i kept by the node
           //printf("G_%d\n",actual->index);
           i = actual->index;

           //printf("OK2\n");

           // Keeps te velocity difference between the halo and the galaxy
           float* deltaV = malloc(3*sizeof(float));
           deltaV[0] = (Hv[j][0] - Gv[i][0]) + Hubb*( deltaPBC(Hp[j][0] , Gp[i][0] ) );
           deltaV[1] = (Hv[j][1] - Gv[i][1]) + Hubb*( deltaPBC(Hp[j][1] , Gp[i][1] ) );
           deltaV[2] = (Hv[j][2] - Gv[i][2]) + Hubb*( deltaPBC(Hp[j][2] , Gp[i][2] ) );

           // Calculates the vector from the point p to the halo (and its norm)
           float* uv = malloc(3*sizeof(float));
           uv[0] = deltaPBC(Hp[j][0] , p[0]);
           uv[1] = deltaPBC(Hp[j][1] , p[1]);
           uv[2] = deltaPBC(Hp[j][2] , p[2]);
           float uvnorm = sqrt( pow(uv[0],2) + pow(uv[1],2) + pow(uv[2],2) );

           // Calculates the relative radial velocity
           float rrv = abs(( uv[0]*deltaV[0] + uv[1]*deltaV[1] + uv[2]*deltaV[2] )/uvnorm );
           free(deltaV);

           //printf("OK3\n");

           // If galaxy fulfill the condition, its environment is saved
           if( rrv - 500 <= 0 ){

             float tnorm = sqrt( pow(deltaPBC(Gp[i][0] , p[0]),2) + pow(deltaPBC(Gp[i][1] , p[1]),2) +
             pow(deltaPBC(Gp[i][2] , p[2]),2) );

             // Saves the distance to the galaxy projected in the sky
             float x =  ( uv[0]*(deltaPBC(Gp[i][0],p[0])) + uv[1]*(deltaPBC(Gp[i][1],p[1])) +
             uv[2]*(deltaPBC(Gp[i][2],p[2])) )/( uvnorm*tnorm );
             Genv[num] = pow( uvnorm/x, 2 )*( 1 - pow(x,2) );

             //printf("%f\n",Genv[num]);
             num ++;

           }

           free(uv);

           // Step
           actual = actual->next;

         }while(actual != 0);

         //printf("_%d_%d_%d_\n",ix,iy,iz);
         iz = rPBC(iz+1);
       }
       iy = rPBC(iy+1);
     }
     ix = rPBC(ix+1);
   }

   //printf("%d___\n",num);
   // Gets the nnth nearest neighbor (nn)
   float maxdist;
   int maxindx;
   for( i = 0; i < nn; i++){
     maxdist = 9e+12;
     maxindx = 0;
     for( l = 0; l < num; l++){
       if( Genv[l] <= maxdist ){
         maxdist = Genv[l];
         maxindx = l;
       }
     }
     // Sets a big number to delete the minimum and find the second one
     if( i != nn-1 ){
       Genv[maxindx] = 9e+12;
     }
   }

   env[j] = maxdist;
   //printf("%f\n",maxdist);
   free(Genv);
}

/*
 * Writes the file of environments
 */

void write_Env( char* nameW ){
  FILE* toWrite = fopen( nameW, "w");
  time_t start = time(NULL);
  printf("Writing Environments...\n");
  unsigned int i;
  for( i = 0; i < nH; i++){
    fprintf(toWrite, "%f\n", env[i]);
  }
  fclose(toWrite);
  printf("Time elapsed: %f\n", (float)(time(NULL) - start));
}

//------------------------------------------------------------------------------
// PBC methods
//------------------------------------------------------------------------------

/*
 * Distance in periodic boundary conditions
 */
float deltaPBC(float x, float y){
  float dx = x - y;
  if(dx >   L * 0.5){
    dx = dx - L;
  }
  else if (dx <= -L * 0.5){
    dx = dx + L;
  }
  return dx;
}

/*
 * Distance in periodic boundary conditions (integer version)
 */
int dPBC(int x, int y){
  int dx = x - y;
  if(dx >  (res*0.5)){
    dx = (int)(dx - res);
  }
  else if (dx <= -(res*0.5)){
    dx = (int)(dx + res);
  }
  return dx;
}

/*
 * Integer position in periodic boundary conditions
 */
int rPBC(int x){
  if(x <  0){
    x = (int)(x + res);
  }
  else if(x >=  res){
    x = (int)(x - res);
  }
  return x;
}


//------------------------------------------------------------------------------
// Struct methods for the list
//------------------------------------------------------------------------------

/*
 * Inicializa una lista vacía dado un índice;
 */
List* iniList( int index ){
  List* li = malloc(sizeof(List));
  li->next = 0;
  li->index = index;
  li->last = li;
  return li;
}

/*
 * Añade un nuevo elemento a la lista
 */
void add(List* lista, unsigned int index){
  // Aumenta el numero de elementos en la lista primera.
  List *li = malloc(sizeof(List));
  li->next = 0;
  //li->n = 1;
  li->index = index;
  li->last = li;
  //lista->n = lista->n +1;
  // Crea el elemento a la lista
  // Añade el nuevo elemento en la cola y cambia la referencia al elemento final en el primer nodo
  (lista->last)->next = li;
  //printf("%d____OK\n",(lista->last)->index);
  lista->last = li;
}
