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
#define man "./environment.x Galaxyn.csv env.data"
// Define the name for the halo and galaxy files
//define Gdata "Galaxy3.csv"
// Simulation constants
#define Hubb 100.0*0.001      // Hubble constant (km/kpc)
#define lh 0.7           // Little h (confirm)
#define L 106500.0          // The size of the box in the simulation in kpc
#define L_h 75000.0         // L/lh
// Method parameters
#define nn 3             // Nearest neighbor for environment definition


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
void get_Env(int j);
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
float** Gp;
float** Gv;

// Vectors of masses of gas and dark matter
float* Ggas;
float* Gdm;

// File to read (Galaxies)
char* Gdata;

// The environment array
float* env;

// The point of view
float* p;

// The number of halos and galaxies
unsigned int nG = 0;

// The grid of indexes for galaxies and halos
List* gridG;


//------------------------------------------------------------------------------
//   Main
//------------------------------------------------------------------------------

int main(int argc, char **argv){

  // The name of the file to write environments
  char* nameW = argv[2];
  Gdata = argv[1];

  // Allocates memory
  allocate_All();

  // Reads the files Halo and Galaxy
  read_File();

  // Gets the environment for all halos
  time_t start = time(NULL);
  printf("Getting Environments ...\n");

  List* actual = gridG->next;

  do{

    printf("Halo %d\n", actual->index+1);
    get_Env( actual->index);

    actual = actual->next;

  }while(actual != 0);

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
  gridG = iniList(0);

  // Galaxies files to read
  FILE *Galaxy;

  // Read the files Halo and Galaxy
  Galaxy = fopen( Gdata, "r" );

  // Read files to obtain number of halos and subhalos
  int test2;
  float tmp;
  do{
    test2 = fscanf(Galaxy, "%f,%f,%f,%f,%f,%f,%f,%f\n", &tmp, &tmp, &tmp, &tmp, &tmp, &tmp ,&tmp, &tmp);
    if(test2 == EOF){
      break;
    }
    nG ++;
  }while( test2!=EOF );

  printf("Number of Galaxies: %d\n",nG);

  // Initialises the point of view
  p = malloc(3*sizeof(float));
  p[0] = 0;
  p[1] = 0;
  p[2] = 0;

  // Environment
  env = malloc(nG*sizeof(float));

  // Allocate galaxies
  int i;
  Gp = malloc(nG*sizeof(float*));
  Gv = malloc(nG*sizeof(float*));
  Ggas = malloc(nG*sizeof(float));
  Gdm = malloc(nG*sizeof(float));

  for( i = 0; i < nG; i++){
    Gp[i] = malloc(3*sizeof(float));
    Gv[i] = malloc(3*sizeof(float));
  }

  // Close the files
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

  // Galaxies files to read
  FILE *Galaxy;

  // Read the files Halo and Galaxy
  Galaxy = fopen( Gdata, "r" );

  // Temporal variables to read
  int test2;
  int numG = 0;
  float x,y,z,vx,vy,vz,mg,mdm;

  do{

    test2 = fscanf(Galaxy, "%f,%f,%f,%f,%f,%f,%f,%f\n", &x, &y, &z, &vx, &vy, &vz, &mg, &mdm );
    if(test2 == EOF){
      break;
    }
    Gp[numG][0] = x/lh;
    Gp[numG][1] = y/lh;
    Gp[numG][2] = z/lh;
    Gv[numG][0] = vx;
    Gv[numG][1] = vy;
    Gv[numG][2] = vz;
    Ggas[numG] = mg;
    Gdm[numG] = mdm;
    //printf("%d__%f\n", numG, x/lh);
    add(gridG,numG);
    numG ++;

  }while( test2!=EOF );
  printf("Number of Galaxies read: %d\n",numG);

  // Close the files
  fclose(Galaxy);

  printf("Time elapsed: %f\n", (float)(time(NULL) - start));

}

/*
 * Given an index, gets the environment of the corresponding halo and saves it to env
 */
void get_Env(int j){

   // Saves the environment of the candidate galaxies
   int num = 0;
   float* Genv = malloc((nG*0.25)*sizeof(float));

   // Indexes for the loops
   int i,l;

   // Actual node of the list
   //printf("OK1.1\n");

   List* actual = gridG->next;

   //printf("OK1.2_____%d\n",(gridG[ix][iy][iz])->next != 0);

   // Loop over the nodes
   do{

     // index i kept by the node
     //printf("G_%d\n",actual->index);
     i = actual->index;
     if( i!= j ){

       //printf("OK2\n");

       // Keeps te velocity difference between the halo and the galaxy
       float* deltaV = malloc(3*sizeof(float));
       deltaV[0] = -(Gv[j][0] - Gv[i][0]) - Hubb*( deltaPBC(Gp[j][0] , Gp[i][0] ) );
       deltaV[1] = -(Gv[j][1] - Gv[i][1]) - Hubb*( deltaPBC(Gp[j][1] , Gp[i][1] ) );
       deltaV[2] = -(Gv[j][2] - Gv[i][2]) - Hubb*( deltaPBC(Gp[j][2] , Gp[i][2] ) );

       // Calculates the vector from the point p to the halo (and its norm)
       float* ph = malloc(3*sizeof(float));
       ph[0] = deltaPBC(Gp[j][0] , p[0]);
       ph[1] = deltaPBC(Gp[j][1] , p[1]);
       ph[2] = deltaPBC(Gp[j][2] , p[2]);
       float phnorm = sqrt( pow(ph[0],2) + pow(ph[1],2) + pow(ph[2],2) );

       // Calculates the relative radial velocity
       float rrv = abs(( ph[0]*deltaV[0] + ph[1]*deltaV[1] + ph[2]*deltaV[2] )/phnorm );
       free(deltaV);

       // If galaxy fulfill the condition, its environment is saved
       if( rrv - 500 <= 0 ){

         // Calculates the vector from the point p to the galaxy (from p to halo + form halo to galaxy)
         float* pg = malloc(3*sizeof(float));
         pg[0] = ph[0] + deltaPBC(Gp[i][0] , Gp[j][0] ); // checked
         pg[1] = ph[1] + deltaPBC(Gp[i][1] , Gp[j][1] );
         pg[2] = ph[2] + deltaPBC(Gp[i][2] , Gp[j][2] );
         float pgnorm = sqrt( pow(pg[0],2) + pow(pg[1],2) + pow(pg[2],2) );

         // Saves the distance^2 to the galaxy projected in the sky
         float x =  ( ph[0]*pg[0] + ph[1]*pg[1] +
         ph[2]*pg[2] )/( phnorm*pgnorm );
         Genv[num] = pow( phnorm/x, 2 )*( 1 - pow(x,2) );
         //printf("%f\n",Genv[num]);
         num ++;
         free(pg);
       }
       free(ph);
     }
     // Step
     actual = actual->next;

   }while(actual != 0);


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
       Genv[maxindx] = 999999999999;
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
  for( i = 0; i < nG; i++){
    fprintf(toWrite, "%f,%f,%f\n", env[i],Ggas[i],Gdm[i]);
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
