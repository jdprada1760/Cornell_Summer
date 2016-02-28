#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define x_size 10
float deltaPBC1(float x, float y);
float deltaPBC2(float x, float y);

int main(int argc, char **argv){
  /*
  int i;
  float x = atof(argv[1]);
  float y = atof(argv[2]);
  printf("dx = %f, %f \n", deltaPBC1(x,y), deltaPBC2(x,y));
  */
  FILE* data = fopen( "data.data", "r" );
  int test,count;
  float x;
  count = 0;
  int c = 0;
  do{

    test = fscanf(data, "%f\n", &x);
    if( x != 0){
      //printf("________%f\n",sqrt(x));
      count++;
    }
    if( x == 9e+12){
      c++;
    }
    //printf("%d_%d_%d_\n",(int)floorf(x/(L/res)),(int)floorf(y/(L/res)),(int)floorf(z/(L/res)));

  }while( test!=EOF );
  printf("FIN__%d__%d\n",count, c);
  return 0;
}

float deltaPBC1(float x, float y){
  float dx = x - y;
  if(dx >   x_size * 0.5){
    dx = dx - x_size;
  }
  else if (dx <= -x_size * 0.5){
    dx = dx + x_size;
  }
  return dx;
}

float deltaPBC2(float x, float y){
  float dx = x - y;
  dx = dx - truncf(dx / x_size) * x_size;
  return dx;
}
