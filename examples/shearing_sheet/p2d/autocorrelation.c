#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ftw.h>
#include <getopt.h>
#include "particle.h"
#include "../../../src/input_params.h"

char *datafilename, *directoryname;
double boxsize_x,boxsize_y,boxsize_z,G,dt,softening,OMEGA;
double resolution_x,resolution_y;
double resolution_theta,resolution_radial;

int N;
double input_interval;

double input_interval;
int serialnumber = 0;
double **autocorrelation;	/* autocorrelation in the xy space */
double **polar;			/* autocorrelation in the theta-radial space */

int arraysize_x,arraysize_y;
int arraysize_theta,arraysize_radial;

void readfile(char *filename);

void read_args(int argc, char *argv[]);

struct particle pair(struct particle p1, struct particle p2){

  struct particle q;

  q.x = p1.x - p2.x;
  q.y = p1.y - p2.y;
  q.z = p1.z - p2.z;  
  q.vx = p1.vx - p2.vx;
  q.vy = p1.vy - p2.vy;
  q.vz = p1.vz - p2.vz;

  q.r = 0.;
  q.m = 0.;

  return q;
}


void compute_autocorrelation(){

  int i,j;

  arraysize_x = (int)(boxsize_x / resolution_x) + 2;
  arraysize_y = (int)(boxsize_y / resolution_y) + 2;

  /* allocate memory for autocorrelation array. */
  autocorrelation = (double **)calloc(arraysize_x,sizeof(double*));
  for(i=0;i<arraysize_y;i++){
    autocorrelation[i] = (double *)calloc(arraysize_x,sizeof(double));
  }

  /* allocate memory for polar autocorrelation array. */
  polar = (double **)calloc(arraysize_theta,sizeof(double*));
  for(i=0;i<arraysize_radial;i++){
    polar[i] = (double *)calloc(arraysize_theta,sizeof(double));
  }

  double resolution_area = resolution_x*resolution_y;

  for(i=0;i<100;i++){
    for(j=0;j<100;j++){

      double _x,_y; double _theta,_radial;
      struct particle p_proj = pair(particles[i],particles[j]);

      if(p_proj.x<-boxsize_x/2.) p_proj.x += boxsize_x;    if(p_proj.y<-boxsize_y/2.) p_proj.y += boxsize_y;
      if(p_proj.x> boxsize_x/2.) p_proj.x -= boxsize_x;    if(p_proj.y> boxsize_y/2.) p_proj.y -= boxsize_y;

      /* calculate autocorrelation in cartesian coordinate. */
      _x	= floor((p_proj.x+boxsize_x/2.)/resolution_x);
      _y	= floor((p_proj.y+boxsize_y/2.)/resolution_y);

      if(0<=(int)_x && (int)_x <=arraysize_x && 0<=(int)_y && (int)_y<=arraysize_y){
	autocorrelation[(int)_x][(int)_y] += p_proj.m/resolution_area;
      }else{
	printf("particle %d out of range. exits at %s %d %s.\n",i,__FILE__,__LINE__,__func__);  exit(1);
      }

      /* calculate autocorrelation in polar coordinate. */
      _theta	= floor(atan2(p_proj.x,p_proj.y)/resolution_theta);
      _radial	= floor(sqrt(p_proj.x*p_proj.x+p_proj.y*p_proj.y)/resolution_theta);

      if(0<=(int)_theta && (int)_theta<=arraysize_theta && 0<=(int)_radial && (int)_radial<=arraysize_radial){
	polar[(int)_theta][(int)_radial] += p_proj.m/(_radial*resolution_radial*resolution_theta);
      }else{
	printf("particle %d out of range. exits at %s %d %s.\n",i,__FILE__,__LINE__,__func__);  exit(1);
      }

    } /* j */
  }   /* i */

}

void write_autocorrelation(){

  char o[1024],o2[1024];
  /* prepare data directory */
  sprintf(o,"%s/autocorrelation[%dx%d]",directoryname,arraysize_x,arraysize_y);
  mkdir(o, S_IRWXU | S_IRWXG | S_IRWXO);
  sprintf(o2, "%s/%010.2f[orb].autocorrelation_xy",o,input_interval/(2.*M_PI/OMEGA)*serialnumber);

  FILE *of; 
  if((of=fopen(o2, "w"))==NULL){
    printf("cannot open file %s from function %s.\n", o2,__func__); exit(1);
  }

  int i,j;
  for(int i=0;i<arraysize_x;i++){
    for(int j=0;j<arraysize_y;j++){
      fprintf(of, "%lf %lf %lf\n", i*resolution_x,j*resolution_y,autocorrelation[i][j]);      
    }
    fprintf(of, "\n");
  }

  fclose(of);
}


int main(int argc, char *argv[])
{
  resolution_x = 5.;   resolution_y = 5.; /* set as default */
  resolution_theta = 1.;   resolution_radial = 5.; /* set as default */

  read_args(argc,argv);

  printf("data file:%s\n", datafilename);

  read_params(datafilename);

  printf("resolution:theta:%lf %lf", resolution_theta,resolution_radial); exit(1);
  /* @TODO set arraysize_theta,radial */

  compute_autocorrelation();

  /* write_autocorrelation(); */

  sprintf(datafilename, "%s/snapshots/%010.2lf[orb].binall",
	  directoryname,input_interval/(2.*M_PI/OMEGA)*(double)serialnumber);

  printf("reading %s resolution:%lf %lf\n", datafilename, resolution_x,resolution_y);

  /* while((result=getopt(argc, argv, "ab:cd:e:"))!=1){ */
  /*   switch(result){ */
  /*   case 'a': */
  /*     printf("option %c applied with %s\n", result, optarg); */
  /*   case 'b': */
  /*     printf("option %c applied with %s\n", result, optarg); */
  /*   case 'c': */
  /*     printf("option %c applied with %s\n", result, optarg); */
  /*   case 'd': */
  /*     printf("option %c applied with %s\n", result, optarg); */
  /*   case 'e': */
  /*     printf("option %c applied with %s\n", result, optarg); */
  /*     break; */
  /*   case ':':   // no value applied */
  /*   case '?':   // invalid option */
  /*     exit(1); */
  /*   } */
  /* } */

  free(autocorrelation); free(polar);

  return 0;
}


