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
#include "../../../src/main.h"

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
double **cumul_polar;

int arraysize_x,arraysize_y;
int arraysize_theta,arraysize_radial;

int once = 1;

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


void reset_autocorrelation(){

  int i,j;

  for(i=0;i<arraysize_x;i++){
    for(j=0;j<arraysize_y;j++){
      autocorrelation[i][j] = 0.;
    }}

  for(i=0;i<arraysize_theta;i++){
    for(j=0;j<arraysize_radial;j++){
      polar[i][j] = 0.;
    }}

}

void compute_autocorrelation(){

  int i,j;
  double scalex = 0.8; double scaley = 0.8;

  arraysize_x = (int)(boxsize_x*scalex / resolution_x) + 2;
  arraysize_y = (int)(boxsize_y*scaley / resolution_y) + 2;

  arraysize_theta = (int)(2.*M_PI/resolution_theta) + 2;
  arraysize_radial = (int)((boxsize_y+boxsize_x)/resolution_radial) + 2;
  
  /* printf("arraysize_theta:%d, resolution_theta:%lf\n", arraysize_theta,resolution_theta); */
  /* printf("arraysize_radial:%d, resolution_radial:%lf\n", arraysize_radial,resolution_radial); exit(1); */

  if(once){
    /* allocate memory for autocorrelation array. */
    autocorrelation = (double **)calloc(arraysize_x,sizeof(double*));
    for(i=0;i<arraysize_y;i++){
      autocorrelation[i] = (double *)calloc(arraysize_x,sizeof(double));
    }

    /* allocate memory for polar autocorrelation array. */
    polar = (double **)calloc(arraysize_theta,sizeof(double*));
    for(i=0;i<arraysize_theta;i++){
      polar[i] = (double *)calloc(arraysize_radial,sizeof(double));
    }
    reset_autocorrelation();

    /* allocate memory for cumulative polar autocorrelation array. */
    cumul_polar = (double **)calloc(arraysize_theta,sizeof(double*));
    for(int i=0;i<arraysize_theta;i++){
      cumul_polar[i] = (double *)calloc(arraysize_radial,sizeof(double));
    }
  }

  double resolution_area = resolution_x*resolution_y;

  for(i=0;i<20;i++){
    up(1);printf("i=%d\n", i);
    for(j=0;j<N;j++){

      double _x,_y; double _theta,_radial;
      struct particle p_proj = pair(particles[i],particles[j]);

      if(p_proj.x<-boxsize_x/2.) p_proj.x += boxsize_x;    if(p_proj.y<-boxsize_y/2.) p_proj.y += boxsize_y;
      if(p_proj.x> boxsize_x/2.) p_proj.x -= boxsize_x;    if(p_proj.y> boxsize_y/2.) p_proj.y -= boxsize_y;

      /* calculate autocorrelation in cartesian coordinate. */
      _x	= floor((p_proj.x+boxsize_x/2.*scalex)/resolution_x);
      _y	= floor((p_proj.y+boxsize_y/2.*scaley)/resolution_y);

      /* printf("_x:%lf _y:%lf\n", _x,_y); */
      if(0<=(int)_x && (int)_x <arraysize_x && 0<=(int)_y && (int)_y<arraysize_y){
	/* autocorrelation[(int)_x][(int)_y] += p_proj.m/resolution_area; */
	autocorrelation[(int)_x][(int)_y] ++;
      }else{
	/* printf("particle %d (%lf,%lf) out of range(0:%d)(0:%d).\nexit from %s %d %s.\n", */
	/*        i,_x,_y,arraysize_x,arraysize_y, */
	/*        __FILE__,__LINE__,__func__);  exit(1); */
	/* printf("particle %d out of range. exits at %s %d %s.\n",i,__FILE__,__LINE__,__func__);  exit(1); */
      }

      /* continue; 		/\* @TODO boxsize issue for theta & radial array *\/ */

      double scale_radial = 10.;

      /* calculate autocorrelation in polar coordinate. */
      _theta	= floor((atan2(p_proj.y,p_proj.x)+M_PI)/resolution_theta);
      _radial	= floor(sqrt(p_proj.x*p_proj.x+p_proj.y*p_proj.y)*scale_radial/resolution_radial);

      /* printf("position:%lf %lf _theta,_radial:%lf %lf\n",p_proj.x,p_proj.y,_theta,_radial); */

      if(0<=(int)_theta && (int)_theta<arraysize_theta && 0<=(int)_radial && (int)_radial<arraysize_radial){
	/* polar[(int)_theta][(int)_radial] += p_proj.m/(_radial*resolution_radial*resolution_theta); */
	polar[(int)_theta][(int)_radial] += 1./_radial;
      }else{
	/* printf("particle %d :%lf %lf %lf\n", i,particles[i].x,particles[i].y,particles[i].z); */
	/* printf("particle %d (%lf,%lf) out of range(0:%d)(0:%d).\nexit from %s %d %s.\n", */
	/*        i,_theta,_radial,arraysize_theta,arraysize_radial, */
	/*        __FILE__,__LINE__,__func__);  exit(1); */
      }

    } /* j */
  }   /* i */

}


void compute_cumulative_polar();

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
      fprintf(of, "%lf %lf %lf\n", (i-.5*arraysize_x)*resolution_x,(j-.5*arraysize_y)*resolution_y,autocorrelation[i][j]);      
    }
    fprintf(of, "\n");
  }

  fclose(of);

  /* write polar autocorrelation. */
  sprintf(o2, "%s/%010.2f[orb].polarauorrelation_tr",o,input_interval/(2.*M_PI/OMEGA)*serialnumber);

  if((of=fopen(o2, "w"))==NULL){
    printf("cannot open file %s from function %s.\n", o2,__func__); exit(1);
  }

  for(int i2=arraysize_theta*1/4;i2<arraysize_theta*5/4;i2++){
    i = i2; 
    if(i2>arraysize_theta-1) i -= arraysize_theta;
    for(int j=0;j<arraysize_radial;j++){
      double o1 = i*resolution_theta+M_PI/2.;
      if(o1>2.*M_PI) o1 -= 2.*M_PI;
      if(o1>M_PI) o1 -= 2.*M_PI;

      fprintf(of, "%lf %lf %lf\n",o1*180./M_PI,j*resolution_radial,polar[i][j]);
    }
    fprintf(of, "\n");
  }

  fclose(of);

  compute_cumulative_polar();

  sprintf(o2, "%s/%010.2f[orb].cumulative_tr",o,input_interval/(2.*M_PI/OMEGA)*serialnumber);
  if((of=fopen(o2, "w"))==NULL){
    printf("cannot open file %s from function %s.\n", o2,__func__); exit(1);
  }

  for(int i2=arraysize_theta*1/4;i2<arraysize_theta*5/4;i2++){
    i = i2; 
    if(i2>arraysize_theta-1) i -= arraysize_theta;
    for(int j=0;j<arraysize_radial;j++){
      double o1 = i*resolution_theta+M_PI/2.;
      if(o1>2.*M_PI) o1 -= 2.*M_PI;
      if(o1>M_PI) o1 -= 2.*M_PI;

      fprintf(of, "%lf %lf %lf\n",o1*180./M_PI,j*resolution_radial,cumul_polar[i][j]);
    }
    fprintf(of, "\n");
  }

  fclose(of);
}

void compute_cumulative_polar(){

  int i,j;

  for(i=0;i<arraysize_theta;i++){
    for(j=0;j<arraysize_radial;j++){
      cumul_polar[i][j] = 0.;
    }}

  for(i=0;i<arraysize_theta;i++){
    for(j=0;j<arraysize_radial;j++){
      for(int k=1;k<j;k++){
	cumul_polar[i][j] += polar[i][j];
      }
    }}

  /* find maximum angle. */
  for(j=0;j<arraysize_radial;j++){
   double max = -1.; int argument = -23;

    for(i=0.5*arraysize_theta;i<arraysize_theta;i++){
      if(cumul_polar[i][j]>max){max=cumul_polar[i][j];argument=i;}
    }

    /* printf("argument:%d\n",argument); */
    cumul_polar[argument][j] = 10000.;
  }

}

int main(int argc, char *argv[]){

  resolution_x = 1.;   resolution_y = 1.; /* set as default */

  resolution_theta = M_PI/180.*0.2;   resolution_radial = 5.; /* set as default */

  read_args(argc,argv);

  printf("data file:%s\n", datafilename);

  read_params(datafilename);

  input_interval = 0.1*2.*M_PI/OMEGA;

  while(1){

    sprintf(datafilename, "%s/snapshots/%010.2lf[orb].binall",
	    directoryname,input_interval/(2.*M_PI/OMEGA)*(double)serialnumber);

    read_params(datafilename);

    printf("particle[%d]:%lf %lf %lf\n",0,particles[0].x,particles[0].y,particles[0].z);

    printf("reading %s resolution:%lf %lf\n",datafilename,resolution_x,resolution_y);

    compute_autocorrelation();

    serialnumber++;

    if(!(serialnumber%20)){
  
      write_autocorrelation();

      reset_autocorrelation();

    }

    if(serialnumber>1000) break;

    /* free(densitymap); */
    once = 0;
  }

  free(autocorrelation); free(polar);

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

  return 0;
}


