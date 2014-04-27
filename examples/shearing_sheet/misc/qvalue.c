#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <ftw.h>
#include <getopt.h>
#include <dirent.h>
#include <string.h>
#include "particle.h"
#include "readfile.h"
#include "../../../src/input_params.h"

char *datafilename, *directoryname;
double boxsize_x,boxsize_y,boxsize_z,G,dt,softening,OMEGA;
double resolution_x,resolution_y;
int N;

int once=1;
double input_interval;
int serialnumber = 0;

double t1 = 0.; double tt1 = 0.1;
double mean=0.; double var=0.; double stdd=0.; /* average, variance & standard deviation for vx */

double c_critical;			/* 3.36G*Sigma/kappa. here Omega=kappa. */

void compute_velocity_dispersion(){

  if(once){
    c_critical = 3.36*surfacedensity()*G/OMEGA;
    printf("surfacedensity:%lf c_critical:%lf", surfacedensity(),c_critical);
  }

  int i;
  mean=0.; var=0.; stdd=0.;

  for(i=0;i<N;i++){
    mean += particles[i].vx;    
    var += particles[i].vx*particles[i].vx;    
  }
  mean = mean/(double)N;
  var = var/(double)N;

  stdd = sqrt(var-mean*mean);
}

/**
   write to file
*/
void out_velocity_dispersion(){

  char o[1024],o2[1024];
  /* prepare data directory */
  sprintf(o,"%s/vdisp",directoryname);
  mkdir(o, S_IRWXU | S_IRWXG | S_IRWXO);
  sprintf(o2, "%s/velocitydispersion",o);
  
  FILE *of; 

  char mode[4]; if(once){sprintf(mode,"w");}else{sprintf(mode,"a");}

  if((of=fopen(o2,mode))==NULL){
    printf("cannot open file %s from function %s.\n", o2,__func__); exit(1);
  }

  fprintf(of, "%lf %.9lf %.9lf %.9lf %lf\n",t1,mean,var,stdd,stdd/c_critical);

  fclose(of);
}

int main(int argc, char *argv[])
{
  read_args(argc,argv);

  strtok(datafilename,"0");

  printf("newfilename:%s\n", datafilename);

  sprintf(datafilename, "%s%010.2lf[orb].binall", datafilename,t1);

  while((read_params(datafilename))==0){

    t1 += tt1;
    input_interval = tt1*2.*M_PI/OMEGA;
    serialnumber++;

    /* struct particle p = particles[0]; */
    /* printf("particle[0].v:%lf %lf %lf\n", p.vx,p.vy,p.vz); */
    compute_velocity_dispersion();
    
    out_velocity_dispersion();

    printf("time:%lf stdd:%lf\n", t1,stdd);

    strtok(datafilename,"0");
    sprintf(datafilename, "%s%010.2lf[orb].binall", datafilename,t1);
    printf("newfilename:%s\n", datafilename);

    once = 0;
  }

  return 0;
}
