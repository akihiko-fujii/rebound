#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include "main.h"
#include "particle.h"

#ifdef INTEGRATOR_SEI 	// Shearing sheet
extern double OMEGA;
#endif
enum collisions_restitution_model collisions_restitution_model_info; /* collision model */
double constant_coefficient_of_restitution;
extern int exit_simulation;

double crosssection(){
  double s = 0.;
  for(int i=0;i<N;i++){
    s += particles[i].r * particles[i].r;
  }
  s *= atan(1.)*4.;
  return s;
}

double surfacedensity(){

  double sigma = 0.;
  for(int i=0;i<N;i++){
    sigma += particles[i].m;
  }
  return sigma/(boxsize_x*boxsize_y);
}

void boxinfo_expand(){
  
  boxsize_x				= boxinfo.boxsize_x;
  boxsize_y				= boxinfo.boxsize_y;
  boxsize_z				= boxinfo.boxsize_z;
  G					= boxinfo.G;
  dt					= boxinfo.dt;
  N					= boxinfo.N;
#ifdef INTEGRATOR_SEI 	// Shearing sheet
  OMEGA					= boxinfo.OMEGA;
#endif 
  collisions_restitution_model_info	= boxinfo.collisions_restitution_model_info;
  constant_coefficient_of_restitution = boxinfo.constant_coefficient_of_restitution;
}


int read_params(char *filename){

  int filesize = 0;
  FILE* inf = fopen(filename,"rb"); 

  if(inf==NULL){
    printf("cannot read file %s. ends.(from function %s)\n", filename,__func__); return 1;
  }

  // Get binary file size
  fseek(inf, 0, SEEK_END); 
  filesize = ftell(inf); 
  fseek(inf, 0, SEEK_SET); 

  long objects = 0;
  int _N;
  double t;
  
  objects += fread(&_N,sizeof(int),1,inf);
  objects += fread(&t,sizeof(double),1,inf);

  // Show warnings when file of binary file is imconpatible. 
  if(filesize != sizeof(int)+sizeof(double)+_N*sizeof(struct particle)+sizeof(struct boxinfo)){
    printf("sizeof(file):%d. sizeof(int)+sizeof(double)+%d*sizeof(particle)+sizeof(boxinfo)=%ld\n",
	   filesize,_N,sizeof(int)+sizeof(double)+_N*sizeof(struct particle)+sizeof(struct boxinfo));
  }

  N = _N;
  particles = calloc(N, sizeof(struct particle));
  for (int i=0;i<_N;i++){
    struct particle p;

    objects += fread(&p,sizeof(struct particle),1,inf);
    particles[i] = p;    
  }

  fread(&boxinfo, sizeof(struct boxinfo), 1, inf);

  fclose(inf);

  boxinfo_expand();  

#define PARAM
#ifdef PARAM
  clr; up(7);
  printf(ANSI_COLOR_BLUE
  	 "boxinfo:\nboxsize:%lf %lf %lf(m)\n"
  	 "G:%lfx1e-11\ndt:%lf(s)\nN:%d\nsoftening:%lf(m)\nOMEGA:%lf(/s)\n"
  	 ANSI_COLOR_RESET,
  	 boxinfo.boxsize_x,boxinfo.boxsize_y,boxinfo.boxsize_z,boxinfo.G*1e11,boxinfo.dt,
  	 boxinfo.N,boxinfo.softening,
#ifdef INTEGRATOR_SEI 	// Shearing sheet
	 boxinfo.OMEGA
#else
	 0.0
#endif
  	 );
  up(5);
  printf(ANSI_COLOR_CYAN
  	 "optical depth: %lf\n" "Toomre wavelength (assuming Omega=kappa): %lf(m)\n"
  	 "Hill sphere:%lf (m)\n" "Normalized with 2R_phys:%lf\n"
	 "coefficient_of_restitution law:%d if constant equals to:%lf\n"
  	 ANSI_COLOR_RESET,
  	 crosssection()/(boxsize_x*boxsize_y),
#ifdef INTEGRATOR_SEI 	// Shearing sheet
  	 4.*M_PI*M_PI*G*surfacedensity()/OMEGA/OMEGA,
  	 pow(2.*particles[0].m*G/(3.*OMEGA*OMEGA), 1./3.),
  	 pow(2.*particles[0].m*G/(3.*OMEGA*OMEGA), 1./3.)/(2.*particles[0].r),
#else
	 0.,0.,0.,
#endif
	 collisions_restitution_model_info,
	 constant_coefficient_of_restitution
  	 );
#endif

  return 0;
}
 



