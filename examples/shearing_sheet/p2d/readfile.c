#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include "particle.h"
#include "collision_resolve.h"
#include "boundaries.h"
#include "collisions.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

int N; int dim1,dim2;		/* autocorrelation resolution */
extern double boxsize_x,boxsize_y,boxsize_z,G,dt,softening,OMEGA;
extern double resolution_x,resolution_y;
extern enum collisions_restitution_model collisions_restitution_model_info; /* collision model */
extern double constant_coefficient_of_restitution;
extern double input_interval;

struct particle* particles;

extern char *datafilename;
extern char *directoryname;

double p2d_crosssection(){
  double s = 0.;
  for(int i=0;i<N;i++){
    s += particles[i].r * particles[i].r;
  }
  s *= atan(1.)*4.;
  return s;
}

double p2d_surfacedensity(){

  double sigma = 0.;
  for(int i=0;i<N;i++){
    sigma += particles[i].m;
  }
  return sigma/(boxsize_x*boxsize_y);
}

/* void boxinfo_expand(){ */
  
/*   boxsize_x				= boxinfo.boxsize_x; */
/*   boxsize_y				= boxinfo.boxsize_y; */
/*   boxsize_z				= boxinfo.boxsize_z; */
/*   G					= boxinfo.G; */
/*   dt					= boxinfo.dt; */
/*   N					= boxinfo.N; */
/* #ifdef INTEGRATOR_SEI 	// Shearing sheet */
/*   OMEGA					= boxinfo.OMEGA; */
/* #endif  */
/*   collisions_restitution_model_info	= boxinfo.collisions_restitution_model_info; */
  
/* } */

void help_show(){

  printf("this is the help for command.\n");

}

int result = 0;
char name[1024];

static struct option long_opt[] = {
  {"help", no_argument, NULL, 'h'},
  /* {"filename", required_argument, NULL, 'f'}, */
    {"id", required_argument, NULL, 'i'},
      {"resolution_x", required_argument, NULL, 'x'},
  {"resolution_y", required_argument, NULL, 'y'},
  {"input_interval", required_argument, NULL, 't'},
    {0,0,0,0}
};

void read_args(int argc, char *argv[]){

  if(argc==1){
    printf(
	   ANSI_COLOR_RED
	   "need at least one option. see -h for help.\n" ANSI_COLOR_RESET); exit(0);
  }
  static char *datadir_path;
  datadir_path = getenv("REBOUND_DATA");    

  int option_index = -11;   

  while((result = getopt_long(argc, argv, "hi:x:y:", long_opt, &option_index))!=-1){

    /* printf("result = %d\n", result); */

    switch(result){
    case 'h':
      help_show(); exit(0);
      break;
      /* case 'f': */
      /*     printf("option %s applied with %s\n", long_opt[option_index].name, optarg); */
      /*     sprintf(name, "%s", optarg); */
      /*     datafilename = calloc(1024, sizeof(char)); */
      /*     sprintf(datafilename, "%s", name); */
      /*     break; */
    case 'x':
      resolution_x = atof(optarg);
      break;
    case 'y':
      resolution_y = atof(optarg);   
      break;
    case 'i':
      sprintf(name, "%s", optarg);
      datafilename = calloc(1024, sizeof(char));
      directoryname = calloc(1024, sizeof(char));
      sprintf(datafilename, "%s/%s/snapshots/0000000.00[orb].binall",datadir_path,name);
      sprintf(directoryname, "%s/%s",datadir_path,name);
      /* sprintf(datafilename, "../restarting_simulation/data/%s/snapshots/0000000.00[orb].binall",name); */
      /* sprintf(directoryname, "../restarting_simulation/data/%s", name); */
      break;
    default:
      printf(
	     ANSI_COLOR_RED "an unknown option: program ends. see -h for help.\n" ANSI_COLOR_RESET);
      exit(1);
    }

  }
}

/* void read_params(char *filename){ */

/*   int filesize = 0; */
/*   FILE* inf = fopen(filename,"rb");  */

/*   if(inf==NULL){ */
/*   printf("cannot read file %s\n.end.(from function %s)\n", filename,__func__); exit(1); */
/*   } */

/*   // Get binary file size */
/*   fseek(inf, 0, SEEK_END);  */
/*   filesize = ftell(inf);  */
/*   fseek(inf, 0, SEEK_SET);  */

/*   long objects = 0; */
/*   int _N; */
/*   double t; */
  
/*   objects += fread(&_N,sizeof(int),1,inf); */
/*   objects += fread(&t,sizeof(double),1,inf); */

/*   // Show warnings when file of binary file is imconpatible.  */
/*   if(filesize != sizeof(int)+sizeof(double)+_N*sizeof(struct particle)+sizeof(struct boxinfo)){ */
/*     printf("sizeof(file):%d. sizeof(int)+sizeof(double)+%d*sizeof(particle)+sizeof(boxinfo)=%ld\n", */
/* 	   filesize,_N,sizeof(int)+sizeof(double)+_N*sizeof(struct particle)+sizeof(struct boxinfo)); */
/*   } */

/*   N = _N; */
/*   particles = calloc(N, sizeof(struct particle)); */
/*   for (int i=0;i<_N;i++){ */
/*     struct particle p; */

/*     objects += fread(&p,sizeof(struct particle),1,inf); */
/*     particles[i] = p;     */
/*   } */

/*   fread(&boxinfo, sizeof(struct boxinfo), 1, inf); */

/*   fclose(inf); */

/*   boxinfo_expand();   */

/*   printf(ANSI_COLOR_BLUE  */
/*     "boxinfo:\nboxsize:%lf %lf %lf(m)\n" */
/*     "G:%lfx1e-11\ndt:%lf(s)\nN:%d\nsoftening:%lf(m)\nOMEGA:%lf(/s)\n" */
/*     ANSI_COLOR_RESET, */
/*     boxinfo.boxsize_x,boxinfo.boxsize_y,boxinfo.boxsize_z,boxinfo.G*1e11,boxinfo.dt, */
/*     boxinfo.N,boxinfo.softening,boxinfo.OMEGA */
/*     ); */

/*   printf(ANSI_COLOR_CYAN */
/*     "optical depth: %lf\n" "Toomre wavelength (assuming Omega=kappa): %lf(m)\n" */
/*     "Hill sphere:%lf (m)\n" "Normalized with 2R_phys:%lf\n" */
/*     ANSI_COLOR_RESET, */
/*     p2d_crosssection()/(boxsize_x*boxsize_y), */
/*     4.*M_PI*M_PI*G*p2d_surfacedensity()/OMEGA/OMEGA, */
/*     pow(2.*particles[0].m*G/(3.*OMEGA*OMEGA), 1./3.), */
/*     pow(2.*particles[0].m*G/(3.*OMEGA*OMEGA), 1./3.)/(2.*particles[0].r) */
/*     ); */
/* } */
 

