#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <ftw.h>
#include <getopt.h>
#include <dirent.h>
#include "particle.h"
#include "readfile.h"
#include "../../../src/input_params.h"
#include "../../../src/main.h"

static unsigned int total = 0;

char *datafilename, *directoryname;
double boxsize_x,boxsize_y,boxsize_z,G,dt,softening,OMEGA;
double resolution_x,resolution_y;

double input_interval;
int serialnumber = 0;
double **densitymap;
int arraysize_x,arraysize_y;

int sum(const char *fpath, const struct stat *sb, int typeflag) {

  /* printf("loading file %s\n", fpath); */
  total += sb->st_size;
  return 0;
}

struct particle projection(struct particle p){
  
  struct particle q;
  
  q.x = p.x;
  q.y = p.y;
  q.z = p.z;
  q.vx = p.vx;
  q.vy = p.vy;
  q.vz = p.vz;
  q.r = p.r;
  q.m = p.m;

  return q;
}

void compute_density(){

  arraysize_x = boxsize_x / resolution_x + 1;
  arraysize_y = boxsize_y / resolution_y + 1;

  densitymap = (double **)calloc(arraysize_x,sizeof(double*));

  int i;
  for(i=0;i<arraysize_x;i++){
    densitymap[i] = (double *)calloc(arraysize_y,sizeof(double));
  }

  double resolution_area = resolution_x*resolution_y;

  for(i=0;i<N;i++){

    double _x,_y;
    struct particle p_proj = projection(particles[i]);
    /* struct particle p_proj = particles[i]; */

    if(p_proj.x<-boxsize_x/2.) p_proj.x += boxsize_x;    if(p_proj.y<-boxsize_y/2.) p_proj.y += boxsize_y;
    if(p_proj.x> boxsize_x/2.) p_proj.x -= boxsize_x;    if(p_proj.y> boxsize_y/2.) p_proj.y -= boxsize_y;

    _x = floor((p_proj.x+boxsize_x/2.)/resolution_x);
    _y = floor((p_proj.y+boxsize_y/2.)/resolution_y);

    if(0<=(int)_x && (int)_x <=arraysize_x && 0<=(int)_y && (int)_y<=arraysize_y){
      densitymap[(int)_x][(int)_y] += p_proj.m/resolution_area;
    }else{
      /* if(0>_x || 0>_y ){ */
      /* 	printf("x,y:%lf %lf arraysize:%d %d\nm", _x,_y,arraysize_x,arraysize_y); exit(1); */
      /* } */
      /* if((int)_x > arraysize_x || (int)_y>arraysize_y){ */
      /* 	printf("x,y:%lf %lf arraysize:%d %d\n", _x,_y,arraysize_x,arraysize_y);exit(1); */
      /* } */
      printf("particle %d out of range. at %s %d %s. end.\n",i,__FILE__,__LINE__,__func__);  exit(1);
    }
  }
}

void write_density(){

  char o[1024],o2[1024];
  /* prepare data directory */
  sprintf(o,"%s/densitymap[%dx%d]",directoryname,arraysize_x,arraysize_y);
  mkdir(o, S_IRWXU | S_IRWXG | S_IRWXO);
  sprintf(o2, "%s/%010.2f[orb].surfacedensity",o,input_interval/(2.*M_PI/OMEGA)*serialnumber);
  
  FILE *of; 
  if((of=fopen(o2, "w"))==NULL){
    printf("cannot open file %s from function %s.\n", o2,__func__); exit(1);
  }

  int i,j;
  for(int i=0;i<arraysize_x;i++){
    for(int j=0;j<arraysize_y;j++){
      fprintf(of, "%lf %lf %lf\n", i*resolution_x,j*resolution_y,densitymap[i][j]);      
    }
    fprintf(of, "\n");
  }

  fclose(of);

}

double density_average(){

  int i,j;
  double sigma = 0.;

  for(int i=0;i<arraysize_x;i++){
    for(int j=0;j<arraysize_y;j++){
      sigma += densitymap[i][j];
    }
  }
  sigma = sigma / (double)((arraysize_x-1)*(arraysize_y-1));

  return sigma;
}


int main(int argc, char **argv) {

  resolution_x = 5.;   resolution_y = 5.; /* set as default */

  read_args(argc,argv);

  read_params(datafilename);

  input_interval = 0.01*2.*M_PI/OMEGA;

  while(1){

    sprintf(datafilename, "%s/snapshots/%010.2lf[orb].binall",
	    directoryname,input_interval/(2.*M_PI/OMEGA)*(double)serialnumber);

    printf("reading %s resolution:%lf %lf\n",datafilename,resolution_x,resolution_y);

    up(1);

    read_params(datafilename);

    /* printf("dirfilename:%s\n", directoryname); */

    compute_density();

    write_density();

    serialnumber++;

    if(serialnumber>1000) break;

    /* free(densitymap); */
  }

  printf("average density:%lf\n", density_average());

  char o[1024],dir[1024];
  sprintf(dir,"%s/densitymap[%dx%d]",directoryname,arraysize_x,arraysize_y);
  sprintf(o, "echo | awk -f ./p2d/density_gnuplot.awk > %s/gnuplot.plt",dir);
  system(o);

  return 0;
}
