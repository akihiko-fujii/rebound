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

double t1 = 0.; double tt1 = 0.1;

int main(int argc, char *argv[])
{
  read_args(argc,argv);

  strtok(datafilename,"0");

  printf("newfilename:%s\n", datafilename);

  sprintf(datafilename, "%s%010.2lf[orb].binall", datafilename,t1);

  while((read_params(datafilename))==0){
    t1 += tt1;
    strtok(datafilename,"0");
    sprintf(datafilename, "%s%010.2lf[orb].binall", datafilename,t1);
    printf("newfilename:%s\n", datafilename);
  }

  /* while((read_params(datafilename))==0){ */
  /*   t1 += tt1*0.1; */
  /* } */


  /* DIR *dp; */
  /* struct dirent *ep;      */

  /* char str11[1024]; */
  /* sprintf(str11,"%s/snapshots",directoryname); */
  /* dp = opendir(str11); */

  /* if (dp != NULL) */
  /*   { */
  /*     struct dirent *ep2; */
  /*     while (ep = readdir (dp)){ */
  /* 	ep2 = ep; */
  /*     } */
  /*     puts (ep2->d_name); */

  /*     (void) closedir (dp); */
  /*   } */
  /* else */
  /*   perror ("Couldn't open the directory"); */

  return 0;
}

