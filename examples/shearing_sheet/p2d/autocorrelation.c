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

void readfile(char *filename);

void read_args(int argc, char *argv[]);

int main(int argc, char *argv[])
{
  printf("Hello,world\n");

  read_args(argc,argv);

  printf("data file:%s\n", datafilename);

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


