/**
 * @file 	problem.c
 * @brief 	Example problem: shearing sheet.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This problem uses shearing sheet boundary
 * conditions. Particle properties resemble those found in 
 * Saturn's rings. 
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"
#include "display.h"

extern double OMEGA;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double); 

double coefficient_of_restitution_bridges(double v); 
double collisions_constant_coefficient_of_restitution_for_velocity(double v);

extern double constant_coefficient_of_restitution;
extern double opening_angle2;

// http://stackoverflow.com/questions/15767691/whats-the-c-library-function-to-generate-random-string

void generate_simulation_id(char *dest, size_t length) {

  srand((unsigned)time(NULL));

  char charset[] = 
    "23456789"
    "abcdefghijkmnprstuvwxy"
    "ABCDEFGHJKLMNPQRSTUVWXY";

  while (length-- > 0) {
    size_t index = (double) rand() / RAND_MAX * (sizeof charset - 1);
    *dest++ = charset[index];
  }
  *dest = '\0';
}

void problem_init(int argc, char* argv[]){

  // Setup constants
#ifdef GRAVITY_TREE
  opening_angle2	= .5;
#endif // GRAVITY_TREE
  OMEGA 			= 0.0001947529221860360;	// 1/s
  G 				= 6.67428e-11;		// N / (1e-5 kg)^2 m^2
  softening 			= 0.1;			// m
  dt 				= 1e-3*2.*M_PI/OMEGA;	// s
#ifdef OPENGL
  // Rotate the box by (display_rotate_z) around the z axis, then 
  display_rotate_z		= 0;			
  // rotate the box by (display_rotate_x) degrees around the x axis	
  display_rotate_x		= 0;			
#ifdef LIBPNG
  system("mkdir png");
#endif // LIBPNG
#endif // OPENGL
  root_nx = 2; root_ny = 2; root_nz = 1;
  nghostx = 2; nghosty = 2; nghostz = 0; 			// Use two ghost rings

  double surfacedensity 	= 580.; 			// kg/m^2
  double particle_density	= 900.;			// kg/m^3
  double particle_radius_min 	= 1.;			// m
  double particle_radius_max 	= 1.;			// m
  double particle_radius_slope 	= -3.;	
  boxsize 			= 100.;			// m
  if (argc>1){						// Try to read boxsize from command line
    boxsize = atof(argv[1]);
  }
  init_box();
	
  // Initial conditions
  /* printf("Toomre wavelength: %f\n",2.*M_PI*M_PI*surfacedensity/OMEGA/OMEGA*G); */


  /**
     Use Bridges et al coefficient of restitution.  
  */
  /* coefficient_of_restitution_for_velocity = coefficient_of_restitution_bridges; */
  /* boxinfo.coefficient_of_restitution_for_velocity = COLLISIONS_CONSTANT_COEFFICIENT_OF_RESTITUTION_FOR_VELOCITY; */
  /* end. */

  /**
     Use constant coefficient of restitution.
  */
  coefficient_of_restitution_for_velocity	= collisions_constant_coefficient_of_restitution_for_velocity;
  boxinfo.collisions_restitution_model_info	= ENUMCOLLISIONS_CONSTANT_COEFFICIENT_OF_RESTITUTION_FOR_VELOCITY;
  constant_coefficient_of_restitution = 0.1;
  /* end. */

  boxinfo.boxsize_x = boxsize_x;
  boxinfo.boxsize_y = boxsize_y;
  boxinfo.boxsize_z = boxsize_z;
  boxinfo.G = G;
  boxinfo.N = N;
  boxinfo.softening = softening;
  boxinfo.OMEGA = OMEGA;

  minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  // small fraction of the shear 

  double total_mass = surfacedensity*boxsize_x*boxsize_y;
  double mass = 0;
  while(mass<total_mass){
    struct particle pt;
    pt.x 		= tools_uniform(-boxsize_x/2.,boxsize_x/2.);
    pt.y 		= tools_uniform(-boxsize_y/2.,boxsize_y/2.);
    pt.z 		= tools_normal(1.);					// m
    pt.vx 		= 0;
    pt.vy 		= -1.5*pt.x*OMEGA;
    pt.vz 		= 0;
    pt.ax 		= 0;
    pt.ay 		= 0;
    pt.az 		= 0;
    double radius 	= tools_powerlaw(particle_radius_min,particle_radius_max,particle_radius_slope);
#ifndef COLLISIONS_NONE
    pt.r 		= radius;						// m
#endif
    double		particle_mass = particle_density*4./3.*M_PI*radius*radius*radius;
    pt.m 		= particle_mass; 	// kg
    particles_add(pt);
    mass += particle_mass;
  }

  // Optical depth
  printf("Optical depth: %f\n", tools_crosssection()/(boxsize_x*boxsize_y));

  // generate simulation id 
  simulation_id = calloc(7, sizeof(char));
  generate_simulation_id(simulation_id,7);
  printf("simulation_id=%s\n", simulation_id); 

  /* prepare data directory */
  char *dirname = calloc(1024, sizeof(char));
  sprintf(dirname, "../restarting_simulation/data/%s",simulation_id);

  struct stat st = {0};

  if(stat(dirname, &st) == -1){
    mkdir(dirname, S_IRWXU | S_IRWXG | S_IRWXO);

    /* prepare data directory */
    sprintf(dirname,"%s/snapshots",dirname);
    mkdir(dirname, S_IRWXU | S_IRWXG | S_IRWXO);

  }else{			/* rare event */
    printf("could not create directory.\n"); exit(1);
  }

  free(dirname);
}

double coefficient_of_restitution_bridges(double v){
	// assumes v in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}

void problem_inloop(){
}

void problem_output(){

#ifdef LIBPNG
  if (output_check(1e-3*2.*M_PI/OMEGA)){
    output_png("png/");
  }
#endif //LIBPNG

  output_timing();

#ifdef OUTPUT

  double output_interval = 0.01*2.*M_PI/OMEGA;

  if (output_check(output_interval)){

    /* printf("cofficient of restitution:%d %lf", ENUMCOLLISIONS_CONSTANT_COEFFICIENT_OF_RESTITUTION_FOR_VELOCITY,constant_coefficient_of_restitution); */
    /* printf("boxinfo:%lf %lf %lf\n", boxinfo.boxsize_x,boxinfo.boxsize_y,boxinfo.boxsize_z); */

    char *o = calloc(1024, sizeof(char));

    sprintf(o,"../restarting_simulation/data/%s/snapshots/%010.2f[orb].binall",
	    simulation_id,t/(2.*M_PI/OMEGA));
    output_binary_all(o);

    sprintf(o,"../restarting_simulation/data/%s/snapshots/%010.2f[orb].bin",
	    simulation_id,t/(2.*M_PI/OMEGA));
    output_binary(o);

    sprintf(o,"../restarting_simulation/data/%s/snapshots/%010.2f[orb].ascii",
	    simulation_id,t/(2.*M_PI/OMEGA));
    output_ascii(o);

    free(o);
  }
#endif

}

  void problem_finish(){
  }
