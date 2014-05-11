/**
 * @file 	problem.c
 * @brief 	Example problem: bouncing balls.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example tests collision detection methods.
 * To change the collision detection algorithm, change the module
 * collisions_direct.c to either collisions_tree.c or 
 * collisions_sweep.c in the Makefile. 
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
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "collisions.h"
#include "collision_resolve.h"

// extern double coefficient_of_restitution; 
enum collisions_restitution_model collisions_restitution_model_info; /* collision model */
extern double constant_coefficient_of_restitution;
extern int exit_simulation;

void problem_init(int argc, char* argv[]){

  // Setup constants
  dt = 1e-4;
  tmax = 10000;
  boxsize = 3;
  constant_coefficient_of_restitution = 1; // elastic collisions
  // Do not use any ghost boxes
  nghostx = 0; nghosty = 0; nghostz = 0;

  // Setup particle structures
  init_box();
  // Initial conditions
  struct particle p;
  p.x  = 0.5; p.y  = 0.0; p.z  = 0.;
  p.vx = 0.0; p.vy = 0.0; p.vz = 0;
  p.ax = 0; p.ay = 0; p.az = 0;
  p.m  = 4.;
  p.r  = 1.;

  particles_add(p);

  p.x  = 0.0; p.y  = 0.; p.z  = 0.;
  p.vx = 1e-5; p.vy =  1.; p.vz =  0;
  p.ax =  0; p.ay =  0; p.az =  0;
  p.m  = 1.;
  p.r  = 1.;

  particles_add(p);

  collisions_search();
  collisions_resolve();
  exit(1);
}

void problem_inloop(){
}

void problem_output(){
}

void problem_finish(){
}
