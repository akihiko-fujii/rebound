/**
 * @file 	collision_resolve.c
 * @brief 	Resolve a single collision. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @section LICENSE
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
#include "particle.h"
#include "collision_resolve.h"
#include "main.h"
#include "boundaries.h"
#include "communication_mpi.h"

double constant_coefficient_of_restitution;
double minimum_collision_velocity = 0;
double collisions_constant_coefficient_of_restitution_for_velocity(double v);
double (*coefficient_of_restitution_for_velocity) (double) = collisions_constant_coefficient_of_restitution_for_velocity;
double 	collisions_plog =0;	/**< Keep track of momentum exchange (used to calculate collisional viscosity in ring systems. */
long	collisions_Nlog =0;	/**< Keep track of Number of collisions. */

/* void (*collision_resolve) (struct collision) = collision_resolve_hardsphere_original; */
void (*collision_resolve) (struct collision) = collision_resolve_hardsphere_revised;

/* #define COLLISION_DEBUG */

void collision_resolve_hardsphere_revised(struct collision c){

#ifdef COLLISION_DEBUG
  printf(ANSI_COLOR_RED "before p1x:%lf %lf %lf\n", particles[c.p1].x,particles[c.p1].y,particles[c.p1].z);
  printf("before p2x:%lf %lf %lf\n", particles[c.p2].x,particles[c.p2].y,particles[c.p2].z);

  printf("before p1v:%lf %lf %lf\n", particles[c.p1].vx,particles[c.p1].vy,particles[c.p1].vz);
  printf("before p2v:%lf %lf %lf\n" ANSI_COLOR_RESET, particles[c.p2].vx,particles[c.p2].vy,particles[c.p2].vz);
#endif

#ifndef COLLISIONS_NONE
  struct particle p1 = particles[c.p1];
  struct particle p2;
#ifdef MPI
  int isloc = communication_mpi_rootbox_is_local(c.ri);
  if (isloc==1){
#endif // MPI
    p2 = particles[c.p2];
#ifdef MPI
  }else{
    int root_n_per_node = root_n/mpi_num;
    int proc_id = c.ri/root_n_per_node;
    p2 = particles_recv[proc_id][c.p2];
  }
#endif // MPI
  //	if (p1.lastcollision==t || p2.lastcollision==t) return;
  struct ghostbox gb = c.gb;
  double x21  = p1.x + gb.shiftx  - p2.x; 
  double y21  = p1.y + gb.shifty  - p2.y; 
  double z21  = p1.z + gb.shiftz  - p2.z; 
  double r = sqrt(x21*x21 + y21*y21 + z21*z21);
  /* double r21 = sqrt(x21*x21 + y21*y21 + z21*z21); */
  double rp   = p1.r+p2.r;
  double oldvyouter;

#ifdef COLLISION_DEBUG
  /* printf("before p*x:%lf %lf %lf\n", p1.vx*p1.m + p2.vx*p2.m, p1.vy*p1.m + p2.vy*p2.m, p1.vz*p1.m + p2.vz*p2.m); */
  printf("p1.m,p2.m:%lf %lf\n", p1.m,p2.m);
  printf("x21,y21,z21:%lf %lf %lf\n", x21,y21,z21);
#endif

  if (x21>0){
    oldvyouter = p1.vy;
  }else{
    oldvyouter = p2.vy;
  }

  if (rp*rp < x21*x21 + y21*y21 + z21*z21) return;

  double vx21 = p1.vx + gb.shiftvx - p2.vx; 
  double vy21 = p1.vy + gb.shiftvy - p2.vy; 
  double vz21 = p1.vz + gb.shiftvz - p2.vz; 

#ifdef COLLISION_DEBUG
  printf("p1.vx,vy,vz:%.3lf %.3lf %.3lf\n", p1.vx,p1.vy,p1.vz);
  printf("p2.vx,vy,vz:%.3lf %.3lf %.3lf\n", p2.vx,p2.vy,p2.vz);
  printf("gb.shiftvx,y,z:%.3lf,%.3lf,%.3lf\n", gb.shiftvx,gb.shiftvy,gb.shiftvz);
  printf("vx21,vy21,vz21:%lf %lf %lf\n", vx21,vy21,vz21);
#endif

  if (vx21*x21 + vy21*y21 + vz21*z21 >0) return; // not approaching

  // Bring the to balls in the xy plane.
  // NOTE: this could probabely be an atan (which is faster than atan2) 
  double theta = atan2(z21,y21);
  double stheta = sin(theta);
  double ctheta = cos(theta);
  double vy21n = ctheta * vy21 + stheta * vz21;
  double y21n = ctheta * y21 + stheta * z21;

  // Bring the two balls onto the positive x axis.
  double phi = atan2(y21n,x21);
  double cphi = cos(phi);
  double sphi = sin(phi);
  double vx21nn = cphi * vx21  + sphi * vy21n;
  double vy21nn = -sphi* vx21  + cphi * vy21n;

#ifdef COLLISION_DEBUG
  printf("phi,cphi:%lf %lf\n", phi,cphi);
  printf("vx21nn, vy21nn:%lf %lf\n", vx21nn,vy21nn);
  printf("r=%lf,rp=%lf,r*vy21nn=%lf\n", r,rp,r*vy21nn);
#endif

  // Coefficient of restitution
  double eps= coefficient_of_restitution_for_velocity(vx21nn);
  double dvx2 = -(1.0+eps)*vx21nn;
  double dvy2 = (r/rp-1.)*vy21nn;

#ifdef COLLISION_DEBUG
  printf("dvy2=%lf,rp*(vy21nn+dvy2)=%lf\n", dvy2,rp*(vy21nn+dvy2));
#endif

#ifdef COLLISION_DEBUG
  printf(ANSI_COLOR_BLUE "dvx2,dvy2:%lf,%lf\n" ANSI_COLOR_RESET, dvx2,dvy2);
#endif
  double minr = (p1.r>p2.r)?p2.r:p1.r;
  double maxr = (p1.r<p2.r)?p2.r:p1.r;
  double mindv= minr*minimum_collision_velocity;
  mindv *= 1.-(r - maxr)/minr;
  if (mindv>maxr*minimum_collision_velocity)mindv = maxr*minimum_collision_velocity;
  if (dvx2<mindv) dvx2 = mindv;

  // added
  double dxx2 = rp-r;
  double dxx2n = cphi * dxx2;
  double dxy2n = sphi * dxx2;
  double dxy2nn = ctheta * dxy2n;
  double dxz2nn = stheta * dxy2n;

  // Now we are rotating backwards
  /* double dvx2n = cphi * dvx2;	 */
  /* double dvy2n = sphi * dvx2; */

  // updated
  double dvx2n = cphi * dvx2 - sphi * dvy2;
  double dvy2n = sphi * dvx2 + cphi * dvy2;

  double dvy2nn = ctheta * dvy2n;
  double dvz2nn = stheta * dvy2n;

  // Applying the changes to the particles.
#ifdef MPI
  if (isloc==1){
#endif // MPI

    const double p1pf = p1.m/(p1.m+p2.m);
    const double p2pf = p2.m/(p1.m+p2.m);
    particles[c.p2].vx -=	p1pf*dvx2n;
    particles[c.p2].vy -=	p1pf*dvy2nn;
    particles[c.p2].vz -=	p1pf*dvz2nn;
    particles[c.p2].lastcollision = t;

    // added
    particles[c.p2].x -=	p1pf*dxx2n;
    particles[c.p2].y -=	p1pf*dxy2nn;
    particles[c.p2].z -=	p1pf*dxz2nn;
#ifdef MPI
  }
#endif // MPI
  particles[c.p1].vx +=	p2pf*dvx2n;
  particles[c.p1].vy +=	p2pf*dvy2nn;
  particles[c.p1].vz +=	p2pf*dvz2nn;

  // added
  particles[c.p1].x +=	p2pf*dxx2n; 
  particles[c.p1].y +=	p2pf*dxy2nn; 
  particles[c.p1].z +=	p2pf*dxz2nn; 

  particles[c.p1].lastcollision = t;

#ifdef COLLISION_DEBUG
  /* printf("after p*x:%lf %lf %lf\n", */
  /*        particles[c.p1].vx*particles[c.p1].m + particles[c.p2].vx*particles[c.p2].m, */
  /*        particles[c.p1].vy*particles[c.p1].m + particles[c.p2].vy*particles[c.p2].m, */
  /*        particles[c.p1].vz*particles[c.p1].m + particles[c.p2].vz*particles[c.p2].m); */
  printf(ANSI_COLOR_RED "after p1x:%lf %lf %lf\n", particles[c.p1].x,particles[c.p1].y,particles[c.p1].z);
  printf("after p2x:%lf %lf %lf\n", particles[c.p2].x,particles[c.p2].y,particles[c.p2].z);
  printf("after p1v:%lf %lf %lf\n", particles[c.p1].vx,particles[c.p1].vy,particles[c.p1].vz);
  printf("after p2v:%lf %lf %lf\n" ANSI_COLOR_RESET, particles[c.p2].vx,particles[c.p2].vy,particles[c.p2].vz);
#endif

  // Return y-momentum change
  if (x21>0){
    collisions_plog += -fabs(x21)*(oldvyouter-particles[c.p1].vy) * p1.m;
    collisions_Nlog ++;
  }else{
    collisions_plog += -fabs(x21)*(oldvyouter-particles[c.p2].vy) * p2.m;
    collisions_Nlog ++;
  }

#endif // COLLISIONS_NONE
}

void collision_resolve_hardsphere_original(struct collision c){
#ifndef COLLISIONS_NONE
	struct particle p1 = particles[c.p1];
	struct particle p2;
#ifdef MPI
	int isloc = communication_mpi_rootbox_is_local(c.ri);
	if (isloc==1){
#endif // MPI
		p2 = particles[c.p2];
#ifdef MPI
	}else{
		int root_n_per_node = root_n/mpi_num;
		int proc_id = c.ri/root_n_per_node;
		p2 = particles_recv[proc_id][c.p2];
	}
#endif // MPI
//	if (p1.lastcollision==t || p2.lastcollision==t) return;
	struct ghostbox gb = c.gb;
	double m21  = p1.m  /  p2.m; 
	double x21  = p1.x + gb.shiftx  - p2.x; 
	double y21  = p1.y + gb.shifty  - p2.y; 
	double z21  = p1.z + gb.shiftz  - p2.z; 
	double rp   = p1.r+p2.r;
	if (rp*rp < x21*x21 + y21*y21 + z21*z21) return;
	double vx21 = p1.vx + gb.shiftvx - p2.vx; 
	double vy21 = p1.vy + gb.shiftvy - p2.vy; 
	double vz21 = p1.vz + gb.shiftvz - p2.vz; 
	if (vx21*x21 + vy21*y21 + vz21*z21 >0) return; // not approaching
	// Bring the to balls in the xy plane.
	// NOTE: this could probabely be an atan (which is faster than atan2)
	double theta = atan2(z21,y21);
	double stheta = sin(theta);
	double ctheta = cos(theta);
	double vy21n = ctheta * vy21 + stheta * vz21;	
	double y21n = ctheta * y21 + stheta * z21;	
	
	// Bring the two balls onto the positive x axis.
	double phi = atan2(y21n,x21);
	double cphi = cos(phi);
	double sphi = sin(phi);
	double vx21nn = cphi * vx21  + sphi * vy21n;		

	// Coefficient of restitution
	double eps= coefficient_of_restitution_for_velocity(vx21nn);
	double dvx2 = -(1.0+eps)*vx21nn/(1.0+m21) ;
	if (dvx2<minimum_collision_velocity){
		dvx2 = minimum_collision_velocity;
	}

	// Now we are rotating backwards
	double dvx2n = cphi * dvx2;		
	double dvy2n = sphi * dvx2;		
	double dvy2nn = ctheta * dvy2n;	
	double dvz2nn = stheta * dvy2n;	

	// Log y-momentum change
	collisions_plog += fabs(dvy2nn*p1.m*x21);
	collisions_Nlog++;

	// Applying the changes to the particles.
#ifdef MPI
	if (isloc==1){
#endif // MPI
		particles[c.p2].vx -=	m21*dvx2n;
		particles[c.p2].vy -=	m21*dvy2nn;
		particles[c.p2].vz -=	m21*dvz2nn;
		particles[c.p2].lastcollision = t;
#ifdef MPI
	}
#endif // MPI
	particles[c.p1].vx +=	dvx2n; 
	particles[c.p1].vy +=	dvy2nn; 
	particles[c.p1].vz +=	dvz2nn; 
	particles[c.p1].lastcollision = t;
#endif // COLLISIONS_NONE
}

double collisions_constant_coefficient_of_restitution_for_velocity(double v){
	return constant_coefficient_of_restitution;
}
