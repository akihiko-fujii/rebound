/**
 * @file 	problem.c
 * @brief 	Example problem: Radiation forces
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2013 Hanno Rein, Dave Spiegel
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
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "problem.h"

void additional_forces();
extern double integrator_epsilon;
double betaparticles = 0.01; 	// Beta parameter. 
				// Defined as the ratio of radiation pressure over gravity.

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 			= 1e-3;	// Initial timestep.
	integrator_epsilon 	= 1e-4;	// Accuracy parameter.
	boxsize 		= 10;	
	tmax			= 1e6;
	N_active		= 1; 	// Only the star is massive.
	problem_additional_forces 	= additional_forces;
	init_box();
	
	
	// Star
	struct particle star;
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = 1.;
	particles_add(star);


	// Planet 
	struct particle planet;
	planet.m  = 1e-3;
	planet.x  = 1; planet.y  = 0.; planet.z  = 0.;
	planet.ax = 0; planet.ay = 0; planet.az = 0;
	planet.vx = 0;
	planet.vy = sqrt(G*(star.m+planet.m)/planet.x);
	planet.vz = 0;
	particles_add(planet);
	
	
	tools_move_to_center_of_momentum();

	// Dust particles
	while(N<2){
		struct particle p; 
		p.m  = 0;	// massless
		double r = 1.;// tools_uniform(1.4,2);
		double v = sqrt(G*(star.m*(1.-betaparticles))/r);
		double phi = tools_uniform(0,2.*M_PI);
		p.x  = r*sin(phi);  p.y  = r*cos(phi); p.z  = 0; 
		p.vx = -v*cos(phi); p.vy = v*sin(phi); p.vz = 0;
		p.ax = 0; p.ay = 0; p.az = 0;
		particles_add(p); 
	}

	system("rm -v r.txt");	
}

void force_radiation(){
	const struct particle star = particles[0];				// cache
#pragma omp parallel for
	for (int i=0;i<N;i++){
		const struct particle p = particles[i]; 			// cache
		if (p.m!=0.) continue; 						// Only dust particles feel radiation forces
		const double prx  = p.x-star.x;
		const double pry  = p.y-star.y;
		const double prz  = p.z-star.z;
		const double pr   = sqrt(prx*prx + pry*pry + prz*prz); 	// distance relative to star
		const double prvx = p.vx-star.vx;
		const double prvy = p.vy-star.vy;
		const double prvz = p.vz-star.vz;

		const double c 		= 1.006491504759635e+04; 		// speed of light.
		const double rdot 	= (prvx*prx + prvy*pry + prvz*prz)/pr; 	// radial velocity relative to star
		const double F_r 	= betaparticles*G*star.m/(pr*pr);

		// Equation (5) of Burns, Lamy, Soter (1979)
		particles[i].ax += F_r*((1.-rdot/c)*prx/pr - prvx/c);
		particles[i].ay += F_r*((1.-rdot/c)*pry/pr - prvy/c);
		particles[i].az += F_r*((1.-rdot/c)*prz/pr - prvz/c);
	}
}

void additional_forces(){
	force_radiation();
}

void problem_inloop(){
	if(output_check(4000.*dt)){
		output_timing();
	}
	if(output_check(M_PI*2000.)){ // output every 1000 years
		FILE* f = fopen("r.txt","a");
		const struct particle star = particles[0];
		for (int i=1;i<N;i++){
			const struct particle p = particles[i]; 
			const double prx  = p.x-star.x;
			const double pry  = p.y-star.y;
			const double prz  = p.z-star.z;
			const double pr   = sqrt(prx*prx + pry*pry + prz*prz); 	// distance relative to star
			fprintf(f,"%e\t%e\n",t,pr);
		}
		fclose(f);
	}
}

void problem_output(){
}

void problem_finish(){
}