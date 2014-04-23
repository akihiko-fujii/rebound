/**
 * @file 	main.h
 * @brief 	Main header file with widely used global variables.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
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
#ifndef _MAIN_H
#define _MAIN_H

#ifndef M_PI
// Make sure M_PI is defined. 
#define M_PI           3.14159265358979323846
#endif

// http://stackoverflow.com/questions/3219393/stdlib-and-colored-output-in-c
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


// http://7ujm.net/etc/esc.html
#define clr		printf("\033[2J") // Clear display
#define clr_right		printf("\033[0K") // カーソル位置からその行の右端までをクリア
#define clr_left		printf("\033[1K") //カーソル位置からその行の左端までをクリア
#define clr_line		printf("\033[2K") //カーソル位置の行をクリア
#define location(x,y)	printf("\033[%d;%dH" ,x,y) //カーソル位置を移動
#define right(x)		printf("\033[%dC" ,x) //カーソルを指定数だけ右に移動
#define left(x)		printf("\033[%dD" ,x) //カーソルを指定数だけ左に移動
#define down(x)		printf("\033[%dB" ,x) //カーソルを指定数だけ下に移動
#define up(x)		printf("\033[%dA" ,x) //カーソルを指定数だけ上に移動

extern double 	softening;	/**< Gravitational softening parameter. Default: 0. */
extern double 	G;		/**< Gravitational constant. Default: 1. */
extern double 	t;		/**< Current simulation time. */
extern double 	tmax;		/**< Maximum simulation time. Simulation stops if t>=tmax. Simulation runs forever if t==0.*/
extern double 	dt;		/**< Current timestep. */
extern int 	exit_simulation;/**< Set to 1 to exit the simulation at the end of the timestep. */
extern int 	N;		/**< Current number of particles on this node. */
extern int 	N_active;	/**< Number of massive particles included in force calculation. Default: N.*/
extern int 	root_nx;	/**< Number of root boxes in x direction. Default: 1. */
extern int 	root_ny;	/**< Number of root boxes in y direction. Default: 1. */
extern int 	root_nz;	/**< Number of root boxes in z direction. Default: 1. */
extern int 	root_n;		/**< Total number of root boxes in all directions, root_nx*root_ny*root_nz. Default: 1. Set in box_init().*/
extern double 	boxsize;	/**< Size of a root box. Needs to be set in problem_init(). */
extern double 	boxsize_x;	/**< Size of the entire box in the x direction, root_nx*boxsize. Set in box_init().*/
extern double 	boxsize_y;	/**< Size of the entire box in the y direction, root_ny*boxsize. Set in box_init().*/
extern double 	boxsize_z;	/**< Size of the entire box in the z direction, root_nz*boxsize. Set in box_init().*/
extern double 	boxsize_max;	/**< Maximum size of the entire box in any direction. Set in box_init().*/

extern double 	timing_initial;	/**< System time at start. Used to meassure total cpu time. */

extern char* 	simulation_id;	/**< Used to identify simulation models. */


/**
 * Initializes the box. 
 * This function needs to be called from problem_init() before any particles are added.
 */
void init_box();

/**
 * Initializes one rootbox with width _boxwidth. 
 * This function needs to be called from problem_init() before any particles are added.
 * @param _boxwidth Size of the box.
 */
void init_boxwidth(double _boxwidth);

/**
 * Main iteration loop.
 * All the work is done within this function.
 * When OpenGL is not used, this function is called by a loop in main(). 
 * When OpenGL is used, this function is called by OpenGL directly. 
 */
void iterate();
#endif
