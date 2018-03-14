/*
**************************************************
* including file for parameters and constants used
* in all "sph.c" files.
**************************************************
*/
#ifndef param_h
#define _param_h_

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define maxn 12000
#define max_int  120000

bool virt_pres; /* virtual pressure to prevent tensile instability */
bool art_visc; /* artificial viscosity for viscosity simul */
bool avg_vel; /* average velocity for xSPH variant */
bool input_file; /* input fluid particles from file */
bool vp_file ; /* input boundary particles from file */

#endif
