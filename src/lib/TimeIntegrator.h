#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "System.h"
#include "input.h"
#ifdef __cplusplus
extern "C" {
#endif

  void TimeIntegrator(System *sys, double t, double dt, const int itime);

#ifdef __cplusplus
}
#endif


#endif
