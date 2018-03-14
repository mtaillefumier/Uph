#ifndef _derivatives_h_
#define _derivatives_h_
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "System.h"

#ifdef __cplusplus
extern "C" {
#endif

void derivatives(System *sys, double t);

#ifdef __cplusplus
}
#endif

/* if using the gcc compiler, add -lm at end of compile command*/
#endif
