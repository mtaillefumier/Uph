#ifndef _NEIGHBORS_H
#define _NEIGHBORS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mkl.h>

typedef struct _domain {
  int *index;
  int NumberOfParticles;
  double xmin, xmax, ymin, ymax;
  int indexSize;
} domain;

typedef struct _grid {
  int dim[2];
  int Length;
  int MaxPoints;
  double dx, dy;
  double complex xhigh, xlow;
  domain *Domain;
  int *WhereIsParticle;
} Grid;

Grid *InitializeGrid(const double complex xlow, const double complex xmax, const double dx, const double dy, const int MaxPoints);
void SortParticlesInToDomain(Grid *g, const double complex *Particles, const int MaxNumberOfParticles);
void FindNeighborsWithinDistance(const Grid *g,
				 const double complex *Particles,
				 const double *hsml,
				 const int ntotal,
				 const int MaxNumberOfParticles,
				 int *Neighbors,
				 int **j_pair,
				 double **rij2);
void DestroyGrid(Grid *g);
#endif
