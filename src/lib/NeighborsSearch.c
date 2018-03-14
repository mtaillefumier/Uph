#include "Neighbors.h"

#include <stdio.h>
#include <string.h>
#include <fenv.h>   
#include <math.h>
#include <omp.h>
Grid *InitializeGrid(const double complex xlow, const double complex xmax, const double dx, const double dy, const int MaxPoints)
{
  Grid *g = (Grid *)malloc(sizeof(Grid));
  double complex diff = xmax-xlow;
  g->dim[0] = lrint(creal(diff)/dx)+1;
  g->dim[1] = lrint(cimag(diff)/dy)+1;
  g->Length = g->dim[0]*g->dim[1];
  g->Domain = (domain *) malloc(sizeof(domain)*g->Length);
  g->MaxPoints = MaxPoints;
  g->WhereIsParticle = (int *)malloc(sizeof(int)*MaxPoints);

  for(int x=0;x<g->dim[0];x++)
    for(int y=0;y<g->dim[1];y++) {
      
      g->Domain[x*g->dim[1]+y].xmin = creal(xlow) + x*dx;
      g->Domain[x*g->dim[1]+y].xmax = g->Domain[x*g->dim[1]+y].xmin + dx;
      g->Domain[x*g->dim[1]+y].ymin = cimag(xlow) + y*dy;
      g->Domain[x*g->dim[1]+y].ymax = g->Domain[x*g->dim[1]+y].ymin + dy;
      g->Domain[x*g->dim[1]+y].index = NULL;
      g->Domain[x*g->dim[1]+y].indexSize = 0;
      g->Domain[x*g->dim[1]+y].NumberOfParticles = 0;
    }
  g->dx=dx;
  g->dy=dy;
  g->xlow = xlow;
  g->xhigh=xmax;
  return g;
}

void SortParticlesInToDomain(Grid *g, const double complex *Particles, const int MaxNumberOfParticles)
{

  if(g->MaxPoints < MaxNumberOfParticles)
    g->WhereIsParticle = realloc(g->WhereIsParticle, sizeof(int)*MaxNumberOfParticles);
  
  for(int d=0;d<g->Length;d++) {
    g->Domain[d].NumberOfParticles = 0;
  }

  fesetround(FE_DOWNWARD);
  for(int particles=0;particles<MaxNumberOfParticles;particles++) {
    double complex tmp = Particles[particles] - g->xlow;
    int xi = lrint(creal(tmp)/g->dx);
    int yi = lrint(cimag(tmp)/g->dy);
    int dom = xi*g->dim[1]+yi;
    if((xi<0)||(xi>=g->dim[0])||(yi<0)||(yi>=g->dim[1])) {
      printf(" Error particle outside the grid %.10lf %.10lf\n", 
	     creal(Particles[particles]), 
	     cimag(Particles[particles]));
      exit(1);
    }

    if((dom >= g->dim[0]*g->dim[1])||(dom<0)) {
      printf("We have a problem\n");
    }
    
    if(g->Domain[dom].indexSize == 0) {
      g->Domain[dom].index = (int *)malloc(sizeof(int)*32);
      g->Domain[dom].indexSize = 32;
    }
    
    if(g->Domain[dom].indexSize <= g->Domain[dom].NumberOfParticles) {
      g->Domain[dom].indexSize += 32;      
      g->Domain[dom].index = realloc(g->Domain[dom].index, sizeof(int)*g->Domain[dom].indexSize);
      printf("%d %d %d %d\n", dom, g->Domain[dom].indexSize, g->Domain[dom].NumberOfParticles, particles);
    }
    
    g->Domain[dom].index[g->Domain[dom].NumberOfParticles] = particles;

    g->WhereIsParticle[particles]=dom;
    g->Domain[dom].NumberOfParticles++;
    //    printf("%.15lf %.15lf %d %d\n", creal(Particles[particles]), cimag(Particles[particles]), xi, yi);
  }
}

void FindNeighborsWithinDistance(const Grid *g,
				 const double complex *Particles,
				 const double *hsml,
				 const int ntotal,
				 const int MaxNumberOfParticles,
				 int *Neighbors,
				 int **j_pair,
				 double **rij2)
{
  memset(Neighbors, 0, sizeof(int)*ntotal);

#pragma omp parallel  
  {
    int threadid = omp_get_thread_num();
    int *tmpi = NULL;
    posix_memalign((void **)&tmpi, 16, sizeof(int)*MaxNumberOfParticles);
    double complex *tvi = NULL;
    posix_memalign((void **)&tvi, 16, sizeof(double complex)*MaxNumberOfParticles);
    double *Norm = NULL;
    posix_memalign((void **)&Norm, 16, sizeof(double)*MaxNumberOfParticles);
    
#pragma omp for 
    for(int particles=0;particles<ntotal;particles++) {
      int x = g->WhereIsParticle[particles]/g->dim[1];
      int y = g->WhereIsParticle[particles]%g->dim[1];
      int dom[9];
      dom[0] = g->WhereIsParticle[particles];
      int dx[9] = {0,1,1,0,-1,-1,-1,0,1};
      int dy[9] = {0,0,1,1,1,0,-1,-1,-1};
      int dd = 1;
      for(int d=1;d<9;d++)
        if(((x+dx[d])>=0)&&((x+dx[d])<g->dim[0]) && ((y+dy[d])>=0)&&((y+dy[d])<g->dim[1])) {
          
          dom[dd] = (x+dx[d])*g->dim[1] + y + dy[d];
          if(g->Domain[dom[dd]].NumberOfParticles)
            dd++;
        }
      /* dom[1] = (x+1)*g->dim[1]+y; */
      
      /* dom[2] = (x+1)*g->dim[1]+y+1; */
      /* dom[3] = (x)*g->dim[1]+y+1; */
      /* dom[4] = (x-1)*g->dim[1]+y+1; */
      /* dom[5] = (x-1)*g->dim[1]+y; */
      /* dom[6] = (x-1)*g->dim[1]+y-1; */
      /* dom[7] = (x)*g->dim[1]+y-1; */
      /* dom[8] = (x+1)*g->dim[1]+y-1; */
      int NumberOfNeighbors = 0;
      for(int d=0;d<dd;d++) {
	if(g->Domain[dom[d]].NumberOfParticles !=0) {
	  // Copy all particles contained in the domain
	  memcpy(&tmpi[NumberOfNeighbors],
		 g->Domain[dom[d]].index,
		 sizeof(int)*g->Domain[dom[d]].NumberOfParticles);
	  
	  // How many particles do we have
	  NumberOfNeighbors += g->Domain[dom[d]].NumberOfParticles;
        }
      }
      // Extract the coordinates
      Neighbors[particles] = 0;
      if(NumberOfNeighbors) {
	vzPackV(NumberOfNeighbors, (MKL_Complex16 *)Particles,
                tmpi,
                (MKL_Complex16 *)tvi);
	
	// Set the second vector with the coordinates of particle i
	for(int ps=0;ps<NumberOfNeighbors;ps++) {
	  tvi[ps] -= Particles[particles];
        }
        
	// Calculate the norm
	vzAbs(NumberOfNeighbors, (MKL_Complex16 *)tvi, Norm);
	// Now check which particle is within 2 h of particle i
	const double hi= hsml[particles];
	for(int ns=0;NumberOfNeighbors-ns;ns++) {
	  const double h = hi + hsml[tmpi[ns]];
          if(particles==0) 
            printf("%d %.10lf %.10lf\n", ns, Norm[ns], h); 
	  if((Norm[ns]<=h)&&(particles!=tmpi[ns])) {
	    // We found one
            j_pair[particles][Neighbors[particles]] = tmpi[ns];
	    rij2[particles][Neighbors[particles]] = Norm[ns]*Norm[ns];
	    Neighbors[particles]++;	  
	  }
	}
      }
    }
    
    free(tmpi);
    free(tvi);
    free(Norm);
  }

  for(int ns=0;ns<ntotal;ns++)
    if(Neighbors[ns] > 20)
      printf("%d\n", Neighbors[ns]);
}

void DestroyGrid(Grid *g)
{
  for(int d=0;d<g->Length;d++)
    if(g->Domain[d].index)
      free(g->Domain[d].index);
  free(g->Domain);
  free(g->WhereIsParticle);
  free(g);
}
