#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include <malloc.h>
#include <stdlib.h>

//#include "Neighbors.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct __system {
    int dim; // Dimension of the system
    int *itype;
    int nx;
    int ny;
    double Lx, Ly;
    double dx, dy;
    int NumberOfInteractingParticles;
  int MaxNumberOfParticles;
  int NumberOfFluidParticles; /* total number of fluid particles */
  int NumberOfBoundaryParticles;  /* number of boundary particles */
  int NumberOfFixedBoundaryParticles;
  int NumberOfMovingBoundaryParticles;
  int NumberOfFrozenParticles;
  int NumberOfSolidParticles;
  int NumberOfVirtualParticles;   /* number of boundary particles */
  int nearNeighbor; /* average number of particles inside support of kernel*/

  double LeftBoundaryPosition, LeftBoundaryVelocity, RightBoundaryPosition, LeftBoundaryAcceleration;
  double **tmp;
  double **Position; /* initialize the position vector */
  double **Velocity; /* initialize the velocity vector */
  double **Velocity_xsph; /* move particles with this average velocity */
  double *mass; /* initialize particle mass */
  double *Pressure; /* initialize particle pressure */
  double *Energy; /* initialize specific energy */
  double *hsml; /* initialize smoothing length h */
  double *density; /* initialize particle density */
  double **Velocity_min;
  double *ddensitydt;
  double **dvdt;
  double **av;
  double *density_min;
  double density0;
  double *viscosity;
  double *ViscousStressxx,*ViscousStressxy,*ViscousStressyx,*ViscousStressyy;
  double coreCoefficient;
  double machNumber;
  double *DynamicPressure, *DynamicPressure_min, *dpdt;
  double *localDensity_min;
  double *meanCompressibility, *meanSpecDensity;
  double pressureKernelNorm;


  double *densityFluctuation, *pressureFluctuation;
  double **velocityFluctuation;
  double *normalizationFactor;// integrate the test function over space.
  double damping_time, tstep;
  int dampingIterations;

  //new variables 25.8.2016
  double *localDensity, *localPressure, **localVelocity;
  double *wc, **dwcdx;
  double timeStep, pecletNumber, viscosity0;
  int virtualType;
  double *tempPosition, *tempPositionGate, gateVelocity;
  int printStep, screenStep;
  double lennardSelect, lidVelocity, wallViscosity;
  double totalKineticEnergy;
  int pointNeighbors;
  double pointPressure, *tempPointPressure;
  double initialEnergy, potentialEnergy;
  //temps
  double *wp;
  int selectBoundaryForce;
  double springConstant0, dampingConstant0;
  double *springConstant, *dampingConstant;
  int *pj_pair;
  double *tempPointVolume;
  double **gradNormFactor;
  int flagInitSimul;
  double *localCompressibility;
  int gravityFlag;
  double *thermalDiffusivity, thermalDiffusivity0;
  int movingBoundary;
  int NumberOfInitializationSteps;
  double consistencyDensityOne, consistencyDensityTwo, consistencyDensityThree;
  double densityUncertainty;
  double *specificVolume;
  double **Position_min;
  double pressureNorm;// integrate kernel about pressure probe point.
  double *Volume;
  double *vorticity, *approxVorticity;
  double pressureProbeX, pressureProbeY;
  double pressureProbeXX, pressureProbeYY;
  double *tempPointPressureTwo;
  int pointNeighborsTwo;
  double pointPressureTwo;
  double *xSpongeVelocity, *ySpongeVelocity; // not used anywhere
  double *spongeFactor, spongeFactor0, spongePosition0;
  int waveAbsorption; // 1... to activate the sponge layer for absorbing reflected waves.


  double enstrophyPower;
  int stressType;
  double *specificDissipationPower, dissipationPower;

  double *specificBoundaryDissipationRate, boundaryDissipationRate;
  double *specificDissipationRate, dissipationRate, *specificDissipationRateTemp;
  int disspationControlVolume;
  double *specificTurbulentDissipationRate, * turbulentKineticEnergy;
  double *controlKineticEnergy, *controlDissipationRate, *controlMaxVelocity, *controlVolumeX;
  int controlVolumeNumber;

  int *damGateFlag;

  double *ddensitydt_approx;
  double *dpdt_approx;
  double **dvdt_approx;

  // temporary variables used for the calculations
  double alpha_cspline, alpha_kernelFour, alpha_deconvolutionFour, alpha_kernelFive, alpha_deconvolutionFive;
  double alpha_wendlandFour, alpha_deconvolutionSix;
  double alpha;
  double SpeedOfSound;
  double slope;
  double stroke, omega, wavePeriod; /* stroke and frequency of wavemaker*/
  double eta; /* smoothing scale*/
  double Gravity; // acceleration due to gravity
  double LennardCoefficient;
  double LocalDepth, FlatBottomLength; /* deep water depth*/
  double epsilon; /* xSPH dispersion factor*/
  int *i_pair, *j_pair;
  int *NeighborsSize;
  int *i_inflow;
  double *rij2, *w;
  double **dwdx;
    //int **Neighbors;
    ///double **DistanceNeighbors;
    //Grid *g;
  //  matrix r;
  //  int numReps;
  //  rep *ri;
  //int NumberOfRepresentatives;
  int inflowParticles;
  double CompressionFactor;
  double AdiabaticConstant;
  int niac;
  int HowManyNeighbors;
  //  matrix x,q;
} System;

double **CreateMatrix(int x, int y);
void DestroyMatrix(double **a);
System *CreateSystem(const int HowManyNeighbors,
		     const int numReps);
void DestroySystem(System *a);

#ifdef __cplusplus
}
#endif

#endif
