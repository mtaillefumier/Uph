#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>
#include <string.h>
#include "param.h"
#include "System.h"

double **CreateMatrix(int x, int y)
{
  double **a = (double **)malloc(sizeof(double *)*x);
  posix_memalign(a, 16, sizeof(double)*x*y);

  for(int i=1;i<x;i++)
    a[i] = a[i-1] + y;
  return a;
}


int **CreateMatrixInt(int x, int y)
{
  int **a = (int **)malloc(sizeof(int *)*x);
  posix_memalign(a, 16, sizeof(int)*x*y);

  for(int i=1;i<x;i++)
    a[i] = a[i-1] + y;
  return a;
}

void DestroyMatrix(double **a)
{
  free(a[0]);
  free(a);
}


System *CreateSystem(const int HowManyNeighbors,
		     const int numReps)
{
  char filename_param[128] = "parameters.dat";
  char filename_fp[128] = "FluidParticle.dat";
  char filename_ffbp[128] = "FixedBoundaryParticle.dat";
  char filename_fmbp[128] = "MovingBoundaryParticle.dat";
  char filename_fsp[128] = "SolidParticle.dat";
  char filename_ffrp[128] = "FrozenParticle.dat";

  FILE *fparam, *ffp, *ffbp, *fmbp, *fsp, *ffrp;
  bool generation_mode = true;

  System *a = (System *)malloc(sizeof(System));
  memset(a, 0, sizeof(System));

  fparam = fopen(filename_param, "r");
  if(fparam == NULL)
  {
    printf("cannot find parameter file parameters.dat\n");
    exit(1);
  }
  fscanf(fparam, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %d %d %d %lf %lf %d %d %lf %lf %lf %d %d %lf %d %d %lf %lf %d %d",
	 &a->stroke,
         &a->wavePeriod,
         &a->LocalDepth,
         &a->FlatBottomLength,
         &a->density0,
         &a->AdiabaticConstant,
         &a->Gravity,
         &a->dx,
         &a->dy,
         &a->dim,
         &a->viscosity0,
         &a->machNumber,
         &a->pecletNumber,
	 &a->epsilon,
         &a->alpha,
	 &a->virtualType,
	 &a->printStep,
	 &a->screenStep,
	 &a->timeStep,
	 &a->lidVelocity,
	 &a->nearNeighbor,
	 &a->selectBoundaryForce,
	 &a->lennardSelect,
	 &a->springConstant0, // 24
	 &a->dampingConstant0,
	 &a->flagInitSimul,
	 &a->gravityFlag,
	 &a->thermalDiffusivity0,
	 &a->movingBoundary,
	 &a->NumberOfInitializationSteps,
	 &a->pressureProbeX,
	 &a->pressureProbeY,
	 &a->stressType,
	 &a->disspationControlVolume);
  fclose(fparam);

//   printf("%lf \n%lf \n%lf \n%lf \n%lf \n%lf \n%lf \n%lf \n%lf \n%d \n%lf \n%lf \n%lf \n%lf \n%lf \n%d \n%d \n%d \n%lf \n%lf \n%d \n%d \n%lf \n%lf \n%lf \n%d \n%d \n%lf \n%d \n%d \n%lf \n%lf \n%lf \n%lf \n%lf \n%lf \n%d",
// 	 a->stroke,
//          a->wavePeriod,
//          a->LocalDepth,
//          a->FlatBottomLength,
//          a->density0,
//          a->AdiabaticConstant,
//          a->Gravity,
//          a->dx,
//          a->dy,
//          a->dim,
//          a->viscosity0,
//          a->machNumber,
//          a->pecletNumber,
// 	 a->epsilon,
//          a->alpha,
// 	 a->virtualType,
// 	 a->printStep,
// 	 a->screenStep,
// 	 a->timeStep,
// 	 a->lidVelocity,
// 	 a->nearNeighbor,
// 	 a->selectBoundaryForce,
// 	 a->lennardSelect,
// 	 a->springConstant, // 24
// 	 a->dampingConstant0,
// 	 a->flagInitSimul,
// 	 a->gravityFlag,
// 	 a->thermalDiffusivity0,
// 	 a->movingBoundary,
// 	 a->NumberOfInitializationSteps,
// 	 a->pressureProbeX,
// 	 a->pressureProbeY,
// 	 a->pressureProbeXX,
// 	 a->pressureProbeYY,
// 	 a->spongeFactor0,
// 	 a->spongePosition0,
// 	 a->waveAbsorption);
//   exit(1);
  // read parameters from files .dat

  // number of fluid particles
  ffp = fopen(filename_fp, "r");
  if(ffp == NULL)
  {
    printf("cannot find fluid particle file: \n");
    exit(1);
  }
  fscanf(ffp, "%d", &a->NumberOfFluidParticles);
  fclose(ffp);

  // number of boundary particles
  ffbp = fopen(filename_ffbp, "r");
  if(ffp == NULL)
  {
    printf("cannot find fixed boundary particle file:\n");
    exit(1);
  }
  fscanf(ffbp, "%d", &a->NumberOfFixedBoundaryParticles);
  fclose(ffbp);

  //moving boundary
  fmbp = fopen(filename_fmbp, "r");
  if(fmbp == NULL)
  {
    printf("cannot find moving boundary file:\n");
    exit(1);
  }
  fscanf(fmbp, "%d", &a->NumberOfMovingBoundaryParticles);
  fclose(fmbp);

  // solid particles
  fsp = fopen(filename_fsp, "r");
  if(fsp == NULL)
  {
    printf("cannot find solid particle file");
    exit(1);
  }
  fscanf(fsp, "%d", &a->NumberOfSolidParticles);
  fclose(fsp);

  if(a->flagInitSimul!=1)
  {
    ffrp = fopen(filename_ffrp, "r");
    if(ffrp==NULL)
    {
      printf("cannot find frozen particle file:\n");
      exit(1);
    }

    fscanf(ffrp, "%d", &a->NumberOfFrozenParticles);
    fclose(ffrp);
  }
  else
  {
    a->NumberOfFrozenParticles = 0;
  }
  a->NumberOfBoundaryParticles = a->NumberOfFixedBoundaryParticles + a->NumberOfMovingBoundaryParticles;
  a->MaxNumberOfParticles =  a->NumberOfSolidParticles + a->NumberOfBoundaryParticles + 100*a->NumberOfFluidParticles + a->NumberOfFrozenParticles;


  a->Position = CreateMatrix(a->MaxNumberOfParticles, a->dim);
  a->Position_min = CreateMatrix(a->MaxNumberOfParticles, a->dim);
  a->Velocity = CreateMatrix(a->MaxNumberOfParticles, a->dim);
  a->Velocity_xsph = CreateMatrix(a->MaxNumberOfParticles, a->dim);
  a->Velocity_min = CreateMatrix(a->MaxNumberOfParticles, a->dim);
  a->dvdt = CreateMatrix(a->MaxNumberOfParticles, a->dim);
  a->av = CreateMatrix(a->MaxNumberOfParticles, a->dim);
  a->velocityFluctuation = CreateMatrix(a->MaxNumberOfParticles, a->dim);
  a->localVelocity = CreateMatrix(a->MaxNumberOfParticles, a->dim);
  a->dwcdx = CreateMatrix(a->MaxNumberOfParticles, a->dim);
  a->gradNormFactor = CreateMatrix(a->MaxNumberOfParticles, a->dim);
  a->tmp = CreateMatrix(64, a->dim*a->MaxNumberOfParticles);
  a->dwdx = CreateMatrix(a->MaxNumberOfParticles, a->dim);

  a->dvdt_approx = CreateMatrix(a->MaxNumberOfParticles, a->dim);


  posix_memalign((void **)&a->Energy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->hsml, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->Pressure, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->itype, 16, sizeof(int)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->density, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->density_min, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->ddensitydt, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->mass, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->i_pair, 16, sizeof(int)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->j_pair, 16, sizeof(int)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->rij2, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->w, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->i_inflow, 16, sizeof(int)*a->MaxNumberOfParticles);

  posix_memalign((void **)&a->viscosity, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->ViscousStressxx, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->ViscousStressxy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->ViscousStressyx, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->ViscousStressyy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->DynamicPressure, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->DynamicPressure_min, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->dpdt, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->normalizationFactor, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->localDensity, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->localPressure, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->wc, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->wp, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->tempPosition, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->tempPositionGate, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->tempPointPressure, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->tempPointVolume, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->pj_pair, 16, sizeof(int)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->localCompressibility, 16, sizeof(int)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->thermalDiffusivity, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->specificVolume, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->NeighborsSize, 16, sizeof(int)*a->NumberOfFluidParticles);
  posix_memalign((void **)&a->Volume, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->densityFluctuation, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->pressureFluctuation, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->localDensity_min, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->meanCompressibility, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->vorticity, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->approxVorticity, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->tempPointPressureTwo, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->xSpongeVelocity, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->ySpongeVelocity, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->spongeFactor, 16, sizeof(double)*a->MaxNumberOfParticles);


  posix_memalign((void **)&a->specificDissipationPower, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->springConstant, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->dampingConstant, 16, sizeof(double)*a->MaxNumberOfParticles);

  posix_memalign((void **)&a->specificBoundaryDissipationRate, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->specificDissipationRate, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->specificDissipationRateTemp, 16, sizeof(double)*a->MaxNumberOfParticles);

  posix_memalign((void **)&a->specificTurbulentDissipationRate, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->turbulentKineticEnergy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->controlKineticEnergy, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->controlDissipationRate, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->controlVolumeX, 16, sizeof(double)*a->MaxNumberOfParticles);

  posix_memalign((void **)&a->damGateFlag, 16, sizeof(int)*a->MaxNumberOfParticles);

  posix_memalign((void **)&a->ddensitydt_approx, 16, sizeof(double)*a->MaxNumberOfParticles);
  posix_memalign((void **)&a->dpdt_approx, 16, sizeof(double)*a->MaxNumberOfParticles);


  a->alpha_cspline = 0;
  a->alpha_wendlandFour = 0;
  a->alpha_kernelFour = 0;
  a->alpha_deconvolutionFour = 0;
  a->alpha_kernelFive = 0;
  a->alpha_deconvolutionFive = 0;
  a->alpha_deconvolutionSix = 0;

  if(a->movingBoundary==1)
  {
    a->SpeedOfSound = sqrt(a->Gravity*a->LocalDepth)/a->machNumber;
  }
  else if(a->movingBoundary==2)
  {
    a->SpeedOfSound = sqrt(2.0*a->Gravity*a->LocalDepth)/a->machNumber;
  }
  else
  {
    a->SpeedOfSound = sqrt(a->Gravity*a->LocalDepth)/a->machNumber;
  }
  a->CompressionFactor = (1./a->AdiabaticConstant) * a->density0 * a->SpeedOfSound * a->SpeedOfSound;

  /***********************************************/
  /* wavemaker */
  a->omega = 2.0*M_PI/a->wavePeriod;   //2.445;
  a->slope = 1.0/15.0;

  a->eta = sqrt(a->nearNeighbor/(4.0*M_PI));
  a->LennardCoefficient = a->Gravity*a->LocalDepth;
  a->coreCoefficient = a->Gravity*a->LocalDepth;

  a->NumberOfInteractingParticles = 0;
  a->pointNeighbors = 0;
  a->totalKineticEnergy = 0.0;
  a->pointPressure = 0.0;
  a->pressureNorm = 1.0;
  a->consistencyDensityOne = a->density0;
  a->consistencyDensityTwo = a->density0;
  a->consistencyDensityThree = a->density0;

  a->enstrophyPower = 0.0;
  a->dissipationPower = 0.0;

  a->dissipationRate = 0.0;
  a->boundaryDissipationRate = 0.0;

  a->controlVolumeNumber = 5; // number of control volumes

  // for sponge layer of wave tank
  a->spongeFactor0 = 0.2;
  a->spongePosition0 = 7.0;
  a->waveAbsorption = 1;
 /*
   printf("finished reading here:\n");
    exit(1);
  */

  // We only need to calculate the neighbors of the particles in the box
  a->HowManyNeighbors = HowManyNeighbors;

  return a;
}

void DestroySystem(System *a)
{
  DestroyMatrix(a->Position);
  DestroyMatrix(a->Velocity);
  DestroyMatrix(a->Velocity_xsph);
  DestroyMatrix(a->Velocity_min);
  DestroyMatrix(a->dvdt);
  DestroyMatrix(a->av);
  DestroyMatrix(a->velocityFluctuation);
  DestroyMatrix(a->localVelocity);
  DestroyMatrix(a->dwcdx);
  DestroyMatrix(a->gradNormFactor);
  DestroyMatrix(a->tmp);
  DestroyMatrix(a->Position_min);
  DestroyMatrix(a->dwdx);

  DestroyMatrix(a->dvdt_approx);

  free(a->Energy);
  free(a->hsml);
  free(a->Pressure);
  free(a->itype);
  free(a->density);
  free(a->density_min);
  free(a->ddensitydt);
  free(a->mass);
  free(a->i_pair);
  free(a->j_pair);
  free(a->rij2);
  free(a->w);
  //  free(a->q.mat);
  //  free(a->x.mat);
  free(a->viscosity);
  free(a->ViscousStressxx);
  free(a->ViscousStressxy);
  free(a->ViscousStressyx);
  free(a->ViscousStressyy);
  free(a->DynamicPressure);
  free(a->DynamicPressure_min);
  free(a->dpdt);
  free(a->normalizationFactor);
  free(a->localDensity);
  free(a->localPressure);
  free(a->wc);
  free(a->wp);
  free(a->tempPosition);
  free(a->tempPositionGate);
  free(a->tempPointPressure);
  free(a->pj_pair);
  free(a->tempPointVolume);
  free(a->localCompressibility);
  free(a->thermalDiffusivity);
  free(a->specificVolume);
  free(a->NeighborsSize);
  free(a->Volume);
  free(a->densityFluctuation);
  free(a->pressureFluctuation);
  free(a->localDensity_min);
  free(a->meanCompressibility);
  free(a->i_inflow);
  free(a->vorticity);
  free(a->approxVorticity);
  free(a->tempPointPressureTwo);
  free(a->xSpongeVelocity);
  free(a->ySpongeVelocity);
  free(a->spongeFactor);

  free(a->specificDissipationPower);
  free(a->dampingConstant);
  free(a->springConstant);

  free(a->specificBoundaryDissipationRate);
  free(a->specificDissipationRate);
  free(a->specificDissipationRateTemp);

  free(a->specificTurbulentDissipationRate);
  free(a->turbulentKineticEnergy);

  free(a->controlDissipationRate);
  free(a->controlKineticEnergy);
  free(a->controlVolumeX);

  free(a->damGateFlag);

  free(a->ddensitydt_approx);
  free(a->dpdt_approx);

  free(a);
}
