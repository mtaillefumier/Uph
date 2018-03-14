/*
***********************************************************
* Input Boundary Particles
* This subprogram generates boundary particles.
***********************************************************
*/
#include "inputbp.h"
#include <omp.h>
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif
void inputbp(System *sys, double t)
{
  FILE *fip, *fop, *fin, *ffp, *ffbp, *fmbp, *fout, *fsp, *ffrp;

  bool boundary_file = true;
  
  char filename_ffp[128] = "FluidParticle.dat";
  char filename_ffbp[128] = "FixedBoundaryParticle.dat";
  char filename_fmbp[128] = "MovingBoundaryParticle.dat";
  char filename_fsp[128] = "SolidParticle.dat";
  char filename_chk[128] = "iocheck.dat";
  char filename_ffrp[128] = "FrozenParticle.dat";
  
  
  /* load boundary particles from file */
  //fixed boundary particle
  ffbp = fopen(filename_ffbp, "r");
  if(ffbp == NULL)
  {
    printf("cannot find fixed boundary file:\n");
    exit(1);
  }
  // moving boundary partilces
  fmbp = fopen(filename_fmbp, "r");
  if(fmbp == NULL)
  {
    printf("cannot find moving boundary file:\n");
    exit(1);
  }
  //fluid particle
  ffp = fopen(filename_ffp, "r");
  if(ffp == NULL)
  {
    printf("cannot find fluid particle file:\n");
    exit(1);
  }
  // solid Particles
  fsp = fopen(filename_fsp, "r");
  if(fsp == NULL)
  {
    printf("cannot find solid particle file:\n");
    exit(1);
  }
  ffrp = fopen(filename_ffrp, "r");
//     // frozen Particles
//     if(sys->flagInitSimul!=1)
//     {
//       ffrp = fopen(filename_ffrp, "r");
//       if(ffrp == NULL)
//       {
// 	printf("cannot find frozen particle file:\n");
//         exit(1);
//       }
//     }
//     else
//     {
//       sys->NumberOfFrozenParticles = 0;
//     }
    
  fscanf(ffp, "%d", &sys->NumberOfFluidParticles);
  fscanf(fsp, "%d", &sys->NumberOfSolidParticles);
  fscanf(ffbp, "%d", &sys->NumberOfFixedBoundaryParticles);
  fscanf(fmbp, "%d", &sys->NumberOfMovingBoundaryParticles);
  if(sys->flagInitSimul!=1)
  {
    fscanf(ffrp, "%d", &sys->NumberOfFrozenParticles); // include frozen particles outside solid boundary 
  }
  else
  {
      sys->NumberOfFrozenParticles = 0;
  }
  sys->NumberOfBoundaryParticles = sys->NumberOfFixedBoundaryParticles + sys->NumberOfMovingBoundaryParticles;
    
  fout = fopen(filename_chk, "w+");
  fprintf(fout, "%d\n", sys->NumberOfSolidParticles);
  // read solid particle 
  for(int i = 0; i < sys->NumberOfSolidParticles; i++)
  {
    fscanf(fsp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf\n", 
	   &sys->Position[sys->NumberOfFluidParticles+i][0], 
	   &sys->Position[sys->NumberOfFluidParticles+i][1], 
	   &sys->Velocity[sys->NumberOfFluidParticles+i][0],
	   &sys->Velocity[sys->NumberOfFluidParticles+i][1],
	   &sys->density[sys->NumberOfFluidParticles+i],
	   &sys->DynamicPressure[sys->NumberOfFluidParticles+i],
	   &sys->mass[sys->NumberOfFluidParticles+i],
           &sys->tempPosition[sys->NumberOfFluidParticles+i],
	   &sys->tempPositionGate[sys->NumberOfFluidParticles+i],
	   &sys->itype[sys->NumberOfFluidParticles+i],
	   &sys->localVelocity[sys->NumberOfFluidParticles+i][0],
	   &sys->localVelocity[sys->NumberOfFluidParticles+i][1]);
    fprintf(fout, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf\n", 
	    sys->Position[sys->NumberOfFluidParticles+i][0], 
	    sys->Position[sys->NumberOfFluidParticles+i][1], 
	    sys->Velocity[sys->NumberOfFluidParticles+i][0],
	    sys->Velocity[sys->NumberOfFluidParticles+i][1],
	    sys->density[sys->NumberOfFluidParticles+i],
	    sys->DynamicPressure[sys->NumberOfFluidParticles+i],
	    sys->mass[sys->NumberOfFluidParticles+i],
            sys->tempPosition[sys->NumberOfFluidParticles+i],
	    sys->tempPositionGate[sys->NumberOfFluidParticles+i],
	    sys->itype[sys->NumberOfFluidParticles+i],
	    sys->localVelocity[sys->NumberOfFluidParticles+i][0],
	    sys->localVelocity[sys->NumberOfFluidParticles+i][1]);
  }
  for(int i = 0; i < sys->NumberOfFixedBoundaryParticles; i++)
  {
    fscanf(ffbp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf\n", 
	   &sys->Position[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+i][0], 
	   &sys->Position[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+i][1], 
	   &sys->Velocity[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+i][0],
	   &sys->Velocity[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+i][1],
	   &sys->density[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+i],
	   &sys->DynamicPressure[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+i],
	   &sys->mass[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+i],
	   &sys->tempPosition[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+i],
	   &sys->tempPositionGate[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+i],
	   &sys->itype[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+i],
	   &sys->localVelocity[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+i][0],
	   &sys->localVelocity[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+i][1]);
  }
  // read moving boundary particle info
  int particleCount = sys->NumberOfFluidParticles + sys->NumberOfSolidParticles + sys->NumberOfFixedBoundaryParticles;
  for(int i = 0; i < sys->NumberOfMovingBoundaryParticles; i++)
  {
    fscanf(fmbp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf\n", 
	   &sys->Position[particleCount+i][0], 
	   &sys->Position[particleCount+i][1], 
	   &sys->Velocity[particleCount+i][0],
	   &sys->Velocity[particleCount+i][1],
	   &sys->density[particleCount+i],
	   &sys->DynamicPressure[particleCount+i],
	   &sys->mass[particleCount+i],
	   &sys->tempPosition[particleCount+i],
	   &sys->tempPositionGate[particleCount+i],
	   &sys->itype[particleCount+i],
	   &sys->localVelocity[particleCount+i][0],
	   &sys->localVelocity[particleCount+i][1]);
  }
  // read frozen particles
  if(sys->flagInitSimul!=1)
  {
    int particleCountFrozen = sys->NumberOfFluidParticles + sys->NumberOfSolidParticles + sys->NumberOfFixedBoundaryParticles + sys->NumberOfMovingBoundaryParticles;
    for(int i = 0; i < sys->NumberOfFrozenParticles; i++)
    {
      fscanf(ffrp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf\n", 
	     &sys->Position[particleCountFrozen+i][0], 
	     &sys->Position[particleCountFrozen+i][1], 
	     &sys->Velocity[particleCountFrozen+i][0],
	     &sys->Velocity[particleCountFrozen+i][1],
	     &sys->density[particleCountFrozen+i],
	     &sys->DynamicPressure[particleCountFrozen+i],
	     &sys->mass[particleCountFrozen+i],
	     &sys->tempPosition[particleCountFrozen+i],
	     &sys->tempPositionGate[particleCountFrozen+i],
	     &sys->itype[particleCountFrozen+i],
	     &sys->localVelocity[particleCountFrozen+i][0],
	     &sys->localVelocity[particleCountFrozen+i][1]);
      fprintf(fout, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf\n", 
	      sys->Position[particleCountFrozen+i][0], 
	      sys->Position[particleCountFrozen+i][1], 
	      sys->Velocity[particleCountFrozen+i][0],
	      sys->Velocity[particleCountFrozen+i][1],
	      sys->density[particleCountFrozen+i],
	      sys->DynamicPressure[particleCountFrozen+i],
	      sys->mass[particleCountFrozen+i],
	      sys->tempPosition[particleCountFrozen+i],
	      sys->tempPositionGate[particleCountFrozen+i],
	      sys->itype[particleCountFrozen+i],
	      sys->localVelocity[particleCountFrozen+i][0],
	      sys->localVelocity[particleCountFrozen+i][1]);
    }
  }
 
 fclose(ffp);
 fclose(ffbp);
 fclose(fmbp);
 fclose(fsp);
 fclose(fout);
 fclose(ffrp);
 
 // implement moving left boundary
 if((sys->movingBoundary==1)&&(sys->flagInitSimul==1))
 {
   sys->LeftBoundaryPosition = 0.5*sys->stroke*(1.0 - cos(sys->omega*t));
   sys->LeftBoundaryVelocity = 0.5*sys->stroke*sys->omega*sin(sys->omega*t);
   sys->LeftBoundaryAcceleration = 0.5*sys->stroke*sys->omega*sys->omega*cos(sys->omega*t);
      
   int sumCount = sys->NumberOfFluidParticles + sys->NumberOfSolidParticles + sys->NumberOfFixedBoundaryParticles;
   for(int i = 0; i < sys->NumberOfMovingBoundaryParticles; i++)
   {
     sys->Position[sumCount+i][0] = sys->LeftBoundaryPosition;
     sys->Velocity[sumCount+i][0] = sys->LeftBoundaryVelocity;
     sys->localVelocity[sumCount+i][0] = sys->LeftBoundaryVelocity; // assign local average velocity to moving boundary particles as well.
   }
 }
 // implement moving top boundary
 else if((sys->movingBoundary==3)&&(sys->flagInitSimul==1))
 {
   int sCounts = sys->NumberOfFluidParticles + sys->NumberOfSolidParticles + sys->NumberOfFixedBoundaryParticles;
   for(int i = 0; i < sys->NumberOfMovingBoundaryParticles; i++)
   {
     sys->tempPosition[sCounts+i] = sys->Position[sCounts+i][0];
   }
   int sumCounts = sys->NumberOfFluidParticles + sys->NumberOfSolidParticles + sys->NumberOfFixedBoundaryParticles;
   for(int i = 0; i < sys->NumberOfMovingBoundaryParticles; i++)
   {
     sys->Position[sumCounts+i][0] = sys->lidVelocity*t + sys->tempPosition[sumCounts+i];
     sys->Velocity[sumCounts+i][0] = sys->lidVelocity;
     sys->localVelocity[sumCounts+i][0] = sys->lidVelocity;
   }
 }
 // implement moving gate as a moving boundary
 else if((sys->movingBoundary==2)&&(sys->flagInitSimul==1))
 {
   double stopTime = 1.5*(sys->LocalDepth-15.0*sys->dx)/sys->lidVelocity; 
   double startTime = 0.0;
   int sumCountGate = sys->NumberOfFluidParticles + sys->NumberOfSolidParticles + sys->NumberOfFixedBoundaryParticles;
   for(int i = 0; i < sys->NumberOfMovingBoundaryParticles; i++)
   {
     if( (t > startTime)&&(t < stopTime))
     {
       sys->Position[sumCountGate+i][1] = sys->lidVelocity*t + sys->tempPositionGate[sumCountGate+i];
       sys->Velocity[sumCountGate+i][1] = sys->lidVelocity;
     }
     else
     {
       sys->Position[sumCountGate + i][1] = sys->lidVelocity*stopTime + sys->tempPositionGate[sumCountGate + i];
       sys->Velocity[sumCountGate + i][1] = 0.0;
     }
   }
 }
 else if(sys->virtualType==5) // sloashing problem
 {
   sys->LeftBoundaryPosition = 0.5*sys->stroke*(1.0 - cos(sys->omega*t));
   sys->LeftBoundaryVelocity = 0.5*sys->stroke*sys->omega*sin(sys->omega*t);
   sys->RightBoundaryPosition = sys->FlatBottomLength + sys->LeftBoundaryPosition;
   sys->LeftBoundaryAcceleration = 0.5*sys->stroke*sys->omega*sys->omega*cos(sys->omega*t);
   int sumCountss = sys->NumberOfFluidParticles + sys->NumberOfSolidParticles + sys->NumberOfFixedBoundaryParticles;

   for(int i = 0; i < sys->NumberOfMovingBoundaryParticles; i++)
   {
     if(sys->itype[sumCountss+i]==2)
     {
       sys->Position[sumCountss+i][0] = sys->LeftBoundaryPosition;
       sys->Velocity[sumCountss+i][0] = sys->LeftBoundaryVelocity;
     }
     else if(sys->itype[sumCountss+i]==3)
     {
       sys->Position[sumCountss+i][0] = sys->RightBoundaryPosition;
       sys->Velocity[sumCountss+i][0] = sys->LeftBoundaryVelocity;
     }
     else if(sys->itype[sumCountss+i]==1)
     {
       sys->Position[sumCountss+i][0] = i*sys->dx*0.5 + sys->LeftBoundaryPosition;
       sys->Velocity[sumCountss+i][0] = sys->LeftBoundaryVelocity;
     }
     else
     {
       sys->LeftBoundaryPosition = 0.0;
       sys->LeftBoundaryVelocity = 0.0;
       sys->LeftBoundaryAcceleration = 0.0;
     }
   }
 }
 else
 {
   // no motion
   sys->LeftBoundaryPosition = 0.0;
   sys->LeftBoundaryVelocity = 0.0;
   sys->LeftBoundaryAcceleration = 0.0;
 }
 // initialize boundary + solid particles
 for(int i = 0; i <(sys->NumberOfBoundaryParticles+sys->NumberOfSolidParticles); i++)
 {
   sys->hsml[sys->NumberOfFluidParticles+i] = sys->eta*sys->dx;
   sys->viscosity[sys->NumberOfFluidParticles+i] = sys->viscosity0;
   sys->localDensity[sys->NumberOfFluidParticles+i] = sys->density0;
   sys->localPressure[sys->NumberOfFluidParticles+i] = 0.0;
//    sys->localVelocity[sys->NumberOfFluidParticles+i][0] = 0.0;
//    sys->localVelocity[sys->NumberOfFluidParticles+i][1] = 0.0;
   sys->normalizationFactor[sys->NumberOfFluidParticles+i] = 1.0;
   sys->tempPointPressure[sys->NumberOfFluidParticles+i] = 0.0;
   sys->tempPointVolume[sys->NumberOfFluidParticles+i] = sys->mass[sys->NumberOfFluidParticles+i]/sys->density0;
   sys->gradNormFactor[sys->NumberOfFluidParticles+i][0] = 0.0;
   sys->gradNormFactor[sys->NumberOfFluidParticles+i][1] = 0.0;
   sys->localCompressibility[sys->NumberOfFluidParticles+i] = sys->DynamicPressure[sys->NumberOfFluidParticles+i] + sys->SpeedOfSound*sys->SpeedOfSound*sys->density0;
   sys->thermalDiffusivity[sys->NumberOfFluidParticles+i] = sys->thermalDiffusivity0;
   sys->specificVolume[sys->NumberOfFluidParticles+i] = 1.0/sys->density0;
   sys->dampingConstant[sys->NumberOfFluidParticles+i] = sys->dampingConstant0;
   sys->vorticity[sys->NumberOfFluidParticles+i] = 0.0;
   sys->approxVorticity[sys->NumberOfFluidParticles+i] = 0.0;
  
   sys->specificDissipationPower[sys->NumberOfFluidParticles+i] = 0.0;
   sys->springConstant[sys->NumberOfFluidParticles+i] = sys->springConstant0;
   sys->specificBoundaryDissipationRate[sys->NumberOfFluidParticles+i] = 0.0;
   sys->specificDissipationRate[sys->NumberOfFluidParticles+i] = 0.0;
   sys->specificTurbulentDissipationRate[sys->NumberOfFluidParticles+i] = 0.0;
   sys->turbulentKineticEnergy[sys->NumberOfFluidParticles+i] = 0.0;
   sys->damGateFlag[sys->NumberOfFluidParticles+i] = 0;
 }
 // initialize frozen particles
 int countFrozen = sys->NumberOfFluidParticles+sys->NumberOfFixedBoundaryParticles+sys->NumberOfSolidParticles+sys->NumberOfMovingBoundaryParticles;
 for(int i = 0; i < sys->NumberOfFrozenParticles; i++)
 {
   sys->hsml[countFrozen+i] = sys->eta*sys->dx;
   sys->viscosity[countFrozen+i] = sys->viscosity0;
   sys->localDensity[countFrozen+i] = sys->density0;
   sys->localPressure[countFrozen+i] = 0.0;
//    sys->localVelocity[countFrozen+i][0] = 0.0;
//    sys->localVelocity[countFrozen+i][1] = 0.0;
   sys->normalizationFactor[countFrozen+i] = 1.0;
   sys->tempPointPressure[countFrozen+i] = 0.0;
   sys->tempPointVolume[countFrozen+i] = sys->mass[countFrozen+i]/sys->density0;
   sys->gradNormFactor[countFrozen+i][0] = 0.0;
   sys->gradNormFactor[countFrozen+i][1] = 0.0;
   sys->localCompressibility[countFrozen+i] = sys->DynamicPressure[countFrozen+i] + sys->SpeedOfSound*sys->SpeedOfSound*sys->density0;
   sys->thermalDiffusivity[countFrozen+i] = sys->thermalDiffusivity0;
   sys->specificVolume[countFrozen+i] = 1.0/sys->density[countFrozen+i];
   sys->dampingConstant[countFrozen+i] = sys->dampingConstant0;
   sys->vorticity[countFrozen+i] = 0.0;
   sys->approxVorticity[countFrozen+i] = 0.0;
  
   sys->specificDissipationPower[countFrozen+i] = 0.0;
   sys->springConstant[countFrozen+i] = sys->springConstant0;
   sys->specificBoundaryDissipationRate[countFrozen+i] = 0.0;
   sys->specificDissipationRate[countFrozen+i] = 0.0;
   sys->specificTurbulentDissipationRate[countFrozen+i] = 0.0;
   sys->turbulentKineticEnergy[countFrozen+i] = 0.0;
   sys->damGateFlag[countFrozen+i] = 0;
 }
 // external boundary force
//  for(int i = sys->NumberOfFluidParticles; i < (sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+sys->NumberOfFixedBoundaryParticles+sys->NumberOfFrozenParticles); i++)
//  {
//    sys->acceleration[i][0] = 0.0;
//    sys->acceleration[i][1] = 0.0;
//  }
}

