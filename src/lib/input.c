/*
*********************************************************************
******* Input file *********
* This is the subprogram for the initial input
* Generates or reads from files all the physical 
* properties of fluid particles only					 
*********************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <malloc.h>
#include "param.h"
#include "System.h"
#include <omp.h>
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif
void input(System *sys)
{
  FILE *fip, *fop, *fin, *ffp, *fout, *fsp;
  double space;
  int tmp = 0;
  bool input_file = true;
   char filename[128] = "FluidParticle.dat";
   char filename_chk[128] = "iocheckfluid.dat";
   char filename_fsp[128] = "SolidParticle.dat";
  /* load initial fluid particles from file */
  
  if(!filename)
    exit(1);
  ffp = fopen(filename, "r");
  if(ffp == NULL)
  {
    printf("cannot find fluid particle input file:\n");
    exit(1);
  }
  fsp = fopen(filename_fsp, "r");
  if(fsp == NULL)
  {
    printf("cannot find solid particle file :\n");
    exit(1);
  }
  sys->NumberOfSolidParticles = 0;
  sys->NumberOfFluidParticles = 0;
  fscanf(ffp, "%d", &sys->NumberOfFluidParticles);
  fscanf(fsp, "%d", &sys->NumberOfSolidParticles);
   
  fout = fopen(filename_chk, "w+");
  fprintf(fout, "%d\n", sys->NumberOfSolidParticles);
//   printf("number of solid = %d\n", sys->NumberOfSolidParticles);
//   exit(1);  
  
  for(int i = 0; i<sys->NumberOfFluidParticles; i++)
  {
    fscanf(ffp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf", 
	   &sys->Position[i][0],
	   &sys->Position[i][1],
	   &sys->Velocity[i][0],
	   &sys->Velocity[i][1],
	   &sys->density[i],
	   &sys->DynamicPressure[i],
	   &sys->mass[i],
	   &sys->tempPosition[i],
	   &sys->tempPositionGate[i],
           &sys->itype[i],
	   &sys->localVelocity[i][0],
	   &sys->localVelocity[i][1]);
  }
//     // read solid particle 
//   for(int i = 0; i < sys->NumberOfSolidParticles; i++)
//   {
//     fscanf(fsp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf\n", 
// 	   &sys->Position[sys->NumberOfFluidParticles+i][0], 
// 	   &sys->Position[sys->NumberOfFluidParticles+i][1], 
// 	   &sys->Velocity[sys->NumberOfFluidParticles+i][0],
// 	   &sys->Velocity[sys->NumberOfFluidParticles+i][1],
// 	   &sys->density[sys->NumberOfFluidParticles+i],
// 	   &sys->DynamicPressure[sys->NumberOfFluidParticles+i],
// 	   &sys->mass[sys->NumberOfFluidParticles+i],
//            &sys->tempPosition[sys->NumberOfFluidParticles+i],
// 	   &sys->tempPositionGate[sys->NumberOfFluidParticles+i],
// 	   &sys->itype[sys->NumberOfFluidParticles+i],
// 	   &sys->localVelocity[sys->NumberOfFluidParticles+i][0],
// 	   &sys->localVelocity[sys->NumberOfFluidParticles+i][1]);
//     fprintf(fout, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf\n", 
// 	    sys->Position[sys->NumberOfFluidParticles+i][0], 
// 	    sys->Position[sys->NumberOfFluidParticles+i][1], 
// 	    sys->Velocity[sys->NumberOfFluidParticles+i][0],
// 	    sys->Velocity[sys->NumberOfFluidParticles+i][1],
// 	    sys->density[sys->NumberOfFluidParticles+i],
// 	    sys->DynamicPressure[sys->NumberOfFluidParticles+i],
// 	    sys->mass[sys->NumberOfFluidParticles+i],
//             sys->tempPosition[sys->NumberOfFluidParticles+i],
// 	    sys->tempPositionGate[sys->NumberOfFluidParticles+i],
// 	    sys->itype[sys->NumberOfFluidParticles+i],
// 	    sys->localVelocity[sys->NumberOfFluidParticles+i][0],
// 	    sys->localVelocity[sys->NumberOfFluidParticles+i][1]);
//   }
//   printf("finished reading solidparticles.dat:\n");
//   exit(1);
  for(int i = 0; i < (sys->NumberOfFluidParticles+sys->NumberOfSolidParticles); i++)
  {
    sys->hsml[i] = sys->eta*sys->dx;
    sys->viscosity[i] = sys->viscosity0;
    sys->localDensity[i] = sys->density0;
    sys->localPressure[i] = 0.0;
//     sys->localVelocity[i][0] = 0.0;
//     sys->localVelocity[i][1] = 0.0;
    sys->normalizationFactor[i] = 1.0;
    sys->tempPointPressure[i] = 0.0;
    sys->tempPointVolume[i] = sys->mass[i]/sys->density0;
    sys->gradNormFactor[i][0] = 0.0;
    sys->gradNormFactor[i][1] = 0.0;
    sys->localCompressibility[i] = sys->DynamicPressure[i] + sys->SpeedOfSound*sys->SpeedOfSound*sys->density0;
    sys->thermalDiffusivity[i] = sys->thermalDiffusivity0;
    sys->specificVolume[i] = 1.0/sys->density0;
    sys->dampingConstant[i] = 0.0;
    sys->vorticity[i] = 0.0;
    sys->approxVorticity[i] = 0.0;
    
    sys->specificDissipationPower[i] = 0.0;
    sys->springConstant[i] = 0.0;
    sys->specificBoundaryDissipationRate[i] = 0.0;
    sys->specificDissipationRate[i] = 0.0;
    sys->specificTurbulentDissipationRate[i] = 0.0;
    sys->turbulentKineticEnergy[i] = 0.0;
  }

  for(int i = 0; i < sys->NumberOfFluidParticles; i++)
  {
    sys->initialEnergy += sys->mass[i]*sys->Gravity*sys->Position[i][1];
    //assign damping coeffienent for wave tank sponge layer
    if(sys->Position[i][0]>=sys->spongePosition0)
    {
      sys->spongeFactor[i] = sys->spongeFactor0;
    }
    else
    {
      sys->spongeFactor[i] = 0.0;
    }
  }
  for(int i = 0; i < sys->NumberOfFluidParticles; i++)
  {
    double widthDam = 0.38;
    double heightColumn = 0.15;
    //assign tag on fluid particles in two compartments of gated dam;
    // for particles in dam column tag=1 and for particles on wet bed tag=-1
    if(sys->Position[i][0]<=widthDam)
    {
      sys->damGateFlag[i] = 1;
    }
    else
    {
      sys->damGateFlag[i] = -1;
    }
  }
  fclose(ffp);
  fclose(fsp);
  fclose(fout);
}