/* 
*****************************************************************************
* subprogram to write data to output file
* writes output in text file with 8 columns
* containing particle number i, 
* position vector components x, y
* velocity field components u, v, 
* particle density density,
* pressure p, and 
* mass m.  
*****************************************************************************
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "param.h"
#include "output.h"

void output(System *sys, int itime)
{
  FILE *fip, *fke, *fpp, *fpo, *finit, *fppu, *fpwr, *fpcv, *fkcv;
  char filename_xv[128];
  // verify this fprintf statement
  memset(filename_xv, 0, sizeof(char)*128);
//   sprintf(filename_xv, "out%09d.dat", itime);
  sprintf(filename_xv, "out%d.dat", itime);
  
  int nsum = sys->NumberOfFluidParticles + sys->NumberOfBoundaryParticles + sys->NumberOfVirtualParticles+sys->NumberOfSolidParticles + sys->NumberOfFrozenParticles; /* fluid + boundary + virtual particles */
  fip = fopen(filename_xv,"w+");
  
  // open and append to file for Kinetic Energy-time 
  fke = fopen("TotalKineticEnergy.dat", "a");
  // open and append to file for PointPressure-time 
  fpp = fopen("PointPressure.dat", "a");
  // open and append to file for FluidParticlesInit-time 
  finit = fopen("FluidParticleInit.dat", "a");
  // open and append to file for PointPressureTwo-time 
  fppu = fopen("PointPressureTwo.dat", "a");
  // open and append to file for Power-time 
  fpwr = fopen("PowerCalc.dat", "a");
   // open and append to file for control volume Power-time 
  fpcv = fopen("ControlVolPower.dat", "a");
   // open and append to file for control volume Power-time 
  fkcv = fopen("ControlVolKinetic.dat", "a");
  
  for(int i = 0; i < nsum; i++)
    {
      
      fprintf(fip, "%d ",i);
      for(int d = 0; d < sys->dim; d++)
	{
	  fprintf(fip, "\t%.10lf", sys->Position[i][d]);
	  
	}
      for(int d = 0; d < sys->dim; d++)
	{
	  fprintf(fip, "\t%.10lf", sys->Velocity[i][d]);
	}
      fprintf(fip, "\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%d\t%.10lf\t%d\t%.10lf\t%.10lf\n", sys->density[i], sys->mass[i], sys->hsml[i], sys->DynamicPressure[i], sys->damGateFlag[i], sys->approxVorticity[i], sys->itype[i], sys->specificDissipationRate[i], sys->localDensity[i]);
    }
    if(!fke)
    {
      perror("error opening file");
      exit(1);
    }
    else
    {
      fprintf(fke, "%.10lf\t%.10lf %.10lf\n", sys->timeStep*itime, sys->totalKineticEnergy, sys->potentialEnergy);
    }
    if(!fpp)
    {
      perror("error opening file");
      exit(1);
    }
    else
    {
      fprintf(fpp, "%.10lf\t%.10lf\t%d\t%.10lf\t%.10lf\n", sys->timeStep*itime, sys->pointPressure, sys->pointNeighbors, sys->pressureNorm, sys->pointPressureTwo);
    }
    if(!fppu)
    {
      perror("error opening file: PointPressureTwo.dat");
      exit(1);
    }
    else
    {
      fprintf(fppu, "%.10lf\t%d\t%.10lf\n", sys->timeStep*itime, sys->pointNeighborsTwo, sys->pointPressureTwo);
    }
    if(!fpwr)
    {
      perror("error opening file: PowerCalc.dat");
      exit(1);
    }
    else
    {
      fprintf(fpwr, "%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", sys->timeStep*itime, sys->enstrophyPower, sys->controlDissipationRate[1], sys->boundaryDissipationRate, sys->dissipationRate);
    }
    if(!fpcv)
    {
      perror("error opening file: ControlVolPower.dat");
      exit(1);
    }
    else
    {
      fprintf(fpcv, "%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", sys->timeStep*itime, sys->controlDissipationRate[1], sys->controlDissipationRate[2], sys->controlDissipationRate[3], sys->controlDissipationRate[4]);
    }
    if(!fkcv)
    {
      perror("error opening file: ControlVolPower.dat");
      exit(1);
    }
    else
    {
      fprintf(fkcv, "%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", sys->timeStep*itime, sys->controlKineticEnergy[1], sys->controlKineticEnergy[2], sys->controlKineticEnergy[3], sys->controlKineticEnergy[4]);
    }
    if(itime==sys->NumberOfInitializationSteps)
    {
      fprintf(finit, "%d\n", sys->NumberOfFluidParticles);
      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
      {
	for(int d = 0; d<sys->dim; d++)
        {
	  fprintf(finit, "%.10lf\t",sys->Position[i][d]);
        }
        for(int d = 0; d<sys->dim; d++)
        {
	  fprintf(finit, "%.10lf\t",sys->Velocity[i][d]);
	}
	fprintf(finit, "%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%d\t%.10lf\t%.10lf\n", sys->density[i], sys->DynamicPressure[i],sys->mass[i], sys->tempPosition[i],sys->tempPositionGate[i], sys->itype[i],sys->localVelocity[i][0],sys->localVelocity[i][1]);
      }
    }
     
  fclose(fip);	
  fclose(fke);
  fclose(fpp);
  fclose(finit);
  fclose(fppu);
  fclose(fpwr);
  fclose(fpcv);
  fclose(fkcv);
}
