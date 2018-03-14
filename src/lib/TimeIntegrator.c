#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "System.h"
#include "derivatives.h"
#include "input.h"
#include "TimeIntegrator.h"

void TimeIntegrator(System *sys, double t, double dt, int itime)
{
  memset(sys->dvdt[0], 0, sizeof(double)*sys->MaxNumberOfParticles*sys->dim);
  memcpy(sys->Velocity_min[0], sys->Velocity[0], sizeof(double)*sys->NumberOfFluidParticles*sys->dim);
  memset(sys->ddensitydt, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memcpy(sys->density_min, sys->density, sizeof(double)*sys->NumberOfFluidParticles);
  memset(sys->dpdt, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memcpy(sys->DynamicPressure_min, sys->DynamicPressure, sizeof(double)*sys->NumberOfFluidParticles);
  memcpy(sys->Position_min[0], sys->Position[0], sizeof(double)*sys->NumberOfFluidParticles*sys->dim);
  
  if(itime!=1)
  {
    memcpy(sys->Velocity_min[0], sys->Velocity[0], sizeof(double)*sys->NumberOfFluidParticles*sys->dim);
    memcpy(sys->density_min, sys->density, sizeof(double)*sys->NumberOfFluidParticles);
    memcpy(sys->DynamicPressure_min, sys->DynamicPressure, sizeof(double)*sys->NumberOfFluidParticles);
    memcpy(sys->Position_min[0], sys->Position[0], sizeof(double)*sys->NumberOfFluidParticles*sys->dim);
      
    for(int i = 0; i < sys->NumberOfFluidParticles; i++)
    {
      sys->density[i] += 0.5*sys->ddensitydt[i]*dt;
      sys->DynamicPressure[i] += 0.5*sys->dpdt[i]*dt;
      for(int d = 0; d < sys->dim; d++)
      {
	sys->Velocity[i][d] += 0.5*sys->dvdt[i][d]*dt;
// 	sys->Position[i][d] += 0.5*sys->Velocity[i][d]*dt;
      }
    }
    
    derivatives(sys, t);
    
    if(itime==1)
    {
      memcpy(sys->Velocity_min[0], sys->Velocity[0], sizeof(double)*sys->NumberOfFluidParticles*sys->dim);
      memcpy(sys->density_min, sys->density, sizeof(double)*sys->NumberOfFluidParticles);
      memcpy(sys->DynamicPressure_min, sys->DynamicPressure, sizeof(double)*sys->NumberOfFluidParticles);
      memcpy(sys->Position_min[0], sys->Position[0], sizeof(double)*sys->NumberOfFluidParticles*sys->dim);
      
      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
      {
	sys->density[i] += 0.5*sys->ddensitydt[i]*dt;
        sys->DynamicPressure[i] += 0.5*sys->dpdt[i]*dt;
	for(int d = 0; d < sys->dim; d++)
	{
	  sys->Velocity[i][d] += 0.5*sys->dvdt[i][d]*dt;
	  sys->Position[i][d] += sys->Velocity[i][d]*dt;
	}
      }
    }
    
    else
    {
      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
      {
	sys->density[i] = sys->density_min[i] + sys->ddensitydt[i]*dt;
	sys->DynamicPressure[i] = sys->DynamicPressure_min[i] + sys->dpdt[i]*dt;
	for(int d = 0; d < sys->dim; d++)
	{
	  sys->Velocity[i][d]  = sys->Velocity_min[i][d] + sys->dvdt[i][d]*dt;
	  sys->Position[i][d] += sys->Velocity[i][d]*dt;
	}
      }
    }
  }
}


