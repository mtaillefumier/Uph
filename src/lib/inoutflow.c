/*
****************************************************************
* This subprogram handles inflow and outflow of particles
* inflow particles are generated at inflow boundary and 
* added to the existing arrays. At the outflow boundary, 
* particles are removed. Only horizontal velocity component
* is possible at inflow.
****************************************************************
*/
#include <stdio.h>
#include "param.h"
#include "System.h"

void inoutflow(System *sys)
{

  int ii;
  double inputmass, inputvol;

  /*Left inflow*/
  ii = sys->NumberOfFluidParticles;
  inputvol = 0.0;
  inputmass = 0.0;
  for(int k = 0; k < sys->inflowParticles; k++)
    {
      int i = sys->i_inflow[k];
      if(sys->Position[i][0] > 5000.1)
	{
	  sys->Position[ii][0] = sys->Position[i][0] - 0.1;
	  sys->Position[ii][1] = sys->Position[i][1];
	  sys->Velocity[ii][0] = sys->Velocity[i][0];
	  sys->Velocity[ii][1] = 0.0;
	  sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
	  sys->Velocity_xsph[ii][1] = 0.0;
	  sys->density[ii] = sys->density[i];
	  sys->mass[ii] = sys->mass[i];
	  sys->Pressure[ii] = sys->Pressure[i];
	  sys->Energy[ii] = sys->Energy[i];
	  sys->hsml[ii] = sys->hsml[i];
	  sys->density_min[ii] = sys->density[i];
	  sys->Velocity_min[ii][0] = sys->Velocity[i][0];
	  sys->Velocity_min[ii][1] = 0.0;
	  /*Track input discharge*/
	  inputmass += sys->mass[ii];
	  inputvol += sys->mass[ii]/sys->density[ii];
	  ++ii;
	}
    }

  /* if(inputmass != 0) */
  /*   { */
      printf("input mass: %lf kg/m\t and input vol: %lf m^2\n", inputmass, inputvol);
    /* } */

  /*Right boundary*/
  /* review this because it fails badly */ 

  int NumberOfFluidParticles = ii;
  for(int i = 0; i < NumberOfFluidParticles; i++)
    {
      if(sys->Position[i][0] > 500000.0)
	{
	  sys->Position[i][0] = sys->Position[NumberOfFluidParticles-1][0];
	  sys->Position[i][1] = sys->Position[NumberOfFluidParticles-1][1];
	  sys->Velocity[i][0] = sys->Velocity[NumberOfFluidParticles-1][0];
	  sys->Velocity[i][1] = sys->Velocity[NumberOfFluidParticles-1][1];
	  sys->Velocity_xsph[i][0] = sys->Velocity_xsph[NumberOfFluidParticles-1][0];
	  sys->Velocity_xsph[i][1] = sys->Velocity_xsph[NumberOfFluidParticles-1][1];
	  sys->mass[i] = sys->mass[NumberOfFluidParticles-1];
	  sys->density[i] = sys->density[NumberOfFluidParticles-1];
	  sys->Pressure[i] = sys->Pressure[NumberOfFluidParticles-1];
	  sys->Energy[i] = sys->Energy[NumberOfFluidParticles-1];
	  sys->hsml[i] = sys->hsml[NumberOfFluidParticles-1];
	  sys->density_min[i] = sys->density_min[NumberOfFluidParticles-1];
	  sys->Velocity_min[i][0] = sys->Velocity_min[NumberOfFluidParticles-1][0];
	  sys->Velocity_min[i][1] = sys->Velocity_min[NumberOfFluidParticles-1][1];
	  NumberOfFluidParticles--;
	}
    }
}	
