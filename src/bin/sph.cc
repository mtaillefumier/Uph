/* ******************************************************
 * This is the main program.
 * It calls " input.c" once, starts time loop with
 * "derivatives.c" and leap frog time integration every
 * time step. Calls "output.c" if requested and then
 * generates output on the screen.
 * ****************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <inttypes.h>
#include "param.h"
#include "input.h"
#include "derivatives.h"
#include "output.h"
#include "System.h"

int main(int argc, char**argv)
{
  System *sys = CreateSystem(128, 128);
  double dt = sys->timeStep;//1e-5; /* set timestep */
  unsigned int maxtimestep = 10; /* set maximum timestep */

  int itime; /* iteration time */

   char filename_ke[128] = "TotalKineticEnergy.dat";
   char filename_pp[128] = "PointPressure.dat";
   char filename_ppu[128] = "PointPressureTwo.dat";
   char filename_init[128] = "FluidParticleInit.dat";
   char filename_pwr[128] = "PowerCalc.dat";
   char filename_pcv[128] = "ControlVolPower.dat";
   char filename_kcv[128] = "ControlVolKinetic.dat";
   FILE *fke, *fpp, *finit, *fppu, *fpwr, *fpcv, *fkcv;

  if(argc == 1)
    {
      printf("Give a maximum number of time step\n");
      exit(1);
    }
  maxtimestep= atoi(argv[1]);

  // parameters for water
//   System *sys = CreateSystem(dim, 0.0, 1000, 7, dt, 32, 128);

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  input(sys);

  // delete file contents for TotalKineticEnergy.dat
  fke = fopen(filename_ke, "w");
  fclose(fke);
  // delete file contents for PointPressure.dat
  fpp = fopen(filename_pp, "w");
  fclose(fpp);

   // delete file contents for FluidParticleInit.dat
  finit = fopen(filename_init, "w");
  fclose(finit);

  // delete file contents for PointPressureTwo.dat
  fppu = fopen(filename_ppu, "w");
  fclose(fppu);

  // delete file contents for PowerCalc.dat
  fpwr = fopen(filename_pwr, "w");
  fclose(fpwr);
   // delete file contents for ControlVolPower.dat
  fpcv = fopen(filename_pcv, "w");
  fclose(fpcv);
  // delete file contents for ControlVolPower.dat
  fkcv = fopen(filename_kcv, "w");
  fclose(fkcv);

  memset(sys->dvdt[0], 0, sizeof(double)*sys->MaxNumberOfParticles*sys->dim);
  memset(sys->ddensitydt, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memcpy(sys->density_min, sys->density, sizeof(double)*sys->NumberOfFluidParticles);
  memset(sys->dpdt, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memcpy(sys->DynamicPressure_min, sys->DynamicPressure, sizeof(double)*sys->NumberOfFluidParticles);

  memset(sys->dvdt_approx[0], 0, sizeof(double)*sys->MaxNumberOfParticles*sys->dim);
  memset(sys->ddensitydt_approx, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->dpdt_approx, 0, sizeof(double)*sys->MaxNumberOfParticles);


  for(itime = 1; itime <= maxtimestep; itime++)
    {
      double t = itime * dt;
      if(itime!=1) {
      	memcpy(sys->density_min, sys->density, sizeof(double)*sys->NumberOfFluidParticles);
	memcpy(sys->DynamicPressure_min, sys->DynamicPressure, sizeof(double)*sys->NumberOfFluidParticles);
	memcpy(sys->localDensity, sys->density, sizeof(double)*sys->NumberOfFluidParticles);
	memcpy(sys->localPressure, sys->DynamicPressure, sizeof(double)*sys->NumberOfFluidParticles);
	memcpy(sys->localVelocity[0], sys->Velocity[0], sizeof(double)*sys->dim*sys->MaxNumberOfParticles);
	memcpy(sys->Velocity_min[0], sys->Velocity[0], sizeof(double)*sys->dim*sys->MaxNumberOfParticles);

// 	memcpy(sys->Position_min[0], sys->Position[0], sizeof(double)*sys->dim*sys->MaxNumberOfParticles);

	/* cblas_daxpy(sys->NumberOfFluidParticles, 0.5*dt, sys->ddensitydt, 1, sys->density, 1); */
	/* memcpy(sys->Velocity_min[0], sys->Velocity[0], sizeof(double)*2*sys->NumberOfFluidParticles); */
	/* cblas_daxpy(sys->NumberOfFluidParticles, 0.5*dt, sys->dvdt[0], 1, sys->Velocity_xsph[0], 1); */
	/* cblas_daxpy(sys->NumberOfFluidParticles, 0.5*dt, sys->dvdt[0], 1, sys->Velocity[0], 1); */
// #pragma omp parallel for
      	for(int i=0;i< (sys->NumberOfFluidParticles) ;i++)
	{
      	  sys->density[i] += 0.5*dt*sys->ddensitydt_approx[i];
	  sys->localDensity[i] = sys->localDensity[i] + 0.5*dt*sys->ddensitydt_approx[i] + sys->densityFluctuation[i];
	  sys->DynamicPressure[i] += 0.5*dt*sys->dpdt_approx[i];
	  sys->localPressure[i] = sys->localPressure[i] + 0.5*dt*sys->dpdt_approx[i] + sys->pressureFluctuation[i];
	}

	for(int i=0;i<sys->NumberOfFluidParticles ;i++)
	{
      	  for(int d=0;d<sys->dim;d++)
	  {
	    sys->Velocity_min[i][d] = sys->Velocity[i][d];
      	    sys->Velocity_xsph[i][d] += 0.5 * dt * sys->dvdt[i][d];
      	    sys->Velocity[i][d] +=  0.5 * dt * sys->dvdt_approx[i][d];
	    sys->localVelocity[i][d] = sys->Velocity[i][d] + 0.5 * dt * sys->dvdt_approx[i][d] + sys->epsilon*sys->velocityFluctuation[i][d];
      	  }
      	}
      }

      if(itime==1) {
	memcpy(sys->Velocity_xsph[0], sys->Velocity[0], sizeof(double)*sys->dim*sys->MaxNumberOfParticles);
// 	memcpy(sys->Velocity_min[0], sys->Velocity[0], sizeof(double)*sys->dim*sys->MaxNumberOfParticles);
// 	memcpy(sys->localPressure, sys->DynamicPressure, sizeof(double)*sys->NumberOfFluidParticles);
// 	memcpy(sys->localVelocity[0], sys->Velocity[0], sizeof(double)*sys->dim*sys->MaxNumberOfParticles);
      }

      derivatives(sys, t);

      if(itime==1) {
      	for(int i=0;i< (sys->NumberOfFluidParticles) ;i++) {
	  sys->localDensity[i] = sys->density[i] + 0.5*dt*sys->ddensitydt[i] + sys->densityFluctuation[i];
      	  sys->density[i] += 0.5*dt*sys->ddensitydt[i];
	  sys->localPressure[i] = sys->DynamicPressure[i] + 0.5*dt*sys->dpdt[i] + sys->pressureFluctuation[i];
	  sys->DynamicPressure[i] += 0.5*dt*sys->dpdt[i];
	}
	for(int i=0;i<sys->NumberOfFluidParticles ;i++) {
      	  for(int d=0;d<sys->dim;d++) {
      	    sys->Velocity_xsph[i][d] = sys->Velocity[i][d] + 0.5 * dt * sys->dvdt[i][d] + sys->epsilon*sys->velocityFluctuation[i][d];
	    sys->localVelocity[i][d] = sys->Velocity[i][d] + 0.5 * dt * sys->dvdt[i][d] + sys->epsilon*sys->velocityFluctuation[i][d];
      	    sys->Velocity[i][d] +=  0.5 * dt * sys->dvdt[i][d];
      	    sys->Position[i][d] +=  dt * sys->Velocity[i][d];//sys->Velocity_xsph[i][d];
	  }
	}
      }
      else {
      	for(int i=0;i< (sys->NumberOfFluidParticles) ;i++) {
	  sys->localDensity[i] = sys->density_min[i] + dt*sys->ddensitydt[i] + sys->densityFluctuation[i];
      	  sys->density[i] = sys->density_min[i] + dt*sys->ddensitydt[i];
	  sys->localPressure[i] = sys->DynamicPressure_min[i] + dt*sys->dpdt[i] + sys->pressureFluctuation[i];
	  sys->DynamicPressure[i] = sys->DynamicPressure_min[i] + dt*sys->dpdt[i];
	}
	for(int i=0;i<sys->NumberOfFluidParticles ;i++) {
      	  for(int d=0;d<sys->dim;d++) {
      	    sys->Velocity_xsph[i][d] = sys->Velocity_min[i][d] + dt * sys->dvdt[i][d] + sys->epsilon*sys->velocityFluctuation[i][d];
	    sys->localVelocity[i][d] =  sys->Velocity_min[i][d] + dt * sys->dvdt[i][d] + sys->epsilon*sys->velocityFluctuation[i][d];
      	    sys->Velocity[i][d] =  sys->Velocity_min[i][d] + dt * sys->dvdt[i][d];
      	    sys->Position[i][d] += dt * sys->Velocity[i][d];//sys->Velocity_xsph[i][d];
	  }
	}
      }
      if( (itime % sys->printStep) == 0)
	{
	  printf("%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
		 sys->Velocity[0][0],
		 sys->Velocity[0][1],
		 sys->dvdt[0][0],
		 sys->dvdt[0][1],
		 sys->Pressure[0],
		 sys->density[0],
		 sys->ddensitydt[0]);
	  output(sys, itime);
  	}
      if( (itime % sys->screenStep) == 0)
  	{
  	  /* write useful info to the console */
  	  printf("******************************************\n");

  	  /* verify %lld using with mingw's gcc */
  	  printf("\nTimestep %d of  %d\n", itime, maxtimestep);
  	  printf("Interaction pairs: %d\n", sys->NumberOfInteractingParticles);
          printf("m %f\t %f \t%f\t%lf\n", sys->mass[0], sys->hsml[0], sys->initialEnergy, sys->SpeedOfSound);
  	  printf("******************************************\n");
  	}
    }

  end = clock();
  cpu_time_used = (double)(end - start)/CLOCKS_PER_SEC;
  printf("CPU time %lf s\n", cpu_time_used);

  DestroySystem(sys);
  return 0;
}
