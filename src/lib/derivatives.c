/* ****************************************************
 * derivative solver
 * the hydrodynamic equation are solved every time step
 * virtual particles are also generated within 2h away
 * from the boundary.
 * makes function calls to inoutflow.h and inputbp.h
 * interaction pairs, kernel W & grad(W) and ddensity/dt
 * are found.
 * use equation of state to compute pressure p(density).
 * compute pressure term, gravity, artificial viscous
 * force and boundary forces.
 * compute dv/dt.
 * influence of average velocity calculated with xSPH.
 * at first timestep a function call to parameterfile.h
 * is made.
 *******************************************************
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "param.h"
#include "inoutflow.h"
#include "inputbp.h"
#include "parameterfile.h"
#include "derivatives.h"
#include <mkl.h>
#include <omp.h>
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

extern   int SearchNeighbors(const double **Position,
                             int *i_pair,
                             int *j_pair,
                             double *rij2,
                             const int NumberOfFluidParticles,
                             const int NumberOfParticles,
                             int ParticleBoundary,
                             const int NumberOfNeighbors,
                             const double dist);

__inline__ void Kernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  sys->alpha_cspline = 10.0/(7.*M_PI*h*h);
  const double q = r/h;
  // testing this again. Verifying how this thing works.
  /*define cubic spline*/
  if(q >= 0.0 && q < 1.0)
    {
      double q2 = q*q;
      double q3 = q*q*q;
      sys->w[niac] = sys->alpha_cspline*(1- 1.5*q2 + 0.75*q3);
      for(int d = 0; d < sys->dim; d++)
        sys->dwdx[niac][d] = sys->alpha_cspline*(-3 + 2.25*q)*dx[d]/(h*h);
    }
  else if(q >= 1.0 && q < 2.0)
    {
      double qp2 = (2-q)*(2-q);
      double qp3 = qp2 * (2-q);
      sys->w[niac] = sys->alpha_cspline*0.25*qp3;
      for(int d = 0; d < sys->dim; d++)
        sys->dwdx[niac][d] = -sys->alpha_cspline*0.75*qp2*dx[d]/(h*r);
    }
}


__inline__ void WendlandFourKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  sys->alpha_wendlandFour = 7.0/(4.0*M_PI*h*h);
  const double q = r/h;
  if(q >= 0.0 && q <= 2.0)
    {
      double q3 = (1-0.5*q)*(1-0.5*q)*(1-0.5*q);
      double q4 = (1-0.5*q)*q3;
      double q2 = q*q;
      sys->w[niac] = sys->alpha_wendlandFour*q4*(1 + 2*q);
      for(int d = 0; d < sys->dim; d++)
        sys->dwdx[niac][d] = -sys->alpha_wendlandFour*5.0*q3*dx[d]/(h*h);
    }
}


__inline__ void SpikyKernel(System *sys, const double *__restrict__ dx, const double r, const double hb, const int niac)
{
  const double alpha_0 = 4.0/(M_PI*hb*hb);
  const double q = r/hb;
  if(q >= 0.0 && q <= 1.0)
    {
      double q2 = q*q;
      double q3 = (1-q2)*(1-q2)*(1-q2);
      sys->w[niac] = alpha_0*q3;
    }
}

__inline__ void NewKernelFour(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  sys->alpha_kernelFour = 6144.0/(2.0*M_PI*h*h*(428.0+75.0*log(5)));//384.0/(2.0*M_PI*h*h*(68.0-15.0*log(5)));
  const double q = r/h;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double pp4 = (1-0.25*q2)*(1-0.25*q2)*(1-0.25*q2)*(1-0.25*q2);
      double pp5 = (1-0.25*q2)*pp4;
      double qq4 = (1+q2)*(1+q2)*(1+q2)*(1+q2);
      double qq5 = (1+q2)*qq4;

      sys->w[niac] = sys->alpha_kernelFour*pp5/qq5;
    }
}
__inline__ void DeconvolutionKernelFour(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  sys->alpha_deconvolutionFour = 384.0/(2.0*M_PI*h*h*(68.0-15.0*log(5)));
  const double q = r/h;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q4 = q2*q2;
      double q6 = q2*q4;
      double q8 = q4*q4;
      double q10 = q4*q6;
      double q12 = q2*q10;
      double q14 = q2*q12;

      double pp1 = (1-0.25*q2);
      double pp2 = (1-0.25*q2)*(1-0.25*q2);
      double pp3 = (1-0.25*q2)*(1-0.25*q2)*(1-0.25*q2);
      double pp4 = (1-0.25*q2)*(1-0.25*q2)*(1-0.25*q2)*(1-0.25*q2);

      double qq4 = (1+q2)*(1+q2)*(1+q2)*(1+q2);
      double qq5 = (1+q2)*(1+q2)*(1+q2)*(1+q2)*(1+q2);
      double qq6 = (1+q2)*(1+q2)*(1+q2)*(1+q2)*(1+q2)*(1+q2);
      double qq7 = (1.0+q2)*qq6;
      double qq8 = qq4*qq4;
      double qq9 = qq4*qq5;

      sys->wc[niac] = sys->alpha_deconvolutionFour*(  pp4/qq4 + 5*pp2*(4-20*q2+q4)/qq6  - 5*(2*q10 - 137*q8+1784*q6-6506*q4+4328*q2-368)/(8*qq8));
      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = -sys->alpha_deconvolutionFour*5.0*pp1*(62.0*q8 - 373.0*q6 + 102.0*q4 + 1246.0*q2 + 1084.0)*dx[d]/(496.0*qq7*h*h);
    }
}

__inline__ void NewKernelFive(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_kernelFive = 63.0/(2.0*M_PI*h*h*4864.0);
  const double q = r/h;
  if(q >= 0.0 && q < 2.0)
    {
      double q2 = q*q;
      double q4 = q2*q2;
      double q6 = q2*q4;
      double q8 = q4*q4;

      double pp5 = (1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2);
      double pp6 = pp5*(1.0-0.25*q2);
      double qq5 = (1+q2)*(1+q2)*(1+q2)*(1+q2)*(1+q2);
      double qq = (1+q2);
      double rr = 1.0/(1.0-0.25*q2);
      double pp = (2.0-q)*(2.0-q)*(2.0-q)*(2.0-q)*(2.0-q);
      double pp4 = (2.0-q)*(2.0-q)*(2.0-q)*(2.0-q);

      sys->w[niac] = alpha_kernelFive*pp*(2.0+q)*(5.0*q2+8.0*q+4.0);
    }
}
// __inline__ void DeconvolutionKernelFive(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
// {
//   const double a = 666409.0/14400000.0;
//   const double b = 151.0*log(5.0)/6144.0;
//
//   sys->alpha_deconvolutionFive = 1.0/(2.0*M_PI*h*h*a*b);
//   const double q = r/h;
//   if(q >= 0.0 && q < 2.0)
//   {
//     double q2 = q*q;
//     double q4 = q2*q2;
//     double q6 = q2*q4;
//     double q8 = q4*q4;
//     double q10 = q4*q6;
//     double q12 = q2*q10;
//     double q14 = q2*q12;
//
//
//     double pp2 = (4.0-q2);
//
//     double qq4 = (1+q2)*(1+q2)*(1+q2)*(1+q2);
//     double qq5 = (1+q2)*(1+q2)*(1+q2)*(1+q2)*(1+q2);
//     double qq6 = (1+q2)*(1+q2)*(1+q2)*(1+q2)*(1+q2)*(1+q2);
//     double qq8 = qq4*qq4;
//     double qq10 = qq4*qq6;
//     for(int d = 0; d <2; d++)
//       sys->dwcdx[niac][d] = -sys->alpha_deconvolutionFive*pp2*(151.0*q14 - 2000.0*q12 + 59802.0*q10 -1083588.0*q8 + 13303143.0*q6 - 60418212.0*q4 + 49681328.0*q2 -6995776.0)*dx[d]/(1536.0*qq10*h*h);
//   }
// }

__inline__ void DeconvolutionKernelFive(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{

  const double alpha_deconvolutionFive = 63.0/(2.0*M_PI*h*h*15200.0);
  const double q = r/h;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q3 = q*q2;
      double q4 = q2*q2;
      double q5 = q*q4;
      double q6 = q2*q4;
      double q8 = q4*q4;
      double q10 = q4*q6;
      double q12 = q2*q10;
      double q14 = q2*q12;
      double q16 = q2*q14;
      double q18 = q2*q16;

      double pp = (2.0-q)*(2.0-q);
      double pp2 = (1.0-0.25*q2)*(1.0-0.25*q2);
      double pp4 = pp2*pp2;

      double qq4 = (1+q2)*(1+q2)*(1+q2)*(1+q2);
      double qq5 = (1+q2)*(1+q2)*(1+q2)*(1+q2)*(1+q2);
      double qq6 = (1+q2)*(1+q2)*(1+q2)*(1+q2)*(1+q2)*(1+q2);
      double qq7 = (1+q2)*qq6;
      double qq11 = qq4*qq7;
      sys->wc[niac] = alpha_deconvolutionFive*pp2*(q8-6.0*q6+81.0*q4-1576.0*q2+336)/qq6;
      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = -alpha_deconvolutionFive*pp*(125.0*q4 - 200.0*q3 -562.0*q2 + 601.0*q + 616.0)*dx[d]/(h*h);


    }
}

__inline__ void DeconvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double a = 51.0/(20.0*M_PI);
  const double b = -51.0/(16.0*M_PI);
  const double c = (153.0*sin(log(2.0)) - 51.0*cos(log(2.0)) )/(160.0*M_PI);
  const double d = (51.0*sin(log(2.0)) + 153.0*cos(log(2.0)) )/(160.0*M_PI);


  const double q = r/h;
  const double eta = 0.01;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;

      double pp = (b/sqrt(q2+eta) + (2.0*c -d)*sin(log(sqrt(q2+eta))) +  (c + 2.0*d)*cos(log(sqrt(q2+eta))) )/(h*h);
      //     double pp = (-1.01461/sqrt(q2+eta) - 0.0660907*sin(log(sqrt(q2+eta))) +  0.714389*cos(log(sqrt(q2+eta))) )/(h*h);
      double pp1 = a + b*sqrt(q2+eta) + ( c*sin(log(sqrt(q2+eta))) +  d*cos(log(sqrt(q2+eta))) )*q2;
      sys->wc[niac] = pp1/(h*h);
      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = pp*dx[d]/(h*h);


    }
}

__inline__ void ConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double a = 51.0/(20.0*M_PI);
  const double b = -51.0/(16.0*M_PI);
  const double c = (153.0*sin(log(2.0)) - 51.0*cos(log(2.0)) )/(160.0*M_PI);
  const double d = (51.0*sin(log(2.0)) + 153.0*cos(log(2.0)) )/(160.0*M_PI);

  const double q = r/h;
  const double eta = 0.0001;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;

      double pp1 = a + b*sqrt(q2+eta) + ( c*sin(log(sqrt(q2+eta))) +  d*cos(log(sqrt(q2+eta))) )*q2;
      double pp2 = b/sqrt(q2+eta) + (3.0*c-4.0*d)*sin(log(sqrt(q2+eta))) + (4.0*c+3.0*d)*cos(log(sqrt(q2+eta)));

      sys->w[niac] = (pp1 + 0.0001*pp2)/(h*h);
      //     double pp1 = -1.01461/sqrt(q2+eta) - 0.84657*sin(log(sqrt(q2+eta))) +  1.36269*cos(log(sqrt(q2+eta)));
      //     double pp2 = (0.1164415*sin(log(sqrt(q2+eta))) +  0.29897370*cos(log(sqrt(q2+eta))))*q2;

      //     for(int d = 0; d <2; d++)
      //       sys->dwcdx[niac][d] = (pp1 + 0.002*pp2)/(h*h);


    }
}

__inline__ void DeconvolutionBiharmonicKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 13.0/(4.0*M_PI*h*h);
  const double c = (153.0*sin(log(2.0)) - 51.0*cos(log(2.0)) )/(160.0*M_PI);
  const double d = (51.0*sin(log(2.0)) + 153.0*cos(log(2.0)) )/(160.0*M_PI);


  const double q = r/h;
  const double eta = 0.01;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      const double varAngle = 0.5*sqrt(3.0)*log(0.5*sqrt(q2+eta));

      double pp = 1.0 - 1.5*sqrt(q2+eta) + (0.5*sqrt(2))*pow( (q2+eta),0.75)*cos(varAngle);
      double pp1 = 0.25*(6.0/sqrt(q2+eta) + sqrt(6.0)*sin(varAngle)/sqrt(sqrt(q2+eta)) - 3.0*sqrt(2.0)*cos(varAngle)/sqrt(sqrt(q2+eta)) );
      sys->wc[niac] = alpha_d*pp;
      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = -alpha_d*pp1*dx[d]/(h*h);


    }
}

__inline__ void ConvolutionBiharmonicKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 13.0/(4.0*M_PI*h*h);
  const double q = r/h;
  const double eta = 0.01;
  const double sensitivity = 0.0001;
  const double coef = 0.75*sqrt(2.0);
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      const double varAngle = 0.5*sqrt(3.0)*log(0.5*sqrt(q2+eta));

      double pp = 1.0 - 1.5*sqrt(q2+eta) + (0.5*sqrt(2))*pow( (q2+eta),0.75)*cos(varAngle);
      double pp2 = coef*(2.0*sin(M_PI/6.0 - varAngle)/sqrt(sqrt(q2+eta)) - sqrt(2.0)/sqrt(q2+eta) );
      double pp1 = 0.25*(6.0/sqrt(q2+eta) + sqrt(6.0)*sin(varAngle)/sqrt(sqrt(q2+eta)) - 3.0*sqrt(2.0)*cos(varAngle)/sqrt(sqrt(q2+eta)) );
      sys->w[niac] = alpha_d*(pp + 0.0*pp2);
      for(int d = 0; d <2; d++)
        sys->dwdx[niac][d] = -alpha_d*pp1*dx[d]/(h*h); // calculate this explicitly


    }
}

__inline__ void LaguerreConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double a0 = 1789404352395344836431623090496.0;
  const double a1 = 572016843554228538028893635008.0;
  const double a2 = 169492307036913903133165008832.0;
  const double a3 = 45205536946879158116113072000.0;
  const double a4 = 10343263090792678269949005100.0;
  const double a5 = 1873321557584287305505054428.0;
  const double a6 = 241247762612017318320324154.0;
  const double a7 = 18798675063995431100502216.0;
  const double a8 = 1059051972436766222690475.0;


  const double normFactor = 24988331289645499024381757625.0/(2*M_PI*84312023993772304141135719107584.0*1397948603616531413951427000.0*h*h);
  const double q = r/h;
  const double eta = 0.0001;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q3 = q*q2;
      double q4 = q*q3;
      double q5 = q*q4;
      double q6 = q*q5;
      double q7 = q*q6;
      double q8 = q*q7;
      double qq4 = (2.0-q)*(2.0-q)*(2.0-q)*(2.0-q);

      sys->w[niac] = normFactor*qq4*(a0 + a1*q + a2*q2 + a3*q3 + a4*q4 + a5*q5 + a6*q6 + a7*q7 + a8*q8);
    }
}

__inline__ void LaguerreDeconvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double a0 = -13312013392139164637528161773686.0;
  const double a1 = 2486209096774867189750717977700.0;
  const double a2 = 1242337862582363668825595512179.0;
  const double a3 = 202839840974935303660981669115.0;
  const double a4 = 57706008503803233528483178229.0;
  const double a5 = 40488652063638442885361523720.0;
  const double a6 = 11231919301882542142240551529.0;
  const double a7 = 1580750718516424206331455589.0;
  const double a8 =164423346806479093197905376.0;
  const double a9 = 12708623669241194672285700.0;

  const double normFactor = 24988331289645499024381757625.0/(2*M_PI*84312023993772304141135719107584.0*1397948603616531413951427000.0*h*h);
  const double q = r/h;
  const double eta = 0.01;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q3 = q*q2;
      double q4 = q*q3;
      double q5 = q*q4;
      double q6 = q*q5;
      double q7 = q*q6;
      double q8 = q*q7;
      double q9 = q*q8;
      double qq2 = (2.0-q)*(2.0-q);

      //     sys->wc[niac] = alpha_d*pp;
      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = normFactor*qq2*(a0 + a1*q + a2*q2 + a3*q3 + a4*q4 + a5*q5 + a6*q6 + a7*q7 + a8*q8 + a9*q9)*dx[d]/(h*h);


    }
}

// __inline__ void ApproximationGaussianConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
// {
//   const double normFactor =  4096/((-3*(-108 + 5*log(5)))*2*M_PI*h*h);
//   const double q = r/h;
//   const double eta = 0.0001;
//   if(q >= 0.0 && q <= 2.0)
//   {
//     double q2 = q*q;
//     double q3 = q*q2;
//     double q4 = q*q3;
//     double qq6 = (1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2);
//     double pp6 = (1.0+q2)*(1.0+q2)*(1.0+q2)*(1.0+q2)*(1.0+q2)*(1.0+q2);
//
//     sys->w[niac] = normFactor*qq6/pp6;
//   }
// }

// __inline__ void ApproximateGaussianDeconvolutionKernelGradient(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
// {
//   const double normFactor =  4096/((-3*(-108 + 5*log(5)))*2*M_PI*h*h);
//   const double q = r/h;
//   const double eta = 0.0001;
//   if(q >= 0.0 && q <= 2.0)
//   {
//     double q2 = q*q;
//     double q3 = q*q2;
//     double q4 = q*q3;
//     double q5 = q*q4;
//     double q6 = q*q5;
//     double q8 = q2*q6;
//     double q10 = q2*q8;
//     double q12 = q2*q10;
//     double q14 = q2*q12;
//     double q16 = q2*q14;
//
//     double num1 = (1 - 0.25*q2);
//     double num3 = (1 - 0.25*q2)*(1 - 0.25*q2)* (1 - 0.25*q2);
//     double num5 = (1 - 0.25*q2)*(1 - 0.25*q2)* (1 - 0.25*q2)* (1 - 0.25*q2)* (1 - 0.25*q2);
//     double mult3 = q6 - 60*q4 + 507*q2 - 132;
//     double mult1 = q12 -147*q10 + 4245*q8 - 3980*q6 + 126630*q4 - 82104*q2 + 9568;
//     double denom7 = (1 + q2)*(1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2);
//     double denom9 = (1 + q2)*(1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2);
//     double denom3 = (1 + q2)*(1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2);
//     double denom11 = (1 + q2)*(1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2)* (1 + q2);
//     double sumNum = 15*num5/denom7 -15*num3*mult3/((2*2*M_PI*normFactor)*denom9) + 45*num1*mult1/(2*4*(2*M_PI*normFactor)*(2*M_PI*normFactor)*denom11);
//
//     double qq4 = (-4 + q2)*(-4 + q2)* (-4 + q2)* (-4 + q2);
//     double pp5 = (1 + q2)*(1 + q2)* (1 + q2)* (1 + q2)* (1 + q2);
//     double pp11 = pp5*pp5*(1 + q2);
//     double sums = 15*(-4 + q2)*(49312 - 207048 *q2 + 396650 *q4 - 128381 *q6 + 2835 *q8 + 127 *q10 + 2567 *q12 - 770 *q14 + 64 *q16)/(65536*pp11) ;
//     for(int d = 0; d<sys->dim; d++)
// //       sys->dwcdx[niac][d] = normFactor*sums*dx[d]/(h*h);
//       sys->dwcdx[niac][d] = -normFactor*sumNum*dx[d]/(h*h);
//
//
//   }
// }

// __inline__ void OrderFiveApproximateGaussianDeconvolutionKernelGradient(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
// {
//   const double normFactor =  1/( (log(1048576)-41/3)*2*M_PI*h*h);   // order five kernel is very stable.
//   const double q = r/h;
//   const double eta = 0.0001;
//   const double coeff1 = 1/(2*M_PI*normFactor);
//   if(q >= 0.0 && q <= 2.0)
//   {
//     double q2 = q*q;
//     double q3 = q*q2;
//     double q4 = q*q3;
//     double q5 = q*q4;
//     double q6 = q*q5;
//     double q8 = q2*q6;
//     double q10 = q2*q8;
//     double q12 = q2*q10;
//     double q14 = q2*q12;
//     double q16 = q2*q14;
//
//     double num1 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)* (1 - 0.25*q2);
//     double num2 = (1 - 0.25*q2)*(1 - 0.25*q2)*(q6-80*q4 +848*q2-640);
//
//     double denom1 = (1 + 0.25*q2)*(1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2);
//     double denom2 = (1 + 0.25*q2)*(1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2);
//
//     double sumNum = 5*num1/denom1 -5*num2*coeff1/(32*denom2);
//
//     for(int d = 0; d<sys->dim; d++)
// //       sys->dwcdx[niac][d] = normFactor*sums*dx[d]/(h*h);
//       sys->dwcdx[niac][d] = -normFactor*sumNum*dx[d]/(h*h);
//
//
//   }
// }

// __inline__ void OrderFiveApproximationGaussianConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
// {
//   const double normFactor =  1/( (log(1048576)-41/3)*2*M_PI*h*h);
//   const double q = r/h;
//   const double eta = 0.0001;
//   if(q >= 0.0 && q <= 2.0)
//   {
//     double q2 = q*q;
//     double qq5 = (1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2);
//     double pp5 = (1.0+0.25*q2)*(1.0+0.25*q2)*(1.0+0.25*q2)*(1.0+0.25*q2)*(1.0+0.25*q2);
//
//     sys->w[niac] = normFactor*qq5/pp5;
//   }
// }

// __inline__ void OrderSixApproximationGaussianConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
// {
//   const double normFactor =  1/( (log(268435456)-289/15)*2*M_PI*h*h);
//   const double q = r/h;
//   const double eta = 0.0001;
//   if(q >= 0.0 && q <= 2.0)
//   {
//     double q2 = q*q;
//     double qq7 = (1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2)*(1.0-0.25*q2);
//     double pp7 = (1.0+0.25*q2)*(1.0+0.25*q2)*(1.0+0.25*q2)*(1.0+0.25*q2)*(1.0+0.25*q2)*(1.0+0.25*q2)*(1.0+0.25*q2);
//
//     sys->w[niac] = normFactor*qq7/pp7;
//   }
// }

// __inline__ void OrderSixApproximateGaussianDeconvolutionKernelGradient(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
// {
//   const double normFactor =  1/( (log(268435456)-289/15)*2*M_PI*h*h);
//   const double q = r/h;
//   const double eta = 0.0001;
//   const double coeff1 = 1/(2*M_PI*normFactor);
//   const double coeff2 = coeff1*coeff1;
//   if(q >= 0.0 && q <= 2.0)
//   {
//     double q2 = q*q;
//     double q3 = q*q2;
//     double q4 = q*q3;
//     double q5 = q*q4;
//     double q6 = q*q5;
//     double q8 = q2*q6;
//     double q10 = q2*q8;
//     double q12 = q2*q10;
//     double q14 = q2*q12;
//     double q16 = q2*q14;
//
//     double num1 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)* (1 - 0.25*q2)* (1 - 0.25*q2)* (1 - 0.25*q2);
//     double num2 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(q6-112*q4 +1616*q2-896);
//     double num3 = (1 - 0.25*q2)*(1 - 0.25*q2)*(q12-280*q10+14352*q8-224000*q6+1157888*q4-1505280*q2+405504);
//
//     double denom1 = (1 + 0.25*q2)*(1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2);
//     double denom2 = (1 + 0.25*q2)*(1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2);
//     double denom3 = (1 + 0.25*q2)*(1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2)* (1 + 0.25*q2);
//
//     double sumNum = 7*num1/denom1 -7*num2*coeff1/(32*denom2) + 21*num3*coeff2/(1024*denom3);
//
//     for(int d = 0; d<sys->dim; d++)
// //       sys->dwcdx[niac][d] = normFactor*sums*dx[d]/(h*h);
//       sys->dwcdx[niac][d] = -normFactor*sumNum*dx[d]/(h*h);
//
//
//   }
// }


__inline__ void BiharmonicDeconvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 1.0/(2.0*M_PI*h*h);


  const double q = r/h;
  const double eta = 0.01;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      const double sums = -2*log(sqrt(q2+eta)) + 1+log(4.0)-4.0/(q2+eta);

      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = alpha_d*sums*dx[d]/(h*h);


    }
}

__inline__ void ConsistentDeconvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 15.0/(72685184*M_PI*h*h);

  const double nord = 1/(2*M_PI*h*h*(45*log(2)-123/4));
  const double q = r/h;
  const double eta = 0.01;
  if(q >= 0.0 && q <= 3.0)
    {
      double q2 = q*q;
      double q4 = q2*q2;
      double q6 = q2*q4;
      double qq2 = (1 - q2/9)*(1 - q2/9);
      double qq4 = (1 - q2/9)*(1 - q2/9)*(1 - q2/9)*(1 - q2/9);
      double pp6 = (1 + q2/9)*(1 + q2/9)*(1 + q2/9)*(1 + q2/9)*(1 + q2/9)*(1 + q2/9);
      double pp8 = (1 + q2/9)*(1 + q2/9)*(1 + q2/9)*(1 + q2/9)*(1 + q2/9)*(1 + q2/9)*(1 + q2/9)*(1 + q2/9);
      const double sums = -20*qq4*nord/(9*pp6) + 160*qq2*(q6-180*q4+4293*q2-7290)/(59049*pp8);

      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = sums*dx[d]/(h*h);


    }
}

__inline__ void OrderFiveApproximateGaussianConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 3.0/((3*log(1048576)-41)*2*M_PI*h*h);
  const double q = r/h;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double qq5 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double pp5 = (1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2);
      sys->w[niac] = alpha_d*qq5/pp5;
    }
}

__inline__ void OrderFiveApproximateGausianDeconvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 3.0/((3*log(1048576)-41)*2*M_PI*h*h);

  const double q = r/h;
  const double eta = 0.01;

  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q4 = q2*q2;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double qq2 = (1 + 0.25*q2)*(1 - 0.25*q2);
      double qq4 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double pp6 = (1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2);
      double pp8 = (1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2);
      const double sums = -5*alpha_d*qq4/pp6 + qq2*(5/(64*M_PI))*(q6-80*q4+848*q2-640)/pp8;

      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = sums*dx[d]/(h*h);
    }
}

__inline__ void OrderSevenApproximateGaussianConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 15.0/((15*log(268435457)-289)*2*M_PI*h*h);
  const double q = r/h;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double qq7 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double pp7 = (1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2);
      sys->w[niac] = alpha_d*qq7/pp7;
    }
}

__inline__ void OrderSevenApproximateGaussianDeconvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 15.0/((15*log(268435457)-289)*2*M_PI*h*h);
  const double alpha_1 = 15.0/((15*log(268435457)-289)*2*M_PI*h*h);

  const double mCoef1 = -4*(-7 + 1/(-(289/15) + 28*log(2)));
  const double mCoef2 = 0.23427071691840673;

  const double q = r/h;
  const double eta = 0.01;

  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q4 = q2*q2;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double q10 = q2*q8;
      double q12 = q2*q10;

      double qq2 = (1 + 0.25*q2)*(1 - 0.25*q2);
      double qq4 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double qq6 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);

      double pp8 = (1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2);
      double pp10 = (1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2);
      double pp12 = (1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2);
      const double sums = -7*alpha_d*qq6/pp8 - alpha_d*mCoef1*qq4*(7/32)*(q6-112*q4+1616*q2-896)/pp10 - mCoef2*(7/(24*1024))*alpha_d*qq2*(q12-280*q10+14352*q8-133682*q6+314920*q4-1023584*q2+405504)/pp12;

      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = sums*dx[d]/(h*h);
    }
}

__inline__ void OrderFiveConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 3.0/((3*log(1048576)-41)*2*M_PI*h*h);
  const double q = r/h;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double qq5 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double pp5 = (1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2);
      sys->w[niac] = alpha_d*qq5/pp5;
    }
}

__inline__ void OrderFiveDeconvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 3.0/((3*log(1048576)-41)*2*M_PI*h*h);

  const double q = r/h;
  const double eta = 0.01;
  const double L1 = -16*(52 - 75*log(2))/(M_PI*(-41 + 60*log(2)));
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q4 = q2*q2;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double qq2 = (1 + 0.25*q2)*(1 - 0.25*q2);
      double qq4 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double pp6 = (1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2);
      double pp8 = (1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2)*(1 + 0.25*q2);
      const double sums = -5*alpha_d*qq4/pp6 - qq2*L1*(5/16)*(q6-80*q4+848*q2-640)/pp8;

      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = sums*dx[d]/(h*h);
    }
}

__inline__ void PolyOrderFiveConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 1386.0/(295.0*2*M_PI*h*h);
  const double q = r/h;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q4 = q2*q2;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double qq6 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double pp6 = (128.0 - 128.0*q2 + 88.0*q4- 48.0*q6 + 23.0*q8)/(128.0);
      sys->w[niac] = alpha_d*qq6*pp6;
    }
}

__inline__ void PolyOrderFiveDeconvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 1386.0/(128.0*295.0*2*M_PI*h*h);

  const double q = r/h;
  const double eta = 0.01;
  const double L1 = -16*(52 - 75*log(2))/(M_PI*(-41 + 60*log(2)));
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q4 = q2*q2;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double q10 = q2*q8;
      double qq2 = (1 + 0.25*q2)*(1 - 0.25*q2);
      double qq3 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double qq5 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double pp6 = (128.0 - 160.0*q2 + 128.0*q4 - 80.0*q6 + 23.0*q8);
      double pp8 = (-5120.0 + 15744.0*q2 - 21792.0*q4 + 17504.0*q6 - 7020.0*q8 + 1035*q10);
      const double sums = -5*alpha_d*qq5*pp6 +(5*144/(2*2*295))*qq3*pp8;

      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = sums*dx[d]/(h*h);
    }
}

__inline__ void PolyOrderSevenConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 2145.0/(1347584.0*M_PI*h*h);
  const double q = r/h;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q4 = q2*q2;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double q10 = q2*q8;
      double q12 = q2*q10;

      double qq = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double pp = (2048 - 3072*q2 + 2816*q4 - 1984*q6 + 1184*q8 - 628*q10 + 305*q12);
      sys->w[niac] = alpha_d*qq*pp;
    }
}

__inline__ void PolyOrderSevenDeconvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_0 = 2145.0/(192512.0*M_PI*h*h);
  const double alpha_1 = 2145.0/(385024.0*M_PI*h*h);
  const double alpha_2 = 6435.0/(385024.0*M_PI*h*h);
  const double deconvCoeff0 = 1.0;
  const double deconvCoeff1 = -131/329;
  const double deconvCoeff2 = 150882/1840097;

  const double q = r/h;
  const double eta = 0.01;
  const double L1 = -16*(52 - 75*log(2))/(M_PI*(-41 + 60*log(2)));

  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q4 = q2*q2;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double q10 = q2*q8;
      double q12 = q2*q10;
      double q14 = q2*q12;
      double q16 = q2*q14;

      double qq2 = (1 + 0.25*q2)*(1 - 0.25*q2);
      double qq3 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double pp3 = (1622016 - 10481664*q2 + 27981824*q4 - 46373376*q6 + 52908480*q8 - 39282928*q10 + 17215824*q12 - 3946488*q14 + 360815*q16);
      double qq5 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double pp5 = (-114688 + 464896*q2 - 809472*q4 + 938880*q6 - 846240*q8 + 547200*q10 - 197288*q12 + 27755*q14);
      double qq7 = (1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2)*(1 - 0.25*q2);
      double pp7 = (2048 - 3584*q2 + 3712*q4 - 2912*q6 + 1912*q8 - 1106*q10 + 305*q12);
      double pp8 = (-5120.0 + 15744.0*q2 - 21792.0*q4 + 17504.0*q6 - 7020.0*q8 + 1035*q10);
      const double sums = -deconvCoeff0*alpha_0*qq7*pp7 + deconvCoeff1*alpha_1*qq5*pp5/2.0 - alpha_2*deconvCoeff2*qq3*pp3/24.0;

      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = sums*dx[d]/(h*h);
    }
}


__inline__ void WendlandCSixConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 39.0/(14.0*M_PI*h*h);
  const double q = r/h;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q4 = q2*q2;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double q10 = q2*q8;
      double q12 = q2*q10;

      double qq = (1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q);
      double pp = (1 +4*q + 6.25*q2 + 4*q*q2);
      sys->w[niac] = alpha_d*qq*pp;
    }
}


__inline__ void WendlandCSixDeconvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_0 = 2145.0/(192512.0*M_PI*h*h);
  const double alpha_1 = 2145.0/(385024.0*M_PI*h*h);
  const double alpha_2 = 6435.0/(385024.0*M_PI*h*h);
  const double deconvCoeff0 = 1.0;
  const double deconvCoeff1 = -131/329;
  const double deconvCoeff2 = 150882/1840097;

  const double alpha_d = 429.0/(42649600.0*M_PI*h*h);

  const double q = r/h;
  const double eta = 0.01;
  const double L1 = -16*(52 - 75*log(2))/(M_PI*(-41 + 60*log(2)));

  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q3 = q*q2;
      double q4 = q2*q2;
      double q5 = q2*q3;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double q10 = q2*q8;
      double q12 = q2*q10;
      double q14 = q2*q12;
      double q16 = q2*q14;


      double qq = (-2 + q)*(-2 + q)*(-2 + q);
      double pp = (1699232 - 1524102 *q - 1630497 *q2 + 1241680 *q3 + 417180 *q4 - 339150 *q5 + 47600 *q6);
      const double sums = alpha_d*pp*qq;

      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = sums*dx[d]/(h*h);
    }
}

__inline__ void SimpleConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 39.0/(112.0*M_PI*h*h);
  const double q = r/h;
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q3 = q*q2;
      double q4 = q2*q2;
      double q5 = q*q4;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double q10 = q2*q8;
      double q12 = q2*q10;

      double pp9 = 8 + 36*q + 68*q2 + 66*q3 + 55*q4;
      double qq9 = (1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q);
      sys->w[niac] = alpha_d*pp9*qq9;
    }
}

__inline__ void SimpleDeconvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_0 = 429.0/(224.0*M_PI*h*h);
  const double alpha_1 = 2145.0/(896.0*M_PI*h*h);
  const double alpha_2 = 19305.0/(3584.0*M_PI*h*h);
  const double deconvCoeff1 = -(48717/269059);
  const double deconvCoeff2 = 94317/269059;
  const double deconvCoeff3 = 150882/1840097;

  const double alpha_d = 429.0/(42649600.0*M_PI*h*h);

  const double q = r/h;
  const double eta2 = 0.01*h*h;
  const double L1 = -16*(52 - 75*log(2))/(M_PI*(-41 + 60*log(2)));
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q3 = q*q2;
      double q4 = q2*q2;
      double q5 = q2*q3;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double q10 = q2*q8;
      double q12 = q2*q10;
      double q14 = q2*q12;
      double q16 = q2*q14;

      double pp4 = (1680 - 12000*q + 40120*q2 - 51448*q3 + 20449*q4)/pow((q2 + eta2),0.5);
      double qq4 = (1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q);

      double pp6 = (-256 + 492*q - 1884*q2 + 1859*q3);
      double qq6 = (1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q);

      double pp8 = (8 + 32*q + 32*q2 + 65*q3);
      double qq8 = (1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q)*(1 - 0.5*q);
      const double sums = -alpha_0*qq8*pp8 - alpha_1*(deconvCoeff1/2)*pp6*qq6-alpha_2*(deconvCoeff2/24)*pp4*qq4;

      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = sums*dx[d]/(h*h);
    }
}

__inline__ void SingularConvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_d = 1.0/(48.0*M_PI*h*h);
  const double q = r/h;

  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q3 = q*q2;
      double q4 = q2*q2;
      double q5 = q*q4;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double q10 = q2*q8;
      double q12 = q2*q10;

      double pp = (-704 - 432*q2 + 108*q4 + 11*q6 + 432*q2*log(4) + 108*q4*log(4) + 64*log(64) + q6*log(64) - 384*log(q) - 864*q2*log(q) - 216*q4*log(q) - 6*q6*log(q));
      sys->w[niac] = alpha_d*pp;
    }
}

__inline__ void SingularDeconvolutionKernel(System *sys, const double *__restrict__ dx, const double r, const double h, const int niac)
{
  const double alpha_0 = 1.0/(4.0*M_PI*h*h);
  const double alpha_1 = 9.0/(2.0*M_PI*h*h);
  const double alpha_2 = 72.0/(1.0*M_PI*h*h);
  const double deconvCoeff1 = -(44/1875);
  const double deconvCoeff2 = 116/1875;
  const double deconvCoeff3 = 150882/1840097;

  const double alpha_d = 429.0/(42649600.0*M_PI*h*h);

  const double q = r/h;
  const double eta2 = 0.01*h*h;
  const double L1 = -16*(52 - 75*log(2))/(M_PI*(-41 + 60*log(2)));
  if(q >= 0.0 && q <= 2.0)
    {
      double q2 = q*q;
      double q3 = q*q2;
      double q4 = q2*q2;
      double q5 = q2*q3;
      double q6 = q2*q4;
      double q8 = q2*q6;
      double q10 = q2*q8;
      double q12 = q2*q10;
      double q14 = q2*q12;
      double q16 = q2*q14;

      double pp0 = (-32 - 144*q2 + 18*q4 + 5*q6 + 144*q2*log(2) + q6*log(8) + 18*q4*log(16) - 144*q2*log(q) - 72*q4*log(q) - 3*q6*log(q))/(q2+eta2);
      double pp1 = (-16 - 16*q2 + 5*q4 + 16*q2*log(4) + q4*log(16) - 32*q2*log(q) - 4*q4*log(q))/(q2+eta2);
      double pp2 = (-4 + q2 + q2*log(4) - 2*q2*log(q))/(q2+eta2);
      const double sums = alpha_0*pp0 +alpha_1*(deconvCoeff1/2)*pp1 + alpha_2*(deconvCoeff2/24)*pp2;

      for(int d = 0; d <2; d++)
        sys->dwcdx[niac][d] = sums*dx[d]/(h*h);
    }
}


/* void SearchNeighbors1(System *sys) */
/* { */
/*   int niac = 0; */
/*   const int particles = sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+sys->NumberOfBoundaryParticles+sys->NumberOfVirtualParticles+sys->NumberOfFrozenParticles; */

/*   sys->NumberOfInteractingParticles = SearchNeighbors((const double **)sys->Position, */
/*                                                       sys->i_pair, */
/*                                                       sys->j_pair, */
/*                                                       sys->rij2, */
/*                                                       sys->NumberOfFluidParticles, */
/*                                                       particles, */
/*                                                       sys->NumberOfFluidParticles + */
/*                                                       sys->NumberOfSolidParticles + */
/*                                                       sys->NumberOfBoundaryParticles - 1, */
/*                                                       40, */
/*                                                       sys->hsml[0] * sys->hsml[0] * 4); */



  /* #pragma omp parallel for */
  /*   for(int i = 0; i < sys->NumberOfFluidParticles; i++) // loop over all particles */
  /*     { */
  /*       int test = (sys->Position[i][0]>= sys->dx); */
  /*       int nei = 0; */
  /*       for(int j = i+1; j < particles; j++) */
  /*         { */
  /*           int check = test&&(j>(sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+sys->NumberOfBoundaryParticles-1)); */
  /*           if(!check) { */
  /*             double dx[2]; */
  /*             for(int d = 0; d < 2; d++) */
  /*               dx[d] = sys->Position[i][d] - sys->Position[j][d]; */
  /*             const double r2 = dx[0]*dx[0] + dx[1]*dx[1]; */
  /*             const double h = 0.5*(sys->hsml[i] + sys->hsml[j]); */
  /*             if(r2 <= 4.*h*h) { */
  /*               sys->Neighbors[i][nei] = j; */
  /*               sys->DistanceNeighbors[i][nei] = r2; */
  /*               nei++; */
  /*             } */
  /*           } */
  /*         } */
  /*       sys->NeighborsSize[i] = nei; */
  /*     } */

  /* */
  /* for(int i = 0; i < sys->NumberOfFluidParticles; i++) // loop over all particles */
  /*   for(int nei=0;nei<sys->NeighborsSize[i];nei++) { */
  /*     const int j = sys->Neighbors[i][nei]; */
  /*     sys->i_pair[niac] = i; */
  /*     sys->j_pair[niac] = j; */
  /*     sys->rij2[niac] = sys->DistanceNeighbors[i][nei]; */
  /*     niac++; */
  /*   } */

  /* sys->NumberOfInteractingParticles = niac; */
/* } */


void CalculateKernel(System *sys)
{
  memset(sys->w, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->dwdx[0], 0, sizeof(double)*sys->MaxNumberOfParticles*2);
  memset(sys->wc, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->dwcdx[0], 0, sizeof(double)*sys->MaxNumberOfParticles*2);
  // SearchNeighbors1(sys);
  const int particles = sys->NumberOfFluidParticles +
    sys->NumberOfSolidParticles +
    sys->NumberOfBoundaryParticles +
    sys->NumberOfVirtualParticles +
    sys->NumberOfFrozenParticles;

  sys->NumberOfInteractingParticles = SearchNeighbors((const double **)sys->Position,
                                                      sys->i_pair,
                                                      sys->j_pair,
                                                      sys->rij2,
                                                      sys->NumberOfFluidParticles,
                                                      particles,
                                                      sys->NumberOfFluidParticles +
                                                      sys->NumberOfSolidParticles +
                                                      sys->NumberOfBoundaryParticles - 1,
                                                      40,
                                                      sys->hsml[0] * sys->hsml[0] * 4);

#pragma omp parallel for
#pragma simd
  for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
    {
      int i = sys->i_pair[k];
      int j = sys->j_pair[k];
      double dx[2];
      double h = 0.5*(sys->hsml[i]+sys->hsml[j]); /* average smoothing length*/
      double hb = 1.2*sys->dx; // radius of spiky boundary kernel.
      for(int d = 0; d < 2; d++)
        dx[d] = sys->Position[i][d] - sys->Position[j][d];
      const double r = sqrt(sys->rij2[k]);
      PolyOrderSevenConvolutionKernel(sys, dx, r, h, k);
      PolyOrderSevenDeconvolutionKernel(sys, dx, r, h, k);

      //       OrderSevenApproximateGaussianConvolutionKernel(sys, dx, r, h, k);
      //       OrderSevenApproximateGaussianDeconvolutionKernel(sys, dx, r, h, k);
    }
}


void PressureGradientForce(System *sys)
{

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    //     #pragma simd
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        const int i = sys->i_pair[k];
        const int j = sys->j_pair[k];

        const double rhoi = sys->localDensity[j]/(sys->density[i]*sys->density[i]*sys->density[j]);
        const double rhoj = sys->localDensity[i]/(sys->density[j]*sys->density[j]*sys->density[i]);
        const double pij = sys->DynamicPressure[i]*rhoi + sys->DynamicPressure[j]*rhoj;

        // sph form of pressure gradient
        //     for(int d = 0; 2-d; d++)
        //     {
        //       sys->dvdt[i][d] -= sys->mass[j]*pij*sys->dwcdx[k][d];
        //     }
        //     for(int d = 0; (2-d)&&(j<sys->NumberOfFluidParticles); d++)
        //     {
        //       sys->dvdt[j][d] += sys->mass[i]*pij*sys->dwcdx[k][d];
        //     }
        for(int d = 0; d<2; d++)
          {
            sys->tmp[threadid][2*i+d] -= sys->mass[j]*pij*sys->dwcdx[k][d];
            //             if(j<sys->NumberOfFluidParticles)
            sys->tmp[threadid][2*j+d] += sys->mass[i]*pij*sys->dwcdx[k][d];
          }
      }
  }

  for(int thd=0;thd<MaxThreads;thd++) {
    #pragma omp parallel for
    for(int d=0;d<sys->NumberOfFluidParticles;d++) {
      sys->dvdt[d][0] += sys->tmp[thd][2*d];
      sys->dvdt[d][1] += sys->tmp[thd][2*d+1];
    }
  }
  //   for(int thd=0;thd<MaxThreads;thd++)
  //     vdAdd(sys->NumberOfFluidParticles*2,
  //           sys->tmp[thd],
  //           sys->dvdt[0],
  //           sys->dvdt[0]);
}



void BodyForce(System *sys)
{
  if(sys->gravityFlag==1) // simulations with gravity
    {
#pragma omp parallel for
      for(int i=0;i<sys->NumberOfFluidParticles;i++)
        {
          sys->dvdt[i][0] -= 0.0;
          sys->dvdt[i][1] -= sys->Gravity;
        }
    }
}



/*Leonnard- Jones boundary force*/

void LeonnardJonesBoundaryForces(System *sys)
{
  const double r0 = sys->dx; /*initial particle spacing*/
  double dx[2];

  for(int j = sys->NumberOfFluidParticles + sys->NumberOfSolidParticles; j < sys->NumberOfFluidParticles + sys->NumberOfSolidParticles + sys->NumberOfBoundaryParticles; j++)
    {
      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
        {
          for(int d = 0; d < 2; d++)
            dx[d] = sys->Position[i][d] - sys->Position[j][d];

          double r2 = dx[0]*dx[0] + dx[1]*dx[1];

          if(r2 < r0*r0)
            {
              double r = r0/sqrt(r2);
              r *= r; // r^2
              r *= r; // r^4
              double pb = sys->LennardCoefficient*r*(r*r-1.)/r2;
              for(int d = 0; d < 2; d++)
                {
                  sys->dvdt[i][d] += pb*dx[d];
                }
            }
        }
    }
}


/**************************************************************************************************/
/* algorithm for predicting the density, pressure and velocity */

void PredictingDensityEquation(System *sys)
{
  /* solve the density equation */

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        double dx[2], dv[2];
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        for(int d = 0; d < sys->dim; d++)
          dv[d] = sys->Velocity[i][d] - sys->localVelocity[j][d];
        for(int d = 0; d < sys->dim; d++)
          {
            double vdw = dv[d]*sys->dwcdx[k][d];

            sys->tmp[threadid][i] += sys->mass[j]*vdw;
            //             if(j<sys->NumberOfFluidParticles)
            sys->tmp[threadid][j] += sys->mass[j]*vdw;
          }
      }
  }
  for(int thd=0;thd<MaxThreads;thd++) {
    vdAdd(sys->NumberOfFluidParticles, sys->tmp[thd], sys->ddensitydt_approx, sys->ddensitydt_approx);
  }
}

void PredictingDynamicPressureEquation(System *sys)
{

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    /* solve the pressure equation */
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        double dx[2], dv[2];

        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        for(int d = 0; d < sys->dim; d++)
          {
            dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
            dx[d] = sys->Position[i][d] - sys->Position[j][d];
          }
        double vdw = 0.0;
        double rdwij = 0.0;
        for(int d = 0; d < sys->dim; d++)
          {
            vdw += dv[d]*sys->dwcdx[k][d];
            rdwij += dx[d]*sys->dwcdx[k][d];
          }

        //         double phij = sys->Gravity*sys->density0*(sys->Position[i][1]-sys->Position[j][1]);
        double pi = (sys->density0*sys->SpeedOfSound*sys->SpeedOfSound + sys->AdiabaticConstant*sys->DynamicPressure[i]);
        double pj = (sys->density0*sys->SpeedOfSound*sys->SpeedOfSound + sys->AdiabaticConstant*sys->DynamicPressure[j]);
        double compi = 1.0/pi;
        double compj = 1.0/pj;
        double speci = 1.0/sys->density[i];
        double specj = 1.0/sys->density[j];

        double volumei = sys->mass[j]/sys->density[i];
        double volumej = sys->mass[j]/sys->density[j];

        double vj = sys->mass[j]/sys->density[i];
        double vi = sys->mass[i]/sys->density[j];

        sys->tmp[threadid][i] += vj*pi*vdw;
        sys->tmp[threadid][j] += vi*pj*vdw;


        double h = 0.5*(sys->hsml[i] + sys->hsml[j]);
        double eta2 = 0.01*h*h;

        double nu = rdwij/(sys->rij2[k] + eta2);
        double pij = 1.0 + pi/pj;
        double pji = 1.0 + pj/pi;

        double epj = 4.0*pi/(pi + pj);
        double epi = 4.0*pj/(pi + pj);

        // laplacian of p
        sys->tmp[threadid][i] += sys->AdiabaticConstant*sys->thermalDiffusivity[i]*volumej*(sys->DynamicPressure[i] - sys->DynamicPressure[j])*nu*epj;
        // 	if(j<sys->NumberOfFluidParticles)
        sys->tmp[threadid][j] += sys->AdiabaticConstant*sys->thermalDiffusivity[j]*volumei*(sys->DynamicPressure[j] - sys->DynamicPressure[i])*nu*epi;
        // density laplacian
        double rho = (sys->density[i] + sys->density[j])/(sys->density[i]*sys->density[j]);
        double rhoi = (sys->density[i] - sys->density[j])*rho;
        double rhoj = (sys->density[j] - sys->density[i])*rho;
        sys->tmp[threadid][i] -= sys->thermalDiffusivity[i]*volumej*pi*rhoi*nu;
        // 	if(j<sys->NumberOfFluidParticles)
        sys->tmp[threadid][j] -= sys->thermalDiffusivity[j]*pj*volumei*rhoj*nu;
      }
  }
  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles, sys->tmp[thd], sys->dpdt_approx, sys->dpdt_approx);
}

void PredictingPressureGradientForce(System *sys)
{

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    //     #pragma simd
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        const int i = sys->i_pair[k];
        const int j = sys->j_pair[k];

        const double rhoi = 1.0/(sys->density[i]*sys->density[i]);
        const double rhoj = 1.0/(sys->density[j]*sys->density[j]);
        const double pij = sys->DynamicPressure[i]*rhoi + sys->DynamicPressure[j]*rhoj;
        for(int d = 0; d<2; d++)
          {
            sys->tmp[threadid][2*i+d] -= sys->mass[j]*pij*sys->dwcdx[k][d];
            // 	if(j<sys->NumberOfFluidParticles)
            sys->tmp[threadid][2*j+d] += sys->mass[i]*pij*sys->dwcdx[k][d];
          }
      }
  }

  for(int thd=0;thd<MaxThreads;thd++) {
#pragma omp parallel for
    for(int d=0;d<sys->NumberOfFluidParticles;d++) {
      sys->dvdt_approx[d][0] += sys->tmp[thd][2*d];
      sys->dvdt_approx[d][1] += sys->tmp[thd][2*d+1];
    }
  }
  //   for(int thd=0;thd<MaxThreads;thd++)
  //     vdAdd(sys->NumberOfFluidParticles*2,
  //           sys->tmp[thd],
  //           sys->dvdt[0],
  //           sys->dvdt[0]);
}

void PredictingBodyForce(System *sys)
{
  if(sys->gravityFlag==1) // simulations with gravity
    {
#pragma omp parallel for
      for(int i=0;i<sys->NumberOfFluidParticles;i++)
        {
          sys->dvdt_approx[i][0] -= 0.0;
          sys->dvdt_approx[i][1] -= sys->Gravity;
        }
    }
}

void PredictingSpringBoundaryForces(System *sys)
{
  const double r0 = 3.0*sys->dx; /*initial particle spacing*/

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        double dx[2], dv[2], initX[2];
        const int i = sys->i_pair[k];
        const int j = sys->j_pair[k];
        double h = 0.5*(sys->hsml[i] + sys->hsml[j]);
        double eta2 = 0.01*h*h;

        if((sys->itype[i]==1)&&((sys->itype[j]==2)||(sys->itype[j]==3)||(sys->itype[j]==5))) // fixed and moving boundary particles
          {
            double rr = 0.0;
            double rvij = 0.0;
            for(int d = 0; d < sys->dim; d++)
              {
                dx[d] = sys->Position[i][d] - sys->Position[j][d];
                dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
                rr += dx[d]*dx[d];
                rvij += dx[d]*dv[d];
              }
            rr = sqrt(rr);
            double nu = rvij/rr ;//+ eta2); // replace rijk with rr*rr
            if((rr < r0))
              {
                double sumX[2];
                for(int d = 0; d < sys->dim; d++)
                  {
                    // r-part
                    sys->tmp[threadid][2*i+d] -= sys->springConstant0*(r0-rr)*dx[d]/(rr*sys->mass[i]);
                    // v-part
                    sys->tmp[threadid][2*i+d] -= sys->dampingConstant0*rvij*dx[d]/(rr*rr*sys->mass[i]);
                  }
              }

          }
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles*2, sys->tmp[thd], sys->dvdt_approx[0], sys->dvdt_approx[0]);
}

void PredictingViscousForcesSingleStep(System *sys)
{
  /* This is the deviatoric stress tensor. As such
     /* its trace is zero.*/

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        double dx[2], dv[2];
        for(int d = 0; d < sys->dim; d++)
          {
            dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
            dx[d] = sys->Position[j][d] - sys->Position[i][d];
          }

        double h = 0.5*(sys->hsml[i] + sys->hsml[j]);
        double eta2 = 0.01*h*h;

        double invR2 = 1.0/(sys->rij2[k] + eta2);

        double rhoi = 1.0/(sys->density[i]*sys->density[i]);
        double rhoj = 1.0/(sys->density[j]*sys->density[j]);

        double effViscosityij = 2.0*(sys->density[i]*sys->viscosity[i])*(sys->density[j]*sys->viscosity[j])/(sys->density[i]*sys->viscosity[i] + sys->density[j]*sys->viscosity[j]);
        double iSigma00 = effViscosityij*(dx[0]*dv[0] - dx[1]*dv[1])*invR2;
        double iSigma01 = effViscosityij*(dx[0]*dv[1] + dx[1]*dv[0])*invR2;
        double iSigma10 = effViscosityij*(dx[1]*dv[0] + dx[0]*dv[1])*invR2;
        double iSigma11 = effViscosityij*(dx[1]*dv[1] - dx[0]*dv[0])*invR2;

        double jSigma00 = effViscosityij*(dx[0]*dv[0] - dx[1]*dv[1])*invR2;
        double jSigma01 = effViscosityij*(dx[0]*dv[1] + dx[1]*dv[0])*invR2;
        double jSigma10 = effViscosityij*(dx[1]*dv[0] + dx[0]*dv[1])*invR2;
        double jSigma11 = effViscosityij*(dx[1]*dv[1] - dx[0]*dv[0])*invR2;

        double acclCx  = (iSigma00*rhoi + jSigma00*rhoj )* sys->dwcdx[k][0];
        acclCx += (iSigma01*rhoi + jSigma01*rhoj )* sys->dwcdx[k][1];
        double acclCy  = (iSigma10*rhoi + jSigma10*rhoj )* sys->dwcdx[k][0];
        acclCy += (iSigma11*rhoi + jSigma11*rhoj )* sys->dwcdx[k][1];
        /* x- component of acceleration */
        sys->tmp[threadid][2*i] += sys->mass[j]*acclCx;
        /* y- component of acceleration */
        sys->tmp[threadid][2*i+1] += sys->mass[j]*acclCy;

        //         if(j<sys->NumberOfFluidParticles) {
        sys->tmp[threadid][2*j] -= sys->mass[i]*acclCx;
        sys->tmp[threadid][2*j+1] -= sys->mass[i]*acclCy;
        //         }
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles*2,
          sys->tmp[thd],
          sys->dvdt_approx[0],
          sys->dvdt_approx[0]);
}

/*********************************End Prediction Algorithm for rates of change *******************************************************/


void DensityEquation(System *sys)
{
  /* solve the density equation */

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        double dx[2], dvij[2], dvji[2];
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        for(int d = 0; d < sys->dim; d++)
          {
            dvij[d] = sys->Velocity[i][d] - sys->localVelocity[j][d];
            dvji[d] = sys->Velocity[j][d] - sys->localVelocity[i][d];
          }
        for(int d = 0; d < sys->dim; d++)
          {

            double vdwij = dvij[d]*sys->dwcdx[k][d];
            double vdwji = dvji[d]*sys->dwcdx[k][d];

            double volumei = sys->mass[i]/sys->density[i];
            double volumej = sys->mass[j]/sys->density[j];

            double mj = volumej*sys->localDensity[j];
            double mi = volumei*sys->localDensity[i];


            sys->tmp[threadid][i] += mj*vdwij;
            //             if(j<sys->NumberOfFluidParticles)
            sys->tmp[threadid][j] -= mi*vdwji;
          }
      }
  }
  for(int thd=0;thd<MaxThreads;thd++) {
    vdAdd(sys->NumberOfFluidParticles, sys->tmp[thd], sys->ddensitydt, sys->ddensitydt);
  }
}

void DynamicPressureEquation(System *sys)
{

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    /* solve the pressure equation */
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        double dx[2], dvij[2], dvji[2];

        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        for(int d = 0; d < sys->dim; d++)
          {
            dvij[d] = sys->Velocity[i][d] - sys->localVelocity[j][d];
            dvji[d] = sys->Velocity[j][d] - sys->localVelocity[i][d];
            dx[d] = sys->Position[i][d] - sys->Position[j][d];
          }
        double vdwij = 0.0;
        double vdwji = 0.0;
        double rdwij = 0.0;
        for(int d = 0; d < sys->dim; d++)
          {
            vdwij += dvij[d]*sys->dwcdx[k][d];
            vdwji += dvji[d]*sys->dwcdx[k][d];
            rdwij += dx[d]*sys->dwcdx[k][d];
          }

        //         double phij = sys->Gravity*sys->density0*(sys->Position[i][1]-sys->Position[j][1]);
        double pi = (sys->density0*sys->SpeedOfSound*sys->SpeedOfSound + sys->AdiabaticConstant*sys->DynamicPressure[i]);
        double pj = (sys->density0*sys->SpeedOfSound*sys->SpeedOfSound + sys->AdiabaticConstant*sys->DynamicPressure[j]);
        double compi = 1.0/pi;
        double compj = 1.0/pj;
        double speci = 1.0/sys->density[i];
        double specj = 1.0/sys->density[j];

        double volumei = sys->mass[i]/sys->density[i];
        double volumej = sys->mass[j]/sys->density[j];

        double reducedVolumei = volumei*sys->localDensity[i]/sys->density[j];
        double reducedVolumej = volumej*sys->localDensity[j]/sys->density[i];

        sys->tmp[threadid][i] += reducedVolumej*pi*vdwij;
        sys->tmp[threadid][j] -= reducedVolumei*pj*vdwji;


        double h = 0.5*(sys->hsml[i] + sys->hsml[j]);
        double eta2 = 0.01*h*h;

        double nu = rdwij/(sys->rij2[k] + eta2);
        double pij = 1.0 + pi/pj;
        double pji = 1.0 + pj/pi;

        double epj = 4.0*pi/(pi + pj);
        double epi = 4.0*pj/(pi + pj);

        /***********************************************************************************************************************************************************/
        // symmetrized form
        //laplacian of pressure
        //       double rpij = (compj + sys->localCompressibility[i])*(sys->DynamicPressure[i]-sys->localPressure[j]);
        //       rpij += (compi + sys->localCompressibility[j])*(sys->localPressure[i]-sys->DynamicPressure[j]);
        //       double rpji = (compi + sys->localCompressibility[j])*(sys->DynamicPressure[j]-sys->localPressure[i]);
        //       rpji += (compj + sys->localCompressibility[i])*(sys->localPressure[j]-sys->DynamicPressure[i]);
        //
        //       sys->dpdt[i] += 0.5*sys->thermalDiffusivity[i]*mj*nu*rpij*pi;
        //       sys->dpdt[j] += 0.5*sys->thermalDiffusivity[j]*mi*nu*rpji*pj;
        //
        //       // laplacian of density
        //       double rhoij = (specj + sys->specificVolume[i])*(sys->density[i]-sys->localDensity[j]);
        //       rhoij += (speci + sys->specificVolume[j])*(sys->localDensity[i]-sys->density[j]);
        //       double rhoji = (speci + sys->specificVolume[j])*(sys->density[j]-sys->localDensity[i]);
        //       rhoji += (specj + sys->specificVolume[i])*(sys->localDensity[j]-sys->density[i]);
        //       sys->dpdt[i] -= 0.5*(sys->thermalDiffusivity[i]/sys->AdiabaticConstant)*mj*nu*rhoij*pi;
        //       sys->dpdt[j] -= 0.5*(sys->thermalDiffusivity[j]/sys->AdiabaticConstant)*mi*nu*rhoji*pj;

        // intuitive formulas
        // laplacian of p
        sys->tmp[threadid][i] += sys->AdiabaticConstant *
          sys->thermalDiffusivity[i] * volumej *
          (sys->DynamicPressure[i] - sys->localPressure[j]) * nu * epj;
        // 	if(j<sys->NumberOfFluidParticles)
        sys->tmp[threadid][j] += sys->AdiabaticConstant * sys->thermalDiffusivity[j] *
          volumei * (sys->DynamicPressure[j] - sys->localPressure[i]) * nu * epi;
        // density laplacian
        double rho = (sys->density[i] + sys->density[j])/(sys->density[i]*sys->density[j]);
        double rhoi = (sys->density[i] - sys->localDensity[j])*rho;
        double rhoj = (sys->density[j] - sys->localDensity[i])*rho;
        sys->tmp[threadid][i] -= sys->thermalDiffusivity[i]*volumej*pi*rhoi*nu;
        // 	if(j<sys->NumberOfFluidParticles)
        sys->tmp[threadid][j] -= sys->thermalDiffusivity[j]*pj*volumei*rhoj*nu;
      }
  }
  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles, sys->tmp[thd], sys->dpdt, sys->dpdt);
}

void EquationOfState(System *sys)
{
#pragma omp parallel for
  for(int i = 0; i < (sys->NumberOfFluidParticles + sys->NumberOfSolidParticles + sys->NumberOfBoundaryParticles); i++)
    {
      const double density = sys->density[i]/sys->density0;
      sys->Pressure[i] = sys->CompressionFactor*(pow(density, sys->AdiabaticConstant) - 1.0);
    }
}



void ClassicShearStress(System *sys)
{
  double dvij[2], dvji[2];
  for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
    {
      int i = sys->i_pair[k];
      int j = sys->j_pair[k];
      for(int d = 0; d < 2; d++)
        {
          dvij[d] = sys->Velocity[i][d] - sys->localVelocity[j][d];
          dvji[d] = sys->Velocity[j][d] - sys->localVelocity[i][d];
        }
      double volumei = sys->mass[i]/sys->density[i];
      double volumej = sys->mass[j]/sys->density[j];

      double reducedVolumei = volumei*sys->localDensity[i]/sys->density[j];
      double reducedVolumej = volumej*sys->localDensity[j]/sys->density[i];
      /* xx-component of gradV */
      sys->ViscousStressxx[i] -= sys->viscosity[i]*reducedVolumej*(dvij[0]*sys->dwcdx[k][0]+dvij[0]*sys->dwcdx[k][0]);

      sys->ViscousStressxx[j] += sys->viscosity[j]*reducedVolumei*(dvji[0]*sys->dwcdx[k][0]+dvji[0]*sys->dwcdx[k][0]);

      /* xy-component of gradV */
      sys->ViscousStressxy[i] -= sys->viscosity[i]*reducedVolumej*(dvij[0]*sys->dwcdx[k][1]+dvij[1]*sys->dwcdx[k][0]);

      sys->ViscousStressxy[j] += sys->viscosity[j]*reducedVolumei*(dvji[0]*sys->dwcdx[k][1]+dvji[1]*sys->dwcdx[k][0]);
      /* yx-component of gradV */
      sys->ViscousStressyx[i] -= sys->viscosity[i]*reducedVolumej*(dvij[1]*sys->dwcdx[k][0]+dvij[0]*sys->dwcdx[k][1]);

      sys->ViscousStressyx[j] += sys->viscosity[j]*reducedVolumei*(dvji[1]*sys->dwcdx[k][0]+dvji[0]*sys->dwcdx[k][1]);
      /* yy-component of gradV */
      sys->ViscousStressyy[i] -= sys->viscosity[i]*reducedVolumej*(dvij[1]*sys->dwcdx[k][1]+dvij[1]*sys->dwcdx[k][1]);

      sys->ViscousStressyy[j] += sys->viscosity[j]*reducedVolumei*(dvji[1]*sys->dwcdx[k][1]+dvji[1]*sys->dwcdx[k][1]);
    }
}


void ViscousForcesTwoStep(System *sys)
{

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        double dx[2];

        double rhoi = sys->localDensity[j]/(sys->density[i]*sys->density[i]*sys->density[j]);
        double rhoj = sys->localDensity[i]/(sys->density[j]*sys->density[j]*sys->density[i]);

        double acclCx  = (sys->ViscousStressxx[i]*rhoi + sys->ViscousStressxx[j]*rhoj )* sys->dwcdx[k][0];
        acclCx += (sys->ViscousStressxy[i]*rhoi + sys->ViscousStressxy[j]*rhoj )* sys->dwcdx[k][1];
        double acclCy  = (sys->ViscousStressyx[i]*rhoi + sys->ViscousStressyx[j]*rhoj )* sys->dwcdx[k][0];
        acclCy += (sys->ViscousStressyy[i]*rhoi + sys->ViscousStressyy[j]*rhoj )* sys->dwcdx[k][1];
        /* x- component of acceleration */
        sys->tmp[threadid][2*i] += sys->mass[j]*acclCx;
        /* y- component of acceleration */
        sys->tmp[threadid][2*i+1] += sys->mass[j]*acclCy;

        if(j<sys->NumberOfFluidParticles) {
          sys->tmp[threadid][2*j] -= sys->mass[i]*acclCx;
          sys->tmp[threadid][2*j+1] -= sys->mass[i]*acclCy;
        }
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles*2,
          sys->tmp[thd],
          sys->dvdt[0],
          sys->dvdt[0]);
}

void ViscousForcesSingleStep(System *sys)
{
  /* This is the deviatoric stress tensor. As such
     /* its trace is zero.*/

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        double dx[2], dv[2], dvij[2], dvji[2];
        for(int d = 0; d < sys->dim; d++)
          {
            dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
            dvij[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
            dvji[d] = sys->Velocity[j][d] - sys->Velocity[i][d];
            dx[d] = sys->Position[j][d] - sys->Position[i][d];
          }

        double h = 0.5*(sys->hsml[i] + sys->hsml[j]);
        double eta2 = 0.01*h*h;

        double invR2 = 1.0/(sys->rij2[k] + eta2);

        double rhoi = sys->localDensity[j]/(sys->density[i]*sys->density[i]*sys->density[j]);
        double rhoj = sys->localDensity[i]/(sys->density[j]*sys->density[j]*sys->density[i]);

        double effViscosityij = 2.0*(sys->density[i]*sys->viscosity[i])*(sys->density[j]*sys->viscosity[j])/(sys->density[i]*sys->viscosity[i] + sys->density[j]*sys->viscosity[j]);
        double iSigma00 = effViscosityij*(dx[0]*dvij[0] - dx[1]*dvij[1])*invR2;
        double iSigma01 = effViscosityij*(dx[0]*dvij[1] + dx[1]*dvij[0])*invR2;
        double iSigma10 = effViscosityij*(dx[1]*dvij[0] + dx[0]*dvij[1])*invR2;
        double iSigma11 = effViscosityij*(dx[1]*dvij[1] - dx[0]*dvij[0])*invR2;

        double jSigma00 = effViscosityij*(dx[0]*dvji[0] - dx[1]*dvji[1])*invR2;
        double jSigma01 = effViscosityij*(dx[0]*dvji[1] + dx[1]*dvji[0])*invR2;
        double jSigma10 = effViscosityij*(dx[1]*dvji[0] + dx[0]*dvji[1])*invR2;
        double jSigma11 = effViscosityij*(dx[1]*dvji[1] - dx[0]*dvji[0])*invR2;

        double acclCx  = (iSigma00*rhoi + jSigma00*rhoj )* sys->dwcdx[k][0];
        acclCx += (iSigma01*rhoi + jSigma01*rhoj )* sys->dwcdx[k][1];
        double acclCy  = (iSigma10*rhoi + jSigma10*rhoj )* sys->dwcdx[k][0];
        acclCy += (iSigma11*rhoi + jSigma11*rhoj )* sys->dwcdx[k][1];
        /* x- component of acceleration */
        sys->tmp[threadid][2*i] += sys->mass[j]*acclCx;
        /* y- component of acceleration */
        sys->tmp[threadid][2*i+1] += sys->mass[j]*acclCy;

        //         if(j<sys->NumberOfFluidParticles) {
        sys->tmp[threadid][2*j] -= sys->mass[i]*acclCx;
        sys->tmp[threadid][2*j+1] -= sys->mass[i]*acclCy;
        //         }
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles*2,
          sys->tmp[thd],
          sys->dvdt[0],
          sys->dvdt[0]);
}

void SpecificPowerDissipation(System *sys)
{
  /* This is the deviatoric stress tensor. As such
     /* its trace is zero.*/

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        double dx[2], dv[2], dvij[2], dvji[2];
        double rvij = 0.0;
        for(int d = 0; d < sys->dim; d++)
          {
            dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
            dvij[d] = sys->Velocity[i][d] - sys->localVelocity[j][d];
            dvji[d] = sys->Velocity[j][d] - sys->localVelocity[i][d];
            dx[d] = sys->Position[j][d] - sys->Position[i][d];
            rvij += dx[d]*dv[d];
          }

        double h = 0.5*(sys->hsml[i] + sys->hsml[j]);
        double eta2 = 0.01*h*h;

        double invR2 = 1.0/(sys->rij2[k] + eta2);

        double rhoi = sys->localDensity[j]/(sys->density[i]*sys->density[i]*sys->density[j]);
        double rhoj = sys->localDensity[i]/(sys->density[j]*sys->density[j]*sys->density[i]);

        double effViscosityij = 2.0*(sys->density[i]*sys->viscosity[i])*(sys->density[j]*sys->viscosity[j])/(sys->density[i]*sys->viscosity[i] + sys->density[j]*sys->viscosity[j]);
        double iSigma00 = effViscosityij*(dx[0]*dvij[0] - dx[1]*dvij[1])*invR2;
        double iSigma01 = effViscosityij*(dx[0]*dvij[1] + dx[1]*dvij[0])*invR2;
        double iSigma10 = effViscosityij*(dx[1]*dvij[0] + dx[0]*dvij[1])*invR2;
        double iSigma11 = effViscosityij*(dx[1]*dvij[1] - dx[0]*dvij[0])*invR2;

        double jSigma00 = effViscosityij*(dx[0]*dvji[0] - dx[1]*dvji[1])*invR2;
        double jSigma01 = effViscosityij*(dx[0]*dvji[1] + dx[1]*dvji[0])*invR2;
        double jSigma10 = effViscosityij*(dx[1]*dvji[0] + dx[0]*dvji[1])*invR2;
        double jSigma11 = effViscosityij*(dx[1]*dvji[1] - dx[0]*dvji[0])*invR2;

        double iudw00 = dv[0]*sys->dwcdx[k][0];
        double iudw01 = dv[0]*sys->dwcdx[k][1];
        double iudw10 = dv[1]*sys->dwcdx[k][0];
        double iudw11 = dv[1]*sys->dwcdx[k][1];

        double judw00 = dv[0]*sys->dwcdx[k][0];
        double judw01 = dv[0]*sys->dwcdx[k][1];
        double judw10 = dv[1]*sys->dwcdx[k][0];
        double judw11 = dv[1]*sys->dwcdx[k][1];

        double jE  = (iSigma00*rhoi + jSigma00*rhoj )* iudw00 + (iSigma01*rhoi + jSigma01*rhoj )* iudw01;
        jE = (iSigma10*rhoi + jSigma10*rhoj )* iudw10 + (iSigma11*rhoi + jSigma11*rhoj )* iudw11;
        double iE  = (iSigma00*rhoi + jSigma00*rhoj )* judw00 + (iSigma01*rhoi + jSigma01*rhoj )* judw01;
        iE = (iSigma10*rhoi + jSigma10*rhoj )* judw10 + (iSigma11*rhoi + jSigma11*rhoj )* judw11;

        sys->tmp[threadid][2*i] -= sys->mass[j]*jE;
        if(j<sys->NumberOfFluidParticles)
          sys->tmp[threadid][2*j] += sys->mass[i]*iE;
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles,
          sys->tmp[thd],
          sys->specificDissipationPower,
          sys->specificDissipationPower);
}

void MonaghanBoundaryForce(System *sys)
{
  double dx[2];
  for(int j = sys->NumberOfFluidParticles+sys->NumberOfSolidParticles; j < sys->NumberOfFluidParticles+sys->NumberOfSolidParticles+sys->NumberOfBoundaryParticles; j++)
    {
      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
        {
          for(int d = 0; d < 2; d++)
            dx[d] = sys->Position[i][d] - sys->Position[j][d];

          double r2 = dx[0]*dx[0] + dx[1]*dx[1];
          double r = sqrt(r2);
          double h = sys->hsml[0];
          double q = r/h;
          double phi_ij = 0.0;
          //define boundary kernel
          if(q>=0.0 && q<2.0)
            {
              double q2 = q*q;
              double q3 = (2.0-q)*(2.0-q)*(2.0-q);
              phi_ij = 0.125*(1.0 + 1.5*q)*q3;
            }
          else
            {
              phi_ij = 0.0;
            }
          // calculate force exerted on fluid particle by boundary
          for(int d = 0; d < 2; d++)
            {
              double gama = (sys->machNumber*sys->SpeedOfSound)*(sys->machNumber*sys->SpeedOfSound)/(sys->density0*sys->dx*sys->dy);
              sys->dvdt[i][d] += gama*phi_ij*sys->mass[j]*dx[d]/r2;
            }
        }
    }
}

void ImprovedMonaghan(System *sys)
{
  const double r0 = sys->dx; /*initial particle spacing*/

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        double dx[2];
        const int i = sys->i_pair[k];
        const int j = sys->j_pair[k];
        if((sys->itype[i]==1)&&((sys->itype[j]==2)||(sys->itype[j]==3)))
          {
            double rr = 0.0;
            for(int d = 0; d < sys->dim; d++)
              {
                dx[d] = sys->Position[i][d] - sys->Position[j][d];
                rr += dx[d]*dx[d];
              }
            rr = sqrt(rr);

            if(rr < r0)
              {
                double h = 0.5*(sys->hsml[i]+sys->hsml[j]);
                double normFactor = 7.0/(64.0*3.142*h*h);
                double eta2 = 0.01*h*h;
                double q = rr/h;
                double phi_ij = 0.0;
                //define boundary kernel
                if(q>=0.0 && q<2.0)
                  {
                    double q2 = q*q;
                    double q4 = (2.0-q)*(2.0-q)*(2.0-q)*(2.0-q);
                    phi_ij = normFactor*(1.0 + 2.0*q)*q4;
                  }
                double boundarySpacing = 5.0; // dx_boundary = dx/4.0 where dx is fluid spacing
                double rd = (rr-sys->dx/boundarySpacing);
                double maxSpeed = (sys->machNumber*sys->SpeedOfSound)*(sys->machNumber*sys->SpeedOfSound)/boundarySpacing;
                for(int d = 0; d < sys->dim; d++)
                  {
                    sys->tmp[threadid][2*i+d] += maxSpeed*phi_ij*dx[d]/(rd*rr);
                  }
              }
          }
      }
  }
  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles*2, sys->tmp[thd], sys->dvdt[0], sys->dvdt[0]);
}

void CorePotential(System *sys)
{
  /* the core potential is used to prevent particle inter-penetration and
     clumping. the coefficient cr represents the strength of the repulsion whereas
     the factor "sigma" gives the core size. Here cr = gd^2 and sigma = dx. Modifications
     have been made to the original formula by Hoover. */

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        double dx[2];

        const double r0 = sys->dx; /*initial particle spacing = core size */
        for(int d = 0; d < 2; d++)
          dx[d] = sys->Position[i][d] - sys->Position[j][d];
        double r2 = dx[0]*dx[0] + dx[1]*dx[1];
        if(r2 < r0*r0)
          {
            double r = sqrt(r2)/r0;
            r *= r; // r^2
            double coreAccel = 2.0*(sys->coreCoefficient/(r0*r0))*(1.0 - r)*(1.0 - r)*(1.0 - r);
            for(int d = 0; d < 2; d++)
              sys->tmp[threadid][2*i+d] += coreAccel*dx[d];
            for(int d = 0; (d<2)&&(j<sys->NumberOfFluidParticles); d++)
              sys->tmp[threadid][2*j+d] -= coreAccel*dx[d];
          }
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles*2, sys->tmp[thd], sys->dvdt[0], sys->dvdt[0]);
}


void LocalVelocityFluctuation(System *sys)
{
  for(int i = 0; i< (sys->NumberOfFluidParticles+sys->NumberOfBoundaryParticles); i++)
    {
      sys->velocityFluctuation[i][0] = 0.0;
      sys->velocityFluctuation[i][1] = 0.0;
    }
  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        for(int d = 0; d < sys->dim; d++)
          {
            double avgpart = (sys->Velocity[i][d] - sys->Velocity[j][d])*sys->w[k];
            double rhoj = 0.5*(sys->density[i] + sys->density[j])/sys->density[i];
            double rhoi = 0.5*(sys->density[i] + sys->density[j])/sys->density[j];
            sys->tmp[threadid][2*i+d] -= sys->mass[j]*avgpart/sys->localDensity[i];
            // 	    if(j<sys->NumberOfFluidParticles)
            sys->tmp[threadid][2*i+d] += sys->mass[i]*avgpart/sys->localDensity[j];
          }
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles*2, sys->tmp[thd], sys->velocityFluctuation[0], sys->velocityFluctuation[0]);

}


void LocalDensityFluctuation(System *sys)
{
  const int niac = sys->NumberOfInteractingParticles;
  for(int i = 0; i< (sys->NumberOfFluidParticles+sys->NumberOfBoundaryParticles); i++)
    {
      sys->densityFluctuation[i] = 0.0;
    }
  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];

        double volumei = sys->mass[i]/sys->density[i];
        double volumej = sys->mass[j]/sys->density[j];

        double avgpart = (sys->density[i] - sys->density[j])*sys->w[k];
        sys->tmp[threadid][i] -= volumej*avgpart;
        //         if(j<sys->NumberOfFluidParticles)
        sys->tmp[threadid][j] += volumei*avgpart;
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles, sys->tmp[thd], sys->densityFluctuation, sys->densityFluctuation);
}

void LocalPressureFluctuation(System *sys)
{
  const int niac = sys->NumberOfInteractingParticles;
  for(int i = 0; i< (sys->NumberOfFluidParticles+sys->NumberOfBoundaryParticles); i++)
    {
      sys->pressureFluctuation[i] = 0.0;
    }
  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < niac; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];

        double volumei = sys->mass[i]/sys->density[i];
        double volumej = sys->mass[j]/sys->density[j];
        double avgpart = (sys->DynamicPressure[i] - sys->DynamicPressure[j])*sys->w[k];
        sys->tmp[threadid][i] -= volumej*avgpart;
        //         if(j<sys->NumberOfFluidParticles)
        sys->tmp[threadid][j] += volumei*avgpart;
      }

  }
  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles, sys->tmp[thd], sys->pressureFluctuation, sys->pressureFluctuation);

}
void LocallyAveragedCompressibility(System *sys)
{
  // then calculate locally averaged density over compact space
#pragma simd
  for(int i = 0; i<(sys->NumberOfFluidParticles+sys->NumberOfBoundaryParticles); i++)
    {
      sys->meanCompressibility[i] = 0.0;
    }
  // calculate sum for <density>
  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        double Volumei = sys->mass[i]/sys->density[i];
        double Volumej = sys->mass[j]/sys->density[j];
        double pij = sys->DynamicPressure[i] - sys->DynamicPressure[j];
        double pi = 1.0/(sys->SpeedOfSound*sys->SpeedOfSound*sys->density0 + sys->DynamicPressure[i]);
        double pj = 1.0/(sys->SpeedOfSound*sys->SpeedOfSound*sys->density0 + sys->DynamicPressure[j]);
        sys->tmp[threadid][i] -= Volumej*pij*sys->w[k];
        if(j<sys->NumberOfFluidParticles)
          sys->tmp[threadid][j] += Volumei*pij*sys->w[k];
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles,
          sys->tmp[thd],
          sys->meanCompressibility,
          sys->meanCompressibility);

}

void LocallyAveragedSpecificVolume(System *sys)
{
  // then calculate locally averaged density over compact space
#pragma simd
  for(int i = 0; i<(sys->NumberOfFluidParticles+sys->NumberOfBoundaryParticles); i++)
    {
      sys->specificVolume[i] = 1.0/sys->density0;//sys->alpha_kernelFour*sys->mass[i];
    }

  // calculate sum for <density>
  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*sys->NumberOfFluidParticles);
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];

        sys->specificVolume[i] += sys->mass[j]*sys->w[k]/(sys->density[j]*sys->density[j]);
        if(j<sys->NumberOfFluidParticles)
          sys->specificVolume[j] += sys->mass[i]*sys->w[k]/(sys->density[i]*sys->density[i]);
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles,
          sys->tmp[thd],
          sys->specificVolume,
          sys->specificVolume);

}

void TotalKineticEnergy(System *sys)
{
  sys->totalKineticEnergy = 0.0;
  sys->potentialEnergy = 0.0;
  for(int i = 0; i < sys->NumberOfFluidParticles; i++)
    {
      sys->totalKineticEnergy += 0.5*sys->mass[i]*(sys->Velocity[i][0]*sys->Velocity[i][0]+sys->Velocity[i][1]*sys->Velocity[i][1]);
      sys->potentialEnergy += sys->mass[i]*sys->Gravity*sys->Position[i][1];
    }
}


void EnstrophyPowerCalculate(System *sys)
{
  sys->enstrophyPower = 0.0;
  double tmp = 0.0;
#pragma omp parallel for reduction (+:tmp)
  for(int i = 0; i < sys->NumberOfFluidParticles; i++)
    {
      tmp -= sys->mass[i]*sys->viscosity[i]*sys->approxVorticity[i]*sys->approxVorticity[i];
    }

  sys->enstrophyPower=tmp;
}

void PowerDissipation(System *sys)
{
  sys->dissipationPower = 0.0;
  double tmp = 0.0;
  for(int i = 0; i < sys->NumberOfFluidParticles; i++)
    {
      tmp -= 0.5*sys->mass[i]*sys->specificDissipationPower[i];
    }
  sys->dissipationPower = tmp;
}

void NeighborsToPressurePoint(System *sys)
{
  double dx[2];
  double temps = 0.0;
  sys->pointNeighbors = 0;
  int pniac = 0;
  double coeff0 = exp(-9.0);
  double coeff1 = 0.5 - 5.0*exp(-9.0);
  double h0 = 1.2*sys->dx;
  double tempVar = 0.0;
  double tempW = 0.0;
  sys->wp[0] = 0.0;
  //   memset(sys->wp, 0, sizeof(double)*sys->MaxNumberOfParticles);
  for(int i = 0; i < sys->NumberOfFluidParticles; i++)
    {
      dx[0] = sys->Position[i][0] - sys->pressureProbeX;
      dx[1] = sys->Position[i][1] - sys->pressureProbeY;
      double r2 = dx[0]*dx[0] + dx[1]*dx[1];
      double r = sqrt(r2);
      double h = sys->hsml[i];
      const double q = r/h;
      if((r<=(3.0*h))&&(sys->DynamicPressure[i]>0.0))
        {
          double q2 = q*q;
          double pp = (2.0-q)*(2.0-q)*(2.0-q)*(2.0-q)*(2.0-q);
          sys->pj_pair[pniac] = i;
          sys->tempPointPressure[pniac] = sys->DynamicPressure[i]*(exp(-r*r/(h*h)) - exp(-9.0))*sys->mass[i]/(2.0*M_PI*h*h*coeff1*sys->density[i]);
          //       sys->tempPointPressure[pniac] = sys->DynamicPressure[i];
          sys->tempPointVolume[pniac] = sys->mass[i]/sys->density[i];
          sys->wp[pniac] = (exp(-r*r/(h*h)) - exp(-9.0))*sys->mass[i]/(2.0*M_PI*h*h*coeff1*sys->density[i]);
          /* integrate over region around pressure probe point */
          tempVar += sys->tempPointPressure[pniac];
          tempW += sys->wp[pniac];
          pniac++;
        }
    }
  sys->pointNeighbors = pniac;

  if(sys->pointNeighbors)
    {
      sys->pressureNorm = tempW;
    }
  else
    {
      sys->pressureNorm = 1.0;
    }
  if(sys->pointNeighbors)
    sys->pointPressure = tempVar/sys->pressureNorm;
}

void NeighborsToPressurePointTwo(System *sys)
{
  for(int i = 0; i < sys->MaxNumberOfParticles; i++)
    {
      sys->tempPointPressureTwo[i] = 0.0;
    }

  double dx[2];
  sys->pointNeighborsTwo = 0;
  int pniac = 0;
  memset(sys->wp, 0, sizeof(double)*sys->MaxNumberOfParticles);
  for(int i = 0; i < sys->NumberOfFluidParticles; i++)
    {
      dx[0] = sys->Position[i][0]- sys->pressureProbeXX;
      dx[1] = sys->Position[i][1]-sys->pressureProbeYY;
      double r2 = dx[0]*dx[0] + dx[1]*dx[1];
      double r = sqrt(r2);
      double h = sys->hsml[i];
      const double q = r/h;
      if((r<(2.0*h))&&(sys->DynamicPressure[i]>0.0))
        {
          sys->pj_pair[pniac] = i;
          sys->tempPointPressureTwo[pniac] = sys->DynamicPressure[i];
          pniac++;
        }
    }
  sys->pointNeighborsTwo = pniac;
}

// void EnsembleAveragedPointPressure(System *sys)
// {
//   sys->pointPressure = 0.0;
//   double tmp = 0.0;
// #pragma omp parallel for reduction (+:tmp)
//    for(int i = 0; i < sys->pointNeighbors; i++)
//    {
//       tmp += sys->tempPointPressure[i]/sys->pointNeighbors;
//    }
//
//    sys->pointPressure=tmp;
// }

void EnsembleAveragedPointPressureTwo(System *sys)
{
  sys->pointPressureTwo = 0.0;
  double tmp = 0.0;
#pragma omp parallel for reduction (+:tmp)
  for(int i = 0; i < sys->pointNeighborsTwo; i++)
    {
      tmp += sys->tempPointPressureTwo[i]/sys->pointNeighborsTwo;
    }

  sys->pointPressureTwo=tmp;
}


void SpringBoundaryForces(System *sys)
{
  const double r0 = 3.0*sys->dx; /*initial particle spacing*/

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        double dx[2], dv[2], initX[2];
        const int i = sys->i_pair[k];
        const int j = sys->j_pair[k];
        double h = 0.5*(sys->hsml[i] + sys->hsml[j]);
        double eta2 = 0.01*h*h;

        if((sys->itype[i]==1)&&((sys->itype[j]==2)||(sys->itype[j]==3)||(sys->itype[j]==5))) // fixed and moving boundary particles
          {
            double rr = 0.0;
            double rvij = 0.0;
            for(int d = 0; d < sys->dim; d++)
              {
                dx[d] = sys->Position[i][d] - sys->Position[j][d];
                dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
                rr += dx[d]*dx[d];
                rvij += dx[d]*dv[d];
              }
            rr = sqrt(rr);
            double nu = rvij/rr ;//+ eta2); // replace rijk with rr*rr
            if((rr < r0))
              {
                double sumX[2];
                for(int d = 0; d < sys->dim; d++)
                  {
                    // r-part
                    sys->tmp[threadid][2*i+d] -= sys->springConstant0*(r0-rr)*dx[d]/(rr*sys->mass[i]);
                    // v-part
                    sys->tmp[threadid][2*i+d] -= sys->dampingConstant0*rvij*dx[d]/(rr*rr*sys->mass[i]);
                  }
              }

          }
      }
  }
  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles*2, sys->tmp[thd], sys->dvdt[0], sys->dvdt[0]);
}


void SerialSpringBoundaryForces(System *sys)
{
  const double r0 = sys->dx; /*initial particle spacing*/
  for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
    {
      double dx[2], dv[2];
      const int i = sys->i_pair[k];
      const int j = sys->j_pair[k];
      if((sys->itype[i]==1)&&(sys->itype[j]==2)) // fixed and moving boundary particles
        {
          double rr = 0.0;
          for(int d = 0; d < sys->dim; d++)
            {
              dx[d] = sys->Position[i][d] - sys->Position[j][d];
              dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
              rr += dx[d]*dx[d];
            }
          rr = sqrt(rr);
          if(rr < r0)
            {
              for(int d = 0; d < sys->dim; d++)
                {
                  sys->dvdt[i][d] -= (sys->springConstant0*dx[d] + sys->dampingConstant0*dv[d])/sys->mass[i];
                }
            }
        }
    }
}

void SerialComputeSpecificBoundaryDissipationRate(System *sys)
{
  const double r0 = sys->dx; /*initial particle spacing*/
  for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
    {
      double dx[2], dv[2];
      const int i = sys->i_pair[k];
      const int j = sys->j_pair[k];
      if((sys->itype[i]==1)&&((sys->itype[j]==2)||(sys->itype[j]==3))) // fixed and moving boundary particles
        {
          double rr = 0.0;
          double rvij = 0.0;
          for(int d = 0; d < sys->dim; d++)
            {
              dx[d] = sys->Position[i][d] - sys->Position[j][d];
              dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
              rr += dx[d]*dx[d];
              rvij +=dx[d]*dv[d];
            }
          rr = sqrt(rr);
          if(rr < r0)
            {
              for(int d = 0; d < sys->dim; d++)
                {
                  sys->specificBoundaryDissipationRate[i] -= sys->springConstant0*rvij*rvij/(rr*rr);
                }
            }
        }
    }
}

// void ComputeSpecificBoundaryDissipationRate(System *sys)
// {
//   const double r0 = sys->dx; /*initial particle spacing*/
//
//   int MaxThreads = 1;
//
// #pragma omp parallel
//   {
//
//     // Determinate the number of threads. Used afterwards for the reduction
//     int threadid = omp_get_thread_num();
//     if(threadid == 0)
//       MaxThreads = omp_get_num_threads();
//     memset(sys->tmp[threadid], 0, sizeof(double)*sys->NumberOfFluidParticles);
// #pragma omp for
//       for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
//         {
//           double dx[2], dv[2];
//           const int i = sys->i_pair[k];
//           const int j = sys->j_pair[k];
// 	  double h = 0.5*(sys->hsml[i] + sys->hsml[j]);
//           double eta2 = 0.01*h*h;
//
//           if((sys->itype[i]==1)&&((sys->itype[j]==2)||(sys->itype[j]==3))) // fixed and moving boundary particles
//             {
//               double rr = 0.0;
// 	      double rvij = 0.0;
//               for(int d = 0; d < sys->dim; d++)
//                 {
//                   dx[d] = sys->Position[i][d] - sys->Position[j][d];
//                   dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
//                   rr += dx[d]*dx[d];
// 		  rvij += dx[d]*dv[d];
//                 }
//               rr = sqrt(rr);
//               if((rr < r0))
// 	      {
// 		sys->tmp[threadid][i] -= sys->dampingConstant0*rvij*rvij/(rr*rr);
// 	      }
//
// 	    }
//         }
//   }
//   for(int thd=0;thd<MaxThreads;thd++)
//     vdAdd(sys->NumberOfFluidParticles, sys->tmp[thd], sys->specificBoundaryDissipationRate, sys->specificBoundaryDissipationRate);
// }

void ComputeBoundaryDissipationRate(System *sys)
{
  sys->boundaryDissipationRate = 0.0;
  double tmp = 0.0;
  // #pragma omp parallel for reduction (+:tmp)
  for(int i = 0; i < sys->NumberOfFluidParticles; i++)
    {
      sys->boundaryDissipationRate += sys->specificBoundaryDissipationRate[i];
    }

  //    sys->boundaryDissipationRate=tmp;
}


void ComputeSpecificDissipationRate(System *sys)
{
  /* This is the dissipation rate within the fluid bulk.*/

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        double dx[2], dv[2], dvij[2], dvji[2];
        double rdwij = 0.0;
        double velocity2 = 0.0;
        for(int d = 0; d < sys->dim; d++)
          {
            dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
            dvij[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
            dvji[d] = sys->Velocity[j][d] - sys->Velocity[i][d];
            dx[d] = sys->Position[j][d] - sys->Position[i][d];
            rdwij += sys->dwcdx[k][d]*dx[d];
            velocity2 += dv[d]*dv[d];
          }

        double h = 0.5*(sys->hsml[i] + sys->hsml[j]);
        double eta2 = 0.01*h*h;

        double rdw = rdwij/(sys->rij2[k] + eta2);

        double rhoi = sys->localDensity[j]/(sys->density[i]*sys->density[i]*sys->density[j]);
        double rhoj = sys->localDensity[i]/(sys->density[j]*sys->density[j]*sys->density[i]);

        double effViscosityij = 2.0*(sys->density[i]*sys->viscosity[i])*(sys->density[j]*sys->viscosity[j])/(sys->density[i]*sys->viscosity[i] + sys->density[j]*sys->viscosity[j]);

        double alphaij = sys->mass[i]*sys->mass[j]*effViscosityij*(rhoi+rhoj)*rdw;

        sys->tmp[threadid][i] -= alphaij*velocity2;
        if(j<sys->NumberOfFluidParticles)
          sys->tmp[threadid][j] -= alphaij*velocity2;
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles,
          sys->tmp[thd],
          sys->specificDissipationRate,
          sys->specificDissipationRate);
}

void ComputeDissipationRate(System *sys)
{
  for(int i = 0; i<sys->NumberOfFluidParticles; i++)
    {
      sys->specificDissipationRateTemp[i] = 0.0;
    }
  sys->dissipationRate = 0.0;
  int count =0;
  double tmp = 0.0;
  // #pragma omp parallel for reduction (+:tmp)
  if(sys->disspationControlVolume==3)
    {
      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
        {
          if((sys->Position[i][0]>3.0)&&(sys->Position[i][0]<3.5)&&(sys->Position[i][1]<0.16)&&(sys->Position[i][1]<0.35))
            {
              sys->specificDissipationRateTemp[count] = sys->specificDissipationRate[i];
              sys->dissipationRate += sys->specificDissipationRate[i];
              count++;
            }

        }
    }
  /* dambreak control volume*/
  else if(sys->disspationControlVolume==2)
    {
      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
        {
          if((sys->Position[i][0]>1.0)&&(sys->Position[i][0]<1.5)&&(sys->Position[i][1]>0.15)&&(sys->Position[i][1]<0.5))
            {
              sys->specificDissipationRateTemp[count] = sys->specificDissipationRate[i];
              sys->dissipationRate += sys->specificDissipationRate[i];
              count++;
            }

        }
    }
  /* hydrostatic control volume */
  else if(sys->disspationControlVolume==1)
    {
      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
        {
          if((sys->Position[i][0]>1.0)&&(sys->Position[i][0]<1.5)&&(sys->Position[i][1]<0.06)&&(sys->Position[i][1]<0.5))
            {
              sys->specificDissipationRateTemp[count] = sys->specificDissipationRate[i];
              sys->dissipationRate += sys->specificDissipationRate[i];
              count++;
            }

        }
    }
  else
    {
      printf("Do nothing: do not compute the dissipation rate\n");
    }
}

void ComputeSpecificTurbulentDissipationRate(System *sys)
{
  /* This is the dissipation rate within the fluid bulk.*/

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        double dx[2], dv[2], dvij[2], dvji[2];
        double rdwij = 0.0;
        double velocity2 = 0.0;
        for(int d = 0; d < sys->dim; d++)
          {
            dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
            dvij[d] = sys->velocityFluctuation[i][d] - sys->velocityFluctuation[j][d];
            dvji[d] = sys->Velocity[j][d] - sys->Velocity[i][d];
            dx[d] = sys->Position[j][d] - sys->Position[i][d];
            rdwij += sys->dwcdx[k][d]*dx[d];
            velocity2 += dvij[d]*dvij[d];
          }

        double h = 0.5*(sys->hsml[i] + sys->hsml[j]);
        double eta2 = 0.01*h*h;

        double rdw = rdwij/(sys->rij2[k] + eta2);

        double rhoi = sys->localDensity[j]/(sys->density[i]*sys->density[i]*sys->density[j]);
        double rhoj = sys->localDensity[i]/(sys->density[j]*sys->density[j]*sys->density[i]);

        double effViscosityij = 2.0*(sys->density[i]*sys->viscosity[i])*(sys->density[j]*sys->viscosity[j])/(sys->density[i]*sys->viscosity[i] + sys->density[j]*sys->viscosity[j]);

        double alphaij = 0.5*sys->mass[i]*sys->mass[j]*effViscosityij*(rhoi+rhoj)*rdw;

        sys->tmp[threadid][i] += alphaij*velocity2;
        if(j<sys->NumberOfFluidParticles)
          sys->tmp[threadid][j] += alphaij*velocity2;
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles,
          sys->tmp[thd],
          sys->specificTurbulentDissipationRate,
          sys->specificTurbulentDissipationRate);
}

void ComputeTurbulentKineticEnergy(System *sys)
{
  for(int i = 0; i < sys->NumberOfFluidParticles; i++)
    {
      double rv2 = 0.0;
      for(int d = 0; d < sys->dim; d++)
        rv2 += sys->velocityFluctuation[i][d]*sys->velocityFluctuation[i][d];
      sys->turbulentKineticEnergy[i] = 0.5*rv2;
    }
}


void ControlVolumeProperties(System *sys)
{
  double xMin = 0.0, xMax;
  if(sys->movingBoundary==1) /* suing pecletNumber = 1/slope */
    {
      xMax = sys->FlatBottomLength + sys->LocalDepth/sys->pecletNumber;
    }
  else
    {
      xMax = sys->FlatBottomLength;
    }
  for(int i = 0; i<(sys->controlVolumeNumber-1); i++)
    {
      sys->controlDissipationRate[i+1] = 0.0;
      sys->controlKineticEnergy[i+1] = 0.0;
    }
  int count =0;
  double tmp = 0.0;

  for(int j = 0; j<=(sys->controlVolumeNumber-1); j++)
    sys->controlVolumeX[j+1] = xMin + j*(xMax - xMin)/(sys->controlVolumeNumber-1);

  // #pragma omp parallel for reduction (+:tmp)
  if(sys->disspationControlVolume==2) /* dam break problem*/
    {
      for(int j = 0; j<=(sys->controlVolumeNumber-2); j++)
        {
          for(int i = 0; i < sys->NumberOfFluidParticles; i++)
            {
              if((sys->Position[i][0]>sys->controlVolumeX[sys->controlVolumeNumber-j-1])&&(sys->Position[i][0]<sys->controlVolumeX[sys->controlVolumeNumber-j])&&(sys->Position[i][1]>=(0.5*sys->LocalDepth)))
                {
                  sys->controlDissipationRate[sys->controlVolumeNumber-j-1] += sys->specificDissipationRate[i];
                  sys->controlKineticEnergy[sys->controlVolumeNumber-j-1] += sys->mass[i]*(sys->Velocity[i][0]*sys->Velocity[i][0] + sys->Velocity[i][1]*sys->Velocity[i][1]);
                }
            }

        }
    }
  if(sys->disspationControlVolume==3) /* breaking wave problem*/
    {
      for(int j = 0; j<=(sys->controlVolumeNumber-2); j++)
        {
          for(int i = 0; i < sys->NumberOfFluidParticles; i++)
            {
              if((sys->Position[i][0]>sys->controlVolumeX[sys->controlVolumeNumber-j-1])&&(sys->Position[i][0]<sys->controlVolumeX[sys->controlVolumeNumber-j])&&(sys->Position[i][1]>sys->LocalDepth))
                {
                  sys->controlDissipationRate[sys->controlVolumeNumber-j-1] += sys->specificDissipationRate[i];
                  sys->controlKineticEnergy[sys->controlVolumeNumber-j-1] += sys->mass[i]*(sys->Velocity[i][0]*sys->Velocity[i][0] + sys->Velocity[i][1]*sys->Velocity[i][1]);
                }
            }

        }
    }

  //    else
  //    {
  //      printf("Do nothing: do not compute the dissipation rate\n");
  //    }
}
void PressureKernelNormalizationFactor(System *sys)
{
  sys->pressureKernelNorm = 1.0;
  for(int k = 0; k < sys->pointNeighbors; k++)
    {
      int j = sys->pj_pair[k];

      sys->pressureKernelNorm += sys->wp[k]*sys->tempPointVolume[k];
    }
}


void SpongeLayerWaveReflectionAbsorption(System *sys)
{
  if((sys->movingBoundary==1)&&(sys->waveAbsorption==1)) // simulations with gravity
    {
#pragma omp parallel for
      for(int i=0;i<sys->NumberOfFluidParticles;i++)
        {
          sys->dvdt[i][0] -= sys->spongeFactor[i]*sys->Velocity[i][0];
          sys->dvdt[i][1] -= 0.0;//sys->spongeFactor[i]*sys->Velocity[i][1];
        }
    }
}

void ApproximateVortcity(System *sys)
{
  const int niac = sys->NumberOfInteractingParticles;
  for(int i = 0; i< (sys->NumberOfFluidParticles+sys->NumberOfBoundaryParticles); i++)
    {
      sys->approxVorticity[i] = 0.0;
    }
  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];

        double dx[2], dv[2];
        double h = 0.5*(sys->hsml[i] + sys->hsml[j]);
        double eta2 = 0.01*h*h;

        double invR2 = 1.0/(sys->rij2[k] + eta2);
        for(int d = 0; d < sys->dim; d++)
          {
            dx[d] = sys->Position[i][d] - sys->Position[j][d];
            dv[d] = sys->Velocity[i][d] - sys->Velocity[j][d];
          }
        double volumei = sys->mass[i]/sys->density[i];
        double volumej = sys->mass[j]/sys->density[j];

        double aij = (dv[1]*dx[0] - dv[0]*dx[1])*sys->w[k]*invR2;
        sys->tmp[threadid][i] += volumej*aij;
        if(j<sys->NumberOfFluidParticles)
          sys->tmp[threadid][j] += volumei*aij;
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles, sys->tmp[thd], sys->approxVorticity, sys->approxVorticity);
}

void VorticityCompute(System *sys)
{
  for(int i = 0; i< (sys->NumberOfFluidParticles+sys->NumberOfBoundaryParticles); i++)
    {
      sys->vorticity[i] = 0.0;
    }

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++) // niac = Total number of interacting particles
      {
        double dx[2], dvij[2], dvji[2];
        int i = sys->i_pair[k]; /* interaction pair : particle i and particle j*/
        int j = sys->j_pair[k];
        for(int d = 0; d < sys->dim; d++)
          {
            dvij[d] = sys->Velocity[i][d] - sys->localVelocity[j][d];
            dvji[d] = sys->Velocity[j][d] - sys->localVelocity[i][d];
          }
        for(int d = 0; d < sys->dim; d++)
          {


            double vdwij = dvij[2]*sys->dwcdx[k][1] - dvij[1]*sys->dwcdx[k][2];
            double vdwji = dvji[2]*sys->dwcdx[k][1] - dvji[1]*sys->dwcdx[k][2];

            double volumej = sys->mass[j]/sys->density[j];
            double volumei = sys->mass[i]/sys->density[i];

            double reducedVolumej = volumej*sys->localDensity[j]/sys->density[i];
            double reducedVolumei = volumei*sys->localDensity[i]/sys->density[j];


            sys->tmp[threadid][i] += volumej*vdwij;
            if(j<sys->NumberOfFluidParticles)
              sys->tmp[threadid][j] -= volumei*vdwji;

          }
      }
  }
  for(int thd=0;thd<MaxThreads;thd++) {
    vdAdd(sys->NumberOfFluidParticles, sys->tmp[thd], sys->vorticity, sys->vorticity);
  }
}





/********************************************************************/
/* This part of the code implements the particle parcking algorithm*/

void CalculateNormalizationFactor(System *sys)
{
  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles); // check the 2, no need for this. remove
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        double Volume0 = sys->mass[0]/sys->density0;
        sys->tmp[threadid][i] += sys->w[k]*Volume0;
        if(j<sys->NumberOfFluidParticles)
          sys->tmp[threadid][j] += sys->w[k]*Volume0;
      }

  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles, sys->tmp[thd], sys->normalizationFactor, sys->normalizationFactor);
}

void CalculateGradientOfNormalizationFactor(System *sys)
{
  for(int i = 0; i<(sys->NumberOfFluidParticles+sys->NumberOfBoundaryParticles); i++)
    {
      for(int d = 0; d<2; d++)
        sys->gradNormFactor[i][d] = 0.0;
    }

  int MaxThreads = 1;

#pragma omp parallel
  {

    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        double Volume0 = sys->mass[0]/sys->density0;
        for(int d = 0; d<2; d++)
          {
            sys->tmp[threadid][2*i+d] += sys->dwcdx[k][d]*Volume0;
            if(j<sys->NumberOfFluidParticles)
              sys->tmp[threadid][2*j+d] -= sys->dwcdx[k][d]*Volume0;
          }
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vdAdd(sys->NumberOfFluidParticles*2, sys->tmp[thd], sys->gradNormFactor[0], sys->gradNormFactor[0]);
}

void InitializationMomentumEquation(System *sys)
{

  int MaxThreads = 1;

#pragma omp parallel
  {
    // Determinate the number of threads. Used afterwards for the reduction
    int threadid = omp_get_thread_num();
    if(threadid == 0)
      MaxThreads = omp_get_num_threads();
    memset(sys->tmp[threadid], 0, sizeof(double)*2*sys->NumberOfFluidParticles);
#pragma omp for
    for(int k = 0; k < sys->NumberOfInteractingParticles; k++)
      {
        int i = sys->i_pair[k];
        int j = sys->j_pair[k];
        double Volume0 = sys->mass[0]/sys->density0;
        double alphab = 0.02;
        double betab = 2.0*sys->Gravity*sys->LocalDepth;
        double dampb = alphab*sqrt(betab/Volume0);
        for(int d = 0; d<2; d++)
          {
            // minimization term
            sys->tmp[threadid][2*i+d] -= betab*sys->gradNormFactor[i][d];
            // damping term
            sys->tmp[threadid][2*i+d] -= dampb*sys->Velocity[i][d];
            if(j<sys->NumberOfFluidParticles) {
              sys->tmp[threadid][2*j+d] -= betab*sys->gradNormFactor[j][d];
              sys->tmp[threadid][2*j+d] -= dampb*sys->Velocity[j][d];
            }
          }
      }
  }

  for(int thd=0;thd<MaxThreads;thd++)
    vzAdd(sys->NumberOfFluidParticles,
          (MKL_Complex16 *)sys->tmp[thd],
          (MKL_Complex16 *)sys->dvdt[0],
          (MKL_Complex16 *)sys->dvdt[0]);

}
void InitialDensityAndPressure(System *sys)
{

#pragma omp parallel for
  for(int i = 0; i <sys->NumberOfFluidParticles; i++)
    {
      double p0 = sys->density0*sys->SpeedOfSound*sys->SpeedOfSound/sys->AdiabaticConstant;
      sys->DynamicPressure[i] = sys->density0*sys->Gravity*(sys->LocalDepth-sys->Position[i][1]);
      sys->density[i] = sys->density0;//*(1.0 + pow(sys->DynamicPressure[i]/p0,1.0/sys->AdiabaticConstant));
    }
}

/*******************************************************************/



void TerminateInitialization(System *sys, double t)
{
  double energyThreshold = 1.0e-7;
  if((sys->flagInitSimul==0)&&(sys->totalKineticEnergy<=energyThreshold)&&(t>0.1))
    {
      printf("Quasi-equilibrium state has been attained:\n");
      printf("Simulation has now been termined:\n");
      exit(1);
    }
  else
    {
      printf("Now running simulation of the problem:\n");
    }
}

/*********************************************************/
/* calculate the density consistency function*/

void DensityConsistencyFunctions(System *sys)
{
  sys->consistencyDensityOne = 0.0;//sys->density0;
  sys->consistencyDensityTwo = 0.0;//sys->density0;
  sys->consistencyDensityThree = 0.0;//sys->density0;
  double tmp = 0.0;
#pragma omp parallel for reduction (+:tmp)
  for(int i = 0; i<sys->NumberOfFluidParticles; i++)
    sys->consistencyDensityOne += sys->localDensity[i]*sys->localDensity[i]/(sys->density[i]*sys->NumberOfFluidParticles);

  sys->consistencyDensityOne = tmp;

  tmp = 0.0;
#pragma omp parallel for reduction (+:tmp)
  for(int i = 0; i<sys->NumberOfFluidParticles; i++)
    tmp += sys->localDensity[i]/sys->NumberOfFluidParticles;

  sys->consistencyDensityTwo = tmp;

  tmp = 0.0;
#pragma omp parallel for reduction (+:tmp)
  for(int i = 0; i<sys->NumberOfFluidParticles; i++)
    tmp += sys->density[i]/sys->NumberOfFluidParticles;

  sys->consistencyDensityThree = tmp;

  sys->densityUncertainty = sys->consistencyDensityOne - sys->consistencyDensityTwo*sys->consistencyDensityTwo/sys->consistencyDensityThree;
}

/*********************************************************/



__inline__ void MovingBoundary(System *sys, double t)
{
  sys->LeftBoundaryPosition = 0.5*sys->stroke*(1.0 - cos(sys->omega*t));
  sys->LeftBoundaryVelocity = 0.5*sys->stroke*sys->omega*sin(sys->omega*t);
  for(int i = 0; i <sys->NumberOfMovingBoundaryParticles; i++)
    {
      const int particleCount = sys->NumberOfFluidParticles + sys->NumberOfFixedBoundaryParticles + sys->NumberOfSolidParticles;
      sys->Position[particleCount + i][0] = sys->LeftBoundaryPosition;
      sys->Velocity[particleCount + i][0] = sys->LeftBoundaryVelocity;
      sys->Position[particleCount + i][1] = 0;
      sys->Velocity[particleCount + i][1] = 0;
    }
}

void derivatives(System *sys, double t)
{
  int ii;
  double xl = 1.0, vl = 1.0;


  /* after moving particles inflow and outflow */

  /* figure out why this if-statement gives an error when " if{}"*/
  /* if(t > 0.00000001) */
  memset(sys->ddensitydt, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->dpdt, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->dvdt[0], 0, sizeof(double)*sys->MaxNumberOfParticles*2);
  memset(sys->ViscousStressxx, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->ViscousStressxy, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->ViscousStressyx, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->ViscousStressyy, 0, sizeof(double)*sys->MaxNumberOfParticles);

  memset(sys->ddensitydt_approx, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->dpdt_approx, 0, sizeof(double)*sys->MaxNumberOfParticles);
  memset(sys->dvdt_approx[0], 0, sizeof(double)*sys->MaxNumberOfParticles*2);




  if(t>0) {
    inoutflow(sys);
  }


  /* with inflow/outflow boundaries, call inputbp.h every timestep */
  inputbp(sys, t);
  //   MovingBoundary(sys, t);

  /* track particles within space from inflow boundary */
  sys->inflowParticles = 0;

#pragma omp parallel for
  for(int i = 0; i < sys->NumberOfFluidParticles; i++)
    {
      if(sys->Position[0][i] <= 0.1)
        {
          sys->i_inflow[sys->inflowParticles] = i;
          sys->inflowParticles++;
        }
    }


  /* boundary particles also need vxsph */

#pragma omp parallel for
  for(int i=0;i<sys->NumberOfBoundaryParticles;i++) {
    for(int d=0;d<2;d++)
      sys->Velocity_xsph[sys->NumberOfFluidParticles + sys->NumberOfSolidParticles + i][d] = sys->Velocity[sys->NumberOfFluidParticles+sys->NumberOfSolidParticles +i][d];
  }

  /* create virtual particles within 2h*/
  ii = sys->NumberOfFluidParticles + sys->NumberOfBoundaryParticles + sys->NumberOfSolidParticles;
  if(sys->virtualType==1) // waves on slope
    {
      //  left boundary
      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0]-sys->LeftBoundaryPosition);
        if(dx < 2.0 * sys->hsml[i]) {
          sys->Position[ii][0] = sys->LeftBoundaryPosition-dx;
          sys->Position[ii][1] = sys->Position[i][1];
          sys->Velocity[ii][0] = -sys->Velocity[i][0];
          sys->Velocity[ii][1] = sys->Velocity[i][1];
          sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
          sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
          sys->mass[ii] = sys->mass[i];
          sys->density[ii] = sys->density[i];
          sys->Pressure[ii] = sys->Pressure[i];
          sys->Energy[ii] = sys->Energy[i];
          sys->hsml[ii] = sys->hsml[i];
          sys->normalizationFactor[ii] = sys->normalizationFactor[i];
          sys->DynamicPressure[ii] = sys->DynamicPressure[i];
          sys->localDensity[ii] = sys->localDensity[i];
          sys->localPressure[ii] = sys->localPressure[i];
          sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
          sys->localVelocity[ii][1] = sys->localVelocity[i][1];
          sys->vorticity[ii] = sys->vorticity[i];
          sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

          sys->itype[ii] = 4;
          ii ++;
        }
      }

      //Bottom boundary

      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
        {
          if((sys->Position[i][1] < 2.0*sys->hsml[i]) && (sys->Position[i][0] <= sys->FlatBottomLength))
            {
              sys->Position[ii][0] = sys->Position[i][0];
              sys->Position[ii][1] = -sys->Position[i][1];
              sys->Velocity[ii][0] = sys->Velocity[i][0];
              sys->Velocity[ii][1] = -sys->Velocity[i][1];
              sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
              sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
              sys->mass[ii] = sys->mass[i];
              sys->density[ii] = sys->density[i];
              sys->Pressure[ii] = sys->Pressure[i];
              sys->Energy[ii] = sys->Energy[i];
              sys->hsml[ii] = sys->hsml[i];
              sys->normalizationFactor[ii] = sys->normalizationFactor[i];
              sys->DynamicPressure[ii] = sys->DynamicPressure[i];
              sys->localDensity[ii] = sys->localDensity[i];
              sys->localPressure[ii] = sys->localPressure[i];
              sys->localVelocity[ii][0] = sys->localVelocity[i][0];
              sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
              sys->vorticity[ii] = sys->vorticity[i];
              sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

              sys->itype[ii] = 4;
              ii ++;
            }
        }

      // edges

      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0] - sys->LeftBoundaryPosition);
        double dy = sys->Position[i][1];

        if((dx < 2.0*sys->hsml[i])&&(dy < 2.0*sys->hsml[i]))
          {
            sys->Position[ii][0] = sys->LeftBoundaryPosition-dx;
            sys->Position[ii][1] = -sys->Position[i][1];
            sys->Velocity[ii][0] = -sys->Velocity[i][0];
            sys->Velocity[ii][1] = -sys->Velocity[i][1];
            sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
            sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
            sys->mass[ii] = sys->mass[i];
            sys->density[ii] = sys->density[i];
            sys->Pressure[ii] = sys->Pressure[i];
            sys->Energy[ii] = sys->Energy[i];
            sys->hsml[ii] = sys->hsml[i];
            sys->normalizationFactor[ii] = sys->normalizationFactor[i];
            sys->DynamicPressure[ii] = sys->DynamicPressure[i];
            sys->localDensity[ii] = sys->localDensity[i];
            sys->localPressure[ii] = sys->localPressure[i];
            sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
            sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
            sys->vorticity[ii] = sys->vorticity[i];
            sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

            sys->itype[ii] = 4;
            ii ++;
          }
      }
      //  slopping boundary
      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double denoms = 1 + sys->slope*sys->slope;
        double nums = 1 - sys->slope*sys->slope;
        double yIntercept = sys->slope*sys->FlatBottomLength;
        double dx = (1.0/sqrt(denoms))*fabs(sys->slope*sys->Position[i][0] - sys->Position[i][1] - yIntercept);
        if( (dx < 3.0*sys->hsml[i])&&(sys->Position[i][0]>sys->FlatBottomLength) ) {
          sys->Position[ii][0] = (nums*sys->Position[i][0] + 2.0*sys->slope*sys->Position[i][1] + 2.0*sys->slope*yIntercept)/denoms;
          sys->Position[ii][1] = (2.0*sys->slope*sys->Position[i][0] - nums*sys->Position[i][1] - 2.0*yIntercept)/denoms;
          sys->Velocity[ii][0] = (nums*sys->Velocity[i][0] + 2.0*sys->slope*sys->Velocity[i][1])/denoms;
          sys->Velocity[ii][1] = (2.0*sys->slope*sys->Velocity[i][0] - nums*sys->Velocity[i][1])/denoms;
          sys->Velocity_xsph[ii][0] = (nums*sys->Velocity_xsph[i][0] + 2*sys->slope*sys->Velocity_xsph[i][1])/denoms;
          sys->Velocity_xsph[ii][1] = (2.0*sys->slope*sys->Velocity_xsph[i][0] - nums*sys->Velocity_xsph[i][1])/denoms;
          sys->mass[ii] = sys->mass[i];
          sys->density[ii] = sys->density[i];
          sys->Pressure[ii] = sys->Pressure[i];
          sys->Energy[ii] = sys->Energy[i];
          sys->hsml[ii] = sys->hsml[i];
          sys->normalizationFactor[ii] = sys->normalizationFactor[i];
          sys->DynamicPressure[ii] = sys->DynamicPressure[i];
          sys->localDensity[ii] = sys->localDensity[i];
          sys->localPressure[ii] = sys->localPressure[i];
          sys->localVelocity[ii][0] = (nums*sys->localVelocity[i][0] + 2.0*sys->slope*sys->localVelocity[i][1])/denoms;
          sys->localVelocity[ii][1] = (2.0*sys->slope*sys->localVelocity[i][0] - nums*sys->localVelocity[i][1])/denoms;
          sys->vorticity[ii] = sys->vorticity[i];
          sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

          sys->itype[ii] = 4;
          ii ++;
        }
      }
    }
  else if(sys->virtualType==2)// lid driven cavity
    {
      // //   //  left boundary
      // //   for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
      // //     double dx = fabs(sys->Position[i][0]-0.0*sys->FlatBottomLength);
      // //     if(dx <= 3.0 * sys->hsml[i]) {
      // //       sys->Position[ii][0] = 0.0 - dx;
      // //       sys->Position[ii][1] = sys->Position[i][1];
      // //       sys->Velocity[ii][0] = -sys->Velocity[i][0];
      // //       sys->Velocity[ii][1] = sys->Velocity[i][1];
      // //       sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
      // //       sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
      // //       sys->mass[ii] = sys->mass[i];
      // //       sys->density[ii] = sys->density[i];
      // //       sys->Pressure[ii] = sys->Pressure[i];
      // //       sys->Energy[ii] = sys->Energy[i];
      // //       sys->hsml[ii] = sys->hsml[i];
      // //       sys->normalizationFactor[ii] = sys->normalizationFactor[i];
      // //       sys->DynamicPressure[ii] = sys->DynamicPressure[i];
      // //       sys->velocityFluctuation[ii][0] = -sys->velocityFluctuation[i][0];
      // //       sys->velocityFluctuation[ii][1] = sys->velocityFluctuation[i][1];
      // //       sys->localDensity[ii] = sys->localDensity[i];
      // //       sys->localPressure[ii] = sys->localPressure[i];
      // //       sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
      // //       sys->localVelocity[ii][1] = sys->localVelocity[i][1];
      // //       sys->itype[ii] = 4;
      // //       ii ++;
      // //     }
      // //   }
      //
      //   //Bottom boundary
      //
      //   for(int i = 0; i < sys->NumberOfFluidParticles; i++)
      //     {
      //       if((sys->Position[i][1] < 2.0*sys->hsml[i]) && (sys->Position[i][0] <= sys->FlatBottomLength))
      //      {
      //        sys->Position[ii][0] = sys->Position[i][0];
      //        sys->Position[ii][1] = -sys->Position[i][1];
      //        sys->Velocity[ii][0] = sys->Velocity[i][0];
      //        sys->Velocity[ii][1] = -sys->Velocity[i][1];
      //        sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
      //        sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
      //        sys->mass[ii] = sys->mass[i];
      //        sys->density[ii] = sys->density[i];
      //        sys->Pressure[ii] = sys->Pressure[i];
      //        sys->Energy[ii] = sys->Energy[i];
      //        sys->hsml[ii] = sys->hsml[i];
      //        sys->normalizationFactor[ii] = sys->normalizationFactor[i];
      //        sys->DynamicPressure[ii] = sys->DynamicPressure[i];
      //        sys->velocityFluctuation[ii][0] = sys->velocityFluctuation[i][0];
      //           sys->velocityFluctuation[ii][1] = -sys->velocityFluctuation[i][1];
      //        sys->localDensity[ii] = sys->localDensity[i];
      //        sys->localPressure[ii] = sys->localPressure[i];
      //        sys->localVelocity[ii][0] = sys->localVelocity[i][0];
      //        sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
      //        sys->itype[ii] = 4;
      //        ii ++;
      //      }
      //     }
      //       //top boundary
      //
      //     for(int i = 0; i < sys->NumberOfFluidParticles; i++)
      //     {
      //       double dx = fabs(sys->Position[i][1]-sys->LocalDepth-sys->dx);
      //       if((dx < 2.0*sys->hsml[i]))
      //      {
      //        sys->Position[ii][0] = sys->Position[i][0];
      //        sys->Position[ii][1] = sys->LocalDepth + 2*dx;
      //        sys->Velocity[ii][0] = sys->Velocity[i][0];
      //        sys->Velocity[ii][1] = -sys->Velocity[i][1];
      //        sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
      //        sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
      //        sys->mass[ii] = sys->mass[i];
      //        sys->density[ii] = sys->density[i];
      //        sys->Pressure[ii] = sys->Pressure[i];
      //        sys->Energy[ii] = sys->Energy[i];
      //        sys->hsml[ii] = sys->hsml[i];
      //        sys->normalizationFactor[ii] = sys->normalizationFactor[i];
      //        sys->DynamicPressure[ii] = sys->DynamicPressure[i];
      //        sys->velocityFluctuation[ii][0] = sys->velocityFluctuation[i][0];
      //           sys->velocityFluctuation[ii][1] = -sys->velocityFluctuation[i][1];
      //        sys->localDensity[ii] = sys->localDensity[i];
      //        sys->localPressure[ii] = sys->localPressure[i];
      //        sys->localVelocity[ii][0] = sys->localVelocity[i][0];
      //        sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
      //        sys->itype[ii] = 4;
      //        ii ++;
      //      }
      //     }
      //  // edges
      //
      //   for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
      //     double dx = fabs(sys->Position[i][0] - sys->LeftBoundaryPosition);
      //     double dy = sys->Position[i][1];
      //
      //     if((dx < 2.0*sys->hsml[i])&&(dy < 2.0*sys->hsml[i]))
      //       {
      //      sys->Position[ii][0] = sys->LeftBoundaryPosition-dx;
      //      sys->Position[ii][1] = -sys->Position[i][1];
      //      sys->Velocity[ii][0] = -sys->Velocity[i][0];
      //      sys->Velocity[ii][1] = -sys->Velocity[i][1];
      //      sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
      //      sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
      //      sys->mass[ii] = sys->mass[i];
      //      sys->density[ii] = sys->density[i];
      //      sys->Pressure[ii] = sys->Pressure[i];
      //      sys->Energy[ii] = sys->Energy[i];
      //      sys->hsml[ii] = sys->hsml[i];
      //      sys->normalizationFactor[ii] = sys->normalizationFactor[i];
      //      sys->DynamicPressure[ii] = sys->DynamicPressure[i];
      //      sys->velocityFluctuation[ii][0] = -sys->velocityFluctuation[i][0];
      //         sys->velocityFluctuation[ii][1] = -sys->velocityFluctuation[i][1];
      //      sys->localDensity[ii] = sys->localDensity[i];
      //      sys->localPressure[ii] = sys->localPressure[i];
      //      sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
      //      sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
      //      sys->itype[ii] = 4;
      //      ii ++;
      //       }
      //   }
      //  //  right boundary
      //  for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
      //     double dx = fabs(sys->Position[i][0]-sys->FlatBottomLength);
      //     if(dx < 2.0 * sys->hsml[i]) {
      //       sys->Position[ii][0] = sys->FlatBottomLength + dx;
      //       sys->Position[ii][1] = sys->Position[i][1];
      //       sys->Velocity[ii][0] = -sys->Velocity[i][0];
      //       sys->Velocity[ii][1] = sys->Velocity[i][1];
      //       sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
      //       sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
      //       sys->mass[ii] = sys->mass[i];
      //       sys->density[ii] = sys->density[i];
      //       sys->Pressure[ii] = sys->Pressure[i];
      //       sys->Energy[ii] = sys->Energy[i];
      //       sys->hsml[ii] = sys->hsml[i];
      //       sys->normalizationFactor[ii] = sys->normalizationFactor[i];
      //       sys->DynamicPressure[ii] = sys->DynamicPressure[i];
      //       sys->velocityFluctuation[ii][0] = -sys->velocityFluctuation[i][0];
      //       sys->velocityFluctuation[ii][1] = sys->velocityFluctuation[i][1];
      //       sys->localDensity[ii] = sys->localDensity[i];
      //       sys->localPressure[ii] = sys->localPressure[i];
      //       sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
      //       sys->localVelocity[ii][1] = sys->localVelocity[i][1];
      //       sys->itype[ii] = 4;
      //       ii ++;
      //     }
      //   }

    }

  else if(sys->virtualType==3)// vortex spin down
    {
      //  left boundary
      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0]-0.0*sys->FlatBottomLength);
        if(dx < 2.0 * sys->hsml[i]) {
          sys->Position[ii][0] = 0.0 - dx;
          sys->Position[ii][1] = sys->Position[i][1];
          sys->Velocity[ii][0] = -sys->Velocity[i][0];
          sys->Velocity[ii][1] = sys->Velocity[i][1];
          sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
          sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
          sys->mass[ii] = sys->mass[i];
          sys->density[ii] = sys->density[i];
          sys->Pressure[ii] = sys->Pressure[i];
          sys->Energy[ii] = sys->Energy[i];
          sys->hsml[ii] = sys->hsml[i];
          sys->normalizationFactor[ii] = sys->normalizationFactor[i];
          sys->DynamicPressure[ii] = sys->DynamicPressure[i];
          sys->localDensity[ii] = sys->localDensity[i];
          sys->localPressure[ii] = sys->localPressure[i];
          sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
          sys->localVelocity[ii][1] = sys->localVelocity[i][1];
          sys->vorticity[ii] = sys->vorticity[i];
          sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

          sys->itype[ii] = 4;
          ii ++;
        }
      }

      //Bottom boundary

      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
        {
          if((sys->Position[i][1] < 2.0*sys->hsml[i]) && (sys->Position[i][0] <= sys->FlatBottomLength))
            {
              sys->Position[ii][0] = sys->Position[i][0];
              sys->Position[ii][1] = -sys->Position[i][1];
              sys->Velocity[ii][0] = sys->Velocity[i][0];
              sys->Velocity[ii][1] = -sys->Velocity[i][1];
              sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
              sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
              sys->mass[ii] = sys->mass[i];
              sys->density[ii] = sys->density[i];
              sys->Pressure[ii] = sys->Pressure[i];
              sys->Energy[ii] = sys->Energy[i];
              sys->hsml[ii] = sys->hsml[i];
              sys->normalizationFactor[ii] = sys->normalizationFactor[i];
              sys->DynamicPressure[ii] = sys->DynamicPressure[i];
              sys->localDensity[ii] = sys->localDensity[i];
              sys->localPressure[ii] = sys->localPressure[i];
              sys->localVelocity[ii][0] = sys->localVelocity[i][0];
              sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
              sys->vorticity[ii] = sys->vorticity[i];
              sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

              sys->itype[ii] = 4;
              ii ++;
            }
        }
      //top boundary

      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
        {
          double dx = fabs(sys->Position[i][1]-sys->LocalDepth-sys->dx);
          if((dx < 2.0*sys->hsml[i]))
            {
              sys->Position[ii][0] = sys->Position[i][0];
              sys->Position[ii][1] = sys->LocalDepth + 2*dx;
              sys->Velocity[ii][0] = sys->Velocity[i][0];
              sys->Velocity[ii][1] = -sys->Velocity[i][1];
              sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
              sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
              sys->mass[ii] = sys->mass[i];
              sys->density[ii] = sys->density[i];
              sys->Pressure[ii] = sys->Pressure[i];
              sys->Energy[ii] = sys->Energy[i];
              sys->hsml[ii] = sys->hsml[i];
              sys->normalizationFactor[ii] = sys->normalizationFactor[i];
              sys->DynamicPressure[ii] = sys->DynamicPressure[i];
              sys->localDensity[ii] = sys->localDensity[i];
              sys->localPressure[ii] = sys->localPressure[i];
              sys->localVelocity[ii][0] = sys->localVelocity[i][0];
              sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
              sys->vorticity[ii] = sys->vorticity[i];
              sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

              sys->itype[ii] = 4;
              ii ++;
            }
        }
      // edges

      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0] - sys->LeftBoundaryPosition);
        double dy = sys->Position[i][1];

        if((dx < 2.0*sys->hsml[i])&&(dy < 2.0*sys->hsml[i]))
          {
            sys->Position[ii][0] = sys->LeftBoundaryPosition-dx;
            sys->Position[ii][1] = -sys->Position[i][1];
            sys->Velocity[ii][0] = -sys->Velocity[i][0];
            sys->Velocity[ii][1] = -sys->Velocity[i][1];
            sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
            sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
            sys->mass[ii] = sys->mass[i];
            sys->density[ii] = sys->density[i];
            sys->Pressure[ii] = sys->Pressure[i];
            sys->Energy[ii] = sys->Energy[i];
            sys->hsml[ii] = sys->hsml[i];
            sys->normalizationFactor[ii] = sys->normalizationFactor[i];
            sys->DynamicPressure[ii] = sys->DynamicPressure[i];
            sys->localDensity[ii] = sys->localDensity[i];
            sys->localPressure[ii] = sys->localPressure[i];
            sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
            sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
            sys->vorticity[ii] = sys->vorticity[i];
            sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

            sys->itype[ii] = 4;
            ii ++;
          }
      }
      //  right boundary
      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0]-sys->FlatBottomLength);
        if(dx < 2.0 * sys->hsml[i]) {
          sys->Position[ii][0] = sys->FlatBottomLength + dx;
          sys->Position[ii][1] = sys->Position[i][1];
          sys->Velocity[ii][0] = -sys->Velocity[i][0];
          sys->Velocity[ii][1] = sys->Velocity[i][1];
          sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
          sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
          sys->mass[ii] = sys->mass[i];
          sys->density[ii] = sys->density[i];
          sys->Pressure[ii] = sys->Pressure[i];
          sys->Energy[ii] = sys->Energy[i];
          sys->hsml[ii] = sys->hsml[i];
          sys->normalizationFactor[ii] = sys->normalizationFactor[i];
          sys->DynamicPressure[ii] = sys->DynamicPressure[i];
          sys->localDensity[ii] = sys->localDensity[i];
          sys->localPressure[ii] = sys->localPressure[i];
          sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
          sys->localVelocity[ii][1] = sys->localVelocity[i][1];
          sys->vorticity[ii] = sys->vorticity[i];
          sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

          sys->itype[ii] = 4;
          ii ++;
        }
      }

    }

  else if(sys->virtualType==5)// sloasing problem
    {
      //  left boundary
      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0]-sys->LeftBoundaryPosition);
        if(dx < 2.0 * sys->hsml[i]) {
          sys->Position[ii][0] = sys->LeftBoundaryPosition - dx;
          sys->Position[ii][1] = sys->Position[i][1];
          sys->Velocity[ii][0] = -sys->Velocity[i][0];
          sys->Velocity[ii][1] = sys->Velocity[i][1];
          sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
          sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
          sys->mass[ii] = sys->mass[i];
          sys->density[ii] = sys->density[i];
          sys->Pressure[ii] = sys->Pressure[i];
          sys->Energy[ii] = sys->Energy[i];
          sys->hsml[ii] = sys->hsml[i];
          sys->normalizationFactor[ii] = sys->normalizationFactor[i];
          sys->DynamicPressure[ii] = sys->DynamicPressure[i];
          sys->localDensity[ii] = sys->localDensity[i];
          sys->localPressure[ii] = sys->localPressure[i];
          sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
          sys->localVelocity[ii][1] = sys->localVelocity[i][1];
          sys->vorticity[ii] = sys->vorticity[i];
          sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

          sys->itype[ii] = 4;
          ii ++;
        }
      }

      //   //Bottom boundary
      //
      //   for(int i = 0; i < sys->NumberOfFluidParticles; i++)
      //     {
      //       double dx = fabs(sys->Position[i][0]-sys->LeftBoundaryPosition);
      //       if((sys->Position[i][1] < 3.0*sys->hsml[i]) && (sys->Position[i][0] <= sys->FlatBottomLength))
      //      {
      //        sys->Position[ii][0] = sys->Position[i][0];
      //        sys->Position[ii][1] = -sys->Position[i][1];
      //        sys->Velocity[ii][0] = sys->Velocity[i][0];
      //        sys->Velocity[ii][1] = -sys->Velocity[i][1];
      //        sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
      //        sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
      //        sys->mass[ii] = sys->mass[i];
      //        sys->density[ii] = sys->density[i];
      //        sys->Pressure[ii] = sys->Pressure[i];
      //        sys->Energy[ii] = sys->Energy[i];
      //        sys->hsml[ii] = sys->hsml[i];
      //        sys->normalizationFactor[ii] = sys->normalizationFactor[i];
      //        sys->DynamicPressure[ii] = sys->DynamicPressure[i];
      //        sys->velocityFluctuation[ii][0] = sys->velocityFluctuation[i][0];
      //           sys->velocityFluctuation[ii][1] = -sys->velocityFluctuation[i][1];
      //        sys->localDensity[ii] = sys->localDensity[i];
      //        sys->localPressure[ii] = sys->localPressure[i];
      //        sys->localVelocity[ii][0] = sys->localVelocity[i][0];
      //        sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
      //        sys->itype[ii] = 4;
      //        ii ++;
      //      }
      //     }

      // edges

      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0] - sys->LeftBoundaryPosition);
        double dy = sys->Position[i][1];

        if((dx < 2.0*sys->hsml[i])&&(dy < 2.0*sys->hsml[i]))
          {
            sys->Position[ii][0] = sys->LeftBoundaryPosition-dx;
            sys->Position[ii][1] = -sys->Position[i][1];
            sys->Velocity[ii][0] = -sys->Velocity[i][0];
            sys->Velocity[ii][1] = -sys->Velocity[i][1];
            sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
            sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
            sys->mass[ii] = sys->mass[i];
            sys->density[ii] = sys->density[i];
            sys->Pressure[ii] = sys->Pressure[i];
            sys->Energy[ii] = sys->Energy[i];
            sys->hsml[ii] = sys->hsml[i];
            sys->normalizationFactor[ii] = sys->normalizationFactor[i];
            sys->DynamicPressure[ii] = sys->DynamicPressure[i];
            sys->localDensity[ii] = sys->localDensity[i];
            sys->localPressure[ii] = sys->localPressure[i];
            sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
            sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
            sys->vorticity[ii] = sys->vorticity[i];
            sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

            sys->itype[ii] = 4;
            ii ++;
          }
      }
      // edges

      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0] - sys->RightBoundaryPosition);
        double dy = sys->Position[i][1];

        if((dx < 2.0*sys->hsml[i])&&(dy < 2.0*sys->hsml[i]))
          {
            sys->Position[ii][0] = sys->RightBoundaryPosition+dx;
            sys->Position[ii][1] = -sys->Position[i][1];
            sys->Velocity[ii][0] = -sys->Velocity[i][0];
            sys->Velocity[ii][1] = -sys->Velocity[i][1];
            sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
            sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
            sys->mass[ii] = sys->mass[i];
            sys->density[ii] = sys->density[i];
            sys->Pressure[ii] = sys->Pressure[i];
            sys->Energy[ii] = sys->Energy[i];
            sys->hsml[ii] = sys->hsml[i];
            sys->normalizationFactor[ii] = sys->normalizationFactor[i];
            sys->DynamicPressure[ii] = sys->DynamicPressure[i];
            sys->localDensity[ii] = sys->localDensity[i];
            sys->localPressure[ii] = sys->localPressure[i];
            sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
            sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
            sys->vorticity[ii] = sys->vorticity[i];
            sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

            sys->itype[ii] = 4;
            ii ++;
          }
      }
      //  right boundary
      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0]-sys->RightBoundaryPosition);
        if(dx < 2.0 * sys->hsml[i]) {
          sys->Position[ii][0] = sys->RightBoundaryPosition + dx;
          sys->Position[ii][1] = sys->Position[i][1];
          sys->Velocity[ii][0] = -sys->Velocity[i][0];
          sys->Velocity[ii][1] = sys->Velocity[i][1];
          sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
          sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
          sys->mass[ii] = sys->mass[i];
          sys->density[ii] = sys->density[i];
          sys->Pressure[ii] = sys->Pressure[i];
          sys->Energy[ii] = sys->Energy[i];
          sys->hsml[ii] = sys->hsml[i];
          sys->normalizationFactor[ii] = sys->normalizationFactor[i];
          sys->DynamicPressure[ii] = sys->DynamicPressure[i];
          sys->localDensity[ii] = sys->localDensity[i];
          sys->localPressure[ii] = sys->localPressure[i];
          sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
          sys->localVelocity[ii][1] = sys->localVelocity[i][1];
          sys->vorticity[ii] = sys->vorticity[i];
          sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

          sys->itype[ii] = 4;
          ii ++;
        }
      }

    }

  else if(sys->virtualType==4)// dam break with moving gate
    {
      //  left boundary
      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0]-0.0*sys->FlatBottomLength);
        if(dx < 2.0 * sys->hsml[i]) {
          sys->Position[ii][0] = 0.0 - dx;
          sys->Position[ii][1] = sys->Position[i][1];
          sys->Velocity[ii][0] = -sys->Velocity[i][0];
          sys->Velocity[ii][1] = sys->Velocity[i][1];
          sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
          sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
          sys->mass[ii] = sys->mass[i];
          sys->density[ii] = sys->density[i];
          sys->Pressure[ii] = sys->Pressure[i];
          sys->Energy[ii] = sys->Energy[i];
          sys->hsml[ii] = sys->hsml[i];
          sys->normalizationFactor[ii] = sys->normalizationFactor[i];
          sys->DynamicPressure[ii] = sys->DynamicPressure[i];
          sys->localDensity[ii] = sys->localDensity[i];
          sys->localPressure[ii] = sys->localPressure[i];
          sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
          sys->localVelocity[ii][1] = sys->localVelocity[i][1];
          sys->vorticity[ii] = sys->vorticity[i];
          sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

          sys->itype[ii] = 4;
          ii ++;
        }
      }

      //Bottom boundary

      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
        {
          if((sys->Position[i][1] < 2.0*sys->hsml[i]) && (sys->Position[i][0] <= sys->FlatBottomLength))
            {
              sys->Position[ii][0] = sys->Position[i][0];
              sys->Position[ii][1] = -sys->Position[i][1];
              sys->Velocity[ii][0] = sys->Velocity[i][0];
              sys->Velocity[ii][1] = -sys->Velocity[i][1];
              sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
              sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
              sys->mass[ii] = sys->mass[i];
              sys->density[ii] = sys->density[i];
              sys->Pressure[ii] = sys->Pressure[i];
              sys->Energy[ii] = sys->Energy[i];
              sys->hsml[ii] = sys->hsml[i];
              sys->normalizationFactor[ii] = sys->normalizationFactor[i];
              sys->DynamicPressure[ii] = sys->DynamicPressure[i];
              sys->localDensity[ii] = sys->localDensity[i];
              sys->localPressure[ii] = sys->localPressure[i];
              sys->localVelocity[ii][0] = sys->localVelocity[i][0];
              sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
              sys->vorticity[ii] = sys->vorticity[i];
              sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

              sys->itype[ii] = 4;
              ii ++;
            }
        }

      // edges

      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0] - sys->LeftBoundaryPosition);
        double dy = sys->Position[i][1];

        if((dx < 2.0*sys->hsml[i])&&(dy < 2.0*sys->hsml[i]))
          {
            sys->Position[ii][0] = sys->LeftBoundaryPosition-dx;
            sys->Position[ii][1] = -sys->Position[i][1];
            sys->Velocity[ii][0] = -sys->Velocity[i][0];
            sys->Velocity[ii][1] = -sys->Velocity[i][1];
            sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
            sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
            sys->mass[ii] = sys->mass[i];
            sys->density[ii] = sys->density[i];
            sys->Pressure[ii] = sys->Pressure[i];
            sys->Energy[ii] = sys->Energy[i];
            sys->hsml[ii] = sys->hsml[i];
            sys->normalizationFactor[ii] = sys->normalizationFactor[i];
            sys->DynamicPressure[ii] = sys->DynamicPressure[i];
            sys->localDensity[ii] = sys->localDensity[i];
            sys->localPressure[ii] = sys->localPressure[i];
            sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
            sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
            sys->vorticity[ii] = sys->vorticity[i];
            sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

            sys->itype[ii] = 4;
            ii ++;
          }
      }
      //  right boundary
      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0]-sys->FlatBottomLength);
        if(dx < 2.0 * sys->hsml[i]) {
          sys->Position[ii][0] = sys->FlatBottomLength + dx;
          sys->Position[ii][1] = sys->Position[i][1];
          sys->Velocity[ii][0] = -sys->Velocity[i][0];
          sys->Velocity[ii][1] = sys->Velocity[i][1];
          sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
          sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
          sys->mass[ii] = sys->mass[i];
          sys->density[ii] = sys->density[i];
          sys->Pressure[ii] = sys->Pressure[i];
          sys->Energy[ii] = sys->Energy[i];
          sys->hsml[ii] = sys->hsml[i];
          sys->normalizationFactor[ii] = sys->normalizationFactor[i];
          sys->DynamicPressure[ii] = sys->DynamicPressure[i];
          sys->localDensity[ii] = sys->localDensity[i];
          sys->localPressure[ii] = sys->localPressure[i];
          sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
          sys->localVelocity[ii][1] = sys->localVelocity[i][1];
          sys->vorticity[ii] = sys->vorticity[i];
          sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

          sys->itype[ii] = 4;
          ii ++;
        }
      }

    }
  else if(sys->virtualType==8)// dam break
    {
      //  left boundary
      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0]-0.0*sys->FlatBottomLength);
        if(dx < 2.0 * sys->hsml[i]) {
          sys->Position[ii][0] = 0.0 - dx;
          sys->Position[ii][1] = sys->Position[i][1];
          sys->Velocity[ii][0] = -sys->Velocity[i][0];
          sys->Velocity[ii][1] = sys->Velocity[i][1];
          sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
          sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
          sys->mass[ii] = sys->mass[i];
          sys->density[ii] = sys->density[i];
          sys->Pressure[ii] = sys->Pressure[i];
          sys->Energy[ii] = sys->Energy[i];
          sys->hsml[ii] = sys->hsml[i];
          sys->normalizationFactor[ii] = sys->normalizationFactor[i];
          sys->DynamicPressure[ii] = sys->DynamicPressure[i];
          sys->localDensity[ii] = sys->localDensity[i];
          sys->localPressure[ii] = sys->localPressure[i];
          sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
          sys->localVelocity[ii][1] = sys->localVelocity[i][1];
          sys->vorticity[ii] = sys->vorticity[i];
          sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

          sys->itype[ii] = 4;
          ii ++;
        }
      }

      //Bottom boundary

      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
        {
          if((sys->Position[i][1] < 2.0*sys->hsml[i]) && (sys->Position[i][0] <= sys->FlatBottomLength))
            {
              sys->Position[ii][0] = sys->Position[i][0];
              sys->Position[ii][1] = -sys->Position[i][1];
              sys->Velocity[ii][0] = sys->Velocity[i][0];
              sys->Velocity[ii][1] = -sys->Velocity[i][1];
              sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
              sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
              sys->mass[ii] = sys->mass[i];
              sys->density[ii] = sys->density[i];
              sys->Pressure[ii] = sys->Pressure[i];
              sys->Energy[ii] = sys->Energy[i];
              sys->hsml[ii] = sys->hsml[i];
              sys->normalizationFactor[ii] = sys->normalizationFactor[i];
              sys->DynamicPressure[ii] = sys->DynamicPressure[i];
              sys->localDensity[ii] = sys->localDensity[i];
              sys->localPressure[ii] = sys->localPressure[i];
              sys->localVelocity[ii][0] = sys->localVelocity[i][0];
              sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
              sys->vorticity[ii] = sys->vorticity[i];
              sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

              sys->itype[ii] = 4;
              ii ++;
            }
        }

      // edges

      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0] - sys->LeftBoundaryPosition);
        double dy = sys->Position[i][1];

        if((dx < 2.0*sys->hsml[i])&&(dy < 2.0*sys->hsml[i]))
          {
            sys->Position[ii][0] = sys->LeftBoundaryPosition-dx;
            sys->Position[ii][1] = -sys->Position[i][1];
            sys->Velocity[ii][0] = -sys->Velocity[i][0];
            sys->Velocity[ii][1] = -sys->Velocity[i][1];
            sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
            sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
            sys->mass[ii] = sys->mass[i];
            sys->density[ii] = sys->density[i];
            sys->Pressure[ii] = sys->Pressure[i];
            sys->Energy[ii] = sys->Energy[i];
            sys->hsml[ii] = sys->hsml[i];
            sys->normalizationFactor[ii] = sys->normalizationFactor[i];
            sys->DynamicPressure[ii] = sys->DynamicPressure[i];
            sys->localDensity[ii] = sys->localDensity[i];
            sys->localPressure[ii] = sys->localPressure[i];
            sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
            sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
            sys->vorticity[ii] = sys->vorticity[i];
            sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

            sys->itype[ii] = 4;
            ii ++;
          }
      }
      //  right boundary
      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0]-sys->FlatBottomLength);
        if(dx < 2.0 * sys->hsml[i]) {
          sys->Position[ii][0] = sys->FlatBottomLength + dx;
          sys->Position[ii][1] = sys->Position[i][1];
          sys->Velocity[ii][0] = -sys->Velocity[i][0];
          sys->Velocity[ii][1] = sys->Velocity[i][1];
          sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
          sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
          sys->mass[ii] = sys->mass[i];
          sys->density[ii] = sys->density[i];
          sys->Pressure[ii] = sys->Pressure[i];
          sys->Energy[ii] = sys->Energy[i];
          sys->hsml[ii] = sys->hsml[i];
          sys->normalizationFactor[ii] = sys->normalizationFactor[i];
          sys->DynamicPressure[ii] = sys->DynamicPressure[i];
          sys->localDensity[ii] = sys->localDensity[i];
          sys->localPressure[ii] = sys->localPressure[i];
          sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
          sys->localVelocity[ii][1] = sys->localVelocity[i][1];
          sys->vorticity[ii] = sys->vorticity[i];
          sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

          sys->itype[ii] = 4;
          ii ++;
        }
      }

    }

  else if(sys->virtualType==7)// hydrostatic tank
    {
      //  left boundary
      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0]-0.0*sys->FlatBottomLength);
        if(dx < 2.0 * sys->hsml[i]) {
          sys->Position[ii][0] = 0.0 - dx;
          sys->Position[ii][1] = sys->Position[i][1];
          sys->Velocity[ii][0] = -sys->Velocity[i][0];
          sys->Velocity[ii][1] = sys->Velocity[i][1];
          sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
          sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
          sys->mass[ii] = sys->mass[i];
          sys->density[ii] = sys->density[i];
          sys->Pressure[ii] = sys->Pressure[i];
          sys->Energy[ii] = sys->Energy[i];
          sys->hsml[ii] = sys->hsml[i];
          sys->normalizationFactor[ii] = sys->normalizationFactor[i];
          sys->DynamicPressure[ii] = sys->DynamicPressure[i];
          sys->localDensity[ii] = sys->localDensity[i];
          sys->localPressure[ii] = sys->localPressure[i];
          sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
          sys->localVelocity[ii][1] = sys->localVelocity[i][1];
          sys->vorticity[ii] = sys->vorticity[i];
          sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

          sys->itype[ii] = 4;
          ii ++;
        }
      }

      //Bottom boundary

      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
        {
          if((sys->Position[i][1] < 2.0*sys->hsml[i]) && (sys->Position[i][0] <= sys->FlatBottomLength))
            {
              sys->Position[ii][0] = sys->Position[i][0];
              sys->Position[ii][1] = -sys->Position[i][1];
              sys->Velocity[ii][0] = sys->Velocity[i][0];
              sys->Velocity[ii][1] = -sys->Velocity[i][1];
              sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
              sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
              sys->mass[ii] = sys->mass[i];
              sys->density[ii] = sys->density[i];
              sys->Pressure[ii] = sys->Pressure[i];
              sys->Energy[ii] = sys->Energy[i];
              sys->hsml[ii] = sys->hsml[i];
              sys->normalizationFactor[ii] = sys->normalizationFactor[i];
              sys->DynamicPressure[ii] = sys->DynamicPressure[i];
              sys->localDensity[ii] = sys->localDensity[i];
              sys->localPressure[ii] = sys->localPressure[i];
              sys->localVelocity[ii][0] = sys->localVelocity[i][0];
              sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
              sys->vorticity[ii] = sys->vorticity[i];
              sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

              sys->itype[ii] = 4;
              ii ++;
            }
        }

      // top
      for(int i = 0; i < sys->NumberOfFluidParticles; i++)
        {
          double dx = fabs(sys->Position[i][1]-sys->LocalDepth);
          if(dx < 2.0*sys->hsml[i])// && (sys->Position[i][1] >= sys->LocalDepth))
            {
              sys->Position[ii][0] = sys->Position[i][0];
              sys->Position[ii][1] = sys->LocalDepth + 2*dx;
              sys->Velocity[ii][0] = sys->Velocity[i][0];
              sys->Velocity[ii][1] = -sys->Velocity[i][1];
              sys->Velocity_xsph[ii][0] = sys->Velocity_xsph[i][0];
              sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
              sys->mass[ii] = sys->mass[i];
              sys->density[ii] = sys->density[i];
              sys->Pressure[ii] = sys->Pressure[i];
              sys->Energy[ii] = sys->Energy[i];
              sys->hsml[ii] = sys->hsml[i];
              sys->normalizationFactor[ii] = sys->normalizationFactor[i];
              sys->DynamicPressure[ii] = sys->DynamicPressure[i];
              sys->localDensity[ii] = sys->localDensity[i];
              sys->localPressure[ii] = sys->localPressure[i];
              sys->localVelocity[ii][0] = sys->localVelocity[i][0];
              sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
              sys->vorticity[ii] = sys->vorticity[i];
              sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

              sys->itype[ii] = 4;
              ii ++;
            }
        }

      // edges

      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0] - sys->LeftBoundaryPosition);
        double dy = sys->Position[i][1];

        if((dx < 2.0*sys->hsml[i])&&(dy < 2.0*sys->hsml[i]))
          {
            sys->Position[ii][0] = sys->LeftBoundaryPosition-dx;
            sys->Position[ii][1] = -sys->Position[i][1];
            sys->Velocity[ii][0] = -sys->Velocity[i][0];
            sys->Velocity[ii][1] = -sys->Velocity[i][1];
            sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
            sys->Velocity_xsph[ii][1] = -sys->Velocity_xsph[i][1];
            sys->mass[ii] = sys->mass[i];
            sys->density[ii] = sys->density[i];
            sys->Pressure[ii] = sys->Pressure[i];
            sys->Energy[ii] = sys->Energy[i];
            sys->hsml[ii] = sys->hsml[i];
            sys->normalizationFactor[ii] = sys->normalizationFactor[i];
            sys->DynamicPressure[ii] = sys->DynamicPressure[i];
            sys->localDensity[ii] = sys->localDensity[i];
            sys->localPressure[ii] = sys->localPressure[i];
            sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
            sys->localVelocity[ii][1] = -sys->localVelocity[i][1];
            sys->vorticity[ii] = sys->vorticity[i];
            sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

            sys->itype[ii] = 4;
            ii ++;
          }
      }
      //  right boundary
      for(int i = 0; i < sys->NumberOfFluidParticles; i++) {
        double dx = fabs(sys->Position[i][0]-sys->FlatBottomLength);
        if(dx < 2.0 * sys->hsml[i]) {
          sys->Position[ii][0] = sys->FlatBottomLength + dx;
          sys->Position[ii][1] = sys->Position[i][1];
          sys->Velocity[ii][0] = -sys->Velocity[i][0];
          sys->Velocity[ii][1] = sys->Velocity[i][1];
          sys->Velocity_xsph[ii][0] = -sys->Velocity_xsph[i][0];
          sys->Velocity_xsph[ii][1] = sys->Velocity_xsph[i][1];
          sys->mass[ii] = sys->mass[i];
          sys->density[ii] = sys->density[i];
          sys->Pressure[ii] = sys->Pressure[i];
          sys->Energy[ii] = sys->Energy[i];
          sys->hsml[ii] = sys->hsml[i];
          sys->normalizationFactor[ii] = sys->normalizationFactor[i];
          sys->DynamicPressure[ii] = sys->DynamicPressure[i];
          sys->localDensity[ii] = sys->localDensity[i];
          sys->localPressure[ii] = sys->localPressure[i];
          sys->localVelocity[ii][0] = -sys->localVelocity[i][0];
          sys->localVelocity[ii][1] = sys->localVelocity[i][1];
          sys->vorticity[ii] = sys->vorticity[i];
          sys->thermalDiffusivity[ii] = sys->thermalDiffusivity[i];

          sys->itype[ii] = 4;
          ii ++;
        }
      }

    }
  else
    {
      int a = 0;
    }
  /* Initialize derivatives: must be after inoutflow.h */
  sys->NumberOfVirtualParticles = ii - sys->NumberOfFluidParticles - sys->NumberOfSolidParticles - sys->NumberOfBoundaryParticles;


  /* Interaction plus kernel definition*/
  CalculateKernel(sys);
  if(sys->flagInitSimul==1)
    {
      /*********************************************/
      /* algorithm for predicitng the rates of change */

      PredictingDensityEquation(sys);
      PredictingDynamicPressureEquation(sys);
      PredictingPressureGradientForce(sys);
      /*************End some predicitng ****************/

      /***************************************/
      /* calculate locally averaged quantities first*/

      LocallyAveragedCompressibility(sys);

      LocallyAveragedSpecificVolume(sys);

      LocalDensityFluctuation(sys);

      LocalPressureFluctuation(sys);

      LocalVelocityFluctuation(sys);
      /**************************************/

      DensityEquation(sys);

      EquationOfState(sys);

      ClassicShearStress(sys);

      DynamicPressureEquation(sys);

      PressureGradientForce(sys);

      if(sys->stressType==1)
        {
          ViscousForcesSingleStep(sys);
          PredictingViscousForcesSingleStep(sys);
        }
      else if(sys->stressType==2)
        {
          ViscousForcesTwoStep(sys);
          PredictingViscousForcesSingleStep(sys);
        }
      else
        {
          printf("No viscous force.");
        }


      BodyForce(sys);

      PredictingBodyForce(sys);

      VorticityCompute(sys);

      ApproximateVortcity(sys);


      SpecificPowerDissipation(sys);

      PowerDissipation(sys);

      SpongeLayerWaveReflectionAbsorption(sys);

      //       MovingBoundaryExcitation(sys);

      /* compute dissipation rate in boundary layer*/
      SerialComputeSpecificBoundaryDissipationRate(sys);
      ComputeBoundaryDissipationRate(sys);

      /* compute dissipation rate in fluid bulk*/
      ComputeSpecificDissipationRate(sys);
      ComputeDissipationRate(sys);

      /* compute specific turbulent dissipation rate */
      ComputeSpecificTurbulentDissipationRate(sys);
      ComputeTurbulentKineticEnergy(sys);

      ControlVolumeProperties(sys);// computes various thermodynamic quantities for multiple control volumes

    }
  else
    {
      CalculateNormalizationFactor(sys);

      CalculateGradientOfNormalizationFactor(sys);

      InitializationMomentumEquation(sys);

      InitialDensityAndPressure(sys);
    }


  if(sys->selectBoundaryForce==1)
    {
      LeonnardJonesBoundaryForces(sys);
      PredictingSpringBoundaryForces(sys); // replace with predicitng Lennard Jones
    }
  else if(sys->selectBoundaryForce==3)
    {
      SpringBoundaryForces(sys);
      PredictingSpringBoundaryForces(sys);
    }
  else if(sys->selectBoundaryForce==4)
    {
      ImprovedMonaghan(sys);
      PredictingSpringBoundaryForces(sys); // replace with predicting MonaghanBoundaryForce
    }
  else if(sys->selectBoundaryForce==5)
    {
      SerialSpringBoundaryForces(sys);
      PredictingSpringBoundaryForces(sys); //
    }
  else
    {
      MonaghanBoundaryForce(sys);
      PredictingSpringBoundaryForces(sys); // replace with predicitng MonaghanBoundaryForce
    }

  TotalKineticEnergy(sys);
  NeighborsToPressurePoint(sys);
  NeighborsToPressurePointTwo(sys);
  PressureKernelNormalizationFactor(sys);
  //   EnsembleAveragedPointPressure(sys);
  EnsembleAveragedPointPressureTwo(sys);
  //   TerminateInitialization(sys, t);
  DensityConsistencyFunctions(sys);

  EnstrophyPowerCalculate(sys);
}
