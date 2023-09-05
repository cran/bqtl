/* hkreg.c */
/* Charles C. Berry */
/* Wed Jan 27 14:43:14 PST 1999 */
/* Do Haley-Knott regression and pseudo-posterior */

#include "lapadj.h"
#ifdef USING_R
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>
#include <R_ext/BLAS.h>
#include <Rmath.h>
#define lgamma(x) lgammafn(x)
#else
int F77_SYMBOL(dqrsl1)();
int F77_SYMBOL(dqr)();
#endif

#define UPRITE(x,y,z) ( (x) < (y) ) ? (x) + (y) * (z) : (y) + (x) * (z)

void hkreg(double *y,double **xc, double **zc, double **txc, double **tzc,
           double *ez,double *wt,double *amnt,
	   double *sigma2,double *llik,double *post,double *hess,
	   double *coefs,longint *nparm,
	   double *xvx,double *xvy, longint *pvt,double *qraux,double *wrksp,
	   longint *rank, double *hkapprox, double *ss0, double *casewt)
{
  int i, j, k, n, nrx, nx2uz, nz2uz, nt2uz, nreg, df;
  double tmp, vsum, ss, sigma, detvar,logpost;
  const double rt2pi = 2.5066282746310002416;
  double tol[1] = {QR_TOL};
  longint dq[2], one, hndrd, df2;


  n = nparm[0];
  nrx = nparm[1];
  nx2uz = nparm[5];
  nz2uz = nparm[6];
  nt2uz = nparm[7]; 
  nreg = nparm[8];


  for (i=0;i<nreg*nreg;i++) 
    xvx[i] = 0.0;
  for (i=0;i<nreg;i++) 
    xvy[i] = 0.0;
  for (i=0;i<nreg;i++) 
    xvx[i+i*nreg] = amnt[i];
  
  for (i=0;i<n;i++) {
    ez[0]=1.0;
    if (*casewt<0.0) {
      xvx[0]+=1.0;
      xvy[0] += y[i];
    }
    else {
      xvx[0]+=casewt[i];
      xvy[0] += y[i]*casewt[i];
    }
    
    if (nrx>1) {
      vsum = 0.0;
      for (k=0;k<nrx;k++) 
	vsum += wt[k+i*nrx];
      for (j=0;j<nx2uz;j++) {
	ez[j+1]=0.0;
	for (k=0;k<nrx;k++) 
	  ez[j+1] += wt[i*nrx+k]*xc[j][k]; /*expected values */
	ez[j+1] /= vsum;
      } /* close j*/
    }
    else {
      vsum =1.0;
    }
     for (j=0;j<nz2uz;j++) ez[j+nx2uz+1] = zc[j][i];     /*fixed covariates */
     for (j=0;j<nt2uz;j++) { /*interaction of cov and locus */
       tmp = 0.0;
       for (k=0;k<nrx;k++) 
         tmp += wt[i*nrx+k]*txc[j][k];
       ez[j+nx2uz+nz2uz+1] = tmp  * tzc[j][i] / vsum;
     } /*close j*/
     for (j=1;j<nreg;j++){     
     if (*casewt<0.0) {
	xvx[j] += ez[j];
	xvy[j] += ez[j]*y[i];
	for (k=1;k<=j;k++) {
	  xvx[j+k*nreg] += ez[j]*ez[k];
	}
      }
      else {
	tmp = casewt[i];
	xvx[j] += tmp*ez[j];
	xvy[j] += tmp*ez[j]*y[i];
	for (k=1;k<=j;k++) {
	  xvx[j+k*nreg] += tmp*ez[j]*ez[k];
	}
      }	
    }
  }
  for (j=0;j<nreg;j++)
    for (k=j+1;k<nreg;k++)
      xvx[j+k*nreg] = xvx[k+j*nreg];

  for (j=0;j<nreg;j++) {  
    for (k=0;k<nreg;k++) 
    hess[k+(nreg+1)*j] = ( hess[j+(nreg+1)*k] =  - xvx[k+j*nreg] );
  }

  /* hess[j+(nreg+1)*k] = - xvx[k+j*nreg] ;
     will fill the sub matrix of hess[0:(nreg-1),0:(nreg-1) with -xvx
     Apr 8, 2012 this is harmless but can be cleaned up
  */
  *rank = nreg;
  dq[0] =nreg;
  dq[1] =nreg;
  for (i=0;i<nreg;i++) pvt[i]= (i+1);
#ifndef USING_R
  F77_CALL(dqr)(xvx, dq, pvt, qraux, tol, wrksp, rank);
#else
  F77_CALL(dqrdc2)(xvx,dq,dq,dq,tol,rank,qraux,pvt,wrksp);
#endif
  
  /*** AFTER HERE xvx IS THE QR DECOMP OF X'X MATRIX ***/

  dq[0] =nreg;
  dq[1] =nreg;
  one = 1;
  hndrd = 100;

#ifndef USING_R
  F77_CALL(dqrsl1)(xvx, dq, qraux, rank, xvy, &one, wrksp, coefs, &hndrd, &one);
#else
  F77_CALL(dqrsl)(xvx,dq,dq,rank,qraux,xvy,wrksp,coefs,coefs,
		    wrksp,wrksp,&hndrd, &one);
#endif

  if (*rank<nreg) { /*re-order coefs*/
	for (i=0;i<*rank;i++) wrksp[pvt[i]-1] = coefs[i];
	for (i=*rank;i<nreg;i++) wrksp[pvt[i]-1] = 0.0;
	CPY(wrksp,coefs,&nreg);
      }

  df = n;
  for (i=0;i<nreg;i++) df += (amnt[i]>0) ? 1 : 0 ;
  df2 = df - nreg;

  ss = 0.0;
  if (*casewt<0.0) {
    for (i=0;i<n;i++) 
      ss += y[i]*y[i];
  }
  else {
    for (i=0;i<n;i++) 
      ss += casewt[i]*y[i]*y[i];
  }

  ss0[0] = ss - xvy[0]*xvy[0]/xvx[0];
  for (j=0;j<nreg;j++) 
    ss -= coefs[j]*xvy[j];

  sigma = sqrt(ss/df);
  *sigma2 = ss/df;
  llik[0] = - df * log(sigma) - ((double) df/2.0)*(1.0+log(2.0*M_PI));
  for (i=0;i<nreg;i++) 
    if (amnt[i]>0)
      llik[0] += log(amnt[i])/2.0;
  if (*casewt > 0.0) {
    tmp = 1.0;
    for (i=0;i<n;i++)
      tmp *= casewt[i];
    llik[0] += log(tmp)/2.0;
  }

  logpost = ((double) -df2) * log(rt2pi);
  if (*casewt > 0.0)
    logpost += log(tmp)/2.0;

  detvar = 1.0;
  for (j=0;j<nreg;j++) {
    if (amnt[j] >0) logpost += log(amnt[j])/2.0;/*post[0] *= sqrt(amnt[j]);*/
    detvar  *= sqrt(fabs(xvx[j*nreg+j]));
  }
  
  logpost += ((double) -df2)/2.0*log(ss/2.0) + lgamma(((double) df2)/2.0);
  /*normalize posteriors to prevent machine under/overflow */
  ss0[0] = ((double) -df2)/2.0*log(ss0[0]/2.0) + lgamma(((double) df2)/2.0);
  post[0] = exp(logpost-ss0[0])/2.0/detvar;
  detvar  *= sqrt(2.0 * ss)/pow(sigma,(double) nreg+1);
  hkapprox[0] = exp(llik[0]-ss0[0])/detvar*pow(2.0*M_PI,(double) (nreg+1)/2.0);


  for (i=0;i<nreg;i++)
    for (j=0;j<nreg;j++) {
      hess[i+j*(nreg+1)]/=(sigma*sigma);
    }

  for (i=0;i<nreg;i++) 
    hess[i+nreg*(nreg+1)] = (hess[(nreg)*(i+1)-1] = 0.0);
  hess[(nreg+1)*(nreg+1)-1] = - (2 * ss)/(*sigma2);

}
