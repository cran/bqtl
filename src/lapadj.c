/* lapadj.c  */
/* Fri Jan 29 11:38:08 PST 1999 */
/* C.C.Berry */


/* lapadj.c is called by Splus/R with a data vector and objects needed
   to set up for laplace adjustment of the posterior under normal
   link. It returns the laplace adjusted posterior, the adjustment
   factor (ratio of exact and laplace for linearized likelihood), the
   loglikelihood at the maximum, and the parameter values that achieve
   the maximum. */

#include "lapadj.h"
#ifdef USING_R
#include <R_ext/Applic.h>
#endif

void bc1wt(double *tab, longint *needProd, double *lambda,
	   double *rtab, longint *nc1, longint *n, 
	   double tmpTab[4], double *odds);

void f2wt(double *tab, longint *needProd, double *lambda, 
	  double *rtab, longint *nc1, longint *n,double *tmpTab, 
	  double *maxtol,longint maxIter[1] );


void hkreg(double *y,double **xc, double **zc, double **txc, double **tzc,
           double *ez,double *wt,double *amnt,
	   double *sigma2,double *llik,double *post,double *hess,
	   double *coefs,longint *nparm,
	   double *xvx,double *xvy, longint *pvt,double *qraux,double *wrksp,
	   longint *rank, double *hkapprox, double *ss0, double *casewt);

void llkEm();
void lapWhl();

void lapadj(longint *crsType, longint *nparm,
            double *y, 
            double *x,
            double *z,
	    double *amnt, double *tab, longint *needProd, double *lambda, double
	    *coefs, double *casewt, double *sigma, double *hkExact, double *hkApprox,
	    double *llk, double *hess, double *postApprox, longint *iter, double *tol, 
	    longint *nem)

{

/* package the longint parms in *nparm
   elements are as follows:
*/
  longint *n,  *nrx,  *ncx,  *nloc, *ncz,  *nx2uz,  *nz2uz,  *nt2uz,
        *nreg, *np, *ptx,  *ptz,  *pttx, *pttz;


  longint npnp[1], itStill[1], rank[1],  reset[1], dq[2] ;
  int i,j;
  double newllk[1], oldllk[1], lprior[1], hkSigma,
    tmpTab[9], sigma2[1], df[1], ss[1], ss0[1],detvar, qrtol=QR_TOL;

  /* dynamic storage decl's follow */

  double *xvx, *xvy, *fits, *res, *vsum, *vsum2, *vysum, *qrxvx,
    *qraux, *wrksp, *dgr, *dparm, *newgr, *oldgr, *newprm, *bk, *bks,
    *qrbk, *paradj, *emprm, *newhss, *curprm, *oldprm, *ez, *wt, *odds,
    *hkCoefs, **zc, **xc, **txc, **tzc;

  longint *pvt, nem1=1;
 
  n=nparm;  
  nrx=n+1;  
  ncx=nrx+1;  
  nloc = ncx+1;  
  ncz = nloc+1;  
  nx2uz = ncz+1;  
  nz2uz = nx2uz + 1;  
  nt2uz = nz2uz + 1;  
  nreg = nt2uz + 1;  
  np = nreg + 1;  
  ptx = np + 1;  
  ptz = ptx +  ( *nx2uz  ? *nx2uz: 1 ) ;  
  pttx =ptz +  ( *nz2uz  ? *nz2uz: 1 );  
  pttz =pttx + ( *nt2uz  ? *nt2uz: 1 ); 
  
  *npnp = *np**np;  
  *itStill =  *iter;  
    
  ALLOC_LONG(pvt,*np);
  ALLOC_DBL(xvx,*nreg**nreg);
  ALLOC_DBL(xvy,*nreg);
  ALLOC_DBL(fits,*nrx);
  ALLOC_DBL(res,*n**nrx);
  ALLOC_DBL(vsum,*n);
  ALLOC_DBL(vsum2,*nrx);
  ALLOC_DBL(vysum,*nrx);
  ALLOC_DBL(qrxvx,*npnp);
  ALLOC_DBL(qraux,*np);
  ALLOC_DBL(wrksp,2**np);
  ALLOC_DBL(dgr,*np);
  ALLOC_DBL(dparm,*np);
  ALLOC_DBL(newgr,*np);
  ALLOC_DBL(newprm,*np);
  ALLOC_DBL(bk,*npnp);
  ALLOC_DBL(bks,*np);
  ALLOC_DBL(qrbk,*npnp);
  ALLOC_DBL(paradj,*np);
  ALLOC_DBL(emprm,*np);
  ALLOC_DBL(newhss,*npnp);
  ALLOC_DBL(curprm,*np);
  ALLOC_DBL(oldprm,*np);
  ALLOC_DBL(oldgr,*np);
  ALLOC_DBL(ez,*nreg);
  ALLOC_DBL(hkCoefs,*nreg);
  ALLOC_DBL(wt,*nrx**n);
  ALLOC_DBL(odds,*n);
  ALLOC_DBLPT(zc,(*nz2uz ? *nz2uz : 1));
  ALLOC_DBLPT(xc,(*nx2uz ? *nx2uz : 1));
  ALLOC_DBLPT(tzc,(*nt2uz ? *nt2uz : 1));
  ALLOC_DBLPT(txc,(*nt2uz ? *nt2uz : 1));
  
/*  init ptrs to x and z columns */
  
  for (i=0;i<*nx2uz;i++) xc[i] = x + ptx[i]**nrx;
  for (i=0;i<*nz2uz;i++) zc[i] = z + ptz[i]**n;
  for (i=0;i<*nt2uz;i++){
    txc[i] = x + pttx[i]**nrx;
    tzc[i] = z + pttz[i]**n;
  }
    
  if (*nrx == 1) { /* NO WEIGHTS */
    for (i=0;i<*n;i++) wt[i] = 1.0;
  }
  else {
    if (*crsType==1) {
      if (*nloc==1) {
	/* one column - just copy tab to wt */
	for (i=0;i<*n;i++)
	  {
	    wt[i*2] = tab[i];
	    wt[1+i*2] = tab[i+*n];
	  }
      }
      else  
	bc1wt(tab,needProd,lambda,wt,nloc,n,tmpTab,odds);
    }
    else
      f2wt(tab,needProd,lambda,wt,nloc,n,tmpTab,tol,iter);
  }
  
  hkreg(y,xc,zc,txc,tzc,ez,wt,amnt,
	sigma2,newllk,hkExact,hess,
	hkCoefs,nparm,
	xvx,xvy, pvt,qraux,wrksp,
	rank, hkApprox, ss0,casewt);
 
  
  if (*rank<*nreg) 
    PROBLEM "deficient rank in hkreg\n" WARNING(NULL_ENTRY) ; 
  
  hkSigma = sqrt(*sigma2);

  /* note llkEm() is written to implicitly include constant
     vector. Therefore, x is just NRX by NREG - 1
     */
 


  for (i=0;i<*nreg;i++) coefs[i] = hkCoefs[i];
    
llkEm(nparm, xc, zc, txc, tzc, &hkSigma, amnt, fits, y, coefs, wt, 
	 res, vsum, vsum2, vysum, lprior, oldllk, xvx, xvy,
	 df, curprm, ss, oldgr, hess, &curprm[*np-1], qrxvx, qraux,
	 pvt, wrksp, &nem1,casewt);

    
  *oldllk += *lprior;
  *llk = *oldllk;

    for (i=0;i<*nreg;i++) 
            emprm[i] = oldprm[i] = hkCoefs[i];
    /*oldprm[i] = curprm[i];*/
    *sigma = sqrt(curprm[*np-1]);
      emprm[*np-1] = oldprm[*np-1] = *sigma2;
      /*oldprm[NP-1] = curprm[NP-1];*/
    *reset = 1;  
    
    lapWhl(y, wt, amnt, nparm, xc, zc, txc, tzc, fits, res, vsum, 
	    vsum2, vysum, xvx, xvy, df, qrxvx, qraux, pvt,
	    wrksp, tol, llk, dgr, dparm, newgr, oldgr, curprm,
	    oldprm, reset, bk, newhss, bks, qrbk, paradj,
	    newprm, emprm, itStill,nem, casewt);
    
    for (i=0;i<*nreg;i++)
      coefs[i] = newprm[i];

  /* get Laplace Approx */

  CPY(newhss,qrxvx,npnp);
  rank[0] = *np;
  dq[0] = *np;
  dq[1] = *np;
#ifndef USING_R
  F77_CALL(dqr)(qrxvx, dq, pvt, qraux, &qrtol, wrksp, rank);
#else
  F77_CALL(dqrdc2)(qrxvx, dq, dq, dq, &qrtol, rank, qraux, pvt, wrksp  );
#endif
  detvar = 1.0;
  for (j=0;j<*np;j++) {
    detvar  *= sqrt(fabs(qrxvx[j**np+j]));
  } 
  postApprox[0]=exp( llk[0] -ss0[0] - log(detvar) + 
		     (double)(*np)/2.0*log(2.0*M_PI));
  /*  *sigma=log(*sigma);*/

  *sigma = log(curprm[*np-1])/2.0;
  *iter -= *itStill;
  
  /* free allocated memory */
    
  DEALLOC_LONG(pvt,*np);
  DEALLOC_DBL(xvx,*nreg**nreg);
  DEALLOC_DBL(xvy,*nreg);
  DEALLOC_DBL(fits,*nrx);
  DEALLOC_DBL(res,*n**nrx);
  DEALLOC_DBL(vsum,*n);
  DEALLOC_DBL(vsum2,*nrx);
  DEALLOC_DBL(vysum,*nrx);
  DEALLOC_DBL(qrxvx,*npnp);
  DEALLOC_DBL(qraux,*np);
  DEALLOC_DBL(wrksp,2**np);
  DEALLOC_DBL(dgr,*np);
  DEALLOC_DBL(dparm,*np);
  DEALLOC_DBL(newgr,*np);
  DEALLOC_DBL(newprm,*np);
  DEALLOC_DBL(bk,*npnp);
  DEALLOC_DBL(bks,*np);
  DEALLOC_DBL(qrbk,*npnp);
  DEALLOC_DBL(paradj,*np);
  DEALLOC_DBL(emprm,*np);
  DEALLOC_DBL(newhss,*npnp);
  DEALLOC_DBL(curprm,*np);
  DEALLOC_DBL(oldprm,*np);  
  DEALLOC_DBL(oldgr,*np);
  DEALLOC_DBL(ez,*nreg);
  DEALLOC_DBL(hkCoefs,*nreg);
  DEALLOC_DBL(wt,*nrx**n);
  DEALLOC_DBL(odds,*n);
  DEALLOC_DBLPT(zc,*nz2uz);
  DEALLOC_DBLPT(xc,*nx2uz);
  DEALLOC_DBLPT(tzc,*nt2uz);
  DEALLOC_DBLPT(txc,*nt2uz);
     
}
 
