/* upbqtl.c  */
/* 28 Oct 2000 */
/* C.C.Berry */


/* upbqtl.c is called by Splus with a data vector and objects needed
   to set up for laplace adjustment of the posterior under normal
   link. It also contains a matric of alternative loci to which models
   are to be fit. For now just the loglik. Later: It returns the
   laplace adjusted posterior, the adjustment factor (ratio of exact
   and laplace for linearized likelihood), the loglikelihood at the
   maximum, and the parameter values that achieve the maximum. */

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

void upbqtl(longint *crsType, longint *nparm,
            double *y, 
            double *x,
            double *z,
	    double *amnt, double *tab, longint *needProd, double *lambda, double
	    *coefs, double *casewt, double *sigma, double *hkExact, double *hkApprox,
	    double *llk, double *hess, double *postApprox, longint *iter, double *tol, 
	    longint *nem,
	    /* the above args are for the usual call to lapadj */
	    
	    longint *loc_mat, longint *n_alt, longint *loc_order, double *result, 
	    double *orig_x, longint *loc_right, double *map_lambda,
	    double *state_matrix, longint *n_state_loc) 
     
     
{
  
  /* package the longint parms in *nparm
     elements are as follows:
  */
  longint *n,  *nrx,  *ncx,  *nloc, *ncz,  *nx2uz,  *nz2uz,  *nt2uz,
    *nreg, *np, *ptx,  *ptz,  *pttx, *pttz;
  
  
  longint N, NRX, NCX, NP, NPNP, NLOC, NX2UZ, NZ2UZ, NT2UZ, NREG,
    itStill[1], rank[1],  reset[1], dq[2],
    *cur_loci, *perm_indx, radix, *radixProd,*newRadixProd ;
  int i,j,k, i_alt;
  double newllk[1], oldllk[1], lprior[1], hkSigma,
    tmp, tmpTab[9], sigma2[1], df[1], ss[1], ss0[1], detvar, qrtol=QR_TOL;
  
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
 
  N = *n;
  NRX = *nrx;
  NCX = *ncx;
  NX2UZ = *nx2uz;
  
  NZ2UZ = *nz2uz;
  NT2UZ = *nt2uz;
  NLOC = *nloc;
  NREG = *nreg;
  NP =   *np;
  NPNP = NP*NP;
   
  ALLOC_LONG(pvt,NP);
  ALLOC_DBL(xvx,NREG*NREG);
  ALLOC_DBL(xvy,NREG);
  ALLOC_DBL(fits,NRX);
  ALLOC_DBL(res,N*NRX);
  ALLOC_DBL(vsum,N);
  ALLOC_DBL(vsum2,NRX);
  ALLOC_DBL(vysum,NRX);
  ALLOC_DBL(qrxvx,NP*NP);
  ALLOC_DBL(qraux,NP);
  ALLOC_DBL(wrksp,2*NP);
  ALLOC_DBL(dgr,NP);
  ALLOC_DBL(dparm,NP);
  ALLOC_DBL(newgr,NP);
  ALLOC_DBL(newprm,NP);
  ALLOC_DBL(bk,NP*NP);
  ALLOC_DBL(bks,NP);
  ALLOC_DBL(qrbk,NP*NP);
  ALLOC_DBL(paradj,NP);
  ALLOC_DBL(emprm,NP);
  ALLOC_DBL(newhss,NP*NP);
  ALLOC_DBL(curprm,NP);
  ALLOC_DBL(oldprm,NP);
  ALLOC_DBL(oldgr,NP);
  ALLOC_DBL(ez,NREG);
  ALLOC_DBL(hkCoefs,NREG);
  ALLOC_DBL(wt,NRX*N);
  ALLOC_DBL(odds,N);
  ALLOC_DBLPT(zc,(NZ2UZ ? NZ2UZ : 1));
  ALLOC_DBLPT(xc,(NX2UZ ? NX2UZ : 1));
  ALLOC_DBLPT(tzc,(NT2UZ ? NT2UZ : 1));
  ALLOC_DBLPT(txc,(NT2UZ ? NT2UZ : 1));
  ALLOC_LONG(cur_loci, NLOC);
  ALLOC_LONG(radixProd, NLOC);
  ALLOC_LONG(newRadixProd, NLOC);
  ALLOC_LONG(perm_indx,NRX);      
/*  init ptrs to x and z columns */
  
  for (i=0;i<NX2UZ;i++) xc[i] = x + ptx[i]*NRX;
  for (i=0;i<NZ2UZ;i++) zc[i] = z + ptz[i]*N;
  for (i=0;i<NT2UZ;i++){
    txc[i] = x + pttx[i]*NRX;
    tzc[i] = z + pttz[i]*N;
  }


  radix = *crsType+1;
  *radixProd=1;
  for (i=1;i<NLOC;i++) radixProd[i] = radixProd[i-1]*radix;
  /* loop thru the loci */
  for (i_alt=0; i_alt<*n_alt; i_alt++){ 

    *itStill =  *iter; /*reset iteration count */

    for (i=0; i<NLOC; i++) cur_loci[i] =
			     loc_mat[loc_order[i+NLOC*i_alt]+NLOC*i_alt];
    if (NLOC>1){
      /* set up needProd */
      for (i=0; i<N; i++)
	for (j=0;j< NLOC-1; j++)
	  needProd[i+N*j] = 
	    (longint) (loc_right[i+N*cur_loci[j]] > cur_loci[j+1]);
      /* calc lambda */
      for (j=0;j< (NLOC-1); j++){
	tmp = 1.0;
	for (k=cur_loci[j]; k<cur_loci[j+1]; k++) 
	  tmp *= map_lambda[k];
	lambda[j]=tmp; 
      }
      /* get perm vector for x */
      
      for (i=0;i<NLOC;i++) {
	newRadixProd[loc_order[i+NLOC*i_alt]]=radixProd[i];
      }
      for (i=0;i<NRX;i++){
	k=0;
	for (j=0;j<NLOC;j++) k+= ((i/newRadixProd[j]) % radix) * radixProd[j];
	perm_indx[i]=k;
      }
      /* revise x */
      for (i=0;i<NCX;i++)
	for (j=0;j<NRX;j++){
	  xc[i][j] = orig_x[ i*NRX + perm_indx[j] ];
	}
    }
    /* revise tab */
    for (k=0;k<radix;k++)
      for (i=0;i<NLOC;i++)
	for (j=0;j<N;j++)
	  tab[j + N*i + N*NLOC*k]=
	    state_matrix[j + N*cur_loci[i] + N**n_state_loc*k];

    if (NRX == 1) { /* NO WEIGHTS */
      for (i=0;i<N;i++) wt[i] = 1.0;
    }
    else {
    if (*crsType==1) {
      if (NLOC==1) {
	/* one column - just copy tab to wt */
	for (i=0;i<N;i++)
	  {
	    wt[i*2] = tab[i];
	    wt[1+i*2] = tab[i+N];
	  }
      }
      else  
	bc1wt(tab,needProd,lambda,wt,&NLOC,&N,tmpTab,odds);
    }
    else
      f2wt(tab,needProd,lambda,wt,&NLOC,&N,tmpTab,tol,iter);
  }
   
  hkreg(y,xc,zc,txc,tzc,ez,wt,amnt,
	sigma2,newllk,hkExact,hess,
	hkCoefs,nparm,
	xvx,xvy, pvt,qraux,wrksp,
	rank, hkApprox,ss0,casewt);
    
  if (*rank<NREG) 
    PROBLEM "deficient rank in hkreg\n" WARNING(NULL_ENTRY) ; 
  
  hkSigma = sqrt(*sigma2);

  /* note llkEm() is written to implicitly include constant
     vector. Therefore, x is just NRX by NREG - 1
     */
 


  for (i=0;i<NREG;i++) coefs[i] = hkCoefs[i];
  
llkEm(nparm, xc, zc, txc, tzc, &hkSigma, amnt, fits, y, coefs, wt, 
	 res, vsum, vsum2, vysum, lprior, oldllk, xvx, xvy,
	 df, curprm, ss, oldgr, hess, &curprm[NP-1], qrxvx, qraux,
	 pvt, wrksp, &nem1,casewt);

  
  *oldllk += *lprior;
  *llk = *oldllk;

    for (i=0;i<NREG;i++) 
            emprm[i] = oldprm[i] = hkCoefs[i];
    /*oldprm[i] = curprm[i];*/
    *sigma = sqrt(curprm[NP-1]);
      emprm[NP-1] = oldprm[NP-1] = *sigma2;
      /*oldprm[NP-1] = curprm[NP-1];*/
    *reset = 1;  
  
    lapWhl(y, wt, amnt, nparm, xc, zc, txc, tzc, fits, res, vsum, 
	    vsum2, vysum, xvx, xvy, df, qrxvx, qraux, pvt,
	    wrksp, tol, llk, dgr, dparm, newgr, oldgr, curprm,
	    oldprm, reset, bk, newhss, bks, qrbk, paradj,
	    newprm, emprm, itStill,nem, casewt);

    for (i=0;i<NREG;i++)
      coefs[i] = newprm[i];

  /* get Laplace Approx */

  CPY(newhss,qrxvx,&NPNP);
  rank[0] = NP;
  dq[0] = NP;
  dq[1] = NP;
#ifndef USING_R
  F77_CALL(dqr)(qrxvx, dq, pvt, qraux, &qrtol, wrksp, rank);
#else
  F77_CALL(dqrdc2)(qrxvx, dq, dq, dq, &qrtol, rank, qraux, pvt, wrksp  );
#endif
  detvar = 1.0;
  for (j=0;j<NP;j++) {
    detvar  *= sqrt(fabs(qrxvx[j*NP+j]));
  } 
  postApprox[0]=exp(llk[0] - ss0[0])/detvar*pow(2.0*M_PI,(double) (NP)/2.0);
  /*  *sigma=log(*sigma);*/

  *sigma = log(curprm[NP-1])/2.0;


  /* save the logpost */
  result[i_alt] = *llk;
  }
  *iter -= *itStill;   /* iterations for last run */
  /* free allocated memory */
    
  DEALLOC_LONG(pvt,NP);
  DEALLOC_DBL(xvx,NREG*NREG);
  DEALLOC_DBL(xvy,NREG);
  DEALLOC_DBL(fits,NRX);
  DEALLOC_DBL(res,N*NRX);
  DEALLOC_DBL(vsum,N);
  DEALLOC_DBL(vsum2,NRX);
  DEALLOC_DBL(vysum,NRX);
  DEALLOC_DBL(qrxvx,NP*NP);
  DEALLOC_DBL(qraux,NP);
  DEALLOC_DBL(wrksp,2*NP);
  DEALLOC_DBL(dgr,NP);
  DEALLOC_DBL(dparm,NP);
  DEALLOC_DBL(newgr,NP);
  DEALLOC_DBL(newprm,NP);
  DEALLOC_DBL(bk,NP*NP);
  DEALLOC_DBL(bks,NP);
  DEALLOC_DBL(qrbk,NP*NP);
  DEALLOC_DBL(paradj,NP);
  DEALLOC_DBL(emprm,NP);
  DEALLOC_DBL(newhss,NP*NP);
  DEALLOC_DBL(curprm,NP);
  DEALLOC_DBL(oldprm,NP);
  DEALLOC_DBL(ez,NREG);
  DEALLOC_DBL(hkCoefs,NREG);
  DEALLOC_DBL(wt,NRX*N);
  DEALLOC_DBL(odds,N);
  DEALLOC_DBLPT(zc,NZ2UZ);
  DEALLOC_DBLPT(xc,NX2UZ);
  DEALLOC_DBLPT(tzc,NT2UZ);
  DEALLOC_DBLPT(txc,NT2UZ);
  DEALLOC_LONG(cur_loci, NLOC);
  DEALLOC_LONG(radixProd, NLOC);
  DEALLOC_LONG(newRadixProd, NLOC);
  DEALLOC_LONG(perm_indx,NRX);      

}
 
