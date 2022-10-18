/* normLogLik.c */
/*Wed Feb  3 13:38:23 PST 1999  */
/* C.C.Berry*/
/* get likelihood, first partials, and second partials*/

#include "lapadj.h"

void normLogLik(
		
		 longint *nparm,
		 double **xc, double **zc, double **txc, double **tzc,
		 double *y, double *fits, double * v, double * sigma,
		 double *amnt, double *beta,
		 double *tmp1, double *tmp2, double *tmp3,
		 double *dl_d1, double *dl_d2, double *llik, double *casewt)


{
  longint i, j, k, k_1, NP, NREG;
  longint *n,  *nrx,  *ncx,  *nloc, *ncz,  *nx2uz,  *nz2uz,  *nt2uz,
        *nreg, *np;


  double yi, r, v_sum, v_adj, sigma_sq, tmp_sig, tmp_sig_2, nrmcns, d__1,
    lgrtN, wt_sig2, w, vrw, vr2w, vr2w2, vr3w2, vr4w2;

  const double rt2pi = 2.5066282746310002416;

  n=nparm;
  nrx=n+1;
  ncx=nrx+1;
  nloc = ncx+1;
  ncz = nloc+1;
  nx2uz = ncz+1;
  nz2uz = nx2uz +1;
  nt2uz = nz2uz +1;
  nreg = nt2uz + 1;
  np = nreg + 1;
     

  nrmcns = log(rt2pi * *sigma);
  lgrtN = 0.0;
  if (*casewt>0.0) {
    lgrtN = 1.0;
    for (i=0;i<*n;i++) 
      lgrtN *= casewt[i];
    lgrtN = log(lgrtN) / 2.0;
  }
  NREG=*nreg;
  NP = *np;
  sigma_sq = *sigma**sigma;

  llik[0] =  lgrtN;
  for (j=0;j<NP;j++) dl_d1[j]=0.0;
  for (j=0;j<NP*NP;j++) dl_d2[j]=0.0;


  for (i=0;i<*n;i++) {
    yi = y[i];
    for (j=0;j<*nz2uz;j++) 
        yi -= zc[j][i]*beta[1+*nx2uz+j];      /* subtract fixed covars */
    wt_sig2 = (*casewt > 0.0) ? casewt[i]/sigma_sq : 1.0/sigma_sq;
    v_sum=0.0;
    tmp_sig = tmp_sig_2 = 0.0;
    for (j=0;j<*nreg;j++) {
      tmp1[j] = 0.0;
 
    }
    
    for (j=0;j<*nrx;j++) 
      v_sum += v[i**nrx+j];
  
    llik[0] += log(v_sum);

    for (k_1=0;k_1<*nrx;k_1++) {
      v_adj = v[i*(*nrx)+k_1]/v_sum;
      r = yi - fits[k_1];
      for (k=0;k<*nt2uz;k++) 
        r -= beta[1+*nx2uz+*nz2uz+k]*txc[k][i]*tzc[k][i];

      w = v_adj*wt_sig2;
      vrw = w*r;
      vr2w = vrw*r;
      vr2w2 = vr2w*wt_sig2;
      vr3w2 = vr2w2*r;
      vr4w2 = vr3w2*r;
      tmp_sig += vr2w;
      tmp_sig_2 += vr4w2-2*vr2w;
      
      tmp1[0] += vrw;                                     /* 1st deriv */
      dl_d2[0] += (vr2w2 - w);                              /*const^2 */
      dl_d2[NP-1] += vr3w2 - 2*vrw;                       /* sigma'const */
      for(j=0;j<*nx2uz;j++) { 
        tmp1[j+1] += xc[j][k_1]*vrw;
        dl_d2[j+1] += xc[j][k_1] * (vr2w2 - w);        /* x'1 */
        dl_d2[(NP-1)+(j+1)*NP] += 
            xc[j][k_1] * (vr3w2 - 2*vrw);                 /* sigma'x */ 
      for (k=j;k<*nx2uz;k++)
        dl_d2[(j+1)*NP+k+1] += 
            xc[j][k_1] * xc[k][k_1] *(vr2w2-w);                /* x'x */
      }

      for(j=0;j<*nz2uz;j++) { 
        tmp1[j+1+*nx2uz] += zc[j][i] * vrw;
        dl_d2[1+j+*nx2uz] += zc[j][i] * (vr2w2 - w);          /* z'1 */
        dl_d2[(NP-1)+(j+1+*nx2uz)*NP] += 
            zc[j][i] * (vr3w2 - 2*vrw);                 /* z'sigma */ 
      for (k=j;k<*nz2uz;k++)
        dl_d2[(j+1+*nx2uz)*NP+k+1+*nx2uz] += 
            zc[j][i] * zc[k][i] *(vr2w2-w);                /* z'z */
      }
      for(j=0;j<*nt2uz;j++) { 
        tmp1[j+1+*nx2uz+*nz2uz] += txc[j][k_1]*tzc[j][i] * vrw;
        dl_d2[1+j+*nx2uz+*nz2uz] +=  
            txc[j][k_1] * tzc[j][i] * (vr2w2 - w);        /* t'1 */
        dl_d2[(NP-1)+(j+1+*nx2uz+*nz2uz)*NP] += 
            txc[j][k_1] * tzc[j][i] * (vr3w2 - 2*vrw);     /* t'sigma */ 
      for (k=j;k<*nt2uz;k++)
        dl_d2[(j+1+*nx2uz+*nz2uz)*NP+k+1+*nx2uz+*nz2uz] += 
            txc[j][k_1] * tzc[j][i] *txc[k][k_1] * tzc[k][i] 
                *(vr2w2-w);                                    /* t't */
      }

      for(j=0;j<*nx2uz;j++)  
      for (k=j;k<*nz2uz;k++)
        dl_d2[(j+1)*NP+k+1+*nx2uz] += 
            xc[j][k_1] * zc[k][i] *(vr2w2-w);                /* z'x */


      for(j=0;j<*nx2uz;j++)  
      for (k=j;k<*nt2uz;k++)
        dl_d2[(j+1)*NP+k+1+*nx2uz+*nz2uz] += 
            xc[j][k_1] * txc[k][k_1] * tzc[k][i] *(vr2w2-w);  /* t'x */

      for(j=0;j<*nz2uz;j++) 
      for (k=j;k<*nt2uz;k++)
        dl_d2[(j+1+*nx2uz)*NP+k+1+*nx2uz+*nz2uz] += 
            zc[j][i] * txc[k][k_1] * tzc[k][i] *(vr2w2-w);  /* t'z */

       } /* ends k_l ?? */
 
    for (j=0;j<NREG;j++) {
      dl_d1[j] += tmp1[j];
      for (k=j;k<NREG;k++) {
	dl_d2[j*NP+k] -= tmp1[j]*tmp1[k];
      }
      dl_d2[(NP-1)+NP*j] -= tmp1[j]*tmp_sig;
    }
    dl_d2[NP*NP-1] += tmp_sig_2 - tmp_sig*tmp_sig;
    dl_d1[NP-1] += tmp_sig;
  }
  dl_d1[NP-1] -= *n;

  /*account for the prior */

  k=0;
  for (i=0;i<NREG;i++) k += (amnt[i] > 0.0) ? 1 : 0;
  if (k>0) {
    for (i = 0; i<NREG; i++) 
      if (amnt[i]>0) {
	d__1 = beta[i] / *sigma;
	llik[0] += log(amnt[i]) / 2.0 - nrmcns - d__1 * d__1 * amnt[i] / 2.0 ;
	dl_d1[i] -= amnt[i]*beta[i]/sigma_sq;
	dl_d1[NP-1] += amnt[i]*beta[i]*beta[i]/sigma_sq;
	dl_d2[i*(NP+1)] -= amnt[i]/sigma_sq;
	dl_d2[NP*(NP-1)+i] += 2.0*amnt[i]*beta[i]/sigma_sq;
	dl_d2[NP*NP-1] -= 2.0 *beta[i]*beta[i]*amnt[i]/sigma_sq;
      }
    dl_d1[NP-1] -= (double) k;

  }

  /* assign the upper triangle of the hessian */

  for (j=0;j<NP;j++)
    for (k=j+1;k<NP;k++) 
        dl_d2[k*NP+j] = dl_d2[j*NP+k];

}
  
