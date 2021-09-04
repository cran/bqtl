/* llkEm.c */
/* Wed Feb  3 13:08:18 PST 1999 */
/* C.C.Berry */

/* llkEm.c is called by lapadj.c and by lapWhl.c */

#include "lapadj.h"
#ifdef USING_R
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>
#else
int F77_NAME(dqrsl1)();
int F77_NAME(dqr)();
#endif

void llkEm(longint *nparm, double **xc,double **zc, double **txc, double **tzc,
	   double *sigma, double *amnt, double *fits, double *y, double *beta,
	   double *wt, double *res, double *vsum, double *vsum2, double *vysum,
	   double *lprior, double *llik, double *xvx, double *xvy, double *df,
	   double *coefs, double *ss, double *grad, double *hess, double *sigma2,
	   double *qrxvx, double *qraux, longint *pvt, double *wrksp,
	   longint *n_em, double *casewt)
     
{
    /* Initialized data  */
    
    double rt2pi = 2.5066282746310002416, tol = QR_TOL, d__1, nrmcns,
	tmp, tmp2, tmp_sig, lgrtN, tmp_yj, tmp_yij, *xpt, *xipt, *xjpt,
	*zipt, *zjpt, *wtpt, *wtend, *xvxpt;

    longint *n,  *nrx,  *ncx,  *nloc, *ncz,  *nx2uz,  *nz2uz,  *nt2uz,
        *nreg, *np, x_dim1, xvx_dim1, hess_dim1, rank, hndrd, dq[2], one;


    int pt_1[1] = {1}, i_iter, i,j,k,isub;
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
    
    lgrtN = 1.0;
    if (*casewt > 0.0) for (i=0;i<*n;i++) lgrtN*=casewt[i];
    lgrtN = log(lgrtN) / 2.0 ;
    
    hess_dim1 = *np; /* this was  *np-1 ??? */ 
    xvx_dim1 = *nreg;
    x_dim1 = *nrx;
    tmp_sig = *sigma; /* initialize in case looping is needed */
    
    for (i=0;i<*nreg;i++) coefs[i]=beta[i];
    for (i_iter=0;i_iter<*n_em;i_iter++){
      nrmcns = rt2pi * tmp_sig;
      
      /* fits - since X can be big, access it by columns */
      for (i=0;i<*nrx;i++) fits[i] = coefs[0];
      for (j=0;j<*nx2uz;j++) {
	tmp=coefs[j+1];
        xpt = xc[j]; /* pt to top of j-th column */
	for (i=0;i<*nrx;i++) {
	  fits[i] += tmp* *xpt++;
	}
      }
      
	for (j=0;j<*n;j++){
          /* if fixed covars are needed, subtract the fitted values from y */
          tmp_yj = y[j];
          for (i=0;i<*nz2uz;i++) tmp_yj -= coefs[i+1+*nx2uz]*zc[i][j];
	  /* res holds posterior wts */
	    for (i=0;i<*nrx;i++){
	      tmp_yij = tmp_yj;
	      for (k=0;k<*nt2uz;k++) 
		tmp_yij -= coefs[ k+1 + *nx2uz + *nz2uz ]*txc[k][i]*tzc[k][j];
	      tmp = (tmp_yij - fits[i]) / tmp_sig;
	      res[i + j * x_dim1] =  (*casewt<0) ?
		exp(-0.5*tmp*tmp ) / nrmcns * wt[i + j * x_dim1] :
		exp(-0.5*tmp*tmp*casewt[j] ) / nrmcns * wt[i + j * x_dim1]; 
	    }
	  
	    
	}
    /* from here on, res only shows up inside loops upto *nrx */
      *llik =  lgrtN;
      
	  for (i=0;i<*n;i++){
	    vsum[i] = 0.0;
	    for (j=0;j<*nrx;j++){
	      vsum[i] += res[j + i * x_dim1];
	    }
	    if ( vsum[i] <= 0) {
	      warning("machine zero likelihood");
	      *llik = -HUGE_VAL;
	      return;
	    }
	    *llik += log(vsum[i]);
	  }

      *lprior = 0.0;
      for (i=0;i<*nreg;i++){
	if (amnt[i] > 0.0) {
	  tmp = coefs[i] / tmp_sig;
	  *lprior = *lprior + log(amnt[i])/2.0 - log(nrmcns) - tmp * 
	      tmp * amnt[i] / 2.0;
	}
      }

      for (j=0;j<*nrx;j++) {
	vsum2[j] = 0.0;
	vysum[j] = 0.0;
      }
      if (*casewt<0.0) {
	for (j=0;j<*n;j++){
	  for (i=0;i<*nrx;i++){
            res [i + j * x_dim1] /= vsum[j];
	    vsum2[i] += res[i + j * x_dim1];
	    vysum[i] += res[i + j * x_dim1] * y[j];
	  }
	}
      }
      else {
	for (j=0;j<*n;j++){
	  for (i=0;i<*nrx;i++){
            res [i + j * x_dim1] /= vsum[j];
	    vsum2[i] += res[i + j * x_dim1] * casewt[j];
	    vysum[i] += res[i + j * x_dim1] * y[j] * casewt[j];
	  }
	}
      }

/* xvx is a dispersion matrix whose row|column partitions are:
    constant, locus varying covariates, subject varying covariates,
    both locus and subject varying covariates.
   Some of the last two partitions may be empty, so we start by
    efficiently filling the upper left submatrix, then if we need to do
    other parts we loop thru subjects and do the remainder.

*/

      for (i=0;i<*nreg;i++){
	xvx[i + i * xvx_dim1] = amnt[i];
	for (j=i+1;j<*nreg;j++){
	  xvx[j + i * xvx_dim1] = 0.0 ;
	}
      }

      xvy[0]=0.0;
      for (k=0;k<*nrx;k++) {
	xvx[0] += vsum2[k];
	xvy[0] += vysum[k];
      }
      qrxvx[0]=xvx[0];


      for (i=1 ; i<(*nx2uz+1) ; i++){
	xvx[i] = F77_CALL(ddot)(nrx,xc[i-1],pt_1,vsum2,pt_1);
	xvy[i] = F77_CALL(ddot)(nrx,xc[i-1],pt_1,vysum,pt_1);
	xvx[i*xvx_dim1] = xvx[i];
/*	qrxvx[i*xvx_dim1] = xvx[i];
	qrxvx[i] = xvx[i];
*/
	for (j=i;j<(*nx2uz+1);j++) {
	  tmp=0.0;
          xipt = xc[i-1];
          xjpt = xc[j-1];
          wtpt = vsum2;
	  for (k=0;k<*nrx;k++){
	    tmp += *xipt++ * *xjpt++ * *wtpt++;
	  }
	  xvx[j+i*xvx_dim1] += tmp;

/*	  xvx[i+j*xvx_dim1] = xvx[j+i*xvx_dim1];
	  qrxvx[j+i*xvx_dim1] = xvx[j+i*xvx_dim1];
	  qrxvx[i+j*xvx_dim1] = xvx[j+i*xvx_dim1];
*/

	}
      }



/* z1 block */
        xvxpt= xvx + *nx2uz + 1; 
        for (j=0;j<*nz2uz;j++){
            tmp= 0.0;
            zjpt = zc[j];
            for (k=0;k<*n;k++) tmp+= *zjpt++;
            *xvxpt = tmp; 
            xvxpt ++;
        
	}
/* zz block[j,j]  */
        for (j=0;j<*nz2uz;j++){
           xvxpt = xvx + (j+ 1 + *nx2uz)*(1 + xvx_dim1); 
            for (i=j;i<*nz2uz;i++){
            tmp= (i == j) ? *xvxpt : 0.0;
            zjpt = zc[j];
            zipt = zc[i];
            for (k=0;k<*n;k++) tmp += *zjpt++ * *zipt++;
            *xvxpt = tmp; 
            xvxpt ++;
            }
        }
/* zy */
        for (j=0;j<*nz2uz;j++){
            tmp= 0.0;
            zjpt = zc[j];
            xpt = y;
            for (k=0;k<*n;k++) tmp += *zjpt++ * *xpt++;
            xvy[1+*nx2uz+j] = tmp; 
        }
        
/* zx block */
        xvxpt = xvx + xvx_dim1 + *nx2uz +1; 
        for (j=0;j<*nz2uz;j++){
            for (i=0;i<*nx2uz;i++){
                zjpt=zc[j];
                tmp=0.0;
                wtpt=res;
                for (isub=0;isub<*n;isub++){
                    xipt=xc[i];
                    wtend = wtpt + *nrx;
		    tmp2 = 0.0;
                    while (wtpt < wtend) tmp2 += *xipt++ * *wtpt++;
                    tmp += tmp2 * *zjpt++;
                }
                *(xvxpt + j + i * xvx_dim1) = tmp;
            }
        }

/* t1 block */
        xvxpt = xvx + 1 + *nx2uz + *nz2uz; 
        for (i=0;i<*nt2uz;i++){
            zipt = tzc[i];
            wtpt = res ;
            for (isub=0;isub<*n;isub++){
                xipt = txc[i];
                wtend = wtpt + *nrx;
                tmp=0.0;
                while (wtpt < wtend) tmp += *xipt++ * *wtpt++;
                *xvxpt += *zipt++ * tmp;
            }
            xvxpt++;         
        }

/* ty block */
	if (*nt2uz!=0) {
	  xvxpt = xvy + 1 + *nx2uz + *nz2uz;
	  
	  for (i=0;i<*nt2uz;i++){
	      *xvxpt = 0.0;
	      zipt = tzc[i];
	      xpt = y;
	      wtpt = res ;
	      for (isub=0;isub<*n;isub++){
		  xipt = txc[i];
		  wtend = wtpt + *nrx;
		  tmp=0.0;
		  while (wtpt<wtend) tmp += *xipt++ * *wtpt++;
		  *xvxpt += *zipt++ * *xpt++ * tmp;
	      }
	      xvxpt++; /* &xvy[i+1] */        
	  }
        }
/*  tt  */
        for (i=0;i<*nt2uz;i++){
            xvxpt = xvx + ( i + 1 + *nx2uz + *nz2uz) * 
                ( xvx_dim1 + 1 ); 
            for (j=i;j<*nt2uz;j++){
               wtpt = res;
               zipt = tzc[i];
               zjpt = tzc[j];
               for (isub=0;isub<*n;isub++){
                    xipt = txc[i];
                    xjpt = txc[j];
                    tmp=0.0;
                    wtend = wtpt + *nrx;
                    while(wtpt<wtend) tmp += *xipt++ * *xjpt++ * *wtpt++;
                    *xvxpt += tmp * *zipt++ * *zjpt++;
                }
                xvxpt++;
            }
        }
    
        for (i=0;i<*nx2uz;i++){
/* tx [1,i] */
            xvxpt = xvx +  (i + 1) * xvx_dim1 + 1 + 
                *nx2uz + *nz2uz; 
            for (j=0;j<*nt2uz;j++){
               wtpt = res;
               zjpt = tzc[j];
               for (isub=0;isub<*n;isub++){
                    xipt = xc[i];
                    xjpt = txc[j];
                    tmp=0.0;
                    wtend = wtpt + *nrx;
                    while(wtpt<wtend) tmp +=  *xipt++ * *xjpt++ * *wtpt++;
                    *xvxpt += tmp *  *zjpt++;
                }
                xvxpt++;
            }
        }
        
        for (i=0;i<*nz2uz;i++){
/* tz[1,i] */
            xvxpt = xvx +  (i + 1 + *nx2uz) * xvx_dim1 + 1 +        
                *nx2uz + *nz2uz; 
            for (j=0;j<*nt2uz;j++){
               zipt = zc[i];
               zjpt = tzc[j];
               wtpt = res;
               for (isub=0;isub<*n;isub++){
                    xjpt = txc[j];
                    tmp=0.0;
                    wtend = wtpt + *nrx;
                    while( wtpt < wtend ) tmp +=  *xjpt++ * *wtpt++;
                    *xvxpt += tmp * *zipt++ * *zjpt++;
                }
                xvxpt++;
            }
        }

/* fill in the other triangle and load qrxvx */
          for (i=0;i<xvx_dim1;i++)
            for (j=0;j<i;j++) {
              xvx[j+i*xvx_dim1] = xvx[i+j*xvx_dim1];
        }
	  i = (xvx_dim1*xvx_dim1);
	  CPY(xvx,qrxvx, &i );	

/* restore res */
    wtpt = res;
    for (i=0;i<*n;i++){
        tmp = vsum[i];
        for (wtend = wtpt + *nrx;wtpt < wtend; wtpt++) *wtpt*=tmp;
    }

      *df = (double) (*n);
      for (i=0;i<*nreg;i++) {
	if (amnt[i] > 0.0) {
	  (*df)++;
	}
      }
 
      rank = *nreg;
      dq[0] = *nreg;
      dq[1] = *nreg;
      for (i=0;i<*nreg;i++) pvt[i]= (i+1);
#ifndef USING_R
      F77_CALL(dqr)(qrxvx, dq, pvt, qraux, &tol, wrksp, &rank);
#else
      F77_CALL(dqrdc2)(qrxvx,dq,dq,dq,&tol,&rank, qraux,pvt,wrksp);
#endif
   
      dq[0] = *nreg;
      dq[1] = *nreg;
      one = 1;
      hndrd = 100;

#ifndef USING_R
      F77_CALL(dqrsl1)(qrxvx, dq, qraux, &rank, xvy, 
	      &one, wrksp, coefs, &hndrd, &one);
#else

      F77_CALL(dqrsl)(qrxvx, dq,dq,&rank,qraux,xvy,wrksp,
	    wrksp,coefs,wrksp,wrksp,&hndrd, &one);
#endif


      if (rank<*nreg) { /*re-order coefs*/
	for (i=0;i<rank;i++) wrksp[pvt[i]-1] = coefs[i];
	for (i=rank;i<*nreg;i++) wrksp[pvt[i]-1] = 0.0;
	CPY(wrksp,coefs,nreg);
	warning( "deficient rank in llkEm" ) ;
      }

      *ss = 0.0;
      if (*casewt<0.0) {
	for (i = 0; i <*n;i++) {
	  *ss += y[i] * y[i];
	}
      }
      else {
	for (i = 0; i <*n;i++) {
	  *ss += y[i] * y[i] * casewt[i];
	}
      }
      for (i = 0; i <*nreg; i++) {
	*ss -= coefs[i]*xvy[i];
      }
 
      if (i_iter==0){
	for (i = 0; i <*nreg; i++) {
 	  grad[i] = (xvy[i] - 
	    F77_CALL(ddot)(nreg,&xvx[i * xvx_dim1],pt_1,beta,pt_1))/
	    (*sigma * *sigma);
	}
		
	d__1 = *sigma, d__1 *= d__1;
		
	grad[*np-1] = (*ss / (d__1 * d__1) - *df / d__1) / 2.0;
	
	for (i = 0; i<*nreg; i++) {
	  for (j = i; j <*nreg; j++) {
	    hess[i + j * hess_dim1] = -xvx[i + j * xvx_dim1] / (d__1);
	    hess[j + i * hess_dim1] = hess[i + j * hess_dim1];
	  }
	  hess[i + ( hess_dim1 - 1 ) * hess_dim1] = 0.0;
	  hess[( 1 + i ) * hess_dim1 - 1] = 0.0;
	}
	
	hess[hess_dim1 * hess_dim1 - 1] = -(*df) / (d__1 * d__1 * 2.0);
      }
      *sigma2 = *ss / *df;
      tmp_sig = sqrt(*sigma2);
    }
    
}
