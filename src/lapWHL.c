/* lapWHL.c */
/* Wed Feb  3 13:32:08 PST 1999 */
/* C.C.Berry */
/* Do the computation to create the Laplace Approximation for
   the Gaussian QTL setup. */

#include "lapadj.h"
#ifdef USING_R
#include <R_ext/Applic.h>
#else
void F77_NAME(dqrsl1)();
void F77_NAME(dqr)();
#endif


void llkEm(longint *, double **, double **,double **,double **,double *,
            double *,double *,double *,double *,double *,double *,double *,
            double *,double *,double *,double *,double *,double *,double *,
            double *,double *,double *,double *,double *,double *,double *,
            longint *, double *, longint *, double *);

void normLogLik( longint*, double **,double **,double **,double **,double *,
                double *,double *,double *,double *,double *,
                double *,double *,double *,double *,double *,double *,double *);

int conCk(double *oldllk, double *newllk, double *paradj, longint *NP,
	  double *tol);

void F77_NAME(hessup)( double *, double *, double *, double *, 
                       double *, double *, int *,double *, double *, 
                       double *, double *, double *, int*, double *);

void lapWhl(double *y, double *wt, double *amnt, longint *nparm, double **xc,
	    double **zc, double **txc, double **tzc, double *fits,
	    double *res, double *vsum, double *vsum2, double *vysum,
	    double *xvx, double *xvy, double *df, double *qrxvx, double *qraux,
	    longint *pvt, double *wrksp, double *tol, double *llk, double *dgr,
	    double *dparm, double *newgrd, double *oldgr, double *curprm,
	    double *oldprm, longint *reset, double *bk, double *newhss,
	    double *bks, double *qrbk, double *paradj, double *newprm,
	    double *emparm, longint *iter, longint *nem, double *casewt)
     
{

    longint j,notdn, nem1=1, firstTry, *n,  *nrx,  *ncx,  *nloc, *ncz,
	*nx2uz,  *nz2uz,  *nt2uz, *nreg, *np;  
    double d__1, sigma, sigma2, ss, oldllk, newllk, lprior, tsttol;
    int ireset[1], inp[1], i, need_search=1;
  
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
    *ireset = *reset;
    *inp=*np;
    firstTry=1;
    oldllk = *llk;
    
    for (j=0;j<*np;j++) pvt[j]=j;
    notdn = 1;
    sigma = sqrt(curprm[*np-1]);
    
    while(notdn == 1L && (*iter)-- > 0) {
	if (sigma > 0.0) {
	    llkEm(nparm, xc, zc, txc, tzc, &sigma, amnt, fits, 
		  y, curprm, wt, res,
		  vsum, vsum2, vysum, &lprior, &newllk, xvx,
		  xvy, df, newprm, &ss, newgrd, newhss, &sigma2,
		  qrxvx, qraux, pvt, wrksp,&nem1, casewt);
	    newllk += lprior;
	    newprm[*np-1] = sigma2;
	    /* check if at max on first loop */
	    if (firstTry==1) { 
		firstTry = 0;
		if (conCk(&oldllk,&newllk,newgrd,np,tol) ==1) 
		    break; /* out of while */
		
	    }
	}
	if (notdn == 1L && (sigma <= 0.0 || newllk < oldllk)) {
	    if (need_search==1) {                                /* do line search */
		if (sigma <= 0.0) {                         /* bail out */
		  warning( "convergence failed" ) ;
		    notdn = 0;
		}
		else {                                      /* line search should work */
		    for (i=0;i<*np;i++) curprm[i] = (curprm[i] + oldprm[i])/2.0 ;
		    sigma = sqrt(curprm[*np-1]);
		    /* now return to top of while and try again */
		}
	    }
	    else {                                               /* do em update */
		sigma = sqrt(oldprm[*np-1]);
		*reset = 1;
		if (*nem==0){
		    CPY(emparm, curprm, np);   /* retrieve last EM update */
		}
		else
		{
		    llkEm(nparm,  xc, zc, txc, tzc, &sigma, 
			  amnt, fits, y, oldprm, wt, res,
			  vsum, vsum2, vysum, &lprior, &newllk, xvx,
			  xvy, df, curprm, &ss, newgrd, newhss, &sigma2,
			  qrxvx, qraux, pvt, wrksp,nem, casewt);
		    curprm[*np-1]=sigma2;
		    sigma = (curprm[*np-1] > 0) ? sqrt(curprm[*np-1]): -1.0;
		    newllk += lprior;
		    need_search = 1;
		}
		if ((NTFINITE(newllk)) == 0) 
		  {
		    tsttol = (d__1 = (newllk - oldllk)*2.0/(newllk + oldllk), ABS(d__1));
		    for (j = 0; j <*np; j++) tsttol += (d__1 = newgrd[j], ABS(d__1));
		    notdn = (tsttol > *tol) ? 1 : 0; 
		} 
		else {
		    notdn = 1;
		}
	    }
	}
	else 
	  {  
	    notdn = (fabs(newllk - oldllk) > *tol) ? 1 : 0; /* fabs needed?? */   
	    CPY(curprm, emparm,np); /*retain curprm as emparm*/
	    F77_CALL(hessup)(dgr, dparm, newgrd, oldgr, curprm, 
			     oldprm, ireset, bk, newhss, bks, qrbk, 
			     paradj, inp, newprm);
	    
	    CPY(emparm, oldprm, np); /* oldprm has best, tested values */
	    sigma = (curprm[*np-1] > 0) ? sqrt(curprm[*np-1]): -1.0;
	    CPY(newprm, emparm, np); /* newprm was the EM update */
	    oldllk = newllk;
	    CPY(newgrd,oldgr, np);
	    if (sigma<0.0) notdn = 1;
	    need_search = 0;
	}
    }
    
    *reset = *ireset;
    sigma = sqrt(oldprm[*np-1]);  /* oldprm = best proven choice */
    /* initialize res and fits for normLogLik */
    llkEm(nparm,  xc, zc, txc, tzc, &sigma, 
	  amnt, fits, y, oldprm, wt, res,
	  vsum, vsum2, vysum, &lprior, &newllk, xvx,
	  xvy, df, curprm, &ss, newgrd, newhss, &sigma2,
	  qrxvx, qraux, pvt, wrksp,&nem1, casewt);
    normLogLik(nparm,xc,zc,txc,tzc,y,fits,res,&sigma,amnt,oldprm, 
	       wrksp,qrxvx,qraux, 
	       newgrd,newhss,llk,casewt); 
}

