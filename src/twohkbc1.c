/* twohkbc1.c
   Thu Jul 29 14:49:18 1999
   C.C.Berry */

/* twohkbc1.c - linearized posterior for two locus models
   returns:
   coefs - posterior averages for the coefficients
   coef1 - post.avg.s for the 1 locus coefficients
   marg - averages for the posterior per locus and model
   cond - the posteriors for the one locus models
   */

#include "lapadj.h"
#define NLOCS *nopt

void condreg(double *varx, double *covxy, double *vary, double *df,
	     double *amount, double *optpri, double *inpri, 
	     longint * nrx, longint *invars ,longint *nin ,longint *optvars ,
	     longint *nopt ,longint *optmax , longint *useopt, longint *pvt, 
	     longint *rank,
	     double *wrksp, double *gama, double *bee, 
	     double *xx, double *xy, 
	     double *zz, double *zy, double *xz, double *beta,
	     double *posterior, 
	     double *marg, double *cond, double *coefs, 
	     double *qraux, double *zraux, longint *zrank, 
	     double *tol);

void twohkbc1(double *varx, double *covxy, double *vary, double *df,
	     double *amount, double *optpri, 
	     longint * nrx, longint *optvars ,
	     longint *nopt, longint *optmax, longint *useopt, longint *pvt, 
	     longint *rank,
	     double *wrksp, double *gama, double *bee, 
	     double *xx, double *xy, 
	     double *zz, double *zy, double *xz, double *beta,
	     double *posterior, 
	     double *marg, double *cond, double *coefwk, double *coefs,
	     double *coef1,
	     double *qraux, double *zraux, longint *zrank, 
	     double *tol)
{
  double inpri[1], tmarg[1], tcond[1], totpst, tcpst; 
  longint invars[2] ,nin[1], nmax[1];
  int i,j1, j2, iloc, locindx;

  totpst = 0.0;
  tcpst = 0.0;
  *nin= 1;
  *nmax = 2;

  for (iloc=0;iloc<*nopt;iloc++) {
    locindx = iloc;

    useopt[locindx] =  FALSE;

    *inpri = optpri[iloc];

    invars[0] = optvars[iloc];

    *tmarg = 0.0;

    condreg(varx, covxy, vary, df,
	    amount, optpri, inpri, 
	    nrx, invars ,nin ,optvars ,
	    nopt ,optmax , useopt, pvt, rank,
	    wrksp, gama, bee, xx, xy, 
	    zz, zy, xz, beta, posterior, 
	    tmarg, tcond, coefwk, qraux, zraux, zrank,
	    tol);

    marg[iloc]=*tmarg;
    cond[iloc]=*tcond;
    totpst += *tmarg;
    tcpst += *tcond;
    coef1[invars[0]] += bee[0]**tcond;
    for (i=0;i<*nopt;i++) 
      if (useopt[i] == 1) {
	j1= *nmax*i;
	j2= *optmax*i;
	coefs[invars[0]] += coefwk[j1]*posterior[i];
	j1 += *nin;
	coefs[optvars[j2]]+= coefwk[j1]*posterior[i];
      } 
     
    useopt[locindx] =  TRUE;

  }
 
  for (i=0;i<NLOCS;i++) {
    coefs[i] /= totpst;
    coef1[i] /= tcpst;
  }

}
