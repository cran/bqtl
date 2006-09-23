/* twohkf2.c
   Wed Jul  7 17:05:48 PDT 1999
   C.C.Berry */

/* twohkf2.c - linearized posterior for two locus models
   returns:
   coefs - posterior averages for the coefficients
   coef1 - post. avg.s for the 1 locus coefficients
   marg - averages for the posterior per locus and model
   cond - the posteriors for the one locus models
   */



#include "lapadj.h"
#define NLOCS *nopt/3

void condreg(double *varx, double *covxy, double *vary, double *df,
             double *amount, double *optpri, double *inpri,
             longint *nrx, longint *invars ,longint *nin, longint *optvars,
             longint *nopt ,longint *optmax , longint *useopt, longint *pvt,
             longint *rank,
             double *wrksp, double *gama, double *bee,
             double *xx, double *xy,
             double *zz, double *zy, double *xz, double *beta,
             double *posterior,
             double *marg, double *cond, double *coefs,
             double *qraux, double *zraux, longint *zrank,
             double *tol);


void twohkf2(double *varx, double *covxy, double *vary, double *df,
	     double *amount, double *optpri, 
	     longint * nrx, longint *optvars ,
	     longint *nopt ,longint *optmax , longint *useopt, longint *pvt, 
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

  for (iloc=0;iloc<*nopt;iloc++) {
    locindx = iloc%NLOCS;
    useopt[locindx] = useopt[locindx+NLOCS] = useopt[locindx+2*NLOCS] = FALSE;
   
    *inpri = optpri[iloc];
    for (i=0;i<2;i++) invars[i] = optvars[iloc*2+i];
    *nin= (invars[1]==-1) ? 1 : 2 ;
    *nmax = *nin+2;
    
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
    if (invars[1]!=-1) 
      coef1[invars[1]] += bee[1]**tcond;
    
    for (i=0;i<*nopt;i++) 
      if (useopt[i] == TRUE) {
	j1= *nmax*i;
	j2= *optmax*i;
	coefs[invars[0]] += coefwk[j1]*posterior[i];
	if (invars[1]!=-1) 
	  coefs[invars[1]] += coefwk[j1+1]*posterior[i];
	j1 += *nin;
	coefs[optvars[j2]]+= coefwk[j1]*posterior[i];
	if (optvars[j2+1]!=-1)
	  coefs[optvars[j2+1]]+= coefwk[j1+1]*posterior[i];
      } 
      
    useopt[locindx] = useopt[locindx+NLOCS] = useopt[locindx+2*NLOCS] = TRUE;
    
  }
  for (i=0;i<NLOCS*2;i++) {
    coefs[i] /= totpst;
    coef1[i] /= tcpst;
  }
   
}
