/* swapf2.c   
   July 5,  1999
   C.C.Berry  */

/* swapf2.c 
   control execution of condreg.c to create samples from an f2 cross
   GIVEN the following:
   nreps - number of cycles to run
   nsteps - the number of steps in each cycle
   varx - a variance-covariance matrix
   covxy - vector of covariances
   vary - the variance of y
   df - the initial degrees of freedom in the setup (often n-1)
   nrx - dimension of varx
   ninopt - number of columns of optvars to include in models
   optblock - matrix of indicators of cols to be omitted if this col is used
   optcur - starting set of cols of optvars to be used
   invars - space for vector of variables to initialize the sequence
   locs - vector giving the locus associated with each variable
   nin - length of invars
   optvars - a matrix of additional variables - one column for each model
   nopt - number of columns of optvars
   optmax - number of rows of optvars
   useopt - logical vector, use this columns of optvars?
   amount - prior parameter for coef corresponding to each variable
   optpri - prior on each model
   inpri - prior on the common submodel
   
   RETURN:
   posterior - vector of posterior probs associated with selected models
   marg - the sum over all posteriors at each step
   cond - the posterior for the sub model at each step 
   coefs - matrix of coefficients from the models
   configs - array of variables used in each selected model
   altmarg - vector of sum of model posteriors
   altcoef - vector of sum of marginal means of coefs
   */


#include "lapadj.h"
     /* MAXMOD - number of model types per locus */
#define MAXMOD 3
     /* OPTROW - number of rows in optvars */
#define OPTROW 2
     /* NLOCS - number of loci */
#define NLOCS (*nopt/3)
     /* MAXVAR - maximum number of variables in all */
#define MAXVAR OPTROW**nstep

void swapf2( longint *nreps, longint *nstep,
             double *varx, double *covxy, double *vary, double *df,
	     double *amount, double *optpri, double *inpri, 
	     longint * nrx, longint *noptuse, longint *optblock, longint *optcur,
	     longint *invars ,longint *nin ,longint *optvars ,
	     longint *nopt ,longint *optmax , longint *useopt, longint *locs, 
	     longint *pvt, longint *rank,
	     double *wrksp, double *gama, double *bee, double *xx, 
	     double *xy, 
	     double *zz, double *zy, double *xz, double *beta, 
	     double *posterior, 
	     double *marg, double *cond, double *coefs, longint *configs, 
	     double *qraux, double *zraux,  longint *zrank,
	     double *tol, double *postwk, double *coefwk,
	     double *altmarg, double *altcoef)
{
  longint oc, NOptU, NVar;
  longint i,j,k, dl, ndelete, curpick, irep, istep, curloc, ncoef;
  double unifr, margwk[1], condwk[1], postsum, tmp;
  
  NOptU = *noptuse;

  seed_in((long *)NULL);

  for (irep = 0; irep < *nreps; irep ++)
    for (istep=0; istep <*nstep; istep++) {
   
      /* remove 1st locus */
      *inpri = 1.0;
      NVar = 0;
      dl = optcur[0];
      for (i=1; i<NOptU; i++){
	oc = optcur[i];
	optcur[i-1]=oc;
	invars[NVar++] = optvars[oc*OPTROW];
	*inpri *= optpri[oc];
	if ( (oc=optvars[OPTROW*oc+1])!= -1)
	  invars[NVar++] = oc;
      }
      
      for (i=NVar;i<MAXVAR;i++) invars[i] = -1;
      
      /* mark models using former 1st locus eligible */
      /* assumed order <additive only, dominance only, both> repeatedly */

      useopt[dl]= TRUE;
      useopt[optblock[(MAXMOD-1)*dl]]= TRUE;
      useopt[optblock[(MAXMOD-1)*dl+1]] = TRUE;
      
      ncoef = OPTROW + NVar; /* needed to index coefwk */
      /* eval models */

      condreg(varx, covxy, vary, df,
	      amount, optpri, inpri, 
	      nrx, invars , &NVar ,optvars ,
	      nopt ,optmax , useopt, pvt, rank,
	      wrksp, gama, bee, xx, xy, 
	      zz, zy, xz, beta, postwk, 
	      &postsum, condwk, coefwk, qraux, zraux, zrank,
	      tol);
      /* pick a model */
      unifr = unif_rand();
      /*     postsum = 0.0;
	     for (i=0;i<*nopt;i++) postsum += postwk[i];*/
      unifr *= postsum;
      tmp = 0.0;

      for (i=0;tmp < unifr;i++) tmp +=postwk[i];

      curpick = i-1;
      optcur[NOptU-1]=curpick;

      /* update altmarg and altcoef */

      for (i=0;i<*nopt;i++) {
	tmp = postwk[i]/postsum;
	if (useopt[i]==TRUE) {
	  altmarg[i] += tmp;
	  for (j=NVar;j<*optmax&&optvars[OPTROW*i+j]>=0;j++) 
	    altcoef[optvars[OPTROW*i+j]] += coefwk[i*ncoef+j]*tmp;
	  for (j=0;j<NVar;j++) 
	    altcoef[invars[j]] += coefwk[i*ncoef+j]*tmp;
	}
      }
      /* add new locus */
      for (i=0;i<*optmax&&optvars[OPTROW*curpick+i]>=0;i++){
	invars[NVar++] = optvars[OPTROW*curpick+i];
      }
    

      /* mark models incompatible with new locus ineligible */ 
      useopt[curpick]= FALSE;
      useopt[optblock[(MAXMOD-1)*curpick]]= FALSE;
      useopt[optblock[(MAXMOD-1)*curpick+1]] = FALSE;
    
      /* store results of current step */
      posterior[irep**nstep+istep] = postwk[curpick];
      marg[irep**nstep+istep] = postsum;
     
      cond[irep**nstep+istep] = *condwk;
      for (i=0 ; i<NVar ; i++) {
	configs[(irep**nstep+istep)*MAXVAR+i] = invars[i];
	coefs[(irep**nstep+istep)*MAXVAR+i] = coefwk[curpick*ncoef+i];
      }
      for (i=NVar ;i<MAXVAR; i++) {
	configs[(irep**nstep+istep)*MAXVAR+i] = -1;
	coefs[(irep**nstep+istep)*MAXVAR+i] = 0.0;
      }
    }
  seed_out((long *)NULL);
}

