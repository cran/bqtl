/* condreg.c */
/* Thu Jul  1 13:04:42 PDT 1999 */
/* C.C.Berry */

/* condreg.c - 
   GIVEN the following:
   varx - a variance-covariance matrix
   covxy - vector of covariances
   vary - the variance of y
   df - the initial degrees of freedom in the setup (often n-1)
   nrx - dimension of varx
   invars - vector of variables included in all models
   nin - length of invars
   optvars - a matrix of additional variables - one row for each model
   nopt - number of rows of optvars
   optmax - number of columns of optvars
   amount - prior parameter for coef corresponding to each variable
   optpri - prior on each model
   inpri - prior on the common submodel
   
   RETURN:
   posterior - vector of posterior probs associated with models
   marg - the sum over all posteriors
   cond - the posterior for the sub model that is common to all models
   coefs - matrix of coefficients from the models
   */

#include "lapadj.h"
#ifdef USING_R
#include <R_ext/Applic.h>
#include <Rmath.h>
#else
int F77_SYMBOL(dqrsl1)();
int F77_SYMBOL(dqr)();
#endif

/* vrload  Jun 30 14:20:35 PDT 1999 
   copy selected rows and cols from a large matrix to a
   small, comapct matrix 
*/

static void vrload(double *from, double *to,  longint *rows, longint
		   *cols, longint *nr, longint *nc, longint *nrfrom)
{ 
    int i,j;
    double *tres, *tfrom;
    
    tres = to;
    
    for (j=0;j<*nc;j++){
	tfrom = from + (int) (*nrfrom*cols[j]);
	for (i=0;i<*nr;i++) {
	    *tres = *(tfrom + rows[i]);
	    tres++;
	}
    }
}

/* resmat  Wed Jun 30 16:14:25 PDT 1999
   find the adjusted varcov matrix, given a varcov matrix,
   regression coefs, and covars for other variables 
*/

static void resmat(double *zz, double *coef, double *covar, longint *nrz, longint
		   *ncz, longint *nrx)
{
    int i,j,k;
    longint one[1];
    one[0] = 1;
    for (i=0;i<*nrz;i++)
	for (j=0;j<*ncz;j++){
	    for (k=0;k<*nrx;k++)
		zz[ i + *nrz * j ] -=  coef[*nrx*i+k]*covar[*nrx*j+k];
	}
}


void condreg(double *varx, double *covxy, double *vary, double *df,
	     double *amount, double *optpri, double *inpri, 
	     longint *nrx, longint *invars ,longint *nin ,longint *optvars ,
	     longint *nopt ,longint *optmax , longint *useopt, longint *pvt, 
	     longint *rank,
	     double *wrksp, double *gama, double *bee, 
	     double *xx, double *xy, 
	     double *zz, double *zy, double *xz, double *beta,
	     double *posterior, 
	     double *marg, double *cond, double *coefs, 
	     double *qraux, double *zraux, longint *zrank, 
	     double *tol)
{
  int i, j, imodel;
  longint   dq[2], dq2[2], nvars, *curvars;
  longint zero[1], one[1], hundred[1];


  double varygx, varygz, tmp, det, *curcoefs, df2, df3;
  
  zero[0]=0; one[0]=1;hundred[0]=100;
  *marg = 0.0;

  /* copy xx, get bee, common model posterior */

  vrload(varx,xx,invars,invars,nin,nin,nrx);
  for (i=0;i<*nin;i++) xx[i*(*nin+1)] += amount[invars[i]];

  dq[0] = dq[1] = *nin;
#ifndef USING_R
  F77_CALL(dqr)(xx,dq,pvt,qraux,tol,wrksp,rank);
#else
  F77_CALL(dqrdc2)(xx,&dq[0],&dq[0],&dq[1],tol,rank,qraux,pvt,wrksp);
#endif

  vrload(covxy,xy,invars,zero,nin,one,nrx);

#ifndef USING_R
  F77_CALL(dqrsl1)(xx,dq,qraux,rank,xy,one,wrksp,bee,hundred,one);
#else 
  F77_CALL(dqrsl)(xx, dq, dq, &dq[1], qraux, xy, &tmp, bee, bee, 
		    &tmp, &tmp, hundred , one);
#endif
  *one=1;*hundred=100;
  varygx = *vary;
  for (i=0;i<*nin;i++) varygx -= bee[i]*xy[i];

  det= *inpri;
  df2=*df;
  for (i=0;i<*nin;i++) 
    if ( amount[invars[i]]>0.0 ) {
      df2++;
      det *=  sqrt( amount[invars[i]] / fabs(xx[i*(*nin+1)]) );
    }
    else {
      det /= sqrt(fabs(xx[i*(*nin+1)]));
    }
  df2 -= *nin;
  *cond = det/pow( varygx / *vary, df2/2.0);
 
  /* loop thru all models - getting posteriors, etc */

  for (imodel=0;imodel<*nopt;imodel++) {
    if (useopt[imodel]==FALSE) {
      posterior[imodel]= 0.0;}
    else
      {
	curcoefs = &coefs[imodel*(*optmax+*nin)];/* pt to coefs*/
	curvars = &optvars[imodel**optmax];      /* pt to extra vars */
	
	for (j=0; (j<*optmax) & (curvars[j] >= 0L) ;j++); /* neg vals ==> no entry*/
	nvars = j;
	
	vrload(varx,zz,curvars,curvars,&nvars,&nvars,nrx); /* regressor varmat */
	for (i=0;i<nvars;i++) zz[i*(nvars+1)] += amount[curvars[i]];
	vrload(varx,xz,invars,curvars,nin,&nvars,nrx); /* cov of in/out vars */

#ifndef USING_R
	F77_CALL(dqrsl1)(xx,dq,qraux,rank,xz,&nvars,wrksp,gama,hundred,one);
#else
	for (i=0;i<nvars;i++) {
	  F77_CALL(dqrsl)(xx, dq, dq, &dq[1], qraux, &xz[i*dq[0]], 
			    &tmp, &gama[i*dq[0]], &gama[i*dq[0]], 
			    &tmp, &tmp,hundred , one);
	}
#endif
	*one=1;*hundred=100;
	resmat(zz,gama,xz,&nvars,&nvars,nin); /* adjust varmat */
	vrload(covxy,zy,curvars,zero,&nvars,one,nrx);
	resmat(zy,gama,xy,&nvars,one,nin);
	dq2[0] = dq2[1] = nvars;
#ifndef USING_R
	F77_CALL(dqr)(zz,dq2,pvt,zraux,tol,wrksp,zrank);
	F77_CALL(dqrsl1)(zz,dq2,zraux,zrank,zy,one,wrksp,
	   &curcoefs[*nin],hundred,one);
#else
	F77_CALL(dqrdc2)(zz,&dq2[0],&dq2[0],&dq2[1],tol,zrank,zraux,pvt,wrksp);
	F77_CALL(dqrsl)(zz, dq2, dq2, zrank, zraux, zy, &tmp,
			  &curcoefs[*nin], &curcoefs[*nin], 
			  &tmp, &tmp,hundred , one);
#endif
	*one=1;*hundred=100;
	varygz = varygx;
	resmat(&varygz,&curcoefs[*nin],zy,one,one,&nvars);
	CPY(bee,curcoefs,nin); /* now adjust low order coefs */
	for (i=0; i<*nin;i++)
	  for (j=0;j<nvars;j++)
	    curcoefs[i] -= curcoefs[*nin+j]*gama[i+j**nin];
	
	tmp= det*optpri[imodel];
	df3=df2;
	for (i=0;i<nvars;i++)
	  if (amount[curvars[i]]>0.0){
	    tmp *=  sqrt( amount[curvars[i]] / fabs(zz[i*(nvars+1)]) );
	    df3++;
	  }
	  else {
	    tmp /= sqrt( fabs(zz[i*(nvars+1)]) ); 
	  }
	df3 -= nvars;
	tmp /= pow(varygz/ *vary,df3/2.0);
	posterior[imodel] = tmp;
	*marg += tmp;
      }
  }
}
