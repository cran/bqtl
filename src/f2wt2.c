/* f2wt.c
   March 5, 1999 
   C.C. Berry
   Create proto.wt.object in make.weight.f2
*/

#include "lapadj.h"

/* 
   adj33
   Tue Mar 16 15:49:32 PST 1999
   fit a table given margins and a genetic parameter
   */

#include "lapadj.h"

static void adj33(double a[3], double b[3], double lambda[1], double tab[9], 
	   double tol[1], longint maxIter[1])
{
  double RowAdj, ColAdj, MaxAdj, MinAdj, p, q;
  int i;
  
  q = (1.0 + lambda[0])/2.0;
  p = (1.0 - lambda[0])/2.0;

  tab[0] = a[0]*q*q;
  tab[1] = a[1]*q*p;
  tab[2] = a[2]*p*p;
  tab[3] = 2.0*a[0]*p*q;
  tab[4] = a[1]*(q*q+p*p);
  tab[5] = 2.0*a[2]*p*q;
  tab[6] = a[0]*p*p;
  tab[7] = tab[1];
  tab[8] = a[2]*q*q;

  MaxAdj = 2.0*tol[0]; MinAdj=0.0;

  for (i=0 ; (i<maxIter[0])& (MaxAdj-MinAdj>tol[0]) ; i++) {

    ColAdj = b[0]/ (tab[0]+tab[1]+tab[2]);
    MinAdj=ColAdj;
    MaxAdj=ColAdj;
    tab[0] *= ColAdj;
    tab[1] *= ColAdj;
    tab[2] *= ColAdj;
 
    ColAdj = b[1]/ (tab[3]+tab[4]+tab[5]);
    tab[3] *= ColAdj;
    tab[4] *= ColAdj;
    tab[5] *= ColAdj;
    if (ColAdj>MaxAdj) {
      MaxAdj = ColAdj;
    }
    else {
      if (ColAdj<MinAdj)
	MinAdj = ColAdj;
    }
    ColAdj = b[2]/ (tab[6]+tab[7]+tab[8]);
    tab[6] *= ColAdj;
    tab[7] *= ColAdj;
    tab[8] *= ColAdj;
    if (ColAdj>MaxAdj) {
      MaxAdj = ColAdj;
    }
    else {
      if (ColAdj<MinAdj)
	MinAdj = ColAdj;
    }
    
    /* do row operations */

    RowAdj = a[0]/(tab[0]+tab[3]+tab[6]);
    tab[0] *= RowAdj;
    tab[3] *= RowAdj;
    tab[6] *= RowAdj;

    if (RowAdj>MaxAdj) {
      MaxAdj = RowAdj;
    }
    else {
      if (RowAdj<MinAdj)
	MinAdj = RowAdj;
    }
       RowAdj = a[1]/(tab[1]+tab[4]+tab[7]);
    tab[1] *= RowAdj;
    tab[4] *= RowAdj;
    tab[7] *= RowAdj;

    if (RowAdj>MaxAdj) {
      MaxAdj = RowAdj;
    }
    else {
      if (RowAdj<MinAdj)
	MinAdj = RowAdj;
    }
       RowAdj = a[2]/(tab[2]+tab[5]+tab[8]);
    tab[2] *= RowAdj;
    tab[5] *= RowAdj;
    tab[8] *= RowAdj;

    if (RowAdj>MaxAdj) {
      MaxAdj = RowAdj;
    }
    else {
      if (RowAdj<MinAdj)
	MinAdj = RowAdj;
    }
    
  }
  
  maxIter[0] = i;
  tol[0] = MaxAdj-MinAdj;
}


void f2wt(double *tab, longint *needProd, double *lambda, 
	  double *rtab, longint *nc1, longint *n,double *tmpTab, 
	  double *maxtol,longint maxIter[1] )
{
  int i, j, k1, k2, k3, tabLen, tabDone, stride1, stride2;
  longint  N, NC, NNC, iter[1];
  double atab[3], btab[3], tol[1];
  
  /* rows move fastest in tmp* and tmpTab */
  
  /* initialize the table of results and the strides used in further
     computations */
  
  /* We suppose that tab is n x nc1 x 3, rtab is 3^nc1 x n */
  /* tmpTab is 3 x 3, lambda is nc1 - 1, odds n, and       */
  /* needProd is n x nc1 - 1                               */
  /*  k1 indexes the (i-1)th dimen    */
  /*  k2 indexes the  i-th dimen      */
  /*  k3 indexes the  lower dims      */


  tabLen = 3;
  for (i = 1; i < *nc1; i++)
    tabLen *= 3;
  N = *n; NC = *nc1; NNC = N*NC; 
  
  for (j = 0; j < N; j++)
    for (i = 0; i < 3; i++)
      rtab[j * tabLen + i] = tab[j + i*NNC];
  /* gradually build up the result */
   
  stride1=1;stride2=3;tabDone=1;
  for (i = 1; i < NC; i++) {
    for (j = 0; j < N; j++) {  /* subject j */
      if (needProd[j+(i-1)*N] == 1) {
	/*got to rake the table*/
	for (k1=0;k1<3;k1++) {
	  atab[k1] = tab[j + (i-1)*N + k1*NNC];
	  btab[k1] = tab[j + i*N + k1*NNC];
	}
	iter[0] = maxIter[0];
	tol[0] = maxtol[0];
	adj33(atab,btab, &lambda[i-1],
	      tmpTab, tol, iter);

	for (k2 = 0; k2<3;k2++) {
	  tmpTab[k2*3] /= atab[0];
	  tmpTab[k2*3 + 1 ] /= atab[1];
	  tmpTab[k2*3 + 2 ] /= atab[2];
	}
      }
      else {
	/* tmpTab is trivial */
	for (k1 = 0; k1<3;k1++)
	  tmpTab[k1*3] = (tmpTab[k1*3+1] = (tmpTab[k1*3+2] = tab[j + i*N + k1*NNC]));
      }
      
      for (k2 = 2; k2 >= 0; k2--) 
	for (k1 = 0; k1 < 3; k1++)  
	  for (k3 = 0; k3 < tabDone; k3++) {
	    rtab[ k3 + k1 * stride1 +  k2*stride2 + j * tabLen] =
	      rtab[k1 * stride1 + k3 + j * tabLen] * tmpTab[k1 + k2 * 3];
	  }		      
    }
    stride1 *=3; stride2*=3;tabDone*=3;
  }
}
