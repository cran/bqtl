/* bc1wt.c
   January 23, 1999 
   C.C. Berry
   Create proto.wt.object in make.weight.bc1
*/

#include "lapadj.h"

/* delta.c - for raking 2 by 2 table */

#include <math.h>

static double delta(double amat[2], double bmat[2], double odds)
{
    double          result, t11, t22, abcd, k1, k2, cc, bb, aa;
    
    t11 = amat[0] * bmat[0];
    t22 = amat[1] * bmat[1];
    
    if (odds == 1.0) {
        result = 0.0;
    }
    else {
        abcd = t11 * t22;
        k1 = t11 + t22;
        k2 = 1.0 - k1;
        cc = abcd * (1.0 - odds);
        bb = k1 + odds * k2;
        aa = (1.0 - odds);
        result = (bb - sqrt(bb * bb - 4.0 * aa * cc)) / (-2.0 * aa);
	
    }
    return result;
}


/* adj22 January 23, 1999
   Rake a 2 by 2 table, then condition on the first margin
*/

static void adj22(double *amat, double *bmat, double *odds, double *result)
{
    
    double          t11, t22, x, dx;
    /* assume amat, bmat are 2, delta is scalar, tmpTab is  2 x 2 */
    
    dx = delta(amat, bmat, odds[0]);
    
    t11 = amat[0] * bmat[0];
    t22 = amat[1] * bmat[1];
    
    x = amat[0];
    if (x == 0.0) {
        result[0] = 0.0;
        result[2] = 0.0;
    }
    else {
    	result[0] = (t11 + dx) / x;
        result[2] = (bmat[1] - dx / x);
    }
    
    x = amat[1];
    if (x == 0.0) {
        result[1] = 0.0;
        result[3] = 0.0;
    }
    else {
        result[3] = (t22 + dx) / x;
        result[1] = (bmat[0] - dx / x);
    }
    
}

void bc1wt(double *tab, longint *needProd, double *lambda,
                      double *rtab, longint *nc1, longint *n, double tmpTab[4], double *odds)
{
  int             i, j, k1, k2, k3, tabLen, tabDone, stride1, stride2;
  longint  N, NC, NNC;
  double          oddz, x, atab[2], btab[2];
  
  /* initialize the table of results and the strides used in further
     computations */
  
  /* We suppose that tab is n x nc1 x 2, rtab is 2^nc1 x n */
  /* tmpTab is 2 x 2, lambda is nc1 - 1, odds n, and       */
  /* needProd is n x nc1 - 1                               */
  /*  k1 indexes the (i-1)th dimen    */
  /*  k2 indexes the  i-th dimen      */
  /*  k3 indexes the  lower dims      */

  tabLen = 2;
  for (i = 1; i < *nc1; i++)
    tabLen *= 2;
  N = *n; NC = *nc1; NNC = N*NC; 

  for (j = 0; j < N; j++)
    for (i = 0; i < 2; i++)
      rtab[j * tabLen + i] = tab[j + i*NNC];
  /* gradually build up the result */
   
  stride1=1;stride2=2;tabDone=1;
  for (i = 1; i < NC; i++) {
    x = lambda[i - 1];
    x = (1.0 + x) / (1.0 - x);
    oddz = x * x;
    for (j = 0; j < N; j++) {  /* subject j */
      odds[j] = (needProd[j+(i-1)*N] == 1) ? oddz : 1.0;
      
      atab[0] = tab[j+(i-1)*N];
      atab[1] = tab[j+(i-1)*N + NNC];
      btab[0] = tab[j+(i)*N];
      btab[1] = tab[j+(i)*N + NNC];
      adj22(atab, btab, &odds[j], tmpTab);
      for (k2 = 1; k2 >= 0; k2--) 
	for (k1 = 0; k1 < 2; k1++)  
	  for (k3 = 0; k3 < tabDone; k3++) {
	    rtab[ k3 + k1 * stride1 +  k2*stride2 + j * tabLen] =
	      rtab[k1 * stride1 + k3 + j * tabLen] * tmpTab[k1 + k2 * 2];
	  }		      
    }
    stride1 *=2; stride2*=2;tabDone*=2;
  }
}
