/* conCk.c */
/* Wed Feb  3 15:35:39 PST 1999 */
/* C.C.Berry */

#include "lapadj.h"

int conCk(double *oldllk, double *newllk, double *paradj, longint *NP,
	  double *tol)
{
  int i;
  double x;
  if ( fabs(*newllk-*oldllk) > *tol )
    return 0;
  else { 
    x = 0.0; 
    for (i=0;i<*NP;i++)
      x += fabs(paradj[i]);
    if (x > *tol)
      return 0;
  }
  return 1;
}
  
