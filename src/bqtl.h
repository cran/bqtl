#include <S.h>



#ifndef USING_R 
/* for use with Splus */
#define NTFINITE(x) (is_na((&x),DOUBLE))? 1: is_inf((&x),DOUBLE)
#define longint long
#include <nonlin.h> 
#include <math.h>
#define F77_SYMBOL(x) F77_CALL(x)
extern double F77_NAME(ddot)();/* in lieu of Blas.h */
#else
/* for use with R */

#include <R.h>
#define longint int
#define NTFINITE(x) (R_FINITE(x))?0:1
#endif


/* funky print macros for debugging */

#define PR(prtext) printf(prtext);
#define PRD(pr_d) printf(#pr_d " = %e \n", (pr_d))
#define PRDM(pr_d,nel) \
{long iii; for (iii=0;iii<(nel);iii++) \
    printf(#pr_d " = %e \n", *(pr_d+iii));}
#define PRDMAT(pr_d,nr,nc) \
{long nrow,ncol; for (nrow=0;nrow<nr;nrow++) \
 {for (ncol=0;ncol<nc;ncol++) \
printf("%3.1e ",*(pr_d+nrow+ncol*nr));printf("\n");};};
#define PRLN printf("%d %s\n",__LINE__,__FILE__);fflush(stdout)
#define PRI(pr_int) printf("%d %s integer = %ld \n",__LINE__,__FILE__, pr_int)
#define CKLOG(x)  if ((x) <= 0.0 ) PRLN 

/* end funky macros */

/* for debugging array overrruns
   UNCOMMENT this group of lines 

   #define ZERO_FLT(x,y) for (i=0;i<(y);i++) (x)[i]=0.0
   #define ALLOC_DBL(x,y) (x) = Calloc( (y) + 1, double); \
   ZERO_FLT((x),(y));(x)[(y)] = 9.87654321 
   #define ALLOC_DBLPT(x,y) (x) = Calloc( (y+1), double *)
   #define ZERO_INT(x,y) for (i=0;i<(y);i++) (x)[i]=0
   #define ALLOC_LONG(x,y) (x)  = Calloc( (y) + 1 , longint); \
   ZERO_INT((x),(y)); (x)[(y)] = 987654
   #define DEALLOC_DBL(x,y) PRD( (x)[(y)] ); Free( (x) ) 
   #define DEALLOC_LONG(x,y) printf("%d\n",(x)[(y)]); Free( (x) ) 
   #define DEALLOC_DBLPT(x,y) Free( (x) )
*/



#define ALLOC_DBL(x,y) (x) = Calloc( (y), double)
#define ALLOC_DBLPT(x,y) (x) = Calloc( (y), double *)
#define ALLOC_LONG(x,y) (x)  = Calloc( (y), longint)
#define DEALLOC_LONG(x,y) Free( (x) )
#define DEALLOC_DBL(x,y) Free( (x) )
#define DEALLOC_DBLPT(x,y) Free( (x) )


#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define CPY(oldx, newx, len) memcpy(newx,oldx,(size_t) *(len) * sizeof(*(oldx)))
#define QR_TOL (1.0e-6)
#define TRUE 1
#define FALSE 0
#define TOL 10E-10 ;
