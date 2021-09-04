#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#define longint int


void lapadj(longint *crsType, longint *nparm,
            double *y, 
            double *x,
            double *z,
	    double *amnt, double *tab, longint *needProd, double *lambda, double
	    *coefs, double *casewt, double *sigma, double *hkExact, double *hkApprox,
	    double *llk, double *hess, double *postApprox, longint *iter, double *tol, 
	    longint *nem);

void swapbc1( longint *nreps, longint *nstep,
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
	      double *altmarg, double *altcoef);

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
	     double *altmarg, double *altcoef);

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
	     double *tol);

void upbqtl(longint *crsType, longint *nparm,
            double *y, 
            double *x,
            double *z,
	    double *amnt, double *tab, longint *needProd, double *lambda, double
	    *coefs, double *casewt, double *sigma, double *hkExact, double *hkApprox,
	    double *llk, double *hess, double *postApprox, longint *iter, double *tol, 
	    longint *nem,
	    /* the above args are for the usual call to lapadj */
	    
	    longint *loc_mat, longint *n_alt, longint *loc_order, double *result, 
	    double *orig_x, longint *loc_right, double *map_lambda,
	    double *state_matrix, longint *n_state_loc);


#define C_DEF(name, n)  {#name, (DL_FUNC) &name, n}


static const R_CMethodDef CEntries[]  = {
	  C_DEF(lapadj, 20  ),
	  C_DEF(swapbc1, 44  ),
	  C_DEF(swapf2, 44  ),
	  C_DEF(twohkbc1, 32  ),
	  C_DEF(twohkf2, 32  ),
	  C_DEF(upbqtl, 29  ),
	  {NULL, NULL, 0}
};



void R_init_bqtl(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
