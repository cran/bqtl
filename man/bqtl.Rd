\name{bqtl}
\alias{bqtl}

\title{ Bayesian QTL Model Fitting }
\description{
  Find maximum likelihood estimate(s) or posterior mode(s) for QTL
  model(s). Use Laplace approximation to determine the posterior mass
  associated with the model(s). 
}
\usage{

bqtl(reg.formula, ana.obj, ...)
}
\synopsis{bqtl(reg.formula, ana.obj, scope, expand.specials, ...)}
\arguments{
 \item{reg.formula}{A formula.object like \code{ y ~ add.PVV4 *
     add.H15C12 }.  The names of the independent
   variables on the right hand side of the formula are the names of loci
   or the names of additive and dominance terms associated with loci.  In
   addition, one can use \code{locus} or \code{configs} terms to specify
one or a collection of terms in a shorthand notation.  See \code{
     \link{locus}} for more details.  The left hand side is the name
   of a trait variable stored in the search path, as a column of the
   data frame \code{data}, or \code{y} if the phenotype variable in
   \code{ ana.obj} is used. }  
 \item{ana.obj}{ The result of \code{ \link{make.analysis.obj} }. }
 \item{\dots}{ Arguments to pass to \code{lapadj}, e.g. \code{rparm} and
   \code{return.hess}}
}
\details{
  This function is a wrapper for \code{lapadj}. It does a lot
  of useful packaging through the \code{ \link{configs} } terms. If there
  is no \code{configs} term, then the result is simply the output of
  \code{ \link{lapadj}} with the \code{call} attribute replaced by the
  call to \code{bqtl}  
}
\value{
 The result(s) of calling \code{ \link{lapadj}}.
 If \code{configs} is used in the \code{reg.formula}, then the
 result is a list with one element for each formula. Each element is the
 value returned by \code{lapadj}
 
}
\references{

Tierney L. and Kadane J.B. (1986) Accurate Approximations for Posterior Moments
and Marginal Densities. \emph{JASA}, \bold{81},82--86. 

}
\author{Charles C. Berry \email{cberry@ucsd.edu} }


\seealso{\code{\link{locus}}, \code{\link{configs}}, \code{\link{lapadj}} }

\examples{

data(little.ana.bc )                        # load  BC1 dataset

loglik( bqtl( bc.phenotype ~ 1, little.ana.bc ) ) #null loglikelihood
                                                  #on chr 1 near cM 25
loglik(bqtl(bc.phenotype~locus(chromo=1,cM=25),little.ana.bc))

little.bqtl <-                              # two genes with epistasis
   bqtl(bc.phenotype ~ m.12 * m.24, little.ana.bc)
summary(little.bqtl)

several.epi <-                              # 20 epistatic models
    bqtl( bc.phenotype ~ m.12 * locus(31:50), little.ana.bc)
several.main <-                             # main effects only
    bqtl( bc.phenotype ~ m.12 + locus(31:50), little.ana.bc)

max.loglik <- max( loglik(several.epi) - loglik(several.main) )

round(
     c( Chi.Square=2*max.loglik, df=1, p.value=1-pchisq(2*max.loglik,1))
      ,2)

five.gene <-                                 ## a five gene model
  bqtl( bc.phenotype ~ locus( 12, 32, 44, 22, 76 ), little.ana.bc , return.hess=TRUE )

regr.coef.table <- summary(five.gene)$coefficients

round( regr.coef.table[,"Value"] +  # coefs inside 95\% CI
          qnorm(0.025) * regr.coef.table[,"Std.Err"] \%o\%
            c("Lower CI"=1,"Estimate"=0,"Upper CI"=-1),3)


\testonly{ rm(five.gene,little.ana.bc,little.bqtl,max.loglik,
              regr.coef.table,several.epi,several.main) }
}
\keyword{regression}

