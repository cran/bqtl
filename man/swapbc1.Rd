\name{swapbc1}
\alias{swapbc1}
\title{
  Sample BC1 or Recombinant Inbred loci  via approximate posterior.
}
\description{
  An MCMC sampler for loci using precomputed dispersion matrices,
  various priors, and a pre-selected set of variables.  For use with BC1
  (backcross) designs and recombinant inbred lines.
}
\usage{
swapbc1(varcov, invars, rparm, nreps, locs=seq(ncol(var.x)),
locs.prior=rep(1, ncol(var.x)),tol=1e-10 )

}
\arguments{
  \item{varcov}{ The result of \code{ \link{make.varcov}}}
  \item{rparm}{ Scalar or vector with \code{nrow(varcov$var.x)} elements;
    the 'ridge' parameters for the independent variables - larger values
    imply more shrinkage or a more concentrated prior for the regresion
    coefficients. }
  \item{invars}{
    Which variables to start in the model. The first of these is
    immediately removed, so it is merely a placeholder.  The number of
    genes in the model is therefore \code{k <- length(invars)}}
  \item{locs}{
    The columns of \code{varcov\$var.x} to use.  The default uses all of
    them.} 
  \item{locs.prior}{
    The prior mass to associate with each variable.  Typically, these sum
    to one, but sometimes they might each be set to one (as in computing
    lod scores). }
  \item{tol}{Used in forming QR decomposition.  Let it be.}
}
\details{
  An MCMC sampler for loci using the object of \code{ make.varcov } is
  executed.  This sampler uses the exact posterior probability under the
  assumed correctness of the regression model using expected genotypes
  given marker values.  This amounts to linearizing the likelihood with
  respect to the  (possibly unknown) locus states.  For models where the
  loci are fully informative markers this is the true posterior.

  The chain is implemented as follows: given a set of regressor
variables to start, one variable is removed, all regressor
variables not in the model are examined to determine the effect of each
on the posterior.  One variable is sampled.  The process is repeated until
each variable has been removed and a new one sampled in its place
(possibly the same variable that was removed is sampled).  And this whole
cycle is repeated \code{nreps] times.  
}
\value{
  A list with components:
  \item{config}{ A k by k by nreps array of the locations sampled in
    each iteration. }
  \item{posteriors}{ A vector of length \code{k*nreps} with the
    posteriors of the models.} 
  \item{coefs}{A k by k matrix of the regression coefficients.}
  \item{call}{The call to \code{swapbc1}}
  \item{cond}{The \code{k*nreps} posterior probabilities of the k-1 gene
    models.}
  \item{marg}{The \code{k*nreps} marginal posteriors for all k gene
    models that could be formed  using  the current k-1 gene model}
  \item{alt.marginal}{A vector with \code{length(locs)} elements.  At
    each step, the posterior associated with  each candidate locus is
    added to an element of this vector.  After all steps are finished,
    the result is normalized to sum to one.  This turns out to be an
    exceedingly stable estimate of the marginal posterior. }
  \item{alt.coef}{A vector with \code{length(locs)} elements.  At
    each step, the product of each posterior times the coefficient
    associated with a candidate locus is 
    added to an element of this vector.  After all steps are finished,
    the result is normalized by the total marginal posterior.  This turns
    out to be an exceedingly stable estimate of the marginal (over all
    models) posterior mean of the regression coefficients. }
}
\seealso{
  \code{\link{swapf2}}
}
\references{Berry C.C. (1998) Computationally Efficient Bayesian QTL Mapping in
  Experimental Crosses. \emph{ASA Proceedings of the Biometrics
    Section,} 164-169.
}
\author{Charles C. Berry \email{cberry@ucsd.edu} }
\keywords{utilities}