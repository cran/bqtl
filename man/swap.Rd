\name{swap}
\alias{swap}

\title{ MCMC sampling of multigene models}
\description{
Given a k-gene model as a starting point, one gene is deleted and
another is sampled in its place. This is done using an approximation to
the posterior. Then another gene is deleted and another sampled,\dots

}
\usage{
swap(varcov, invars, rparm, nreps, ana.obj, ...)
}

\arguments{
  \item{varcov}{The result of \code{\link{make.varcov}} }
  \item{invars}{ Vector of
    numerical indexes of \code{ana.obj$reg.names} telling which
    variables to start in the model. The first of these is immediately
    removed, so it is merely a placeholder.  The number of genes in the
    model is therefore \code{k <- length(invars)} (except when
    \code{ana.obj$method=="F2"} when it is \code{k <-
      length(unique(col(ana.obj$reg.names)[invars]))})  }
  \item{rparm}{Scalar or vector with \code{nrow(varcov$var.x)} elements;
    the 'ridge' parameters for the independent variables - larger values
    imply more shrinkage or a more concentrated prior for the regresion
    coefficients.  }
  \item{nreps}{ How many cycles (of \code{k} samples each) to perform.  }
  \item{ana.obj}{ An \code{analysis.object} --- see \link{make.analysis.obj}}
  \item{\dots}{Additional arguments override the default choices of
    candidate loci (\code{locs}), prior for locus (\code{locs.prior}),
    or method specified by \code{ana.obj}. Also, the default prior for
    model (\code{combo.prior} when \code{ana.obj$method=="F2" can be
      overridden. See \link{swapbc1} and \link{swapf2} for details. } 
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
cycle is repeated \code{nreps} times.  
}
\value{
  A list with components:
  \item{config}{ A k by k by nreps array (or, for
    \code{ana.obj$method=="F2"}, a 2k by k by nreps array) of the
    locations (variables) sampled in each iteration.} 
  \item{posteriors}{ A vector of length \code{k*nreps} with the
    posteriors of the models.} 
  \item{coefs}{A k by k matrix of the regression coefficients(or, for
    \code{ana.obj$method=="F2"}, a 2k by nreps matrix).}
  \item{call}{The call to \code{swap}}
  \item{cond}{The \code{k*nreps} posterior probabilities of the k-1 gene
    models.}
  \item{marg}{The \code{k*nreps} marginal posteriors for all k gene
    models that could be formed  using  the current k-1 gene model}
  \item{alt.marginal}{A vector with \code{length(locs)} (or \code{2*length(locs)}) elements.  At
    each step, the posterior associated with  each candidate locus is
    added to an element of this vector.  After all steps are finished,
    the result is normalized to sum to one.  This turns out to be a
    stable estimate of the marginal posterior. }
  \item{alt.coef}{A vector with \code{length(locs)} (or
    \code{2*length(locs)}) elements.  At each step, the product of each
    posterior times the coefficient(s) associated with a candidate locus
    is added to an element of this 
    vector.  After all steps are finished, the result is normalized by
    the total marginal posterior.  This turns out to be a stable
    estimate of the marginal (over all models) posterior mean of the
    regression coefficients. }
}
\references{Berry C.C. (1998) Computationally Efficient Bayesian QTL Mapping in
  Experimental Crosses. \emph{ASA Proceedings of the Biometrics
    Section,} 164-169. 
}
\author{Charles C. Berry \email{cberry@ucsd.edu} }
\examples{
data( little.ana.bc )
little.vc <- varcov( bc.phenotype~locus(all), little.ana.bc)
little.4 <- swap( little.vc, c(1,15,55,75), rparm=1, 50, little.ana.bc )
little.4.smry <- summary( little.4 )
print(c("Bayes Factor (3 vs 4)"=little.4.smry$ratio$mean))
par(mfrow=c(3,2))
plot( little.ana.bc, little.4.smry$loc.posterior, type="h",
 ylab="E(genes)" )
rm(little.4,little.vc,little.ana.bc)
}
\keyword{regression}