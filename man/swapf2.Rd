\name{swapf2}
\alias{swapf2}
\title{
  Sample F2 loci  via approximate posterior 
}
\description{
  An MCMC sampler for loci using precomputed dispersion matrices,
  various priors, and a pre-selected set of variables. For use with F2
  intercross design.

Using precomputed dispersion matrices, various priors, and a pre-selected set of variables, one locus is removed, all other loci are examined to determine the effect of each on the posterior. One locus is sampled. The process is repeated until each locus has been removed and a new one sampled in its place (possibly the same one that was removed is sampled). 
}
\usage{
swapf2(varcov, invars, rparm, nreps, ana.obj, locs = <<see below>>, 
    locs.prior = <<see below>>, combo.prior = <<see below>>, tol = 1e-10)
} 

\arguments{
  \item{varcov}{ The result of \code{ make.varcov }. The columns of
    \code{varcov$var.x} must alternate 'additive' and 'dominance' terms.}
  \item{rparm}{
    The 'ridge' parameters for the independent variables - larger values
    imply more shrinkage or a more concentrated prior for the regresion
    coefficients. }
  \item{nreps}{How many cycles of MCMC to perform}
  \item{ana.obj}{A object produced by \code{\link{make.analysis.obj}}}
  \item{invars}{A vector of variable indexes. This determines which
    variables to start in the model. If both additive and 
    dominance terms are to be used, they should occupy adjacent
    locations in \code{invars}. The variable(s) associated with the
    first locus is (are) immediately removed, serving only as
    placeholder(s). If there are k loci associated with the variables,
    then all subsequent models have k loci, although the nuimber of
    variables may vary according to the selection of one or both of the
    'additive' or 'dominance' terms.  } 
  \item{locs}{
    The pairs of columns of \code{varcov$var.x} to use. The default
    uses all of them.}
  \item{locs.prior}{ Vector whose elements are 
    the prior masses to associate with each locus. Typically, these sum
    to one, but sometimes they might each be set to one (as in computing
    lod scores). The default value sets them all to 1.0.
     }
  \item{combo.prior}{
    The prior probability for each term or combination of terms for the
    phenotypic effect at a locus. Typically, there will be three of
    these - one for the 'additive' term (linear in number of alleles
    from one parent strain), the 'dominance' term (quadratic in allele
    number), or both terms. The default sets them all to 1/3.
  }

  \item{tol}{Used in forming QR decomposition. Let it be.}
}

\value{
  A list with components:
  \item{configs}{ A 2k by k by nreps array of indexes of variables sampled in
    each of the nreps iterations. Models using less than 2k variables
    \code{configs[,i,j]} will contain one or more zeroes in the last
    position(s) }
  \item{posteriors}{ A vector of length \code{k*nreps} with the
    posteriors of the models sampled.} 
  \item{coefs}{A 2k by k by nreps matrix of the regression
    coefficients. Models using less than 2k variables 
    \code{configs[,i,j]} will contain one or more zeroes in the last
    position(s) }
  \item{call}{The call to \code{swapf2}}
  \item{cond}{The \code{k*nreps} posterior probabilities of the k-1 gene
    models.} 
  \item{marg}{The \code{k*nreps} marginal posteriors for all k gene
    models that could be formed using the current k-1 gene model)}
  \item{alt.marginal}{A vector with \code{length(locs)} elements. At
    each step, the posterior associated with  each candidate locus is
    added to an element of this vector. After all steps are finished,
    the result is normalized to sum to one. This turns out to be an
    exceedingly stable estimate of the relative marginal posterior. }
  \item{alt.coef}{A vector with \code{2*length(locs)} elements. At
    each step, the product of each posterior times the coefficient
    associated with a candidate variable is 
    added to an element of this vector. After all steps are finished,
    the result is normalized by the total marginal posterior. This turns
    out to be an rather stable estimate of the marginal (over all
    models) posterior mean of the regression coefficients. }
}

\details{
  A call to \code{swapf2} is used to obtain the results. This function
  is really just a wrapper. 
}
\seealso{
  \code{\link{swapbc1}}
}

\references{Berry C.C. (1998) Computationally Efficient Bayesian QTL Mapping in
  Experimental Crosses. \emph{ASA Proceedings of the Biometrics
    Section,} 164-169.
}
\author{Charles C. Berry \email{cberry@ucsd.edu} }

\keyword{utilities}