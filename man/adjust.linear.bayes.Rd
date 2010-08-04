\name{adjust.linear.bayes}
\alias{adjust.linear.bayes}

\title{Use Laplace Approximations to improve linear approximations to
  the posterior  }
\description{
 The approximation provided by \code{linear.bayes} can be improved by
 performing Laplace approximations. This function is a development
 version of a wrapper to do that for all of the returned by
 \code{linear.bayes}. 
}
\usage{
adjust.linear.bayes(lbo, ana.obj=lbo$call$ana.obj, \dots)
}

\arguments{
  \item{lbo}{ The object returned by \code{linear.bayes} }
  \item{ana.obj}{ The \code{analysis.object} used to create
    \code{lbo}. This need not be given explicitly, iff the original
    version is in the search path. }
  \item{\dots}{ currently unused }
}

\value{A list of class \code{"adjust.linear.bayes"} containing:
  \item{odds}{A vector, typically of length k giving the odds for models
    of size 1, 2, ..., k under a uniform posterior relative to a model
    with no genes.}
  \item{loc.posterior}{The marginal posterior probabilities by locus}
  \item{coefficients}{The marginal posterior means of the coefficients}
  \item{one.gene.adj}{Results of fits for one gene models}
  \item{n.gene.adj}{Results of fits for modles with more than one gene}
  \item{call}{the call to \code{adjust.linear.bayes}}
  
}
\note{ For large \code{linear.bayes} objects invloving many gene models,
  this can require a very long time to run.}
\author{Charles C. Berry \email{cberry@ucsd.edu} }
\seealso{\code{\link{linear.bayes}}}
\keyword{regression}
