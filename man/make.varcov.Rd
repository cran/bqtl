\name{make.varcov}
\alias{make.varcov}
\title{ Create moment matrices}
\description{
  Create a moment matrix of the marker variables and of the regressors by
  the phenotype variable.  For use in regression modelling on the markers.
}
\usage{
make.varcov(regressor.matrix, y, subset=is.finite(y), casewt=NULL)
}

\arguments{
 \item{regressor.matrix}{The object produced by \code{\link{make.regressor.matrix}} }
 \item{y}{ A vector of phenotype information with the same number of
   elements as there are rows in \code{regressor.matrix} }
 \item{subset}{ Logical vector with the same number of
   elements as there are rows in \code{regressor.matrix} to indicate
   which rows to keep.  }
 \item{casewt}{Optional vector of case weights. }
}

\value{
 A list with components
  \item{var.x }{Moment matrix of the marker regressor variables}
  \item{cov.xy }{Moment matrix of the marker regressor variables versus
    the phenotype variable}
  \item{var.y}{The Second central moment of the phenotype variable}
  \item{df}{ \code{ sum(subset==TRUE) - 1}
}}
\author{Charles C. Berry \email{cberry@ucsd.edu} }
\note{ It is generally NOT a good idea to do regressions on
  ill-conditioned designs using the moment matrices like this.  The
  excuse for doing so here is twofold.  First, calculations using this
  method are used to perform importance sampling, so minor numerical
  inaccuracies in computing the probabilites used in sampling get
  straightened out by the importance weights.  Second, it will typically
  be the case that a prior is set on the regression coefficients and
  this results in a positive constant (aka a 'ridge' parameter) being
  added to diagonal of \code{varcov$var.x} and this reduces the
  ill-conditioning.  Of course the rational for using the method is to
  speed the sampling, and it is very effective at doing so.} 

\keyword{utilities}