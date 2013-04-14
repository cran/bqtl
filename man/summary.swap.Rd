\name{summary.swap}
\alias{summary.swap}

\title{Summarize Gibbs samples for a k-gene model }
\description{
Calculate marginal posteriors for location, posterior means for
coefficients, and the Bayes Factor for k vs k-1 genes  
}
\usage{
\method{summary}{swap}(object, method=NULL, ncoef=length(object$alt.coef), nloc=object$nloc,\dots)
}

\arguments{
  \item{object}{ The result of \code{\link{swap}} }
  \item{method}{Optional. One of the supported methods, see \code{\link{make.analysis.obj}} }
  \item{ncoef}{Optional. The number of coefficients in the class of
    models. Typically, \code{2*nloc} for \code{method=="F2"} and
    \code{nloc} for all other methods  }
  \item{nloc}{Optional. The number of loci in the sample space. }
  \item{\dots}{ignored}
}
\value{
  A list with components:
  \item{loc.posterior }{A vector of (marginal) posterior odds for each
    locus compared to a no gene model}
  \item{coefs}{Posterior means of coefficients.}
  \item{ratio}{A list with components \code{mean}, an estimate of the
  Bayes Factor for k versus k-1 gene models, and \code{var}, an estimate
  of its variance.}
}
\author{Charles C. Berry \email{cberry@ucsd.edu} }

\keyword{methods}
