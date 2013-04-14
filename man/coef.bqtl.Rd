\name{coef.bqtl}
\alias{coef.bqtl}
\alias{coef.bqtl.list}
\title{ Extract Coefficients from fitted objects }
\description{
 Return a vector or matrix of coefficients as appropriate
}
\usage{
\method{coef}{bqtl}(object,\dots)
}
\arguments{
  \item{object}{The object returned by \code{bqtl}}
  \item{\dots}{ignored}
}
\value{ A vector (if \code{bqtl} returned a single fit) or matrix (if
 \code{bqtl} returned a list with more than one fit)
}
\author{Charles C. Berry \email{cberry@ucsd.edu} }
\seealso{\code{\link{bqtl}}}
\keyword{methods}