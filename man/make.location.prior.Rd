\name{make.location.prior}
\alias{make.location.prior}

\title{ Provide a default prior}
\description{
 Uses the map distances as a means of assigning a prior for chromosomal
 location.  Basically, this function attempts to assign equal weight
 according to the spacing or markers and 'virtual' markers.
}
\usage{
make.location.prior( x, add.2.end=0, normalize=TRUE )
}

\arguments{
 \item{x}{ \eqn{x = e^{-mgd}}{x = exp(-mgd)} ,where mgd is the map distance in Morgans }
 \item{add.2.end}{ How many Morgans to extend the first and last
   interval on each chromosome }
 \item{normalize}{ If TRUE, let the result sum to 1.0 }
}

\value{
  A vector of \code{length(x)} whose sum is one, if \code{normalize==TRUE}
}

\author{Charles C. Berry, \email{cberry@ucsd.edu}}



\keyword{utilities}