\name{make.loc.right}
\alias{make.loc.right}

\title{ Keep track of fully informative markers or states}
\description{
  Helps speed computations in multigene models by allowing a quick
  assessment of whether two loci are independent given the marker
  information for the individual.
}
\usage{
make.loc.right(marker.frame, marker.distances)
}

\arguments{
 \item{marker.frame}{ A \code{marker.frame.object}  } 
 \item{marker.distances}{Actually a misnomer, this is a vector with a
   zero in the last position of each chromosome.  }
}

\value{
 A matrix of the same dimension as \code{marker.frame} whose elements
 index the column on the next (right) fully informative marker.
}

\author{Charles C. Berry \email{cberry@ucsd.edu} }


\examples{

}

\keyword{manip}