\name{little.mf.5}
\alias{little.mf.5}
\non_function{}
\title{ Package of Simulated Marker Map Information }
\description{
The \code{little.mf.5} data frame has 114 rows and 9 columns
consisting of \code{little.map.frame} plus 64 'virtual' marker loci
}
\format{
This data frame contains the following columns:
\describe{
\item{marker.name}{
The marker names taken from \code{little.map.frame} and those created to
fill virtual markers in between actual markers.
}
\item{dx}{
a vector of locations
}
\item{prior}{
weights to be used in sampling and Bayesian computations
}
\item{pos.type}{
a factor with levels
\code{left} 
\code{right} 
\code{center} 
}
\item{is.marker}{
\code{TRUE} for the 50 markers, \code{FALSE} for the 'virtual' markers
}
\item{pos.plot}{
a vector of plotting positions
}
\item{lambda}{
transformed recombination fractions
}
\item{locus}{
an abbreviated locus name
}
\item{chr.num}{
the chromosome number 1, 2, 3, 4, or 5.
}
}
}

\keyword{datasets}
