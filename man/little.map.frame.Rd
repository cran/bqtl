\name{little.map.frame}
\alias{little.map.frame}
\non_function{}
\title{ Package of Simulated Marker Map Information }
\description{
The \code{little.map.frame} data frame has 50 rows and 9 columns that
describe the marker map of \code{\link{little.map.dx}} in the format
produced by \code{\link{make.map.frame}. \code{little.map.dx} has the
  minimal data needed to construct this.}
}
\format{
This data frame contains the following columns:
\describe{
\item{marker.name}{
a factor with levels
\code{m.1} 
\code{m.2} 
\dots
\code{m.50} 
}
\item{cM}{
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
always \code{TRUE} for these data
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
