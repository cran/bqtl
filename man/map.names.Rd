\name{map.names}
\alias{map.names}
\alias{map.names.default}
\alias{map.names.bqtl}
\alias{map.names.bqtl.list}
\alias{map.names.analysis.object}
\alias{map.names.map.frame}
\title{Look up names of markers or loci}
\description{This is  a generic helper function with methods  that will
  return the names of markers or loci. 

}
\usage{
map.names(x,...)
map.names.map.frame(x,...)
map.names.bqtl(x, ..., ana.obj)
}

\arguments{
  \item{x}{ An object that has marker names in it. Methods for objects
    of the 
    \code{map.frame}, \code{analysis.object},\code{bqtl}, and
    \code{bqtl.list} class.}
  \item{\dots}{For \code{class(x)=="analysis.object"} or
    \code{class(x)=="map.frame"}, arguments \code{chromo} and \code{cM}
    can be used as in \code{map.index}
  }
  \item{ana.obj}{(optional) analysis object in which to look up the names}
}
\value{
A character vector 
}
\author{Charles C. Berry \email{cberry@ucsd.edu} }

\seealso{ \code{\link{map.index}}, \code{\link{map.location}}}

\examples{
data(little.ana.bc)

map.names(little.ana.bc,chromo=1,cM=24)

map.names(little.ana.bc,chromo=c(1,1),cM=c(40,55))

fit <- bqtl( bc.phenotype ~ locus(23,42) , little.ana.bc )

map.names( fit )

\testonly{rm(fit, little.ana.bc)}
}
\keyword{utilities}
