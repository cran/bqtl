\name{map.location}
\alias{map.loc}
\alias{map.location}
\alias{map.location.default}
\alias{map.location.map.frame}
\alias{map.location.analysis.object}
\alias{map.location.bqtl}
\alias{map.location.bqtl.list}
\title{ Report map location  }
\description{
 Report the chromosome number and location of loci in a genetic map.
}
\usage{
map.location(x, y, \dots)
map.loc(x, chromo=3)
map.loc(x, chromo=c(4,6), cM=c(22,13) )
map.loc(x, map.names="m.12")
map.loc(x)
}

\arguments{
  \item{x}{A object of class \code{map.frame},\code{analysis.object},
\code{bqtl}, or \code{bqtl.list}}
  \item{\dots}{Other arguments usage depend on the class of \code{x} }
  \item{y}{A vector of row numbers or \code{map.names} specifying which
    subset of the \code{map.frame} of \code{x} is to be returned}
  \item{chromo}{A vector of chromosome numbers}
  \item{cM }{(Optional) map distance vector. If the same chromosome
    number is used twice in \code{chromo}, \code{cM} must be
    ordered. If \code{cM} is omitted, all loci on each chromosome listed
    in \code{chromo} will be included.)

  \item{map.names}{A vector of \code{map.names}}
}
\details{ It is often helpful to refer to genentic loci by their
locations. The methods of \code{map.location} (alias \code{map.loc})
will extract the row index, chromosome number and location, and the name
for specified loci. For direct lookups of the loci in a \code{map.frame}
or \code{analysis.object}, one must specify \code{y} or \code{chromo} or
\code{map.names}. When \code{class(x)=="bqtl"} \code{map.location}s of
terms used in a call to \code{bqtl} are returned. When \code{cM} is
used, an attempt will be made to match the location; if the match fails,
the nearest locus will be used. When there are two elements in
\code{chromo} and two in \code{cM}, all the map locations in between the
matching loci will be returned.

}
\value{
An object of class \code{map.location} which inherits from
\code{map.frame}. It has columns:

  \item{chr.num }{The chromosome number}
  \item{cM }{The location in centiMorgans on that chromosome.}
  \item{marker.name }{The name by which that marker is known}
  
  \item{attr(,"row.names")}{An index of the locations}

}
\author{Charles C. Berry \email{cberry@ucsd.edu}}

\seealso{\code{\link{make.map.frame}}}

\examples{

data(little.ana.bc)

map.loc(little.ana.bc, c(1,15,45))
map.loc(little.ana.bc,chromo=3,cM=22)
map.loc(little.ana.bc,"m.12")
rm(little.ana.bc)
}
\keyword{\keyword{bqtl}