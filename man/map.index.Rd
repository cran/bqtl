\name{map.index}
\alias{map.index}
\alias{map.index.default}
\alias{map.index.analysis.object}
\title{Look up numerical index(es) of map locations}
\description{
One way to index a locus (loci) in a genetic map is by the numerical index of
its row (their rows). \code{map.index} performs a lookup in a specific
\code{map.frame} given one (or two) chromosome number(s) and one (or
two) map distance(s).
}
\usage{
map.index(x, ... )
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x }{A \code{map.frame} or \code{analysis.object}}
  \item{...}{ For methods that look up a location in a \code{map.frame}
    the following named arguments may be used:
  \code{chromo }{A chromosome number or 2 ordered numbers }
  \code{cM }{(Optional) map distance or two. If the same chromosome
    number is used twice in \code{chromo}, \code{cM} must be
    ordered. If \code{cM} is omitted, all loci on the chromosome will
  be included.}
}}
\details{
  It is often convenient to refer to genetic loci or regions by the
  numerical index(es) in a \code{map.frame}. \code{map.index} allows
  lookups according to the approximate map location.
}
\value{
A numerical vector of one or more row numbers. If only \code{chromo} is
specified, all row numbers on the specified chromosome are returned. If
\code{chromo} has two elements, then all row numbers on those
chromosomes with numbers in \code{range(chromo)} will be returned. If
one of each of \code{chromo} and \code{cM} are specified, then the row
number of the closest locus will be returned. For two of each, row
numbers in the range of the closest matches will be returned.
}
\author{Charles C. Berry \email{cberry@ucsd.edu}}
\seealso{ \code{\link{make.map.frame}} for a description of how map
  information is organized.}
\examples{
data(little.ana.bc)
map.index(little.ana.bc,chromo=1,cM=25) # locus nearest 1,25
index.chr.1 <- map.index(little.ana.bc,chromo=1)
fit.on.1 <- bqtl(bc.phenotype~locus(index.chr.1),little.ana.bc)
summary( loglik( fit.on.1 ) )

}

\keyword{methods}