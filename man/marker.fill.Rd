\name{marker.fill}
\alias{marker.fill}
\title{
Map Positions Between Markers
}
\description{


Given a set of markers, one wants to create a finer map at a given
resolution. marker fill takes a a collection of marker distances and a
desired resolution and finds positions that are intermediate and at
that resolution.


}
\usage{
marker.fill(map.frame, reso, return.nint = FALSE)


}
\arguments{
\item{map.frame}{

A map.frame.object.

}
\item{reso}{
The desired interval between loci in the same metric as \code{map.frame$cM}
}
\item{return.nint}{Whether to output a vector of number of intervals to
  produce in each existing interlocus interval}
}
\value{
If \code{return.nint} is \code{TRUE}, a vector of integers is
returned. It indicates how many intervals to place between this marker
and the next to achive the desired minimum distance.

If \code{return.nint} is \code{FALSE}, a vector of distances is
returned. The names attribute has suffixes 
added to indicate positions filled to the 'right' of existing
markers. Thus if markers 'mark.01' and 'mark.02' are in succession at
a distance of 3 and reso==1, then the value associated with 'mark.01'
(which was 3) becomes 1, a value of 1 is associated with new loci
called 'mark.01.1' and 'mark.01.2' in created with values of 1
each. The returned vector is ordered by chromosome, then marker or
filled locus.
}

\seealso{
  \code{\link{make.map.frame}}
}

\examples{


data( little.map.frame, lib.loc=.lib.loc )
little.nint <- marker.fill( little.map.frame, reso=5, TRUE )
cbind(nint=little.nint,cM=little.map.frame$cM)[1:10,]
rm( little.map.frame, little.nint )
}
% Converted by Sd2Rd version 1.10.
\keyword{utilities}