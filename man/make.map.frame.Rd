\name{make.map.frame}
\alias{make.map.frame}

\title{Create marker map specifications}
\description{

A map.frame.object describes a marker map and additional loci that may
be used in a QTL study.  Each row pertains to one locus.  Names of
markers, abbreviated names, distances, and other necessary and useful
information are bundled.

}
\usage{
make.map.frame(dx,chr.num = NULL, prior = make.location.prior(lambda), 
    morgan = 100, nint = NULL, reso = NULL)

}

\arguments{

  \item{dx}{ An object of class \code{"map.frame"} or class
\code{"data.frame"} or a vector or a
    data.frame with a column named \code{cM}, \code{M}, or \code{dx} or
    whose first column gives location  
    on each chromosome in  centiMorgans (from start of
    chromosome or Morgans if \code{M} was the column name).  It is best
    if \code{names(dx)} (for vector arguments) 
    or \code{row.names(dx)} (for data.frame arguments) give names of
    markers for later reference, but this isn't really necessary.  } 
  \item{chr.num}{ (Optional) Vector of  chromosome  numbers}
  \item{prior}{ (Optional) Vector of Prior probabilities for the loci }
  \item{morgan}{ 100 if centiMorgans, 1 if Morgans }
  \item{nint}{ (Optional) Vector of one plus number of 'virtual' markers to be
    inserted after each locus }
  \item{reso}{Maximum distance between loci. If necessary fill in with 'pseudo-markers'}
}
\details{

  The QTL analysis depends on information about the marker map and on
  specifications of the loci to be studied.  The 'map.frame' contains
  this information. 

}
\value{
  A data frame with components:

  \item{marker.name}{The full text identifier for this marker,
    e.g.  "HH.360L.Col" is a marker on chromosome 1 of arabidopsis
    thaliana, and names like this can be used for reference
    purposes.  'Virtual' markers have a suffix appended to the name of the
    previous marker.}
  \item{cM }{Location on the chromosome.  If this is a marker of a locus
    that was input via \code{dx}, then it is just the value of \code{dx}.}
  \item{pos.type}{"left" if it is the first locus on this
    chromosome,"right" if it is last, or "center" otherwise. }
  \item{is.marker}{TRUE if this was actually a marker, FALSE if it is a
    'virtual' marker}
 \item{pos.plot}{Plotting position for this locus.  Typically the same as
   \code{dx}.}
 \item{lambda}{Twice the recombination fraction minus one.}
 \item{locus}{An abbreviation for the locus of the form
   "C.<\code{chr.num}>.<\code{cM}>"}
 \item{chr.num}{The chromosome number.}

}

\author{ Charles C. Berry \email{cberry@ucsd.edu}
	
 }

 \note{ The idea in having all of this bundled together is to make it
   easier for plot and summary methods to be implemented and to allow
   convenient references in formula based methods.}

\examples{

data( little.map.dx )
little.map.frame <- make.map.frame( little.map.dx )
plot( little.map.frame ) # there is a plot method
# add 'virtual' markers to map
little.mf.5 <-   make.map.frame(little.map.frame,reso=5)
print(little.mf.5[1:10,],digits=1) # show a few rows
plot( little.mf.5 ) # notice the 'virtual' markers added
\testonly{ rm( little.map.dx,little.map.frame ,little.mf.5) }
}

\keyword{manip}