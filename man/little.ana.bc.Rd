\name{little.ana.bc}
\alias{little.ana.bc}
\title{ A simulated dataset }
\description{
A simulation of a BC1 cross of 150 organisms with a genome of around 500
cM consisting of 5 chromosomes. The format is that created by
\code{\link{make.analysis.obj}}

This dataset is built up from several others. The basic data are:

\itemize{
\item{little.bc.pheno}{A vector of phenotype data}

\item{little.bc.markers}{A \code{map.frame} of marker data} and

\item{little.dx}{A data frame with 50 rows and 2 columns that specify
  the map locations of a simulated set of markers }
}
These are used to construct

\itemize{
\item{little.mf.5}{A \code{map.frame} with 'pseudo-markers' at least
  every 5 cM made from
  
  \code{little.mf.5 <-  make.map.frame(little.map.frame,
    nint=marker.fill( little.map.frame, reso=5, TRUE ))}
}
}
Then phenotype, covariate, and marker data are combined
with \code{little.mf.5}
\itemize{
\item{little.bc.pheno}{ A \code{data.frame} with the
  variable \code{bc.phenotype}}
\item{little.bc.markers}{A \code{data.frame} with marker state information}
}
}
\seealso{The examples in \code{\link{make.analysis.obj}}}


\keyword{datasets}