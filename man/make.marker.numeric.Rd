\name{make.marker.numeric}
\alias{make.marker.numeric}

\title{Translate a marker.frame.object to numeric matrix }
\description{
  Not to be called directly by users. This utility function simply
  returns the coded numeric values corresponding to the allele states.
}
\usage{
make.marker.numeric(marker.frame, level.names=NULL)
}

\arguments{
  \item{marker.frame}{ A data.frame.object consisting of factors or
    character vectors that encode the allele states.}
  \item{level.names}{ A vector of length 6 to translate the levels
    attribute or character codes into allele states that
    \code{make.state.matrix} understands. If necessary, dummy codes are
    used to fill the vector.}
}

\value{
  A matrix, for which column i is
  \code{match(as.character(marker.frame[,i]),level.names)}
}

\author{Charles C. Berry \email{cberry@ucsd.edu} }
