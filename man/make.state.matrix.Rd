\name{make.state.matrix}
\alias{make.state.matrix}
\title{
  Create state.matrix.object
}
\description{
  Create a state.matrix.object to be used encode marker information in a
  form in which it can be used in subsequent calculations. 
}
\usage{
make.state.matrix(marker.frame, marker.distances, method="F2")
}
\arguments{
  \item{marker.frame}{
    Actually, this is a misnomer. This is NOT
    a\code{marker.frame.object}. Rather it is obtained by by call like
    \code{make.marker.numeric(marker.frame.object)} (see \code{
      \link{make.marker.numeric} }) and it is coerced to 
    a matrix. It  encodes marker allele states. One column is used for each 
    marker or pseudo-marker (basically a placeholder with all missing
    values). The entries are in 1:6, if NA's are present, they are
    recoded to 6.  The columns are arranged in
    linkage groups with presumed order reflected in the actual order of the
    columns.  
  }
  \item{marker.distances}{
    Distances between the markers in the 'lambda' metric. -log(lambda)/2
    is the Haldance map distance. Linkage groups are separated by values
    of 0.0. 
  }
  \item{method}{
  method = "F2" is the default, and   "BC1", "RI.self", and "RI.sib" are
  other options. The assumed setup is as follows (strains are A and a):

  \tabular{lccc}{
    marker state \tab  F2.code \tab  BC.code\tab RI.code \cr
    \code{"AA"}  \tab    1    \tab     1   \tab    1\cr
    \code{"Aa"} \tab     2  \tab      2 \tab \cr
    \code{"aa"}  \tab    3   \tab        \tab      2 \cr
    \code{"A-"} (not aa)    \tab  4      \tab   \tab     \cr
    \code{"a-"} (not AA) \tab    5    \tab     \tab        \cr
    \code{"--"} (unknown) \tab    6 \tab        6   \tab    6\cr
  }
}}
\value{
  n by k by q array. q is  3 for method="F2" and 2 for others methods.  Each
  element encodes the probability of the allele state conditional on the
  marker states. 
}
\references{
 Lander E.S. and Green P. (1987) Construction of multilocus genetic linkage maps
 in humans. \emph{Proceedings of the National Academy of Sciences of the
 United States of America}, \bold{84}(8), 2363--7.
 
  Jiang C. and Zeng Z-B.  (1997) Mapping quantitative trait loci with dominant and
  missing markers in various crosses from tow inbred lines. \emph{Genetica}
  \bold{101}, 47-58.  
}
\note{
  It might have been better to design this array so that the third
  subscript moves fastest.  In large problems, the current structure may
  involve excessive memory access.}
