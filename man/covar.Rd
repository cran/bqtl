\name{covar}
\alias{covar}

\title{Treat locus as covariate }
\description{
  Sometimes it is helps speed computations to linearize the likelihood
  or at least a part of it w.r.t. the locus allele values. Both
  'Haley-Knott regression' and 'composite interval mapping' use this
  approach. \code{covar} provides a mechanism for creating formula
  objects that specify such linearizations. 
}
\usage{
covar(x,...,scope=<see below>, method=<see below>)
}
\arguments{
  \item{x}{The name of a locus (except for F2 designs, when it is the
    name of an effect like 'add.m.32') or any argument of the sort that 
    \code{\link{locus}} allows.}
  \item{...}{If \code{x} evaluates to a single value, then additional
    atomic elements may be included as with \code{locus}.}
  \item{scope}{Not supplied by the user. see \link{locus}}
  \item{method}{Not supplied by the user. see \link{locus}}

}
\details{
  The function \code{covar} actually only returns \code{x}. The real
  work is done by a \code{covar} function that is hidden inside of
  \code{bqtl}, where the arguments are parsed as for \code{locus}. Each
  of the return values from \code{locus} is prefixed by "covar(" and
  suffixed by ")". If \code{x} is a name of a locus or effect, then
  \code{paste("covar(",deparse(x),")")} is returned. Later, when
  \code{bqtl} calls \code{lapadj}, terms like \code{covar(PVV4.1)} are
  recognized as requiring a linearization w.r.t. effect 'PVV4.1'.}
}
\value{
  a character string or vector 

  \references{

    HALEY, C. S. and S. A. KNOTT, 1992  A simple regression method for
    mapping quantitative trait loci in line crosses using flanking
    markers. Heredity 69:315-324.

    Knapp SJ, Bridges WC, and Birkes D. Mapping quantitative trait loci
    using molecular marker linkage maps. Theoretical and Applied Genetics
    79: 583-592, 1990.
    
    ZENG, Z.-B., 1994 Precision mapping of quantitative trait
    loci. Genetics 136:1457-1468
  }
 \author{Charles C. Berry \email{cberry@ucsd.edu}}

   \seealso{\code{\link{locus}}, \code{\link{add}}, \code{\link{dom}},
     \code{\link{configs}} 
   }


