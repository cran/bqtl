\name{marker.levels}
\alias{bc1.levels}
\alias{f2.levels}
\alias{ri.levels}
\title{ Define marker level codes }
\description{
The coding scheme used to define marker.levels is set up by these
functions. \emph{BQTL} has defaults that these functions can help the
user to redefine.
}
\usage{
bc1.levels( AA="AA", Aa="Aa", miss.val="--")
ri.levels(  AA="AA", aa="aa", miss.val="--")
f2.levels(  AA="AA", Aa="Aa", aa="aa", not.aa="A-", not.AA="a-",
miss.val="--")

make.analysis.obj(data, map.frame, marker.frame,
marker.levels=f2.levels(),method="F2")
}
\arguments{
  \item{AA}{Always used: the code for the homozygous state from one parent line}
  \item{Aa}{F2 and BC1 setups: the code for the heterozygous state }
  \item{aa}{F2 and RI setups: the code for the homozygous state for the
    other parent line }
  \item{not.aa}{F2 only: the code for a dominant marker that rules out  \code{aa} }
  \item{not.AA}{F2 only: the code for a dominant marker that rules out \code{AA} }
  \item{miss.val}{The character string for a missing (unknown) allele
    state. \code{NA}s are automatically detected, so this is only needed
    if string values are used to denote missing values.}
}
\details{It is essential that the codes intended by the user be
clearly understood by \emph{BQTL}. It is hoped that thees functions
provide a bridge between the internals of \emph{BQTL} and the user's
view of the marker codes. Numeric values can be used, but they will be
coerced to character values.
\value{A vector with 6 elements corresponding to the values of 
  \code{AA, \code{Aa}, \code{aa}, \code{not.aa}, \code{not.AA}, and
    \code{miss.val}. For RI and BC1 setups, those that do not apply will
    be unnamed and set to \code{"nil"}
}
\author{ Charles C. Berry \email{cberry@ucsd.edu}}
\seealso{ \code{\link{make.analysis.obj}} }

\examples{

### show the defaults:

f2.levels()
bc1.levels()
ri.levels()

### suppose that 1,2,3 are codes used in F2:

f2.levels(1,2,3)

### show what would happen changing "Aa" to "H"

f2.levels(Aa="H")
bc1.levels(Aa="H")

}
