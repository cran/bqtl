\name{locus}
\alias{locus}
\alias{add}
\alias{dom}
\title{Lookup loci or effects for genetic model formulas }
\description{
Convert numeric indexes to names of regressors for a genetic model.
One or many genetic models can be specified through the use of this
function.  It is used on the right hand side of a formula in the
\code{\link{bqtl}} function. 

}
\synopsis{
locus(x, ..., scope, method, chromo, cM, ana.obj)}
\usage{
locus(x, ..., scope, method, chromo, cM, ana.obj)
add(x)
dom(x)
}

\arguments{

 \item{x}{ 

   Typically an integer, an integer vector, or an array whose elements
   are integers. These index loci described in a \code{map.frame}
   object. 
   
   However, \code{x} can also be a character string, vector, et cetera,
   in which case the elements must belong to \code{names(scope)}. }

\item{\dots}{ 
  Optional arguments (usually integers) to be used when
  \code{is.atomic(x)} is TRUE. 
}
\item{chromo }{A chromosome number or 2 ordered numbers. The loci on the
  chromosome or in the range of chromosome numbers are used. If
\code{chromo} is used, \code{x} must not be used.}
\item{cM }{(Optional) map distance or two giving a location near a locus
  or range of locations from which loci will be included. If the one
  chromosome number is specified in \code{chromo}, \code{cM} must be 
  ordered. If \code{cM} is omitted, all loci on the chromosome(s) will
  be included.}

\item{scope}{ (Optional and)
  Usually not supplied by the user.  Rather \code{bqtl} fills this in
  automatically.  A vector of regressor names, like the \code{reg.names}
  component returned by \code{\link{make.analysis.obj}}.}

\item{method}{ (Optional and) 
  Usually not supplied by the user. Like \code{scope}, \code{bqtl} takes
  care of filling this in with "BC1", "F2", et cetera as appropriate.}
\item{ana.obj}{Usually not specified by the user. This is the
  \code{analysis.object} to be used to lookup loci if a \code{chromo}
  argument is used.}  }
\details{ \code{locus} is used in the model
formula notation of \code{bqtl}, possibly more than once, and possibly
with regressors named in the usual manner.  \code{locus} is intended to
speed up the specification and examination of genetic models by allowing
many models to be specified in a shorthand notation in a single model
formula.  The names of genetic loci can consist of marker names, names
that encode chromosome number and location, or other shorthand
notations.  The names of terms in genetic models will typically include
the names of the locus and may prepend "add." or "dom." or similar
abbreviations for the 'additive' and 'dominance' terms associated with
the locus.

When used as in \code{bqtl( y ~ locus(34), my.analysis.obj )}, it will
look up the term or terms corresponding to the 34$^{th}$ locus. When
this is passed back to \code{bqtl}, it is pasted into a text string that
will become a formula and is subsequently processed to yield the fit for
a one gene model.

When used as in \code{bqtl( y ~ locus(34,75,172), my.analysis.obj)} it
looks up each term and returns a result to \code{bqtl} that results in
fitting a 3 gene model (without interaction terms). 

When \code{x} is a vector or array, the processing typically returns
pieces character strings for many model formulas.  \code{bqtl(y ~
locus(26:75), \dots)} results in a list of 50 different one gene model
fits from \code{bqtl} for the terms corresponding to the 26th through
the 75th variables.  \code{ bqtl(y ~
locus(cbind(c(15,45,192),c(16,46,193))),\dots)} returns two three gene
models.  And more generally, whenever \code{is.array(x)} is TRUE, the
columns (or slices) specify \code{ dim(x)[1]/length(x)} different
models. 

The \code{chromo} argument performs a lookup of loci on the chromosome
via the function \code{\link{map.index}}. If \code{cM} is also given,
the locus nearest that location is used. If two values are given for
\code{cM} all loci in the range are used.

\code{add(x)} and \code{dom(x)} are alternatives that specify that only
the \emph{additive} or \emph{dominance} terms in an F2 intercross. 
}
\value{

A character vector whose element(s) can be parsed as the right hand side
of a model formula(s). 

}
\author{Charles C. Berry \email{cberry@ucsd.edu} }

\seealso{ 

  \code{configs},
  \code{\link{bqtl}}, and the examples there for a sense of how to use
  \code{locus}, 
  \code{\link{make.analysis.obj}} for the setup that encodes the marker
  map and the marker information.
  
}
\keyword{regression}