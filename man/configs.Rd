\name{configs}
\alias{configs}

\title{Lookup loci or effects for genetic model formulas }
\description{
Convert numeric indexes to names of regressors for a genetic model.
One or many genetic models can be specified through the use of this
function.  It is used on the right hand side of a formula in the
\code{\link{bqtl}} function. 

}
\usage{
configs(x, ..., scope=<see below>)

bqtl( y ~ PVV.4 + configs(14,17), my.analysis.object )

bqtl( y ~ configs(14,17) * configs(133,245), my.analysis.object )
}

\arguments{

 \item{x}{ 

Typically an integer, an integer vector, an array, or a list with a
\code{configs} component such as returned by \code{\link{swapbc1} }.
However, it can also be a character string, vector, et cetera, in which
case the elements must belong to \code{names(scope)} }

\item{\dots}{ 
Optional arguments to be used when \code{is.atomic(x)} is TRUE. }
\item{scope}{ (Optional and) 

Usually not supplied by the user.  Rather \code{bqtl} fills this in
automatically.  A vector of regressor names, like the \code{reg.names}
component returned by \code{\link{make.analysis.obj}}.  When
\code{mode(x)} is \code{"character"}, then \code{names(scope)}
\bold{must} be non-\code{NULL} }

}
\details{ \code{configs} is used in the model formula notation of
\code{bqtl}, possibly more than once, and possibly with regressors named
in the usual manner.  \code{configs} is intended to speed up the
specification and examination of genetic models by allowing many models
to be specified in a shorthand notation in a single model formula.  The
names of genetic loci can consist of marker names, names that encode
chromosome number and location, or other shorthand notations.  The names
of terms in genetic models will typically include the names of the locus
and may prepend "add." or "dom." or similar abbreviations for the
'additive' and 'dominance' terms associated with the locus. 

When used as in \code{bqtl( y ~ configs(34), my.analysis.obj )}, it will
look up the term \code{my.analysis.obj$reg.names[34]}.  When
this is passed back to \code{bqtl}, it get pasted into the formula and
is subsequently processed to yield the fit for a one gene model.

When used as in \code{bqtl( y ~ configs(34,75,172), my.analysis.obj)} it
looks up each term and returns a result to \code{bqtl} that results in
fitting a 3 gene model (without interaction terms). 

When \code{x} is a vector, array, or list, the processing typically
returns pieces of many model formulas.  \code{bqtl(y ~ configs(26:75),
\dots)} results in a list of 50 different one gene model fits from
\code{bqtl} for the terms corresponding to the 26th through the 75th
variables.  \code{ bqtl(y ~
configs(cbind(c(15,45,192),c(16,46,193))),\dots)} returns two four gene
models.  And more generally, whenever \code{is.array(x)} is TRUE, the
columns (or slices) specify \code{ dim(x)[1]/length(x)} different
models.  When \code{x$configs} is an array, this also happens.  This turns
out to be useful when the result of running \code{\link{swapbc1}} or
\code{\link{swapf2}} is treated as an importance sample.  In such a case,
\code{bqtl(y ~ configs(my.swap),my.analysis.obj)} will return a list in
which element \code{i} is the ith sample drawn when \code{my.swap <-
swapbc1(\dots)} was run.

}
\value{

A character vector whose element(s) can be parsed as the right hand side
of a model formula. 

}
\author{Charles C. Berry \email{cberry@ucsd.edu} }

\seealso{ 

\code{\link{bqtl}} and the examples there for a sense of how to use
\code{configs}, 
\code{\link{make.analysis.obj}} for the setup that encodes the marker
map and the marker information,
\code{\link{swapbc1}} and 
\code{\link{swapf2}} for generating samples to be screened by \code{bqtl}.

 }

