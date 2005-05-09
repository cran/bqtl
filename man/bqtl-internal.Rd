\name{bqtl-internal}
\alias{\%equiv\%}
\alias{map.dx}
\alias{rhs.bqtl}
\alias{uniq.config}
\alias{zero.dup}
\alias{version.bqtl}
\title{Internal BQTL functions}
\description{
 Internal bqtl functions and objects
}
\usage{
x \%equiv\% y
 map.dx(lambda, theta, min.lambda)
rhs.bqtl(reg.terms, ana.obj, bqtl.specials, local.covar, scope, 
    expand.specials = NULL, method, ...) 
zero.dup(x,dig=6)
uniq.config(swap.obj)

}
\arguments{

  \item{lambda}{ (2*(recomb fraction-1/2) }
  \item{theta}{recomb fraction}
  \item{min.lambda}{smallest map distance to use}
  \item{reg.terms}{a formula}
  \item{ana.obj}{an analysis.object}
  \item{bqtl.specials}{a vector of acceptable special names}
  \item{local.covar}{a function}
  \item{scope}{vector of strings}
  \item{expand.specials}{logical,whether to use expand.grid on the loci}
  \item{method}{e'g' "F2", "BC!", etc}
  \item{...}{not sure}
  \item{swap.obj}{result of swap}
  \item{x}{numeric vector or matrix}
  \item{y}{numeric vector or matrix}
  \item{dig}{how many significant digits to use}

}
\details{
  These are not to be called by the user.
}
\keyword{utilities}
