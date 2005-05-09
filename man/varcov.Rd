\name{varcov}
\alias{varcov}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Create moment matrices}
\description{
  Create a moment matrix of the marker variables and of the regressors by
  the phenotype variable.  For use in regression modelling on the markers.
}
\synopsis{ varcov(x,ana.obj,partial=NULL,scope=ana.obj$reg.names,... ) }
\usage{
varcov(x, ana.obj, partial=NULL, scope=<<see below>>)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A formula to specify the dependent and independent variables
    to be used in subsequent calculations e.g \code{trait ~ locus(.) }}
  \item{ana.obj}{An \code{analysis.object},
    see\code{\link{make.analysis.obj}}  }
  \item{partial}{A formula whose right hand side specifies variables to
    be treated as covariates. }
  \item{scope}{ Usually not explicitly used. Optional vector of variable
    names. } 
}
\details{This is just a wrapper for \code{\link{make.varcov}}.}
\value{
 A list with components
  \item{var.x }{Moment matrix of the marker regressor variables}
  \item{cov.xy }{Moment matrix of the marker regressor variables versus
    the phenotype variable}
  \item{var.y}{The Second central moment of the phenotype variable}
  \item{df}{The degrees of freedom, when no variables are specified in
    \code{partial} it is \code{ sum(subset==TRUE) - 1}
}}
\author{Charles C. Berry \email{cberry@ucsd.edu} }
\note{ It is generally NOT a good idea to do regressions on
  ill-conditioned designs using the moment matrices.  The
  excuse for doing so here is twofold.  First, calculations using this
  method are used to perform importance sampling, so minor numerical
  inaccuracies in computing the probabilites used in sampling get
  straightened out by the importance weights.  Second, it will typically
  be the case that a prior is set on the regression coefficients and
  this results in a positive constant (aka a 'ridge' parameter) being
  added to diagonal of \code{varcov()$var.x} and this reduces the
  ill-conditioning.  Of course the rational for using the method is to
  speed the sampling, and it is very effective at doing so.}
\seealso{The examples in \code{\link{swap}} and \code{\link{twohk}}.}
\keyword{utilities}