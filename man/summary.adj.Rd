\name{summary.adj}
\alias{summary.adj}

\title{Summarize Laplace approximations }
\description{ The linear approximations of \code{swap} are much improved
  by the use a Laplace approximations for loci that are not
  markers. This function combines the results of a call like
  \code{bqtl(y~configs(swap.obj),\dots)} with the data in
  \code{swap.obj} to provide improved posteriors, et cetera}
\usage{

\method{summary}{adj}(object, n.loc, coef.znames, mode.names=c("add",
"dom"), imp.denom=NULL, swap.obj=NULL,\dots)
}
\arguments{
  \item{object}{Typically, this is the result of a call like
  \code{bqtl(y~configs(swap.obj),\dots)}  }
  \item{n.loc}{The number of genes in this model }
  \item{coef.znames}{\code{map.names} for the sample space }
  \item{mode.names}{\code{NULL} except for "F2", in which case it is
    uusally c("add","dom") }
  \item{imp.denom}{Optional, and only used when some sampling scheme
    other than the default MCMC generates \code{object} }
  \item{swap.obj}{The result of a call to \code{swap} }
  \item{\dots}{unused}
}
\details{
  There are a lot of details. This sections nneds to be revised to
  reflect them.
}
\value{
A list with components
  \item{adj }{This multiplier adjusts the posterior odds for k vs k-1
    gene models }
  \item{var }{An estimate of the variance of \code{adj}}
  \item{coef}{Posterior means of coefficients}
  \item{loc}{Marginal Posterior for location for k gene model}
\item{hk.ratio.mean}{argh! I need to look this up}
}
\references{Berry C.C. (1998) Computationally Efficient Bayesian QTL Mapping in
  Experimental Crosses. \emph{ASA Proceedings of the Biometrics
    Section,} 164-169. 
}
\author{Charles C. Berry \email{cberry@ucsd.edu} }
\keyword{methods}
