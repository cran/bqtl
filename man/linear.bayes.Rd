\name{linear.bayes}
\alias{linear.bayes}

\title{ Bayesian QTL mapping via Linearized Likelihood}
\description{
  The Bayesian QTL models via a likelihood that is linearized w.r.t. a
  fixed genetic model. By default, all one and two gene models (without
  epistasis) are fitted and a MCMC sampler is used to fit 3,4, and 5 gene
  and (optionally) larger models. 
}
\usage{
linear.bayes(x, ana.obj, partial=NULL, rparm=<see below>, specs=<see below>, scope=<see below>, subset=<see below>, casewt=<see below>, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ a formula giving the QTL and the candidate loci or a
    \code{varcov} object }
  \item{ana.obj}{An analysis.object, see \code{\link{make.analysis.obj}}
  }
  \item{partial}{a formula giving covariates to be controlled }
  \item{rparm}{ A ridge parameter. A value of 1 is suggested, but the
    default is 0. }
  \item{specs}{ An optional list with components \code{gene.number} (to
  indicate the model sizes), \code{burn.in} (to indicate the number of
  initial MCMC cycles to discard), and \code{n.cycles} (to indicate how
  many MCMC cycles to perform for each model size). If no values are
  supplied, \code{specs} defaults to
  \code{list(gene.number=c(1,2,3,4,5),burn.in=1,n.cycles=c(0,0,200,100,100))    }}
  \item{scope}{Not generally used. If supplied this will be passed to
    \code{varcov}. } 
  \item{subset}{Not generally used. If supplied this will be passed to
    \code{varcov}.  }
  \item{casewt}{Not generally used. If supplied this will be passed to
    \code{varcov}. }
  \item{\dots}{optional arguments to pass to \code{twohk} and \code{swap} }
}
\details{
  This function is a wrapper for
  \code{\link{varcov}}, \code{\link{twohk}}, \code{\link{swap}}, and
  \code{\link{summary.swap}}, and a better understanding of optional
  arguments and the object generated is gained from their
  documentation.
}
\value{
  \item{hk}{The object returned by \code{\link{twohk}}}
  \item{swaps}{A list of objects returned by calls to
    \code{\link{swap}}. The i$^{th}$ element in \code{swaps} is for i gene
    models.}
  \item{smry}{A list of objects returned by calls to
    \code{\link{summary.swap}}. Some elements may be \code{NULL} if no
    samples were requested or if the sampling process yielded degenerate
  results. Usually, this happens if no posterior is specified for the
  regression coefficients, i.e. if \code{rparm=0} was used or implied }
  \item{odds}{A Vector of odds (relative to a no gene setup) for each
    model size evaluated. The odds are computed under a prior that
    places equal weights on models of each size considered (and are,
    therefore,  Bayes Factors). If models of size 1 and 2 are not
    evaluated or if some degenerate results were encountered, this will
    be \code{NULL}}
  \item{coefs}{A vector of posterior means of the regression
    coefficients. If models of size 1 and 2 are not
    evaluated or if some degenerate results were encountered, this will
    be \code{NULL}}
  \item{loc.posterior}{A vector of locus-wise posterior probabilities
    that the interval covered by this locus contains a gene.If models of
    size 1 and 2 are not evaluated or if some degenerate results were
    encountered, this will be \code{NULL}}
  \item{call}{The call that generated this object}
}
\references{ Berry C.C.(1998) Computationally Efficient Bayesian QTL Mapping
  in Experimental Crosses. \emph{ASA Proceedings of the Biometrics
    Section}. 164--169.
  Also available from
  \url{http://hacuna.ucsd.edu/bqtl/} }

\author{Charles C. Berry \email{cberry@ucsd.edu} }
\examples{
data( little.ana.bc )
little.lin <- linear.bayes( bc.phenotype~locus(all), little.ana.bc, rparm=1 )
par(mfrow=c(2,3))
plot( little.ana.bc, little.lin$loc.posterior, type="h" )
little.lin$odds
par(mfrow=c(1,1))
plot(fitted(little.lin), residuals(little.lin))
\testonly{rm(little.lin,little.ana.bc)}
}

