\name{summary.bqtl}
\alias{summary.bqtl}

\title{Summarize bqtl object }
\description{
  Extract coefficients (and related stats), loglikelihood, and residual
  standard error of the trait.  
}
\usage{
summary.bqtl(bqtl.obj)
}

\arguments{
  \item{bqtl.obj}{The result of \code{link{bqtl}} }

}
\value{
  A list containing
  \item{coefficients }{Either a vector of regression coefficents, or if
    {bqtl.obj} was created via \code{bqtl(\dots,return.hess=TRUE)}, then
    a matrix with coefficients, standard errors, t-statistics, and p-values}
  \item{loglik }{
    the loglikelihood or log posterior}
  \item{std.res}{The residual standard deviation of the trait}
  \item{N}{The counts of all observations, the number omitted, and the
    number used in the fit}
}

\author{Charles C. Berry \email{cberry@ucsd.edu}}
\examples{
data(little.ana.bc)
fit <- bqtl( bc.phenotype~locus(4)*locus(45), little.ana.bc,
return.hess=TRUE )
summary(fit)
\testonly{rm(little.ana.bc,fit)}
}
