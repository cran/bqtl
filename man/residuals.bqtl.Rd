\name{residuals.bqtl}
\alias{residuals.bqtl}
\title{ Residuals from QTL models}
\description{
 The phenotype data, estimated coefficients, and expected locus values
are used to find fitted values for the QTL model
}
\usage{
residuals(object)
residuals.bqtl(object)
}
\arguments{
  \item{object}{ An object of class \code{\link{bqtl}}}}
}
\details{
 The estimated coefficients for a specific QTL model fit are used along
 with the  expected locus values (conditionally on the marker values)
 are used to find fitted values for the QTL model; these are subtracted
 from the origianl trait values to get residuals. This is \emph{not} the
 only  way in which such fits could be obtained; one could condition the
 expected marker values on \emph{both} the trait value and the marker
 values. One could also define fitted values for specific genotype
 combinations, e.g. for a backcross with k animals and a two gene model
 4 fitted values could be determined for each animal leading to 2*2*k
 residuals. 
}
\value{ A vector with as many elements trait values used in the original
  fitted model.

 
}
\author{Charles C. Berry \email{cberry@ucsd.edu}}

\seealso{ \code{\link{bqtl}}}

\examples{

data(little.ana.bc)

fit.pheno <- bqtl(bc.phenotype~locus(15)+locus(42),little.ana.bc)

summary(residuals(fit.pheno))

plot( fitted( fit.pheno ), residuals( fit.pheno) )

\testonly{rm(little.ana.bc,fit.pheno)}

}
\keyword{methods}