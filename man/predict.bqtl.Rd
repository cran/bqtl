\name{predict.bqtl}
\alias{predict.bqtl}
\alias{fitted.bqtl}
\title{ fitted values from QTL models}
\description{
 The estimated coefficients and expected locus values are used to find
 fitted values for the QTL model
}
\usage{
predict(object, newdata)
fitted(object, newdata)
}
\arguments{
  \item{object}{ An object of class \code{\link{bqtl}}}}
  \item{newdata}{An optional data.frame for which fitted values are to
    be found. If not specified, the a search for the original data frame
    for the fit will be made. }
}
\details{
 The estimated coefficients for a specific QTL model fit are used along
 with the  expected locus values (conditionally on the marker values)
 are used to find fitted values for the QTL model. This is \emph{not} the only
 way in which such fits could be obtained; one could condition the
 expect marker values on \emph{both} the trait value and the marker
 values. One could also define fitted values for specific genotype
 combinations, e.g. for a backcross with k animals and a two gene model
  4 fitted values could be determined for each animal leading to 2*2*k
  values. In fact, using \code{newdata} one can do this.
}
\value{ A vector with as many elements as rows in newdata (after
  removing missing data) or in the original model.frame.

 
}
\author{Charles C. Berry \email{cberry@ucsd.edu}}

\seealso{ \code{\link{bqtl}}}

\examples{

data(little.ana.bc)

fit.pheno <- bqtl(bc.phenotype~locus(15)+locus(42),little.ana.bc)

summary(predict(fit.pheno))

genotype.grid <- expand.grid( c(-1,1), c(-1,1) )     # set up a grid
names(genotype.grid) <- map.names( fit.pheno )       # use matching names

fit.vals <- predict( fit.pheno, genotype.grid )      # make predictions
cbind( genotype.grid, fit.vals )                     # print them!

}

