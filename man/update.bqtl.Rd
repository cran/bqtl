\name{update.bqtl}
\alias{update.bqtl}

\title{Get loglikelihoods for many models of a common form }
\description{
For a single type of model, this function evaluates multiple models that
differ only in terms of the loci involved. The looping is all done by
internal C functions, so this is faster than simply using \code{bqtl} to
do the same thing.
}
\usage{
update.bqtl(setup, loc.mat, ana.obj)
}

\arguments{
  \item{setup}{The object returned by  \code{bqtl( <...> , setup=TRUE )} }
  \item{loc.mat}{A matrix of locus numbers, s.t. \code{nrow(loc.mat)}
    equals the number of loci in \code{setup} }
  \item{ana.obj}{An \code{analysis.object). Usually the one used in \code{setup) }
}
\details{
In order to avoid the computational overhead of running large loops of
very repetitive operations in R/S, \code{update.bqtl} used after the
\code{setup=TRUE} option in \code{\link{bqtl}} will loop through the
loci specified in \code{loc.mat} using internal C code. This is many
times faster than running the same code via \code{bqtl}.
}
\value{
  For now it only returns the loglikelihood. But it would be trivial to
  build an option that would allow other quantities computed to be
  returned, and this should probably be done. However, some care is
  needed to keep objects from becoming unmanageably large.
}

\author{Charles C. Berry \email{cberry@ucsd.edu} }

\seealso{\code{\link{bqtl}}}

\examples{
data( little.ana.bc )
little.setup <-
bqtl( bc.phenotype~locus(1)*locus(2), little.ana.bc, setup=TRUE )
combos <- t( as.matrix( expand.grid( 1:21, 44:64 ) ) )
little.update <- update.bqtl(little.setup, combos, little.ana.bc)
little.res <- matrix( little.update, nr=21 )
image( 1:21, 44:64, little.res )
rm(little.ana.bc, little.update, little.res  )
}

