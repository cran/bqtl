\name{twohkbc1}
\alias{twohkbc1}
\alias{twohkf2}

\title{ One and Two Gene Models Using Linearized Posterior}
\description{
  Fits all one and two gene models (without interactions aka
  'epistasis') in an intercross, backcross, or recombinant inbred
  line. Uses a linear 
  approximation to the likelihood, i.e. the expected allele states are
  used.
}
\usage{
twohkbc1(varcov, rparm, locs=<<see below>> ), locs.prior=<<see below>> )

twohkf2(varcov, ana.obj, rparm, locs=<<see below>>, locs.prior=<<see below>>,
combo.prior=<<see below>>) 
}
\arguments{
  \item{varcov}{ An object produced by \code{ \link{make.varcov}} }
  \item{ana.obj}{An object produced by \code{ \link{make.analysis.obj}}}
 \item{rparm}{ The 'ridge' parameters for the independent variables -
   larger values imply more shrinkage or a more concentrated prior for
   the regresion coefficients. } 
 \item{locs}{ The columns (or pairs of columns for \code{twohkf2}) of
   \code{varcov\$var.x} to use. The default uses all of them.}
  \item{locs.prior}{ 
    The prior mass to  associate with each locus. Typically, these sum
    to one, but sometimes they might each be set to one (as in computing
    lod scores). }
\item{combo.prior}{Only valid for \code{twohkf2}. 
    The prior probability for each term or combination of terms for the
    phenotypic effect at a locus. Typically, there will be three of
    these - one for the 'additive' term (linear in number of alleles
    from one parent strain), the 'dominance' term (quadratic in allele
    number), or both terms. The default sets them all to 1/3.
  }
}
\details{
  The marginal posterior (integrating over regression parameters and
  dispersion) is calculated for each one and two gene model under the
  assumed correctness of the regression model using expected genotypes
  given marker values. This amounts to linearizing the likelihood with
  respect to the  (possibly unknown) locus states. For models where the
  loci are fully informative markers this is the true posterior.
}
\value{
 A list with components:
 \item{loc.1}{The marginal posterior for each one gene model.  For
   \code{twohkf2} this is a matrix of 3 columns; the first for models
   with additive terms, the second for dominance terms, and the third
   for both. The sum over all three columns yields the marginal
   posterior for the locus.}
 \item{loc.2}{The marginal posterior for each locus - obtained by summing
   over all two gene models that include that locus. For
   \code{twohkf2} this is a matrix of 3 columns; the first for models
   with additive terms, the second for dominance terms, and the third for both. }
 \item{coefs.1}{The regression coefficients for the genetic effect for 
   each locus. For \code{twohkf2}, this is a matrix with two rows; the
   first is for the 'additive effect' and the second is for the
   'dominance' effect.}
 \item{coefs.2}{The marginal posterior mean of regression coefficients
   for the genetic effect for each locus - obtained by averaging over
   all two gene models that include that locus according to the
   posterior masses. For \code{twohkf2}, this is a matrix with two rows; the
   first is for the 'additive effect' and the second is for the
   'dominance' effect.}
}
\references{
  Haley C.S. and Knott S.A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking
  markers. \emph{Heredity} \bold{69},315-324. } 
\author{Charles C. Berry \email{cberry@ucsd.edu} }
\keyword{utilities}