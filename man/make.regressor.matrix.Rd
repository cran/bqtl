\name{make.regressor.matrix}
\alias{make.regressor.matrix}
\title{
  Create regressors using expected marker values
}
\description{
  Create regression variables for markers and loci between or near
  markers by imputation conditional on known marker states.
}
\usage{
make.regressor.matrix(state.matrix, mode.mat=NULL)
}
\arguments{
\item{state.matrix}{
A state.matrix.object - see \code{\link{make.state.matrix}} for more details
}
\item{mode.mat}{
  A matrix which indicates the values of regressor variables
  corresponding to the allele states.  If mode.mat=NULL (the default) a
  mode.mat is inferred from the dimensions of state.matrix.  For the F2
  intercross these are typically additive and dominance codes like
  (-1,0,1) and (1,-1,1).  For BC1 backcross and RI lines, typically the
  values are (-1,1).
}}
\value{
  A matrix with variables suitable for use as regressors.
}
\references{
  Haley C.S. and Knott S.A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking
  markers. \emph{Heredity} \bold{69},315-324.
}
\seealso{
  \code{\link{make.state.matrix}}
}
