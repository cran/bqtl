\name{predict.linear.bayes}
\alias{predict.linear.bayes}
\alias{fitted.linear.bayes}
\alias{residuals.linear.bayes}
\title{ Residuals or Predicted Values for linear.bayes objects }
\description{
 The \code{linear.bayes} object returns fitted coefficients. These are
 used to construct predicted values. Since the fitting process for
 \code{linear.bayes} objects is based on moments of centered variables,
 the 'intercept' is lost; see `Details' below.
}
\usage{
\method{residuals}{linear.bayes}(object, \dots)
\method{predict}{linear.bayes}(object, newdata = lb.call$ana.obj, return.resids =
                 FALSE, ...)
\method{fitted}{linear.bayes}(\dots)
}
\arguments{
  \item{object}{An object returned by \code{linear.bayes} }
  \item{\dots}{possibly the following}
  \item{newdata}{Optional \code{data.frame} in which to do the calculations }
  \item{return.resids}{Not usually set by the user. } }

\details{ Since the linear.bayes \code{object} is based on a moment
matrix, some information is lost thsat must be reconstructed or
assumed. The intercept and possibly the coefficients for control
variates are aong these. Also, when the call to \code{linear.bayes}
supplied the moment matrix rather than formulae with which to create
one, then it is unclear what variable was used as the regressand and
hence which variable to use in ofrming residuals. So, in that case,
\code{residuals} will report an error }

\value{
  A vector of predicted values or residuals
}
\author{ Charles C. Berry \email{cberry@ucsd.edu}}

\seealso{ \code{\link{linear.bayes}}}

\keyword{methods}