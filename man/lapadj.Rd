\name{lapadj}
\alias{lapadj}

\title{Approximate marginal posterior for chosen model}

\usage{ lapadj(reg.formula, loc.right, marker.distances, state.matrix,
marker.names=dimnames(state.matrix)[[2]], rparm=NULL, casewt=NULL,
tol=9.9e-09, return.hess=F, mode.names=NULL, mode.mat=NULL,
nc = 1), method="BC1", maxit=100, nem=1, ...)  }


\arguments{

  \item{reg.formula}{ A formula, like
    \code{y~add.X.3+dom.X.3+add.x.45*add.x.72}}
  \item{loc.right}{ See\code{\link{make.analysis.obj}}, which returns
    objects like this. It is a matrix of pointers to the next marker with a
    known state on the current chromosome (if any). }
  \item{marker.distances}{Distances between the markers in the
    'lambda' metric. \code{-log(lambda)/2} is the Haldane map
    distance. Linkage groups are separated by values of 0.0.  }
  \item{state.matrix}{ See\code{\link{make.analysis.obj}}, which returns
    objects like this. An n by k by q array. q is 2 for method="BC1" and
    3 for method="F2". Each element encodes the probability of the
    allele state conditional on the marker states. see
    \code{\link{make.state.matrix}} for more details.  } 
  \item{marker.names}{The names by which the markers are known.}
  \item{rparm}{

    One of the following:

    A scalar that will be used as the ridge parameter for all
      design terms except for the intercept ridge parameter which is set
      to zero 

    
   
      A vector who named elements  can be matched by the design term
      names returned in \code{$reg.vec}. If no term named
      "intercept" is provided, \code{rparm["intercept"]} will be set to
      zero. 

 
      A vector with \code{(q-1)*k} elements (this works when there are
      no interactions specified). If names are provided, these will be
      used for matching.

    Positive entries are 'ridge' parameters or variance ratios in a
    Bayesian prior for the regression coefficients. Larger values imply
    more shrinkage or a more concentrated prior for the regresion
    coefficients.   
  } 
  \item{tol}{Iteration control parameter  }
  \item{return.hess}{Logical, include the Hessian in the output? }
  \item{mode.names}{names to use as \code{dimnames(mode.mat)[[2]]}  }
  \item{mode.mat}{Not usually set by the user. A matrix which indicates the values of regressor
    variables corresponding to the allele states. If \code{mode.mat} is
    not given by the user,
   \code{ana.obj$mode.mat} is used.}
  \item{method}{ Currently, "BC1", "F2", "RI.self" and "RI.sib" are recognized. }
  \item{maxit}{ Maximum Number of iterations to perform}
  \item{nem}{ Number of EM iterations to use in reinitializing the
    pseudo-Hessian }
  \item{\dots}{ other objects needed in fitting }
}
\description{
  lapadj provides the Laplace approximation to the marginal posterior
  (over coefficients and dispersion parameter) for a given genetical
  model for a quantitative trait. A by-product is the parameter value
  corresponding to the maximum posterior or likelihood.
}
\details{
  The core of this function is a quasi-Newton optimizer due to Minami
  (1993) that has a computational burden that is only a bit more than
  the EM algorithm, but features fast convergence.  This is used to find
  the mode of the posterior. Once this is in hand, one can find the
  Laplace approximation to the marginal likelihood. In addition, some
  useful quantities are provided that help in estimating the marginal
  posterior over groups of models. 
}
\value{
  A list with components to be used in constructing approximations to the
  marginal posterior. These are:
  \item{adj}{The ratio of the laplace approximation to the posterior for
    the correct likelihood to the laplace approximation to the posterior
    for  the linearized likelihood}
  \item{logpost}{The logarithm of the posterior or likelihood at the mode}
  \item{parm}{the location of the mode}
  \item{posterior}{The laplace approximation of the marginal posterior
    for the exact likelihood}
  \item{hk.approx}{ Laplace approximation to the linearized likelihood}
  \item{hk.exact}{ Exact marginal posterior for the linearized likelihood}
  \item{reg.vec}{A vector of the variables used}
  \item{rparm}{
    Values of ridge parameters used in this problem.}
}
\references{ Berry C.C.(1998) Computationally Efficient Bayesian QTL Mapping
  in Experimental Crosses. \emph{ASA Proceedings of the Biometrics
  Section}. 164--169.
  
  Minami M. (1993) Variance estimation for simultaneous response growth curve
models.  Thesis (Ph. D.)--University of California, San Diego,
Department of Mathematics.
  }
  \author{Charles C. Berry \email{cberry@ucsd.edu} }
 \keyword{regression}