\name{A Starting Point}

\title{Some Introductory Comments}
\description{Some pointers to a few key functions in \emph{BQTL}}

\section{New to R?}{
  \itemize{
    \item Be sure to check out all of the free documentation
  that comes with R.

  \item The \code{example} function is very
  helpful in getting familiar with a new function. You type
  \code{example(fun)} and the examples in the documentation for
  \code{fun} are run, then you can read the documentaiton to get a bette
  sense of what is really going on. My personal favorite is to type
  \code{par(ask=T)}, hit the 'enter' key, then \code{example(image)},
    and 'enter' again; after each display you hit the 'enter' key to get
    to the next one.

  \item 
    \code{library(bqtl)} is needed to load the \emph{BQTL} functions and data
    sets.}
  
}
\section{Key Functions}{
  \describe{
    \item{Data Input}{ }
   
      \describe{
	\item{\code{make.map.frame}}{defines the map,}
	
	\item{\code{marker.levels}}{defines the coding scheme for
	  marker levels, }
	
	\item{\code{make.analysis.obj}}{combines marker data, phenotype
	  data,and the \code{map.frame} to create an object that can be
	  used by data analysis functions.}
	
      }
    
      \item{Maximum Likelihood Methods}{ }
      
      \describe{
       \item{\code{bqtl}}{does a host of things from marker regression
	 and interval mapping to full maximum likelihood. The best way to
	 get started is to run \code{example(bqtl)} and take a look at
	 the resulting output.}
       
       \item{\code{locus}}{ is very helpful in specification of runs.} 
     }
     
     \item{Approximate Bayesian Analysis: }{

       Try \code{example(linear.bayes)}, then look at the help page for it.}
   
 }
 }
\author{Charles C. Berry \email{cberry@ucsd.edu}}


\keyword{Introduction }
