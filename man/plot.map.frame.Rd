\name{plot.map.frame}
\alias{plot.map.frame}
\alias{plot.analysis.object}
\title{ plots by chromosome location }
\description{
  Multiple x-y plots are formed using chromosome numbers
  (\code{chr.num}) and positions (\code{pos.plot}) specified in a object
  of the sort created by \code{\link{make.map.frame}}
}
\usage{
plot(x, ...)

plot.map.frame(mf.obj, y, fun=<see below>, type="l",
 include.rug=TRUE, title.string=NULL, ylab=deparse(substitute(y)),
 xlab="Location", ...)
}

\arguments{
 \item{mf.obj}{A \code{map.frame.object} or an \code{analysis.object}  }
 \item{y}{(optional) A vector with as many elements or a
   matrix with as many rows as \code{nrow(mf.obj)} .  If omitted, a plot
   will be drawn in a single frame representing the chromosomes as
   horizontal bars and giving tic marks to show the locations markers
   and virtual markers (if any).  }
 \item{fun}{ A plotting function to be used \emph{after} the plot axes and
   labels have been drawn.  The current default \code{if (y.type ==
     "matrix") matlines else lines} usually is good enough.  But a
   fancier function could be used for a fancier plot. }
 \item{type}{ \code{"l"} for lines, "p" for points, et cetera.  see \code{\link{par}} }
 \item{include.rug}{ if \code{TRUE} place a tick on the x-axis at each marker location }
 \item{title.string}{(optional) label to prepend to each title }
 \item{ylab}{plot label for y-axis, see \code{\link{par}}  }
 \item{xlab}{plot label for x-axis, see \code{\link{par}}  }
 \item{\dots}{(optional args to pass to \code{fun}  }
}
\details{
  This function enables drawing graphs that depend on chromosome and
  chromosome location.  Typically, one will use a command like
  \code{par(mfrow=c(nrows,ncols))} first to set up a page on which
  multiple plots will be drawn.  However, one can draw one plot per page
  on postscript devices by leaving \code{par(mfrow=c(1,1))}
}
\value{
  \code{NULL} - this function is called only for its side effects
}

\author{Charles C. Berry \email{cberry@ucsd.edu} }

\seealso{
  \code{\link{plot}},
\code{\link{lines}}, and 
\code{\link{matlines}} for general information on plotting functions;
\code{\link{par}} for optional arguments to add as arguments; and 
\code{\link{make.map.frame}} for the details on the object the drives this function.

}

\examples{

data( little.ana.bc )
null.llk <- loglik(bqtl(bc.phenotype~1,little.ana.bc))
llk <- loglik( bqtl( bc.phenotype~locus(all), little.ana.bc) ) - null.llk
.old.par <- par(mfrow=c(2,3))
plot.map.frame(little.ana.bc$map.frame,llk)
par(.old.par)
\testonly{ rm(null.llk, llk, little.ana.bc, .old.par ) }

}
\keyword{hplot,methods	}