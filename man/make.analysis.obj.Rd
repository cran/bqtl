\name{make.analysis.obj}
\alias{make.analysis.obj}

\title{Set up data for QTL mapping }
\description{
Create commonly used objects for the analysis of a backcross or
intercross experiment or of recombinant inbred lines.
}
\usage{
make.analysis.obj(data, map.frame, marker.frame, marker.levels=NULL,
method="F2", casewt=NULL,varcov=FALSE,mode.mat=NULL)
}

\arguments{
 \item{data}{ A \code{data.frame} (or vector) of phenotype and
   (optionally) covariate information }
 \item{map.frame}{ A \code{map.frame.object} (see \code{ \link{make.map.frame} })  encoding the map
   information and other details of the study }
 \item{marker.frame}{ A \code{marker.frame.object}.  A matrix or data.frame
   of marker state information.  }
 \item{marker.levels}{A vector of length six or \code{NULL}. If
   \code{NULL} then the defaults for the elements are:

   \tabular{lccc}{
 Element    \tab  F2.default \tab  BC.default\tab RI.default \cr
 1  \tab   \code{"AA"}           \tab     \code{"AA"}
 \tab     \code{"AA"}   \cr
 2  \tab   \code{"Aa"}           \tab     \code{"Aa"}
 \tab      \code{"aa"}  \cr
   3  \tab   \code{"aa"}           \tab    \code{"nil"}
   \tab       \code{"nil"}  \cr
   4  \tab   \code{"A-"}  \tab  \code{"nil"}
   \tab   \code{"nil"}     \cr
   5  \tab  \code{"a-"}   \tab    \code{"nil"}
   \tab   \code{"nil"}      \cr
   6  \tab  \code{"--"}     \tab    \code{"--"}
   \tab   \code{"--"} \cr
  }

   NA's are allowed in \code{marker.frame} as well as the sixth
   element(\code{"--"} by default) to denote missing data.  To use other
   coding schemes replace "AA" and "aa" by codes for homozygous states,
   "Aa" by the code for heterozygotes, "A-" by the code for 'not aa',
   "a-" by the code for 'not AA', and \code{"--"} by the missing code.
   Positions 3:5 are just placeholders if \code{method!="F2"}, but must
   be present.} 
 \item{method}{One of "F2", "BC1", "RI.self", or "RI.sib" }
 \item{casewt}{If there are multiple observations on one genotype (such
   as in recombinant inbreds) this can be used to assign a weight to
   each observation. The wisdom of doing this is debatable. }
 \item{varcov}{If FALSE, don't create a varcov.object. Otherwise give an
   index into data to select a dependent variable. See \code{ \link{varcov}}
 }
 \item{mode.mat}{If \code{NULL} use the default. For
   \code{method=="F2"} ( and the default \code{marker.levels} of
   \code{AA}, \code{Aa}, and \code{aa} ), 
   this is a 3 by 2 matrix:
   \tabular{lrr}{
     Genotype \tab add \tab dom \cr
     AA \tab -1 \tab -1 \cr
     Aa \tab 0 \tab 1 \cr
     AA \tab 1 \tab -1 \cr
   }
   For \code{method=="BC1"} ( and the default \code{marker.levels} of
   \code{AA} and \code{Aa} ),it is
   \tabular{lr}{
     Genotype \tab \cr
     AA \tab -1 \cr
     Aa \tab 1 \cr
   }
   and for RIL methods ( and the default \code{marker.levels} of
   \code{AA} and \code{aa} ),it is
   \tabular{lr}{
     Genotype \tab \cr
     AA \tab -1 \cr
     aa \tab 1 \cr
   }
   Other choices of \code{marker.levels} will relabel the corresponding
rows.
 }
}
\details{
  A lot of stuff is bundled together in one object. The function is
  really just a wrapper calling other \code{make.*} functions.
  
}
\value{
  A list with components
  \item{data }{\code{data.frame} of phenotype, covariate information,
    and regressors created by \code{ make.regressor.matrix }}
  \item{varcov}{A varcov.object. See  \code{ make.varcov } }
  \item{reg.names}{The names of the regressors from \code{
      make.regressor.matrix }}
  \item{method}{The \code{ method } argument in the call.}
  \item{state.matrix}{See \code{ make.state.matrix }}
  \item{loc.right}{See \code{ make.loc.right }}
  \item{map.frame}{See \code{ make.map.frame }}
  \item{casewt}{The \code{ casewt } argument}
  \item{mode.mat}{The \code{mode.mat} used}
  \item{call}{The function call}
}

\author{Charles C. Berry \email{cberry@ucsd.edu}}
\note{ This can be quite a LARGE object.It might be better
  in crosses with lots (say, thousands) of markers, or in which many
  'virtual' markers are used, or on computers with limited RAM to
  store each component separately. Not all components are used in
  every type of analysis. }



\seealso{

  \code{\link{make.map.frame}} for definition of the marker map, The
  internally used functions are: \code{\link{make.loc.right}},
  \code{\link{make.state.matrix}}, \code{\link{make.regressor.matrix}},
  \code{\link{make.varcov}},  and \code{\link{make.marker.numeric}} }

\examples{
data( little.bc.pheno )
data( little.mf.5 )
data( little.bc.markers )
names(little.bc.pheno)
little.ana.bc <- make.analysis.obj(little.bc.pheno$bc.phenotype,
                                   little.mf.5,little.bc.markers,
                                   method="BC1")
summary( little.ana.bc )

\testonly{rm( little.ana.bc, little.bc.pheno, little.mf.5, little.bc.markers )}
}
\keyword{manip}