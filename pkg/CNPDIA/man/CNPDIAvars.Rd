\name{CNPDIAvars}
\alias{CNPDIA0D}
\alias{CNPDIA1D}
\alias{CNPDIAsvar}
\alias{MPBDIA0D}
\alias{MPBDIA1D}
\alias{MPBDIAsvar}
\title{
  Functions to retrieve variables for the CNPDIA and MPBDIA model.
}
\description{

  \code{CNPDIA0D, CNPDIA1D, CNPDIAsvar} retrieve the (0-dimensional or 1-dimensional) output variables or the state variables of CNPDIA model solutions.  When called with \code{out} not specified, will return the names of these variables and their units.

\code{MPBDIA0D, MPBDIA1D, MPBDIAsvar} retrieve the (0-dimensional or 1-dimensional) output variables or the state variables of MPBDIA model solutions.  When called with \code{out} not specified, will return the names of these variables and their units.
}
\usage{
  CNPDIA0D(out, as.vector = FALSE, which = NULL) 
  CNPDIA1D(out, which = NULL) 
  CNPDIAsvar(out, which = NULL) 
  MPBDIA0D(out, as.vector = FALSE, which = NULL) 
  MPBDIA1D(out, which = NULL) 
  MPBDIAsvar(out, which = NULL) 
}
\arguments{

  \item{out }{an output object returned by \link{CNPDIAsolve}, \link{MPBDIAsolve},
   or \link{CNPDIAdyna}, \link{MPBDIAdyna}. 
   If \code{NULL}, \code{CNPDIAparms} or \code{MPBDIAparms} will return the default (parameter) values.
  }
  \item{as.vector }{if \code{TRUE} will return the parameter vector, else a data.frame that also contains the units.
  }
  \item{which }{if not \code{NULL}, a vector with names of the variables/parameters to return.
  }
}
\value{
 \code{CNPDIA0D} and \code{CNPDIA1D} return the output variables of the solution as a vector or data.frame. 
 For dynamic runs, the output is averaged over the mean of the run.
 
 \code{CNPDIA1D} always returns the sediment depth and the porosity as the first two columns.
}

\author{
  Karline Soetaert
}


\examples{

# defaults
  CNPDIAsvar()
  head(CNPDIA0D())
  CNPDIA1D()
  
# some runs to work with  
  defsteady <- CNPDIAsolve()
  defdyn    <- CNPDIAdyna()

# altered steady-state run
  out <- CNPDIAsolve(parms = list(Cflux = 1000))

# 0-D outputs
  cbind(steady = CNPDIA0D(defsteady, as.vector = TRUE), 
        dyna = CNPDIA0D(defdyn, as.vector = TRUE), 
        ref = CNPDIA0D(out))

# 1-D outputs
  head(CNPDIA1D(out))
  head(CNPDIA1D(defdyn, which = c("O2", "TOC")))

}

\references{
  Soetaert K, PMJ Herman and JJ Middelburg, 1996a.
  A model of early diagenetic processes from the shelf to abyssal depths.
  Geochimica Cosmochimica Acta, 60(6):1019-1040.

  Soetaert K, PMJ Herman and JJ Middelburg, 1996b.
  Dynamic response of deep-sea sediments to seasonal variation: a model.
  Limnol. Oceanogr. 41(8): 1651-1668.
}

\keyword{ utilities }
