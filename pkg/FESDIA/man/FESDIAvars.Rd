\name{FESDIAvars}
\alias{FESDIA0D}
\alias{FESDIA1D}
\alias{FESDIAsvar}
\alias{FESDIAextract}
\alias{FESDIApH}
\alias{PHDIA0D}
\alias{PHDIA1D}
\alias{PHDIAsvar}
\title{
  Functions to retrieve variables for the PHDIA and FESDIA model.
}
\description{

  \code{FESDIA0D, FESDIA1D, FESDIAsvar} retrieve the (0-dimensional or 1-dimensional) output variables or state variables of FESDIA model solutions. When called with \code{out} not specified, will return the names of these variables and their units.

  \code{PHDIA0D, PHDIA1D, PHDIAsvar} retrieve the (0-dimensional or 1-dimensional) output variables or state variables of PHDIA model solutions. When called with \code{out} not specified, will return the names of these variables and their units.

  \code{FESDIAextract} extracts model output variables of any dimension.

  \code{FESDIApH} estimates the pH profile(s) from FESDIA model solutions.

}
\usage{
  FESDIA0D(out, as.vector = FALSE, which = NULL) 
  FESDIA1D(out, which = NULL) 
  FESDIAsvar(out, which = NULL) 
  
  PHDIA0D(out, as.vector = FALSE, which = NULL) 
  PHDIA1D(out, which = NULL) 
  PHDIAsvar(out, which = NULL) 
  
  FESDIAextract(out, which = NULL) 
  FESDIApH(out)
}
\arguments{

  \item{out }{an output object returned by \link{FESDIAsolve} or \link{FESDIAdyna}. If \code{NULL}, \code{FESDOIAparms} will return the default (parameter) values.
  }
  \item{as.vector }{if \code{TRUE} will return the parameter vector, else a data.frame that also contains the units.
  }
  \item{which }{if not \code{NULL}, a vector with names of the variables to return.
  }

}
\value{
 \code{FESDIA0D} and \code{FESDIA1D} return the output variables of the solution as a vector or data.frame. 
 For dynamic runs, the output is averaged over the run.
 
 \code{FESDIA1D} always returns the sediment depth and the porosity as the first two columns.
}

\author{
  Karline Soetaert
}


\examples{

# =====================
# defaults
# =====================

  FESDIAsvar()
  head(FESDIA0D())
  FESDIA1D()
  
# =====================
# variables of runs  
# =====================

# defaults
  defsteady <- FESDIAsolve()
  defdyn    <- FESDIAdyna()

# altered steady-state run
  out  <- FESDIAsolve(parms = list(Cflux = 10))
  dyna <- FESDIAdyna(parms = list(Cflux = 10), CfluxForc = list(amp = 1))

# 0-D outputs
  cbind(steady.default = FESDIA0D(defsteady, as.vector = TRUE), 
        dyna.default   = FESDIA0D(defdyn, as.vector = TRUE), 
        out            = FESDIA0D(out,  as.vector = TRUE),
        dyna           = FESDIA0D(dyna, as.vector = TRUE))

# 1-D outputs
  head(FESDIA1D(out))
  head(FESDIA1D(defdyn, which = c("O2", "TOC")))

# and the pH
  mf <- par (mfrow = c(1,2))
  plot(x = FESDIApH(out), y = FESDIAdepth(dyna), ylim = c(10,0), type = "l")

  image2D(FESDIApH(dyna), y = FESDIAdepth(dyna), ylim = c(10,0))


  par(mfrow = mf)
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
