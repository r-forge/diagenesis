\name{FESDIAsolve}
\alias{FESDIAsolve}

\title{
  Steady-state solution of the FESDIA model, calculating C, N, P, O2, Fe and S dynamics in the sediment.
}
\description{

  \code{FESDIAsolve} finds the steady-state solution of the FESDIA model; pH can be calculated afterwards, ie ignoring carbonate dynamics.
    
}
\usage{
  FESDIAsolve (parms = list(), yini = NULL, gridtype = 1, Grid = NULL, 
      porosity = NULL, bioturbation = NULL, irrigation = NULL,   
      surface = NULL, diffusionfactor = NULL, 
      dynamicbottomwater = FALSE, ratefactor = NULL, 
      calcpH = FALSE, verbose = FALSE, 
      method = NULL, times = c(0, 1e+06), ...)
}
\arguments{

  \item{parms}{ A list with parameter values.Available parameters can be listed using function \link{FESDIAparms}.
    
  See details of \link{FESDIAparms}.
  }

  \item{gridtype }{Type of grid: \code{1} for cartesian, \code{2} for cylindrical, \code{3} for spherical.
  }

 \item{Grid}{ If specified: either an object, as returned by \code{setup.grid.1D} from the package \code{ReacTran}, 
   a vector of length 101 with the transport distances (from mid to mid of layers, upper interface = diffusive boundary layer), 
   or one number with the constant layer thickness. 
   If \code{NULL}, it is defined as \code{setup.grid.1D(x.up = 0, dx.1 = 0.01, N = 100, L = 100)}, i.e. the total length is 100 cm, 
   the first layer is 0.01 cm thick and layers are increasing with depth for 100 layers.  
  } 

 \item{porosity}{ If specified, either an object with porosities ([-]) as returned by \code{setup.prop.1D} from the package \code{ReacTran}, 
   a vector of length 101 with the porosities defined at the layer interfaces, or one number with the constant porosity. 
   If \code{NULL}, it is defined by the parameters \code{por0}, \code{pordeep} and \code{porcoeff} as: 
   \code{(pordeep + (por0 - pordeep) * exp(-pmax(0, x.int)/porcoeff))}, where \code{x.int} is the distance, from the surface of the layer interface. 
   Note that the porosity values should be consistent withe the \code{Grid} - and should be inbetween 0 and 1.
  }

 \item{bioturbation}{ If specified, either an object with bioturbation rates (units [cm2/d]) as returned by \code{setup.prop.1D} from the package \code{ReacTran}, 
   a vector of length 101 with the bioturbation defined at the layer interfaces, or one number with the constant bioturbation. 
   If \code{NULL}, it is defined by the parameters \code{biot}, \code{biotdepth} and \code{biotatt} as: 
   \code{biot * exp(-pmax(0, (x.int-biotdpeth))/biotatt)}, where \code{x.int} is the distance, from the surface of the layer interface. 
   Note that the bioturbation values should be consistent with the \code{Grid}.
  }

 \item{irrigation}{ If specified, either an object with irrigation rates (units [/d]) as returned by \code{setup.prop.1D} from the package \code{ReacTran}, 
   a vector of length 100 with the irrigation rates defined at the layer centres, or one number with the constant rates. 
   If \code{NULL}, it is defined by the parameters \code{irr}, \code{irrdepth} and \code{irratt} as: 
   \code{irr * exp(-pmax(0, (x-irrdepth)/irratt))}, where \code{x} is the distance, from the surface, of the layer centres. 
   Note that the irrigation values should be consistent with the \code{Grid}.
  }

 \item{surface}{ If specified, either an object with surface areas (units [cm2]) as returned by \code{setup.prop.1D} from the package \code{ReacTran}, 
   a vector of length 101 with the surface areas defined at the layer interfaces, or one number with the constant surface area. 
   If \code{NULL}, it is defined by the parameter \code{gridtype}, and the \code{Grid} as: 
   \code{surface = 1} for \code{gridtype == 1}, \code{surface = rev(2*pi*Grid$x.int)} for \code{gridtype == 2}  and \code{surface = rev(pi*(Grid$x.int)^2)} for \code{gridtype == 3}.

   Note that the surface values should be consistent with the \code{Grid}.
  }

 \item{diffusionfactor}{ The multiplication factor necessary to go from molecular diffusion to effective sediment diffusion, i.e. that takes into account tortuosity.
   If specified, either an object with these factors ([-]) as returned by \code{setup.prop.1D} from the package \code{ReacTran}, 
   a vector of length 101 with these factors defined at the layer interfaces, or one number with the constant factor. 
   If \code{NULL}, it is set equal to the porosity. Note that the factors should be consistent with the \code{Grid}.
  }

  \item{yini }{Initial guess of the steady-state solution.
  }

 \item{dynamicbottomwater}{ If \code{TRUE}, then the concentrations in the water overlying the sediment will also be dynamically described, and with water height equal to \code{Hwater}. Note that this will slow down the simulation.   
  }

 \item{ratefactor}{ \code{NULL} or a list, detailing the forcing function for the  biogeochemical rate multiplication factor. If not specified (or \code{NULL}), then it is assumed to be 1, if a data series, then the mean of the data will be used. 
    If a \code{list}, it should contain either a data time series (\code{list(data = )}) or parameters determining the periodicity of the seasonal signal (defined as \code{list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)}. see details. 
  } 
  
 \item{method}{ The method to be used for estimating steady-state, will be passed as argument to the solver \link[rootSolve]{steady.1D}. The default is to use the "stode" method. Other options are to use "stodes", "runsteady" or "mix". The last option combines runsteady with stode, i.e. the model is first dynamically run for the times specified in argument \code{times}, and then the final value from this run used as initial guess for the steady-state estimated using \code{stode}
  }
  
 \item{times}{ start and end time of the dynamic run - only if \code{method} is one of \code{"runsteady"} or \code{"mix"}. 
  }
 
 \item{verbose}{ If TRUE, will write progession to the screen .
  } 
 
 \item{calcpH}{ if TRUE, the pH will be calculated - note that this will 
  not include the effects of calcium carbonate precipitation or dissolution
  on pH.
  Note also that the pH can be estimated afterwards by running \link{FESDIApH}. See examples.
  }

  \item{... }{Any argument passed to the steady-state solver.
  }
}

\value{
 \code{FESDIAsolve} returns an object of class \code{FESDIAstd} or \code{PHDIAstd}, and of class \code{steady1D}, as generated by the solvers from R-package \code{rootSolve} (\link{steady.1D}[rootSolve]).

 It contains a.o. the elements:
 \itemize{
   \item \code{y}, with the state variables at steady-state (\code{FDET, SDET},\code{O2, NO3, NH3}, \code{PO4, FeP, CaP}, \code{DIC, Fe, FeOH3}, \code{H2S, SO4, CH4, ALK}).
   \item \code{O2flux, O2deepflux, NO3flux, NO3deepflux, ...}, 
      \code{FDETflux, FDETdeepflux, SDETflux, SDETdeepflux, FePdeepflux, CaPdeepflux, FeOH3deepflux, OrgCflux, OrgNflux, OrgPflux}, the sediment-water and burial fluxes, in nmol/cm2/d.
   \item \code{DINDIPflux, DINDIPmean, DINDIPdeep}, the dissolved nitrogen to phosphorus ratio of flux, sediment concentrations and deep (burial) concentration.
   \item \code{TotMin, TotOxic, TotDenit, TotFered TotBSR, TotMeth}, total mineralisation, total oxic mineralisation, denitrification, iron reduction, biological sulphate reduction, and methanogenesis, in nmol/cm2/d.
   \item \code{PartOxic, PartDenit, PartFered, PartBSR, PartMethano}, the fraction of mineralisation due to oxic, denitrification, iron reduction, biological sulphate reduction, and methanogenesis.
   \item \code{TotNitri, TotFeoxid, TotH2Soxid, TotCH4oxid, TotAOM, TotFePprod, TotCaPprod, TotFePdesorp, TotCaPdiss, TotNprod, TotPprod, TotNH3ads}, integrated rates, nmol/cm2/d.
   \item \code{PartPremoved, PartNremoved},  the total P and N removed, relative to its production.
   \item \code{TOC}, the Total organic carbon concentration profile, \%.
   \item \code{Cprod, Nprod, Pprod, Oxicmin, Denitrific, Feredmin, BSRmin, Methmin, nitri,Feoxid, H2Soxid, CH4oxid, AOM, FeSprod, FePadsorp, CaPprod}, rate profiles, nm/cm3 liquid/d.
   \item \code{FePdesorp, CaPdiss}, rate profiles, nm/cm3 solid/d.
  }
  
 The meaning and units of these columns can be assessed via the R-functions:
 
 \code{FESDIAsvar()}, \code{FESDIA1D()}, \code{FESDIA0D()}. See \link{FESDIA0D}.

  
}

\author{
  Karline Soetaert
}
\details{
  To solve the model, a steady-state solver from package rootSolve
  (here we used \code{\link[rootSolve]{steady.1D}}) is
  used. 

}


\examples{

#===========================================
# Show parameter values
#===========================================

  FESDIAparms()

#===========================================
# Runs with different carbon fluxes
#===========================================

  out  <- FESDIAsolve(parms = list(Cflux = 2*1e5/12/365), calcpH = TRUE)
  out2 <- FESDIAsolve(parms = list(Cflux = 20*1e5/12/365), calcpH = TRUE)
  out3 <- FESDIAsolve(parms = list(Cflux = 200*1e5/12/365), calcpH = TRUE)

  par(mar = c(3,3,3,2))
  plot(out, out2, out3, ylim = c(20, 0), mfrow = c(4, 4), 
    which = c(1:14,16), lwd = 2, lty = 1)

  plot(out, out2, out3, which = "pH", ylim = c(20,0))

#===========================================
# Runs that are difficult to solve
#===========================================

# simple run
  out1 <- FESDIAsolve(parms = c(por0 = 0.5, MPBprod = 0))

# simple run used as initial condition for second run
  out2 <- FESDIAsolve(parms = c(por0 = 0.5, MPBprod = 1e3), yini = out1$y)

# second run used as initial condition for third run
# use mixed method: first dynamic run, then steady-state solver
  out3 <- FESDIAsolve(parms = c(por0 = 0.5, MPBprod = 1e4), yini = out2$y, method = "mixed")

  out4 <- FESDIAsolve(parms = c(por0 = 0.5, MPBprod = 5e4), yini = out3$y, 
     method = "runsteady", times = c(0, 1e6))
  
  FESDIAbudgetO2(out1, out2, out3, out4, which = "Rates") 

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
