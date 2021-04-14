\name{CNPDIAsolve}
\alias{CNPDIAsolve}
\title{
  Steady-state solution of C, N, P and O2 dynamics in the sediment.
}
\description{

  \code{CNPDIAsolve} finds the steady-state solution of the CNPDIA model.
}
\usage{
  CNPDIAsolve (parms = list(), yini = NULL, gridtype = 1, Grid = NULL, 
     porosity = NULL, bioturbation = NULL, irrigation = NULL, surface = NULL, 
     diffusionfactor = NULL, dynamicbottomwater = FALSE, 
     ratefactor = NULL, verbose = FALSE, ...)
  }
\arguments{

  \item{parms}{ A list with parameter values, Available parameters can be listed using function \link{CNPDIAparms}.

    See details.
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
   If \code{NULL}, it is defined by the parameters \code{biot}, \code{mixdepth} and \code{mixatt} as: 
   \code{biot * exp(-pmax(0, (x.int-mixL))/mixatt)}, where \code{x.int} is the distance, from the surface of the layer interface. 
   Note that the bioturbation values should be consistent withe the \code{Grid}.
  }

 \item{irrigation}{ If specified, either an object with irrigation rates (units [/d]) as returned by \code{setup.prop.1D} from the package \code{ReacTran}, 
   a vector of length 100 with the irrigation rates defined at the layer centres, or one number with the constant rates. 
   If \code{NULL}, it is defined by the parameters \code{irr}, \code{mixdepth} and \code{mixatt} as: 
   \code{irr * exp(-pmax(0, (x-mixdepth)/mixxatt))}, where \code{x} is the distance, from the surface, of the layer centres. 
   Note that the irrigation values should be consistent withe the \code{Grid}.
  }

 \item{surface}{ If specified, either an object with surface areas (units [cm2]) as returned by \code{setup.prop.1D} from the package \code{ReacTran}, 
   a vector of length 101 with the surface areas defined at the layer interfaces, or one number with the constant surface area. 
   If \code{NULL}, it is defined by the parameter \code{gridtype}, and the \code{Grid} as: 
   \code{surface = rep(1, 101)} for \code{gridtype == 1}, \code{surface = rev(2*pi*Grid$x.int)} for \code{gridtype == 2}  and \code{surface = rev(pi*(Grid$x.int)^2)} for \code{gridtype == 3}.

   Note that the surface values should be consistent withe the \code{Grid}.
  }

  
 \item{diffusionfactor}{ The multiplication factor necessary to go from molecular diffusion to effective sediment diffusion, i.e. that takes into account tortuosity.
   If specified, either an object with these factors ([-]) as returned by \code{setup.prop.1D} from the package \code{ReacTran}, 
   a vector of length 101 with these factors defined at the layer interfaces, or one number with the constant factor. 
   If \code{NULL}, it is set equal to the porosity. Note that the factors should be consistent withe the \code{Grid}.
  }

  \item{yini }{Initial guess of the steady-state solution.
  }

 \item{dynamicbottomwater}{ If \code{TRUE}, then the concentrations in the water overlying the sediment will also be dynamically described, and with water height equal to \code{Hwater}. Note that this will slow down the simulation.   
  }
 \item{ratefactor}{ \code{NULL} or a list, detailing the forcing function for the  biogeochemical rate multiplication factor. If not specified (or \code{NULL}), then it is assumed to be 1 and constant. 
    If a \code{list}, it should contain either a data time series (\code{list(data = )}) or parameters determining the periodicity of the seasonal signal (defined as \code{list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)}. see details. 
  }  
 \item{verbose}{ If TRUE, will write progession to the screen .
  }
  \item{... }{Any argument passed to the steady-state solver.
  }
}

\value{
 \code{CNPDIAsolve} returns an object of class \code{FESDIAstd}, and of class \code{steady1D}, as generated by the solvers from R-package \code{rootSolve} (\link{steady.1D}[rootSolve]).

 It contains, a.o. the elements:
 \itemize{
   \item \code{y}, with the state variables at steady-state (\code{FDET, SDET, O2, NO3, NO2, NH3, ODU, PO4, FeP, CaP, DIC}).
   \item \code{O2flux, O2deepflux, NO3flux, NO3deepflux, NO2flux, NO2deepflux, NH3flux, NH3deepflux, 
      ODUflux, ODUdeepflux, PO4flux, PO4deepflux, DICflux, DICdeepflux, 
      FDETflux, FDETdeepflux, SDETflux, SDETdeepflux, FePdeepflux, CaPdeepflux, 
      OrgCflux, OrgNflux, OrgPflux}, the sediment-water and burial fluxes, in nmol/cm2/d.
   \item \code{DINDIPflux, DINDIPmean, DINDIPdeep}, the dissolved nitrogen to phosphorus ratio of flux, sediment concentrations and deep (burial) concentration.
   \item \code{TotMin, TotOxic, TotDenit, TotAnoxic}, total mineralisation, total oxic mineralisation, denitrification and anoxic mineralisation, in nmol/cm2/d.
   \item \code{PartOxic, PartDenit, PartAnoxic}, the fraction of mineralisation due to oxic, denitrification and anoxic mineralisation.
   \item \code{TotNitri, TotODUoxid, TotFePprod, TotCaPprod, TotFePdesorp, TotCaPdiss, TotNprod, TotPprod, TotNH3ads}, integrated rates, nmol/cm2/d.
   \item \code{PartPremoved, PartNremoved},  the total P and N removed, relative to its production.
   \item \code{TOC}, the Total organic carbon concentration profile, \%.
   \item \code{Cprod,Nprod,Pprod,Oxicmin,Denitrific,anoxicmin,nitri,oduox,odudepo,FePadsorp,CaPprod}, rate profiles, nm/cm3 liquid/d.
   \item \code{FePdesorp,CaPdiss}, rate profiles, nm/cm3 solid/d.
  }
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

  CNPDIAparms()
  par(mar = c(3,3,3,3))

#===========================================
# Runs with different carbon fluxes
#===========================================

  out  <- CNPDIAsolve()
  out2 <- CNPDIAsolve(parms = list(Cflux = 200*1e5/12/365))
  out3 <- CNPDIAsolve(parms = list(Cflux = 2*1e5/12/365))
  plot(out, out2, out3, xyswap = TRUE, grid = CNPDIAdepth(out), 
    ylim = c(20,0), mfrow = c(3, 4))

  CNPDIAbudgetO2(out)  
  
#===========================================
# long-distance reactions
#===========================================

  out2 <- CNPDIAsolve(parms = list(Cflux = 100*1e5/12/365))
  out2b <- CNPDIAsolve(parms = list(Cflux = 100*1e5/12/365, 
               rSurfODUox = 10, ODUoxdepth = 1), yini = out2$y)
  out2c <- CNPDIAsolve(parms = list(Cflux = 100*1e5/12/365,  
               rSurfODUox = 10, ODUoxdepth = 5), yini = out2b$y)
  plot(out2, out2b, out2c, xyswap = TRUE, grid = CNPDIAdepth(out2), 
     which = c("O2","ODU"), ylim = list(c(1,0), c(10,0)), 
     log = list("","x"))

  plot(out2b, out2c, xyswap = TRUE, grid = CNPDIAdepth(out2), 
     which = c("O2","ODU"), ylim = list(c(1,0), c(10,0)))

#===========================================
# Dynamic bottom water concentrations
#===========================================

  std <- CNPDIAsolve(parms = list(Cflux = 20*1e5/12/365))
  BWs <- CNPDIAparms(std, as.vector = TRUE)[c("O2bw", "NO3bw", 
       "NO2bw", "NH3bw", "ODUbw", "DICbw", "PO4bw")]
  yini <- rbind(c(FDETbw = 0, SDETbw = 0, BWs, FePbw = 0, 
        CaPbw = 0, Padsbw=0), std$y)
  out2 <- CNPDIAdyna(parms = list(Cflux = 20*1e5/12/365, Hwater = 50), 
      dynamicbottomwater = TRUE, times = seq(0, 1, by = 0.01), 
      yini = yini)
  
  plot(out2, which = c("O2bw", "NO3bw", "O2flux", "NO3flux"))
  image(out2, which = c("O2","NO3","NH3","DIC"), legend = TRUE)
\dontrun{
# 
  }
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
