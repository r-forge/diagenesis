\name{FESDIAdyna}
\alias{FESDIAdyna}
\title{
  Dynamic solution of the FESDIA model, calculating C, N, P, O2, Fe and S dynamics in the sediment.
}
\description{
  \code{FESDIAdyna} dynamically runs the FESDIA model; pH can be calculated afterwards, ie ignoring carbonate dynamics.
}

\usage{
  FESDIAdyna (parms = list(), times = 0:365, spinup = NULL, yini = NULL, 
     gridtype = 1, Grid = NULL, porosity = NULL, bioturbation = NULL, 
     irrigation = NULL, surface = NULL, diffusionfactor = NULL, 
     dynamicbottomwater = FALSE,
     CfluxForc = NULL,FeOH3fluxForc = NULL, CaPfluxForc = NULL,  
     O2bwForc   = NULL,  NO3bwForc  = NULL,  NO2bwForc  = NULL, 
     NH3bwForc  = NULL,  FebwForc   = NULL,  H2SbwForc  = NULL,  
     SO4bwForc  = NULL,  CH4bwForc  = NULL,  PO4bwForc  = NULL,  
     DICbwForc  = NULL,  ALKbwForc  = NULL,  wForc      = NULL,  
     biotForc   = NULL,  irrForc    = NULL,  rFastForc  = NULL,  
     rSlowForc  = NULL,  pFastForc  = NULL,  MPBprodForc= NULL,  
     gasfluxForc = NULL, HwaterForc = NULL,  ratefactor = NULL, 
     calcpH = FALSE, verbose = FALSE, ...) 

}
\arguments{

  \item{parms}{ A list with parameter values.Available parameters can be listed using function \link{FESDIAparms}.
    
  See details of \link{FESDIAparms}.
  }
 
 \item{times}{ Output times for the dynamic run
  }

 \item{spinup}{ Spinput times for the dynamic run; not used for output; the outputted simulation starts from the final values of the spinup run.
  }

 \item{CfluxForc, FeOH3fluxForc, CaPfluxForc}{ \code{NULL} or a list, detailing the forcing function for the deposition fluxes of organic carbon, FeP and CaP.     
    If \code{NULL} then the corresponding parameter value (Cflux, FeOH3flux, CaPflux) is used. 
    If a \code{list}, it should contain either a data time series (\code{list(data = )}) or parameters determining the periodicity of the seasonal signal (defined as \code{list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)}. see details. 
  }

 \item{O2bwForc, NO3bwForc, NO2bwForc, NH3bwForc, CH4bwForc, FebwForc, SO4bwForc, H2SbwForc, PO4bwForc, DICbwForc, ALKbwForc}{ \code{NULL} or a list, detailing the forcing function for the boundary concentrations. 
    If \code{NULL} then the corresponding parameter value (O2bw, NO3bw, NO2bw, NH3bw, CH4bw, Febw, SO4bw, H2Sbw, PO4bw, DICbw, ALKbw) is used. 
    If a \code{list}, it should contain either a data time series (\code{list(data = )}) or parameters determining the periodicity of the seasonal signal (defined as \code{list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)}. see details. 
  }  

 \item{wForc, biotForc, irrForc}{ \code{NULL} or a list, detailing the forcing function forthe advection, bioturbation and irrigation rates (units of cm/d, cm2/d and /d).     
     If \code{NULL} then the corresponding parameter value (w, biot, irr) is used. 
    If a \code{list}, it should contain either a data time series (\code{list(data = )}) or parameters determining the periodicity of the seasonal signal (defined as   \code{list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)}. see details. 
 }  

 \item{rFastForc, rSlowForc, pFastForc}{ \code{NULL} or a list, detailing the forcing function for the decay rate of fast and slow detritus, and the fraction of fast detritus in organic flux (units of /d, /d, -). 
    If \code{NULL} then the corresponding parameter value (\code{rFast, rSlow, pFast}) is used. 
    If a \code{list}, it should contain either a data time series (\code{list(data = )}) 
    or parameters determining the periodicity of the seasonal signal (defined as \code{list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)}. see details. 
  } 
  
 \item{MPBprodForc}{ \code{NULL} or a list, detailing the forcing function for the maximal microphytobenthos production rate (units of mmolO2/m3/d). 
    If \code{NULL} then the corresponding parameter value (\code{MPBprod}) is used. 
    If a \code{list}, it should contain either a data time series (\code{list(data = )}) or parameters determining the periodicity of the seasonal signal (defined as \code{list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)}. see details. 
  }  

 \item{gasfluxForc}{ \code{NULL} or a list, detailing the forcing function for the intensity of exchange with the air (units of cm/d). 
    If \code{NULL} then the corresponding parameter value (gasflux) is used. 
    If a \code{list}, it should contain either a data time series (\code{list(data = )}) or parameters determining the periodicity of the seasonal signal (defined as \code{list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)}. see details. 
    A gasflux>0 represents exchange rate for O2 and DIC, not for other dissolved substances. There can still be deposition of C, FeP and CaP.
  }  

 \item{HwaterForc}{ \code{NULL} or a list, detailing the forcing function for the height of the water above the sediment - if \code{dynamicbottomwater} = TRUE. 
    If \code{NULL} then the corresponding parameter value (Hwater) is used. 
    If a \code{list}, it should contain either a data time series (\code{list(data = )}) or parameters determining the periodicity of the seasonal signal (defined as \code{list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)}. see details. 
  }  
 
 \item{ratefactor}{ \code{NULL} or a list, detailing the forcing function for the  biogeochemical rate multiplication factor. If not specified (or \code{NULL}), then it is assumed to be 1 and constant. 
    If a \code{list}, it should contain either a data time series (\code{list(data = )}) or parameters determining the periodicity of the seasonal signal (defined as \code{list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)}. see details. 
  }  

 \item{dynamicbottomwater}{ If \code{TRUE}, then the concentrations in the water overlying the sediment will also be dynamically described, and with water height equal to \code{Hwater}. Note that this will slow down the simulation.   
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
   Note that the porosity values should be consistent with the \code{Grid} - and should be inbetween 0 and 1.
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

 \item{yini}{ The condition at which to inialise the dynamic simulation, i.e. a vector or matrix, with the values of 
   \code{FDET, SDET}, \code{O2, NO3, NO2}, \code{NH3, DIC, Fe}, 
   \code{FeOH3, H2S, SO4}, \code{CH4, PO4, FeP, CaP, Pads, ALK}, 
   each in 100 layers, and in that order. 
   If \code{NULL}, the default, then the steady-state solution based on the mean forcing is used (and obtained with \link{FESDIAsolve}) .
  } 

 \item{calcpH}{ if TRUE, the pH will be calculated - note that this will 
  not include the effects of calcium carbonate precipitation or dissolution
  on pH . 
  Note also that the pH can be estimated afterwards by running \link{FESDIApH}. See examples.
  }

 \item{verbose}{ If TRUE, will write progession to the screen .
  } 

  \item{...}{ Any argument passed to the dynamic solver (\link{ode.1D}[deSolve])
  }
}

\value{
 A matrix of class \code{FESDIAdyn} or \code{PHDIAdyn} and \code{deSolve}, as generated by the solver from R-package \code{deSolve} (\link{ode.1D}).
 
 It contains several output columns, the first is \code{time}. The meaning and units of these columns can be assessed via the R-functions:
 
 \code{FESDIAsvar()}, \code{FESDIA1D()}, \code{FESDIA0D()}. See \link{FESDIA0D}.

}

\author{
  Karline Soetaert
}
\details{
  Several parameters can also be described as forcing functions. They are, for \code{FESDIAdyna}:
  \code{Cflux, FeOH3flux, CaPflux}, \code{O2bw, NO3bw, NO2bw, NH3bw}, \code{Febw, H2Sbw, 
        SO4bw, CH4bw}, \code{PO4bw, DICbw, ALKbw, w, biot, irr}, 
        \code{rFast, rSlow, pFast}, \code{MPBprod, gasflux, Hwater}.
             
  The forcing functions are prescribed as a \code{list} that either contains a data series or specifies a periodic signal.
  The list is defined as: \code{list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)}
  
  \itemize{
  \item Forcing functions imposed by a data series are set with item \code{data}, and contains a time series to replace the parameter. \code{data} should be a two-columned matrix with times (first column) and values (second column). The values should be in the same units as the parameters.

   The time series should embrace both arguments \code{times} and \code{spinup}.  
  
  \item if a periodic signal, the list should contain \code{amp, period, phase, pow} and \code{min}: parameters determining the periodicity of the seasonal signal, and in the same units as the parameters. 
   From this the forcing function time series is estimated, e.g. for \code{CfluxForc} as follows: CfluxForc = \code{max(min, Cflux*(1 + (amp*sin((times-phase)/period*2*pi))^pow)}, where \code{Cflux} is the parameter value. Used only if \code{data} is \code{NULL}.   If \code{amp} = 0, or \code{pow} = 0, then the forcing will be kept constant and equal to the parameter value. 
  }
}

\note{
  The model application starts by estimating the steady-state condition
  of the model. This steady-state condition is then used as a starting
  condition for a dynamic simulation, with sinusoidal forcing
  of organic carbon deposition flux (truncated to be >= 0).

  If only one year is run, the dynamic run is likely not at equilibrium 
  with the transient forcing.

}


\examples{

#===========================================
# Show parameter values
#===========================================

  FESDIAparms()

#===========================================
# Default run
#===========================================

  out <- FESDIAdyna()
  image2D(out, ylim = c(20, 0), which = 1:9)
  FESDIAbudgetO2(out)  

#===========================================
# Altered parameter values
#===========================================

  out2 <- FESDIAdyna(parms = list(pFast = 0.95, Cflux = 500), 
    spinup = 0:365)
  
  image2D(out2, ylim = c(20, 0), which = c("TOC", "O2", "NO3", "PO4"))

#===========================================
# A forcing function flux timeseries, 
# and including pH
#===========================================

  Cflux <- cbind (time = c(0, 100,  150,  175, 200, 250, 365),
                  flux = c(1, 1,    1000, 800, 1200, 10, 1)) 
  
  out3 <- FESDIAdyna(parms = list(pFast = 0.9), 
                     Cflux = list(data = Cflux), times = 0:365,
                     calcpH = TRUE)
                     
  image2D(out3, which = c("TOC", "O2", "NO3", "PO4", "ALK", "pH"), 
       ylim = c(15, 0))


# same, but more elaborate:
  out4 <- FESDIAdyna(parms = list(pFast = 0.9), 
                     Cflux = list(data = Cflux), times = 0:365)

  pH <- FESDIApH(out4)
  
  plot3D::image2D(pH, ylim = c(15,0), y = FESDIAdepth(out4), 
                  x = out4[,"time"], main = "pH")

#===========================================
# a simple forcing function 
# (varying the part of deposition in fast decay C)
#===========================================

  out4 <- FESDIAdyna(parms = list(pFast = 0.9), 
    pFastForc = list(amp = 0.8), times = 0:365)

  image2D(out4, which = 1:9, ylim = c(20, 0))

  plot(out4, which = c("OrgCflux", "O2flux"), mfrow = NULL)  
  head(FESDIA0D(out4))
  FESDIAbudgetC(out4)

#===========================================
# long-distance reactions (LDR)
#===========================================

# first without LDR

  out <- FESDIAsolve(parms = list(Cflux = 100*1e5/12/365), 
                     calcpH = TRUE)
  dyn <- FESDIAdyna(parms = list(Cflux = 100*1e5/12/365), 
                             yini = out$y, calcpH = TRUE,
                             times = seq(0,5, length.out = 30))
  matplot1D(dyn, which = c("O2", "H2S", "ALK", "pH"), 
     ylim = list(c(1,0), c(10,0), c(10,0)))

# then including LDR

  out2 <- FESDIAsolve(parms = list(Cflux = 100*1e5/12/365, 
                               rSurfH2Sox = 5, ODUoxdepth = 5), 
                      yini = out$y, calcpH = TRUE, method = "mixed")
  plot(out, out2, which = c("O2", "H2S", "ALK", "pH"), 
     ylim = list(c(1,0), c(10,0), c(10,0), c(10,0)), mfrow = c(2,4))
     

  dyn <- FESDIAdyna(parms = list(Cflux = 100*1e5/12/365,  
                             rSurfH2Sox = 5, ODUoxdepth = 5), 
                    yini = out$y, calcpH = TRUE, 
                    times = seq(0, 10, length.out = 100))

  matplot1D(dyn, which = c("O2","H2S","ALK", "pH"), 
     ylim = list(c(1,0), c(10,0), c(10,0), c(10,0)))

  matplot1D(dyn, which = "pH",  ylim = c(1,0), col = jet.col(100))

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
