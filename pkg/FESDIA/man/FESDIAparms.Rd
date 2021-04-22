\name{FESDIAparms}
\alias{FESDIAparms}
\alias{PHDIAparms}
\alias{FESDIApor}
\alias{FESDIAbiot}
\alias{FESDIAirr}
\alias{FESDIAdepth}
\alias{FESDIAgrid}
\alias{FESDIAdx}

\title{
  Functions to retrieve parameters, porosity, depth and sediment grid, irrigation and bioturbation, for the FESDIA and PHDIA models}
  
\description{

  \code{PHDIAparms, FESDIAparms} retrieve the parameters PHDIA and FESDIA model solutions.

  \code{FESDIAdepth, FESDIAdx} retrieve the sediment depths and layer thicknesses of PHDIA or FESDIA model solutions.
  
  \code{FESDIAbiot, FESDIApor, FESDIAirr} retrieve the bioturbation, porosity, and irrigation profiles of PHDIA and FESDIA model solutions.
  
}

\usage{
  FESDIAparms(out = NULL, as.vector = FALSE, which = NULL)
  PHDIAparms(out = NULL, as.vector = FALSE, which = NULL)
  
  FESDIAdepth(out = NULL)
  FESDIAgrid(out = NULL)
  FESDIAbiot(out)
  FESDIApor(out)
  FESDIAirr(out)
}

\arguments{

  \item{out }{an output object returned by \link{FESDIAsolve} or \link{FESDIAdyna}. If \code{NULL}, \code{FESDOIAparms} will return the default (parameter) values.
  }
  \item{as.vector }{if \code{TRUE} will return the parameter vector, else a data.frame that also contains the units.
  }
  \item{which }{if not \code{NULL}, a vector with names of the parameters to return.
  }

}

\value{
 \code{FESDIA0D} and \code{FESDIA1D} return the output variables of the solution as a vector or data.frame. 
 For dynamic runs, the output is averaged over the mean of the run.
 
 \code{FESDIA1D} always returns the sediment depth and the porosity as the first two columns.
}

\details{
  The parameters and their meaning are the following (with default values):
  \itemize{
    \item Cflux     = 20*1e5/12/365  ,  nmolC/cm2/d  - Carbon deposition
    \item rFast     = 25/365   , /d         , decay rate fast decay detritus
    \item rSlow     = 0.5/365  , /d         , decay rate slow decay detritus
    \item pFast     = 0.9      , -          , fraction fast detritus in flux
    \item NCrFdet   = 0.16     , molN/molC  , NC ratio fast decay detritus
    \item NCrSdet   = 0.13     , molN/molC  , NC ratio slow decay detritus
    \item PCrFdet   = 9.75e-03 , molP/molC  , PC ratio fast decay det.
    \item PCrSdet   = 9.75e-03 , molP/molC  , PC ratio slow decay det.
    \item FeOH3flux = 1        , nmol/cm2/d  , FeOH3 deposition rate
    \item CaPflux   = 0        , nmolP/cm2/d , deposition rate of CaP

    \item O2bw      = 300      , mmol/m3     Oxygen conc in bottom water
    \item NO3bw     = 10       , mmol/m3     Nitrate
    \item NH3bw     = 1        , mmol/m3     Ammonium
    \item CH4bw     = 0        , mmol/m3     Methane
    \item PO4bw     = 0.5      , mmol/m3     Phoshpate
    \item DICbw     = 2200     , mmol/m3     dissolved inorganic carbon
    \item Febw      = 0        , mmol/m3     dissolved iron
    \item H2Sbw     = 0        , mmol/m3     sulphide
    \item SO4bw     = 30000    , mmol/m3     sulphate

    \item w         = 0.1/365000, cm/D      , advection rate
    \item mixL      = 5        , cm         , the depth of bioturbation layer
    \item biot      = 1        , cm2/yr     , the bioturbation rate
    \item irr       = 0        , /d         , the irrigation rate
    \item gasflux   = 0        , cm/d       , piston velocity for dry flats - exchange of O2 and DIC only+deposition 
  
    \item MPBprod   =     0    , mmol/m3/d  , maximal rate of picrophytobenthos production - range: 5000-5e4
    \item kMPB      =     4    , /cm        , exponential decay 
    \item ksDIN     =  0.01    , mmol/m3    , N limitation constant
    \item ksPO4     = 0.001    , mmol/m3    , P limitation constant
    \item ksDIC     =     1    , mmol/m3    , C limitation constant
  
    \item NH3Ads    = 1.3      , -          , Adsorption coeff ammonium
    \item rnit      = 20.      , /d         , Max nitrification rate
    \item ksO2nitri = 1.       , umolO2/m3  , half-sat O2 in nitrification
    \item ksO2oxic  = 3.       , mmolO2/m3  , half-sat O2 in oxic mineralisation
    \item ksNO3denit= 30.      , mmolNO3/m3 , half-sat NO3 in denitrification
    \item kinO2denit= 1.       , mmolO2/m3  , half-sat O2 inhib denitrification
    \item kinNO3anox= 1.       , mmolNO3/m3 , half-sat NO3 inhib anoxic degr
    \item kinO2anox = 1.       , mmolO2/m3  , half-sat O2 inhib anoxic min
    \item temperature = 10     , dgC  - for estimation of diffusion coefficients
    \item salinity    = 35     ,
    \item TOC0        = 0.5    , the background C concentration, %
    \item rFePadsorp  = 1e-5  , /d, FeP adsorption rate
    \item rCaPprod    = 0      , /d, CaP production rate
    \item rCaPdiss    = 0      , /d, CaP dissolution rate
    \item CPrCaP      = 0.2869565 , Ca:P ratio (mol/mol) - Ca10(PO4)4.6(CO3)1.32F1.87(OH)1.45
    \item ksFeOH3     = 12500.   , mmolFeOH3/m3  half-sat FeOH3 in iron red  
    \item kinFeOH3    = 12500.   , mmolFeOH3/m3  half-sat FeOH3 inhib BSR
    \item ksSO4BSR    = 1600.    , mmolSO4/m3  half-sat SO4 in sulfate reduction
    \item kinSO4Met   = 1000     , mmolSO4/m3, half-sat SO4 inhibition for methanogenesis
    \item rFeox       = 0.3      , /d/nmol/cm3  oxidation constant for iron by O2 (bimolecular rate law)
    \item rH2Sox      = 5e-4     , /d/nmol/cm3  oxidation constant for diss Sulfide by O2 (bimolecular rate law)
    \item rFeS        = 1e-3     , /d/nmol/cm3  oxidation constant for diss Sulfide by O2 (bimolecular rate law)
    \item rCH4ox      = 27       , /d/nmol/cm3  oxidation constant for CH4 by O2 (bimolecular rate law)
    \item rAOM        = 3e-5     , /d/nmol/cm3  oxidation constant for AOM CH4 by SO4 (bimolecular rate law)
    \item por0        = 0.9      , -  surface porosity
    \item pordeep     = 0.5      , -  deep porosity
    \item porcoeff    = 0.3      , cm porosity coefficient
    \item gridtype    = 1        , 1 = cartesian, 2 = cylindrical, 3 = spherical
  }


}

\author{
  Karline Soetaert
}


\examples{

# default parameters
  defparms  <- FESDIAparms(as.vector = TRUE)
  defparms
  
# Some runs to work with  
  defsteady <- FESDIAsolve()
  defdyn    <- FESDIAdyna()

# altered steady-state run
  out <- FESDIAdyna(parms = list(Cflux = 10), CfluxForc = list(amp = 1))
  cbind(default = defparms, out = FESDIAparms(out))
  
# grid used for outputs

  pm <- par(mfrow = c(2, 2))
  plot(FESDIApor(out), FESDIAdepth(out), ylim = c(10,0), 
       type = "l", ylab = "cm", xlab = "-", main = "porosity")
  plot(FESDIAbiot(out), FESDIAdepth(out), ylim = c(10,0), 
       type = "l", ylab = "cm", xlab = "cm2/d", main = "bioturbation")
  image(out, which = "NH3", grid = FESDIAdepth(out), ylim = c(10,0), 
        main = "NH3", mfrow = NULL, legend = TRUE, ylab = "cm")
  matplot.1D(out, which = "NH3", xyswap = TRUE, grid = FESDIAdepth(out), 
        type = "l", col = "grey", ylim = c(10,0), mfrow = NULL, 
        ylab = "cm", xlab = "mmol/m3")      
        
  par(mfrow = pm)

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
