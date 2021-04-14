\name{CNPDIAparms}
\alias{CNPDIAparms}
\alias{MPBDIAparms}
\alias{CNPDIApor}
\alias{CNPDIAbiot}
\alias{CNPDIAirr}
\alias{CNPDIAdx}
\alias{CNPDIAdepth}
\alias{CNPDIAgrid}
\title{
  Functions to retrieve parameters, porosity, depth and sediment grid, irrigation and bioturbaion, for the CNPDIA and MPBDIA model.
}
\description{

  \code{CNPDIAparms, CNPDIAdepth, CNPDIAdx} retrieve the parameters, sediment depths and layer thicknesses of CNPDIA model solutions.
  
  \code{CNPDIAbiot, CNPDIApor, CNPDIAirr} retrieve the bioturbation, porosity, and irrigation profiles of CNPDIA model solutions.

  \code{MPBDIAparms} retrieves the parameters of the MPBDIA model
}
\usage{
  CNPDIAparms(out = NULL, as.vector = FALSE, which = NULL)
  MPBDIAparms(out = NULL, as.vector = FALSE, which = NULL)
  CNPDIAdepth(out = NULL)
  CNPDIAgrid(out = NULL)
  CNPDIAdx(out = NULL)
  CNPDIApor(out)
  CNPDIAbiot(out)
  CNPDIAirr(out)
}
\arguments{

  \item{out }{an output object returned by \link{CNPDIAsolve}, \link{MPBDIAsolve},
   or \link{CNPDIAdyna}, \link{MPBDIAdyna}. 
   If \code{NULL}, \code{CNPDIAparms} or \code{MPBDIAparms} will return the default (parameter) values.
  }
  \item{as.vector }{if \code{TRUE} will return the parameter \code{vector}, else a \code{data.frame} that also contains the units.
  }
  \item{which }{if not \code{NULL}, a vector with names of the variables/parameters to return.
  }
}

\author{
  Karline Soetaert
}

\details{
  For the CNPDIA model, the parameters and their meaning are the following (with default values):
  \itemize{
    \item Cflux       ,  5.00e+02 ,  nmolC/cm2/d    ,                total organic C deposition
    \item pFast       ,  9.00e-01 ,            -    ,                  part FDET in carbon flux
    \item FePflux     ,  0.00e+00 ,  nmolP/cm2/d    ,                    deposition rate of FeP
    \item CaPflux     ,  0.00e+00 ,  nmolP/cm2/d    ,                    deposition rate of CaP
    \item rFast       ,  6.85e-02 ,           /d    ,                           decay rate FDET
    \item rSlow       ,  1.37e-04  ,          /d    ,                           decay rate SDET
    \item NCrFdet     ,  1.51e-01 ,    molN/molC    ,                             NC ratio FDET
    \item NCrSdet     ,  1.51e-01 ,    molN/molC    ,                             NC ratio SDET
    \item PCrFdet     ,  9.43e-03 ,    molP/molC    ,                             PC ratio FDET
    \item PCrSdet     ,  9.43e-03 ,    molP/molC    ,                             PC ratio SDET
    \item BCupLiq     ,  2.00e+00 ,            -  , upper boundary liq. 1:flux, 2:conc, 3:0-grad
    \item BCdownLiq   ,  3.00e+00 ,            -  , lower boundary liq. 1:flux, 2:conc, 3:0-grad
    \item O2bw        ,  3.00e+02 ,      mmol/m3    , upper boundary O2  -if BC=1: flux, 2:conc
    \item NO3bw       ,  1.00e+01 ,      mmol/m3    , upper boundary NO3 -if BC=1: flux, 2:conc
    \item NO2bw       ,  0.00e+00 ,      mmol/m3    , upper boundary NO2 -if BC=1: flux, 2:conc
    \item NH3bw       ,  1.00e+00 ,      mmol/m3    , upper boundary NH3 -if BC=1: flux, 2:conc
    \item ODUbw       ,  0.00e+00 ,      mmol/m3    , upper boundary ODU -if BC=1: flux, 2:conc
    \item PO4bw       ,  5.00e-01 ,      mmol/m3    , upper boundary PO4 -if BC=1: flux, 2:conc
    \item DICbw       ,  2.20e+03 ,      mmol/m3    , upper boundary DIC -if BC=1: flux, 2:conc
    \item O2dw        ,        NA ,      mmol/m3    , lower boundary O2  -if BC=1: flux, 2:conc
    \item NO3dw       ,        NA ,      mmol/m3    , lower boundary NO3 -if BC=1: flux, 2:conc
    \item NO2dw       ,        NA ,      mmol/m3    , lower boundary NO2 -if BC=1: flux, 2:conc
    \item NH3dw       ,        NA ,      mmol/m3    , lower boundary NH3 -if BC=1: flux, 2:conc
    \item ODUdw       ,        NA ,      mmol/m3    , lower boundary ODU -if BC=1: flux, 2:conc
    \item PO4dw       ,        NA ,      mmol/m3    , lower boundary PO4 -if BC=1: flux, 2:conc
    \item DICdw       ,        NA ,      mmol/m3    , lower boundary DIC -if BC=1: flux, 2:conc
    \item w           ,  2.74e-07 ,         cm/d    ,                            advection rate
    \item biot        ,  2.74e-03 ,        cm2/d    ,                  bioturbation coefficient
    \item biotdepth   ,  5.00e+00 ,           cm    ,                      depth of mixed layer
    \item biotatt     ,  1.00e+00 ,           cm    ,         attenuation coeff below biotdepth
    \item irr         ,  0.00e+00 ,           /d    ,                       bio-irrigation rate
    \item irrdepth    ,  5.00e+00 ,           cm    ,                  depth of irrigated layer
    \item irratt      ,  1.00e+00 ,           cm    ,          attenuation coeff below irrdepth
    \item gasflux     ,  0.00e+00 ,         cm/d    ,             piston velocity for dry flats
    \item NH3Ads      ,  1.30e+00 ,            -    ,                 Adsorption coeff ammonium
    \item rnitri1     ,  2.00e+01 ,           /d    ,      Max nitrification rate step1 (NH3ox)
    \item rnitri2     ,  2.00e+01 ,           /d    ,      Max nitrification rate step2 (NO2ox)
    \item ksO2nitri   ,  1.00e+00 ,    mmolO2/m3    ,              half-sat O2 in nitrification
    \item ranammox    ,  1.00e-01 ,           /d    ,                              Anammox rate
    \item ksNO2anammox,  1.00e-01 ,     mmolN/m3    ,                   half-sat NO2 in anammox
    \item rODUox      ,  2.00e+01 ,           /d    ,       Max rate ODU oxidation in one layer
    \item rSurfODUox  ,  0.00e+00 ,           /d    ,         Max rate ODU oxidation with BW O2
    \item ODUoxdepth  ,  5.00e+00 ,           cm    ,        Max depth ODU oxidation with BW O2
    \item ODUoxatt    ,  1.00e+00 ,           cm    ,           depth attenuation ODU oxidation
    \item ksO2oduox   ,  1.00e+00 ,    mmolO2/m3    ,           half-sat O2 in oxidation of ODU
    \item ksO2oxic    ,  3.00e+00 ,    mmolO2/m3    ,        half-sat O2 in oxic mineralisation
    \item ksNO3denit  ,  3.00e+01 ,   mmolNO3/m3    ,           half-sat NO3 in denitrification
    \item kinO2denit  ,  1.00e+00 ,    mmolO2/m3    ,         half-sat O2 inhib denitrification
    \item kinNO3anox  ,  1.00e+00 ,   mmolNO3/m3    ,            half-sat NO3 inhib anoxic degr
    \item kinO2anox   ,  1.00e+00 ,    mmolO2/m3    ,              half-sat O2 inhib anoxic min
    \item pdepo       ,        NA ,            -    ,  part ODU prod lost (NA:estimated from w)
    \item rdepo       ,  0.00e+00 ,           /d    ,                          ODU removal rate
    \item temperature ,  1.00e+01 ,          dgC    ,                               temperature
    \item salinity    ,  3.50e+01 ,          psu    ,                                  salinity
    \item TOC0        ,  5.00e-01 ,            %    ,                    refractory Carbon conc
    \item rFePdesorp  ,  1.00e-02 ,           /d    ,                       rate FeP desorption
    \item rFePadsorp  ,  3.00e-01 ,           /d    ,                       rate FeP adsorption
    \item rCaPprod    ,  0.00e+00 ,           /d    ,                       rate CaP production
    \item rCaPdiss    ,  0.00e+00 ,           /d    ,                      rate CaP dissolution
    \item CPrCaP      ,  2.87e-01 ,    molC/molP    ,                           C:Pratio in CaP
    \item rPads       ,  0.00e+00 ,           /d    ,                       adsorption rate PO4
    \item rPdes       ,  0.00e+00 ,           /d    ,             desorption rate of adsorbed P
    \item maxPads     ,  1.00e+03 ,mmolP/m3solid    ,              Max adsorbed P concentration
    \item por0        ,  9.00e-01 ,            -    ,                          surface porosity
    \item pordeep     ,  5.00e-01 ,            -    ,                             deep porosity
    \item porcoeff    ,  3.00e-01 ,           cm    ,                porosity decay coefficient
    \item formationtype, 1.00e+00 ,            - , formationfactor, 1=sand,2=fine sand,3=general
    \item dilution    ,  0.00e+00 ,           /d    ,       relaxation towards background conc 
    \item Hwater      ,  1.00e+01 ,           cm    ,                 height of water over core
    \item Cfall       ,  1.00e+02 ,         cm/d    ,      fall speed of organic C (FDET, SDET)
    \item FePfall     ,  1.00e+02 ,         cm/d    ,                         fall speed of FeP
    \item CaPfall     ,  1.00e+02 ,         cm/d    ,                         fall speed of CaP
    \item MPBprod     ,  0.00e+00 ,    mmol/m3/d    ,               maximal MPB production rate
    \item kdSed       ,  2.00e+01 ,          /cm    ,              light extinction coefficient
    \item kNO3upt     ,  3.00e+00 ,     mmolN/m3    ,               NO3 limitation constant MPB
    \item kNH3upt     ,  3.00e+00 ,     mmolN/m3    ,               NH3 limitation constant MPB
    \item kPO4upt     ,  1.00e-01 ,     mmolP/m3    ,                 P limitation constant MPB
    \item kDICupt     ,  1.00e+00 ,     mmolC/m3    ,                 C limitation constant MPB
  }

}
\examples{

# defaults
  defparms  <- CNPDIAparms(as.vector = TRUE)
  defparms

# a run to work with  
  defsteady <- CNPDIAsolve()
  defdyn    <- CNPDIAdyna()

  out <- CNPDIAdyna(parms = list(Cflux = 1000), CfluxForc = list(amp = 0.8))

# grid used for outputs
  pm <- par(mfrow = c(2, 2))
  plot(CNPDIApor(out), CNPDIAdepth(out), ylim = c(10,0), 
       type = "l", ylab = "cm", xlab = "-", main = "porosity")
  plot(CNPDIAbiot(out), CNPDIAdepth(out), ylim = c(10,0), 
       type = "l", ylab = "cm", xlab = "cm2/d", main = "bioturbation")
  image(out, which = "NH3", grid = CNPDIAdepth(out), ylim = c(10,0), 
        main = "NH3", mfrow = NULL, legend = TRUE, ylab = "cm",
        clab = "mmol/m3")
  matplot.1D(out, which = "NH3", xyswap = TRUE, grid = CNPDIAdepth(out), 
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
