\name{FESDIA-package}
\alias{FESDIA-package}
\docType{package}
\title{
  The FESDIA early diagenetic model.
}
\description{

  FESDIA is a 1-D model of Carbon, nitrogen, phosphorus and oxygen diagenesis in a marine sediment.
  It is based on the OMEXDIA model, extended with P, Fe, S and CH4 dynamics.

  The model describes fourteen state variables, in \bold{100} layers:
  \itemize{
    \item 2 fractions of organic carbon (FDET,SDET): fast and slow
      decaying, solid substance.
    \item Oxygen (O2), dissolved substance.
    \item Nitrate (NO3), dissolved substance.
    \item Ammonium (NH3), dissolved substance.
    \item Phosphate (PO4), dissolved substance.
    \item Iron-bound P (FeP), P bound to iron oxides, solid substance.
    \item Ca-bound P (CaP), P bound to iron oxides, solid substance.
    \item dissolved inorganic carbon (DIC), dissolved substance.
    \item reduced iron (Fe), dissolved substance.
    \item iron oxides (FeOH3), solid substance.
    \item sulphide (H2S), dissolved substance.
    \item sulphate (SO4), dissolved substance.
    \item methane (CH4), dissolved substance.
    \item alkalinity (ALK), dissolved substance
  }
  
  Concentrations of liquids and solids are expressed in [nmol/cm3 liquid] and  [nmol/cm3 solid] respectively.
  See Soetaert et al., 1996 for further details of the original OMEXDIA model.

  The model is implemented in fortran and linked to R.
}

\author{
  Karline Soetaert
}
\references{
  Soetaert K, PMJ Herman and JJ Middelburg, 1996a.
  A model of early diagenetic processes from the shelf to abyssal depths.
  Geochimica Cosmochimica Acta, 60(6):1019-1040.

  Soetaert K, PMJ Herman and JJ Middelburg, 1996b.
  Dynamic response of deep-sea sediments to seasonal variation: a model.
  Limnol. Oceanogr. 41(8): 1651-1668.
}

\keyword{ package }
