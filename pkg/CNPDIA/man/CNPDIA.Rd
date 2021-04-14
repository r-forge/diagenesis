\name{CNPDIA-package}
\alias{CNPDIA-package}
\alias{CNPDIA}
\docType{package}
\title{
  An early diagenetic model for C, N, P and O2 cycle in marine sediments.
}
\description{

  CNPDIA is a 1-D model of Carbon, nitrogen, phosphorus and oxygen diagenesis in a marine sediment.
  It is based on the OMEXDIA model, extended with simple P dynamics and with 
  microphytobenthos dynamics.

  The model describes twelve state variables, in \bold{100} layers:
  \itemize{
    \item 2 fractions of organic carbon (FDET,SDET): fast and slow
      decaying, solid substance.
    \item Oxygen (O2), dissolved substance.
    \item Nitrate (NO3), dissolved substance.
    \item Nitrite (NO2), dissolved substance.
    \item Ammonium (NH3), dissolved substance.
    \item Oxygen demand unit (ODU), dissolved substance, as a
      lump sum of all reduced substances other than ammonium.
    \item Phosphate (PO4), dissolved substance
    \item Iron-bound P (FeP), P bound to iron oxides, solid substance
    \item Ca-bound P (CaP), P bound to iron oxides, solid substance
    \item dissolved inorganic carbon (DIC), dissolved substance
    \item adsorbed P (Pads), solid substance
  }
  Concentrations of liquids and solids are expressed in [nmol/cm3 liquid] and  [nmol/cm3 solid] respectively.
  See Soetaert et al., 1996 for further details of the original OMEXDIA model.

  MPBDIA extends the CNPDIA model with microphytobenthos dynamics. In addition to the CNPDIA state variables it also describes the following state variables:
    \itemize{
    \item MPBC, microphytobenthos carbon, [mmolC/m3 solid]
    \item MPBN, microphytobenthos nitrogen, [mmolN/m3 solid]
    \item MPBC, microphytobenthos chlorophyll, [mgChl/m3 solid]
    \item EPS, extracellular polymeric substance, [mmolC/m3 solid]
    \item PSIIin, fraction of photosynthesis system II that is in inhibited state [-]

  }

  The models are implemented in fortran and linked to R.

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
