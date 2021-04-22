

FESDIAextract <- function(out, which = NULL) {
  
  if (missing(out)) 
    stop("object 'out' should be given")
  
  if ("steady1D" %in% (class(out))) {
    if (which %in% colnames(out$y))
      W <- out$y[,which]
    else
      W <- out[[which]]
  } else W <- subset(out, which = which)  
  if (any(is.null(W))) W <- 0  
  return(W)
}


FESDIApH <- function (out) {

  P    <- FESDIAparms(out, as.vector = TRUE)
  
  temp    <- P["temperature"]
  S       <- P["salinity"]
  D       <- sw_dens(S = S, t = temp)
  
  SumCO2  <- as.vector(FESDIAextract(out, which = "DIC")/1000/D)  # convert to mol/kg
  SumNH4  <- as.vector(FESDIAextract(out, which = "NH3")/1000/D)
  SumH2S  <- as.vector(FESDIAextract(out, which = "H2S")/1000/D)
  SumH3PO4<- as.vector(FESDIAextract(out, which = "PO4")/1000/D)
  SumH2SO4<- as.vector(FESDIAextract(out, which = "SO4")/1000/D)
  SumHNO3 <- as.vector(FESDIAextract(out, which = "NO3")/1000/D)
  SumHNO2 <- as.vector(FESDIAextract(out, which = "NO2")/1000/D)
  TA      <- FESDIAextract(out, which = "ALK")/1000/D
  
  DIM     <- dim(TA)
  TA      <- as.vector(TA)
  len     <- length(TA)
  
  A <- .Fortran ("EstimatePH", Temp = temp, Sal = S, depth = 0., TA = TA, 
                 SumCO2 = SumCO2, SumH3PO4 = SumH3PO4, SumNH4 = SumNH4, SumH2S =SumH2S,  
                 SumH2SO4 = SumH2SO4, maxIter = as.integer(10), numvals = as.integer(len), 
                 pH  = as.double(rep(8., times = len)), KS = as.double(rep(0., 20)),
                 Cts = as.double(rep(0., times = 12)), PACKAGE = "FESDIA")
  
#    PH <- aquaenv(S = S, t = temp, SumCO2 = SumCO2, SumNH4 = SumNH4, SumH2S = SumH2S, 
#                SumH3PO4 = SumH3PO4, SumH2SO4 = SumH2SO4, SumHNO3 = SumHNO3, TA = TA,
#                speciation = FALSE, skeleton = TRUE)$pH
  PH       <- A$pH
  K        <- A$KS    # dissociation constants
  
  names(K) <- c("K_W","K_HF","K_CO2","K_HCO3","K_BOH3","K_NH4","K_H2S",
                "K_HS","K_H3PO4","K_H2PO4","K_HPO4","K_HSO4","K_H2SO4",
                "K0_CO2", "Ksp_Calcite", "Ksp_Aragonite",   "TotToFree",
                "SWSToFree", "SumBorate", "SumF")
  
  if (length(DIM)) dim(PH) <- DIM
  attr(PH, "Constants")    <- K
  PH
}


dpH.int <- function(out){
  por <- FESDIApor(out)
  dx  <- FESDIAdx(out)
  phr <- dpH.Rates(out)
  Rates <- apply(phr, MARGIN = 2, FUN = function(x) sum(x*por*dx))
  
}

  
dpH.Rates <- function(out, asH = FALSE){

  por <- FESDIApor(out)
  
  getVars <- function(name, solid = FALSE){
    Var <- as.vector(FESDIAextract(out, which = name))     # convert to mol
    if (solid) Var <- Var*(1-por)/por
    Var
  }

  dPH.numeric <- function(dSumCO2   = 0, dSumBOH3 = 0, dSumH2S = 0,     # change in summed conc due to process
                        dSumSiOH4 = 0, dSumNH4  = 0, dSumH3PO4 = 0,   # default = 0
                        dSumHNO3  = 0, dSumHNO2 = 0, dSumHF = 0, dSumH2SO4 = 0,    
                        dTA = 0,  dC = 1e-8                   # change in total alkalinity due to process        
                        ) {                      # pH range for which to estimate dpHdProcess
  
    P       <- FESDIAparms(out, as.vector = TRUE)
    temp    <- P["temperature"]
    S       <- P["salinity"]
    D       <- sw_dens(S = S, t = temp)
  
  # reference settings
    pH      <- FESDIAextract(out, which = "pH")
    SumCO2  <- as.vector(FESDIAextract(out, which = "DIC")/1000/D)  # convert to mol/kg
    SumNH4  <- as.vector(FESDIAextract(out, which = "NH3")/1000/D)
    SumH2S  <- as.vector(FESDIAextract(out, which = "H2S")/1000/D)
    SumH3PO4<- as.vector(FESDIAextract(out, which = "PO4")/1000/D)
    SumH2SO4<- as.vector(FESDIAextract(out, which = "SO4")/1000/D)
    SumHNO3 <- as.vector(FESDIAextract(out, which = "NO3")/1000/D)
    SumHNO2 <- as.vector(FESDIAextract(out, which = "NO2")/1000/D)
    TA      <- as.vector(FESDIAextract(out, which = "ALK")/1000/D)
  
  # perturbed settings
    delt    <- dC/1000/D
  
    SumCO2p   <- SumCO2   + delt*dSumCO2
    SumNH4p   <- SumNH4   + delt*dSumNH4
    SumH2Sp   <- SumH2S   + delt*dSumH2S
    SumH3PO4p <- SumH3PO4 + delt*dSumH3PO4
    SumH2SO4p <- SumH2SO4 + delt*dSumH2SO4
    SumHNO3p  <- SumHNO3  + delt*dSumHNO3
    SumHNO2p  <- SumHNO2  + delt*dSumHNO2
    TAp       <- TA       + delt*dTA

    DIM     <- dim(TA)
    TA      <- as.vector(TA)
    len     <- length(TA)
   # reference pH
    pHr <- .Fortran ("EstimatePH", Temp = temp, Sal = S, depth = 0., TA = TA, 
                 SumCO2 = SumCO2, SumH3PO4 = SumH3PO4, SumNH4 = SumNH4, SumH2S =SumH2S,  
                 SumH2SO4 = SumH2SO4, maxIter = as.integer(10), numvals = as.integer(len), 
                 pH  = as.double(rep(8., times = len)), KS = as.double(rep(0., 20)),
                 Cts = as.double(rep(0., times = 12)), PACKAGE = "FESDIA")$pH
    if (asH) pHr <- 10^-pHr
   # perturbed pH
    pHp  <-.Fortran ("EstimatePH", Temp = temp, Sal = S, depth = 0., TA = TAp, 
                 SumCO2 = SumCO2p, SumH3PO4 = SumH3PO4p, SumNH4 = SumNH4p, SumH2S =SumH2Sp,  
                 SumH2SO4 = SumH2SO4p, maxIter = as.integer(10), numvals = as.integer(len), 
                 pH  = as.double(rep(8., times = len)), KS = as.double(rep(0., 20)),
                 Cts = as.double(rep(0., times = 12)), PACKAGE = "FESDIA")$pH
    if (asH) pHp <- 10^-pHp
  

    # estimate dpHdspec by numerical differencing. Divide by 1e6 so that pH change is per micromol rather than per mol.
    dpH_dSpec <- (pHp-pHr)/(dC*1e6)  # units in per umolC  - DENSITY????
  
    return(dpH_dSpec)
  }

 DICprodMin   <- getVars("DICprodMin")
 NC           <- getVars("DINprodMin")/DICprodMin
 PC           <- getVars("DIPprodMin")/DICprodMin
 NC[is.nan(NC)] <- 0
 PC[is.nan(PC)] <- 0
 P       <- FESDIAparms(out, as.vector = TRUE)
  
 NCf    <- P["NCrFdet"]
 PCf    <- P["PCrFdet"]
 
 pHeffect <- data.frame(
 Oxicmin.pH    = dPH.numeric(dSumCO2=1,   dSumNH4=NC, dSumH3PO4=PC,                dTA=NC-PC    ) *getVars("Oxicmin",    TRUE),
 Denitrific.pH = dPH.numeric(dSumCO2=1,   dSumNH4=NC, dSumH3PO4=PC, dSumHNO3=-0.8, dTA=0.8+NC-PC) *getVars("Denitrific", TRUE),             
 Feredmin.pH   = dPH.numeric(dSumCO2=1,   dSumNH4=NC, dSumH3PO4=PC,                dTA=NC-PC+8   ) *getVars("Feredmin",   TRUE),
 BSRmin.pH     = dPH.numeric(dSumCO2=1,   dSumNH4=NC, dSumH3PO4=PC, dSumH2SO4=-0.5, dSumH2S=0.5,dTA=NC-PC+1) * getVars("BSRmin",     TRUE), 
 Methmin.pH    = dPH.numeric(dSumCO2=0.5, dSumNH4=NC, dSumH3PO4=PC,                             dTA=NC-PC)   * getVars("Methmin",    TRUE), 

 Nitri1.pH      = dPH.numeric(dSumNH4=-1,  dSumHNO2=1,             dTA=-2) * getVars("Nitri1"),
 Nitri2.pH      = dPH.numeric(dSumHNO2=-1, dSumHNO3=1)                     * getVars("Nitri2"),
 Anammox.pH     = dPH.numeric(dSumNH4=-1, dSumHNO2=-1,            dTA=0)* getVars("Anammox"),
 Feoxid.pH      = dPH.numeric(                                    dTA=-2)* getVars("Feoxid"),
 H2Soxid.pH     = dPH.numeric(dSumH2S=-1, dSumH2SO4=1,            dTA=-2)* getVars("H2Soxid"),
 CH4oxid.pH     = dPH.numeric(dSumCO2=1,                          dTA=0)* getVars("CH4oxid"),
 AOM.pH         = dPH.numeric(dSumH2SO4=-1, dSumH2S=1, dSumCO2=1, dTA=2)* getVars("AOM"),

 FeSprod.pH     = dPH.numeric(dSumH2S=-1.5, dTA=0)* getVars("AOM"),
 adsorption.pH  = dPH.numeric(              dTA=1)* getVars("AOM"),                           
 CO2release.pH  = dPH.numeric(dSumCO2=-1,   dTA=0)* getVars("AOM"),
 NH3release.pH  = dPH.numeric(dSumNH4=-1,   dTA=-1)* getVars("AOM"),
 NH4release.pH  = dPH.numeric(dSumNH4=-1,   dTA=0)* getVars("AOM"),

 MPBCprod.pH      = dPH.numeric(dSumCO2=-1)* getVars("MPBCprod"),
 MPBuptakeNH3.pH  = dPH.numeric(dSumNH4=-NCf, dTA=-NCf)* getVars("MPBuptakeNH3"),
 MPBuptakeNO3.pH  = dPH.numeric(dSumHNO3=-NCf, dTA= NCf)* getVars("MPBuptakeNO3"),
 MPBuptakePO4.pH  = dPH.numeric(dSumH3PO4=-PCf, dTA= PCf)* getVars("MPBuptakeNO3"),
 CaCO3prod.pH     = dPH.numeric(dSumCO2=-1, dTA=-2)* getVars("CaCO3prod"),
 CaCO3diss.pH     = dPH.numeric(dSumCO2=1,  dTA=2)* getVars("CaCO3diss"),
 ARAGdiss.pH      = dPH.numeric(dSumCO2=1,  dTA=2)* getVars("ARAGdiss")
)
 pHeffect
}