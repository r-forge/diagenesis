#####################################################################################################
######                           FESDIA: C, N, P, O2 diagenesis                                ######
######                                     BUDGETTING                                          ######
#####################################################################################################

# Utility functions
dVal <- function(out)  {# derivative
  toty <- c("time", "TotalFDET", "TotalSDET", "TotalO2", "TotalNO3",       
   "TotalNO2", "TotalNH3",  "TotalDIC",  "TotalFe", "TotalFeOH3", "TotalH2S",
   "TotalCH4", "TotalPO4", "TotalFeP", "TotalCaP", "TotalPads") 
  if (inherits(out, c("PHDIAstd", "PHDIAdyn")))
    toty <- c(toty, "TotalCaCO3", "TotalARAG")

  OUT <- out[c(1, nrow(out)), toty]
  as.list((OUT[2,]-OUT[1,])/(OUT[2,1]-OUT[1,1]))
}

getMean0D <- function(out){ # mean value, taking into account unequal times
  OUT <- FESDIA0D(out)
  on <- OUT$name
  OUT <- as.list(OUT$value)
  names(OUT) <- on
  return(OUT)
}

#====================#
# budget wrapper     #
#====================#

## -----------------------------------------------------------------------------

FESDIAbudget_all <- function(out, ..., 
   which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat"), 
   func = FESDIAbudgetO2_one, args) {

  which <- match.arg(which, choices = c("All", "Rates", "Fluxes", "Losses", "Fluxmat"))
  NM <- unlist(lapply(args[-1], as.character))

  ALL <- list(out, ...)
  
  budg <- func(out)  

  if (length(ALL) > 1) {

  
  budgFlux    <- unlist(budg$Fluxes)
  budgRates   <- unlist(budg$Rates)
  budgLoss    <- budg$Losses
  budgdC      <- budg$dC
  budgDelta   <- budg$Delta
  budgFluxmat <- as.vector(budg$Fluxmat)
  RES <- unlist(budg)

    for ( i in 2:length(ALL)) {
    budg <- func(ALL[[i]])
    
    budgFlux    <- cbind(unlist(budgFlux),    unlist(budg$Fluxes))
    budgRates   <- cbind(unlist(budgRates),   unlist(budg$Rates))
    budgLoss    <- cbind(unlist(budgLoss),    budg$Losses)
#    budgPerturb <- cbind(unlist(budgPerturb), budg$Perturb)
    budgdC      <- cbind(unlist(budgdC),      budg$dC)
    budgDelta   <- cbind(unlist(budgDelta),   budg$Delta)
    budgFluxmat  <- cbind(unlist(budgFluxmat),    as.vector(budg$Fluxmat))

    } 
  
    cn <- rep(names(budg$Fluxes), each = 4)
    nc <- nchar(cn)
    cn <- paste (cn, c("surf","deep","perturb","net"),sep="")

    rownames(budgFlux) <- cn

        budg <- list(Fluxes = budgFlux, Rates = budgRates, Losses = budgLoss, #Perturb = budgPerturb,
           dC = budgdC, Delta = budgDelta)
   }
  if (which == "Rates")
     budg <- budg$Rates
  else if (which == "Fluxes")
     budg <- budg$Fluxes
  else if (which == "Losses")
     budg <- budg$Losses
  else if (which == "Fluxmat")
     budg <- budgFluxmat
  return(budg)
}    

## --------------------------------------------------------------------------------------------------

FESDIAbudgetO2 <- function(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudget_all(out, ..., which = which, func = FESDIAbudgetO2_one, args = sys.call())  
FESDIAbudgetC <- function(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudget_all(out, ..., which = which, func = FESDIAbudgetC_one, args = sys.call())  
FESDIAbudgetN <- function(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudget_all(out, ..., which = which, func = FESDIAbudgetN_one, args = sys.call())  
FESDIAbudgetP <- function(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudget_all(out, ..., which = which, func = FESDIAbudgetP_one, args = sys.call())  
FESDIAbudgetS <- function(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudget_all(out, ..., which = which, func = FESDIAbudgetS_one, args = sys.call())  
FESDIAbudgetFe <- function(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudget_all(out, ..., which = which, func = FESDIAbudgetFe_one, args = sys.call())  


FESDIAbudgetO2_one <- function(out) {
  
  
  pertFlux <- attributes(out)$perturbFlux
  pF <- NULL
  if (! is.null(pertFlux)){
   pF <- colSums(pertFlux)/diff(range(out[,1]))
  }

  dC <- NULL
  if ("deSolve" %in% class(out)){
    dC  <- dVal(out)
    out <- getMean0D(out)
  }
  
  Cflux <- out$O2flux - out$O2deepflux

  Fluxes <- data.frame(O2 = c(out$O2flux, out$O2deepflux))
  if(! is.null(pF)) 
     Fluxes <- rbind(Fluxes, perturb = unlist(pF[c("O2")]))
  else
     Fluxes <- rbind(Fluxes, perturb = c(0, 0))
  Fluxes           <- rbind(Fluxes, Fluxes[1,] - Fluxes[2,] + Fluxes[3,])                     
  rownames(Fluxes) <- c("surface", "bottom", "perturb", "netin")

  Rates <- data.frame(
    Nitrification = out$TotNitri1*1.5 + out$TotNitri2*0.5,
    FeOxidation  = 0.25*out$TotFeoxid,
    H2Soxidation = 2.*out$TotH2Soxid,
    CH4oxidation = 2.*out$TotCH4oxid,
    H2Soxid.dist = 2.*out$TotH2Soxsurf,
    CH4oxid.dist = 2.*out$TotCH4oxsurf,
    OxicMineralisation =  out$TotOxic,
    MPBO2production = out$TotO2prod,
    MPBO2respiration = 0)
 
  Rates$Total = sum(Rates)   
  rownames(Rates) <- "nmolO2/cm2/d" 
  
  mat <- matrix (nrow = 8, ncol = 8, data = 0)
  rownames(mat) <- colnames(mat) <- c("Ext", "O2", "NO2", "NO3", "DIC", "SO4", "FeOH3", "Burial")
  mat["Ext", "O2"] <- out$O2flux * (out$O2flux > 0)      
  mat["O2", "Ext"] <- -out$O2flux * (out$O2flux < 0)      
  mat["O2", "Burial"] <- out$O2deepflux* (out$O2deepflux > 0)   
  mat["Burial", "O2"] <- -out$O2deepflux* (out$O2deepflux < 0)   
  mat["O2", "DIC"] <- Rates$OxicMineralisation + Rates$MPBO2respiration + Rates$CH4oxidation + Rates$CH4oxid.dist
  mat["O2", "NO2"] <- out$TotNitri1*1.5
  mat["O2", "NO3"] <- out$TotNitri2*0.5
  mat["O2", "SO4"] <- Rates$H2Soxidation + Rates$H2Soxid.dist
  mat["O2", "FeOH3"] <- Rates$FeOxidation
  mat["DIC", "O2"] <- Rates$MPBO2production

  if (is.null(dC))
  dC <- rowSums(mat) - colSums(mat)
  else 
    dC <- c(O2 = dC$TotalO2)
  return(list(
     Fluxes = Fluxes, Rates = Rates, Losses = out$O2deepflux, dC = c(dC, sum = sum(dC)), 
     Delta = Cflux + Rates$Total, Fluxmat = mat))                 
 

}


FESDIAbudgetC_one <- function(out) {

  Parms <- FESDIAparms(out, TRUE)  
  hasPH <- inherits(out, c("PHDIAstd", "PHDIAdyn"))
  
  pertFlux <- attributes(out)$perturbFlux
  pF <- NULL
  if (! is.null(pertFlux)){
   pF <- colSums(pertFlux)/diff(range(out[,1]))
  }
    
  dV <- NULL
  if ("deSolve" %in% class(out)){
    dV  <- dVal(out)
    out <- getMean0D(out)
  }

  Fluxes <- data.frame(FDET    = c(out$FDETflux, out$FDETdeepflux), 
                       SDET    = c(out$SDETflux, out$SDETdeepflux), 
                       DIC     = c(out$DICflux,  out$DICdeepflux),
                       CH4     = c(out$CH4flux, out$CH4deepflux),
                       CinCaP  = c(out$CaPflux,  out$CaPdeepflux)*Parms[["CPrCaP"]])
  if (hasPH) Fluxes <- cbind (Fluxes, 
                       CaCO3   = c(out$CaCO3flux, out$CaCO3deepflux), 
                       ARAG    = c(out$ARAGflux,  out$ARAGdeepflux))
  else  Fluxes <- cbind (Fluxes, 
                       CaCO3   = c(0, 0), 
                       ARAG    = c(0,  0))
  
  if(! is.null(pF)) Fluxes <- rbind(Fluxes, 
    perturb = c(unlist(pF[c("FDET", "SDET", "DIC", "CaP")])*c(1,1,1,Parms[["CPrCaP"]])), 0,0)
  else
     Fluxes <- rbind(Fluxes, perturb = c(0, 0, 0, 0, 0, 0))

  Fluxes           <- rbind(Fluxes, Fluxes[1,] - Fluxes[2,] + Fluxes[3,])                     
  rownames(Fluxes) <- c("surface", "bottom", "perturb", "netin")
  Fluxes$Total <- rowSums(Fluxes)
  influx <- Fluxes[1,"Total"]
  burial <- Fluxes[2,"Total"]
  Cflux   <- influx - burial

  Rates <- data.frame(
    OxicMineralisation = out$TotOxic,
    Denitrification    = out$TotDenit,
    IronReduction      = out$TotFered,
    SulphateReduction  = out$TotBSR,
    Methanogenesis     = out$TotMeth,
    TotalMineralisation   = out$TotMin,  
    CH4oxidation       = out$TotCH4oxid,
    CH4oxid.dist       = out$TotCH4oxsurf,
    CH4oxidAOM         = out$TotAOM,
    MPBDICuptake       = -out$TotO2prod,
    MPBFDETproduction  = out$TotO2prod,
    MPBResp            = 0,
    CaPprecipitation   = out$TotCaPprod*Parms[["CPrCaP"]],
    CaPdissolution     = out$TotCaPdiss*Parms[["CPrCaP"]],
    CaCO3dissolution   = 0,
    ARAGdissolution    = 0,
    CaCO3production    = 0
    )
  if (hasPH){
    Rates$CaCO3dissolution  <- out$TotCaCO3diss
    Rates$ARAGdissolution   <- out$TotARAGdiss
    Rates$CaCO3production   <- out$TotCaCO3prod
    }
    rownames(Rates) <- "nmolC/cm2/d"  

      # derivatives
  mat <- matrix (nrow = 9, ncol = 9, data = 0)
  rownames(mat) <- colnames(mat) <- c("Ext", "DET", "DIC", "CaP", "CH4", "MPB", "CaCO3", "ARAG", "Burial")
  mat["Ext", "DET"]    <- out$FDETflux      + out$SDETflux
  mat["DET", "Burial"] <- out$FDETdeepflux  + out$SDETdeepflux
  mat["Ext", "DIC"]    <-  out$DICflux * (out$DICflux > 0) 
  mat["DIC", "Ext"]    <- -out$DICflux * (out$DICflux < 0)
  mat["Burial", "DIC"] <- - out$DICdeepflux*(out$DICdeepflux < 0)
  mat["DIC", "Burial"] <- + out$DICdeepflux*(out$DICdeepflux > 0) 
  mat["Ext", "CaP"] <- 0
  mat["CaP", "Ext"] <- out$CaPdeepflux*Parms[["CPrCaP"]] 

  mat["MPB", "DET"] <- Rates$MPBFDETproduction 
  mat["DET", "DIC"] <- Rates$TotalMineralisation-0.5*Rates$Methanogenesis
  mat["DIC", "MPB"] <- Rates$MPBDICuptake
  mat["DIC", "CaP"] <- Rates$CaPprecipitation
  mat["CaP", "DIC"] <- Rates$CaPdissolution
  mat["MPB", "DIC"] <- Rates$MPBResp 
  mat["DET", "CH4"] <- 0.5*Rates$Methanogenesis
  mat["CH4", "DIC"] <- Rates$CH4oxidation + Rates$CH4oxid.dist + Rates$CH4oxidAOM

  if (hasPH){
    mat["Ext", "CaCO3"]    <-  out$CaCO3flux * (out$CaCO3flux > 0) 
    mat["CaCO3", "Ext"]    <- -out$CaCO3flux * (out$CaCO3flux < 0)
    mat["Ext", "ARAG"]    <-  out$ARAGflux * (out$ARAGflux > 0) 
    mat["ARAG", "Ext"]    <- -out$ARAGflux * (out$ARAGflux < 0)
    mat["Burial", "CaCO3"] <- - out$CaCO3deepflux*(out$CaCO3deepflux < 0)
    mat["CaCO3", "Burial"] <- + out$CaCO3deepflux*(out$CaCO3deepflux > 0) 
    mat["Burial", "ARAG"] <- - out$ARAGdeepflux*(out$ARAGdeepflux < 0)
    mat["ARAG", "Burial"] <- + out$ARAGdeepflux*(out$ARAGdeepflux > 0) 
    mat["DIC", "CaCO3"]    <- Rates$CaCO3production
    mat["CaCO3", "DIC"]    <- Rates$CaCO3dissolution
    mat["ARAG", "DIC"]    <- Rates$ARAGdissolution
  }
  # derivatives
  if (is.null(dV)) 
    dC <- rowSums(mat) - colSums(mat)
  else
    dC <- c(DET = dV$TotalFDET+dV$TotalSDET, DIC = dV$TotalDIC, 
      CaP = dV$TotalCaP, CH4 = dV$TotalCH4  , CaCO3 = dV$TotalCaCO3, ARAG = dV$TotalARAG)


  return(list(Fluxes = Fluxes, Rates = Rates, Losses = burial, dC = c(dC, sum = sum(dC)), 
     Delta = influx- burial, Fluxmat = mat))                 
 
}

FESDIAbudgetS_one <- function(out) {

  Parms <- FESDIAparms(out, TRUE)  

  pertFlux <- attributes(out)$perturbFlux
  pF <- NULL
  if (! is.null(pertFlux)){
   pF <- colSums(pertFlux)/diff(range(out[,1]))
  }
    
  dV <- NULL
  if ("deSolve" %in% class(out)){
    dV  <- dVal(out)
    out <- getMean0D(out)
  }

  Fluxes <- data.frame(H2S = c(out$H2Sflux, out$H2Sdeepflux), 
                       SO4  = c(out$SO4flux, out$SO4deepflux))

  if(! is.null(pF)) Fluxes <- rbind(Fluxes, 
    perturb = unlist(pF[c("H2S", "SO4")]))
  else
    Fluxes <- rbind(Fluxes, perturb = c(0, 0))

  Fluxes           <- rbind(Fluxes, Fluxes[1,] - Fluxes[2,] + Fluxes[3,])                     
  rownames(Fluxes) <- c("surface", "bottom", "perturb", "netin")
  Fluxes$Total <- rowSums(Fluxes)
  influx <- Fluxes[1,"Total"]
  burial <- Fluxes[2,"Total"]
  Cflux   <- influx - burial

  Rates <- data.frame(
    SulphateReduction = 0.5*out$TotBSR,
    H2Soxidation      = out$TotH2Soxid,
    H2Soxid.dist      = 2.*out$TotH2Soxsurf,
    AOM               = out$TotAOM,
    FeSproduction     = out$TotFeSprod 
  )      
  rownames(Rates) <- "nmolS/cm2/d" 


  mat <- matrix (nrow = 5, ncol = 5, data = 0)
  rownames(mat) <- colnames(mat) <- c("Ext", "H2S", "SO4", "Burial", "FeS")
  mat["Ext", "H2S"]    <- out$H2Sflux*(out$H2Sflux>1)
  mat["H2S", "Ext"]    <- -out$H2Sflux*(out$H2Sflux<1)
  mat["H2S", "Burial"] <- out$H2Sdeepflux*(out$H2Sdeepflux>0)
  mat["Burial", "H2S"] <- -out$H2Sdeepflux*(out$H2Sdeepflux<0)
  mat["Ext", "SO4"]    <-  out$SO4flux * (out$SO4flux > 0) 
  mat["SO4", "Ext"]    <- -out$SO4flux * (out$SO4flux < 0)
  mat["Burial", "SO4"] <- - out$SO4deepflux*(out$SO4deepflux < 0)
  mat["SO4", "Burial"] <- + out$SO4deepflux*(out$SO4deepflux > 0) 

  mat["SO4", "H2S"] <- Rates$SulphateReduction
  mat["H2S", "SO4"] <- Rates$H2Soxidation + Rates$H2Soxid.dist + Rates$AOM 
  mat["H2S", "FeS"] <- Rates$FeSproduction

  # derivatives
  if (is.null(dV)) 
    dC <- rowSums(mat) - colSums(mat)
  else
    dC = c(H2S = dV$TotalH2S, SO4 = dV$TotalSO4, FeS = dV$TotalFeS)

  return(list(Fluxes = Fluxes, Rates = Rates, Losses = burial+dC[["FeS"]], dC = c(dC, sum = sum(dC)), 
     Delta = sum(dC), Fluxmat = mat))                 

}

FESDIAbudgetFe_one <- function(out) {

  Parms <- FESDIAparms(out, TRUE)  

  pertFlux <- attributes(out)$perturbFlux
  pF <- NULL
  if (! is.null(pertFlux)){
   pF <- colSums(pertFlux)/diff(range(out[,1]))
  }
    
  dV <- NULL
  if ("deSolve" %in% class(out)){
    dV  <- dVal(out)
    out <- getMean0D(out)
  }

  Fluxes <- data.frame(Fe    = c(out$Feflux,    out$Fedeepflux), 
                       FeOH3 = c(out$FeOH3flux, out$FeOH3deepflux))

  if(! is.null(pF)) Fluxes <- rbind(Fluxes, 
    perturb = unlist(pF[c("Fe", "FeOH3")]))
  else
     Fluxes <- rbind(Fluxes, perturb = c(0, 0))

  Fluxes           <- rbind(Fluxes, Fluxes[1,] - Fluxes[2,] + Fluxes[3,])                     
  rownames(Fluxes) <- c("surface", "bottom", "perturb", "netin")
  Fluxes$Total <- rowSums(Fluxes)
  influx <- Fluxes[1,"Total"]
  burial <- Fluxes[2,"Total"]
  Cflux   <- influx - burial

   Rates <- data.frame(
    IronReduction = 4.*out$TotFered,
    FeOxidation   = out$TotFeoxid,
    FeSproduction = out$TotFeSprod 
  )      
  rownames(Rates) <- "nmolFe/cm2/d" 

  
  mat <- matrix (nrow = 5, ncol = 5, data = 0)
  rownames(mat) <- colnames(mat) <- c("Ext", "Fe", "FeOH3", "Burial", "FeS")
  mat["Ext", "FeOH3"] <- out$FeOH3flux*(out$FeOH3flux>1)
  mat["FeOH3", "Ext"] <- -out$FeOH3flux*(out$FeOH3flux<1)
  mat["FeOH3", "Burial"] <- out$FeOH3deepflux*(out$FeOH3deepflux>0)
  mat["Burial", "FeOH3" ] <- -out$FeOH3deepflux*(out$FeOH3deepflux<0)
  mat["Ext", "Fe"]    <-  out$Feflux * (out$Feflux > 0) 
  mat["Fe", "Ext"]    <- -out$Feflux * (out$Feflux < 0)
  mat["Burial", "Fe"] <- - out$Fedeepflux*(out$Fedeepflux < 0)
  mat["Fe", "Burial"] <- + out$Fedeepflux*(out$Fedeepflux > 0) 

  mat["FeOH3", "Fe"] <- Rates$IronReduction
  mat["Fe", "FeOH3"] <- Rates$FeOxidation
  mat["Fe", "FeS"] <- Rates$FeSproduction

  # derivatives
  if (is.null(dV)) 
    dC <- rowSums(mat) - colSums(mat)
  else
    dC = c(H2S = dV$TotalFe, SO4 = dV$TotalFeOH3, FeS = dV$TotalFeS)

  return(list(Fluxes = Fluxes, Rates = Rates, loss = burial+dC[["FeS"]], dC = c(dC, sum = sum(dC)), 
     Delta = sum(dC), Fluxmat = mat))                 
}



 

FESDIAbudgetN_one <- function(out) {

  Parms <- FESDIAparms(out, TRUE)  

  pertFlux <- attributes(out)$perturbFlux
  pF <- NULL
  if (! is.null(pertFlux)){
   pF <- colSums(pertFlux)/diff(range(out[,1]))
  }
    
  dV <- NULL
  if ("deSolve" %in% class(out)){
    dV  <- dVal(out)
    out <- getMean0D(out)
  }

  Fluxes <- data.frame(FDET_N = c(out$FDETflux        *Parms[["NCrFdet"]] , 
                                  out$FDETdeepflux    *Parms[["NCrFdet"]] ), 
                       SDET_N = c(out$SDETflux        *Parms[["NCrSdet"]] , 
                                  out$SDETdeepflux)   *Parms[["NCrSdet"]] , 
                       NO3  = c(out$NO3flux, out$NO3deepflux),
                       NO2  = c(out$NO2flux, out$NO2deepflux),
                       NH3  = c(out$NH3flux, out$NH3deepflux))
  if(! is.null(pF)) Fluxes <- rbind(Fluxes, 
    perturb = unlist(pF[c("FDET", "SDET", "NO3", "NO2", "NH3")])*c(Parms[["NCrFdet"]],Parms[["NCrSdet"]],1,1,1))
  else
     Fluxes <- rbind(Fluxes, perturb = c(0, 0, 0, 0))

  Fluxes           <- rbind(Fluxes, Fluxes[1,] - Fluxes[2,] + Fluxes[3,])                     
  rownames(Fluxes) <- c("surface", "bottom", "perturb", "netin")
  Fluxes$Total <- rowSums(Fluxes)
  influx <- Fluxes[1,"Total"]
  burial <- Fluxes[2,"Total"]
  Cflux   <- influx - burial

  Rates <- data.frame(
    NH3production   = out$TotNH3prod,
    Denitrification = out$TotDenit*0.8,
    MPBNO3consumption  = out$TotMPBNO3uptake,
    MPBNH3consumption  = out$TotMPBNH3uptake,
    Nitrification1     = out$TotNitri1,
    Nitrification2     = out$TotNitri2,
    Anammox            = out$TotAnammox,
    NH3adsorption      = out$TotNH3ads,
    N2production       = out$TotDenit*0.8 + 2*out$TotAnammox,
    MPBNdeath          = out$TotMPBNO3uptake + out$TotMPBNH3uptake,
    NH3prodMPBdeath    = 0,
    DETNprodMPBdeath   = out$TotMPBNO3uptake + out$TotMPBNH3uptake
  )
  

  rownames(Rates) <- "nmolN/cm2/d" 

  # derivatives
  mat <- matrix (nrow = 8, ncol = 8, data = 0)
  rownames(mat) <- colnames(mat) <- c("Ext", "DET", "NH3", "NO3", "NO2", "MPB", "N2", "Burial")
  mat["Ext", "DET"] <- out$FDETflux*Parms[["NCrFdet"]]  + out$SDETflux*Parms[["NCrSdet"]] 
  mat["DET", "Burial"] <- out$FDETdeepflux*Parms[["NCrFdet"]] + out$SDETdeepflux *Parms[["NCrSdet"]]
  
  mat["Ext", "NH3"]    <-  out$NH3flux * (out$NH3flux > 0) 
  mat["NH3", "Ext"]    <- -out$NH3flux * (out$NH3flux < 0)
  mat["Burial", "NH3"] <- - out$NH3deepflux*(out$NH3deepflux < 0)
  mat["NH3", "Burial"] <- + out$NH3deepflux*(out$NH3deepflux > 0) 

  mat["Ext", "NO3"]    <-  out$NO3flux * (out$NO3flux > 0) 
  mat["NO3", "Ext"]    <- -out$NO3flux * (out$NO3flux < 0)
  mat["Burial", "NO3"] <- - out$NO3deepflux*(out$NO3deepflux < 0)
  mat["NO3", "Burial"] <- + out$NO3deepflux*(out$NO3deepflux > 0) 

  mat["Ext", "NO2"]    <-  out$NO2flux * (out$NO2flux > 0) 
  mat["NO2", "Ext"]    <- -out$NO2flux * (out$NO2flux < 0)
  mat["Burial", "NO2"] <- - out$NO2deepflux*(out$NO2deepflux < 0)
  mat["NO2", "Burial"] <- + out$NO2deepflux*(out$NO2deepflux > 0) 
  
  mat["MPB", "DET"] <- Rates$DETNprodMPBdeath 
  mat["DET", "NH3"] <- Rates$NH3production
  mat["NH3", "MPB"] <- Rates$MPBNH3consumption
  mat["NH3", "NO2"] <- Rates$Nitrification1
  mat["NH3", "N2"]  <- Rates$Anammox
  mat["MPB", "NH3"] <- Rates$NH3prodMPBdeath
  mat["NO3", "N2"]  <- Rates$Denitrification
  mat["NO3", "MPB"] <- Rates$MPBNO3consumption
  mat["NO2", "NO3"] <- Rates$Nitrification2
  mat["NO2", "N2"]  <- Rates$Anammox
  
  # derivatives
  if (is.null(dV)) 
    dC <- rowSums(mat) - colSums(mat)
  else
    dC = c(DET = dV$TotalFDET*Parms[["NCrFdet"]] + dV$TotalSDET*Parms[["NCrSdet"]], 
      NH3 = dV$TotalNH3*(1+Parms[["NH3Ads"]]), NO3 = dV$TotalNO3, NO2 = dV$TotalNO2)

  return(list(Fluxes = Fluxes, Rates = Rates, Losses = -dC["N2"]+burial, 
     dC = c(dC, sum = sum(dC)), Delta = sum(dC), Fluxmat = mat))                 

}

FESDIAbudgetP_one <- function(out) {

  Parms <- FESDIAparms(out, TRUE)  

  pertFlux <- attributes(out)$perturbFlux
  pF <- NULL
  if (! is.null(pertFlux)){
   pF <- colSums(pertFlux)/diff(range(out[,1]))
  }
    
  dV <- NULL
  if ("deSolve" %in% class(out)){
    dV  <- dVal(out)
    out <- getMean0D(out)
  }

  Fluxes <- data.frame(FDET_P = c(out$FDETflux     *Parms[["PCrFdet"]] , 
                                  out$FDETdeepflux *Parms[["PCrFdet"]] ), 
                       SDET_P = c(out$SDETflux     *Parms[["PCrSdet"]] , 
                                  out$SDETdeepflux *Parms[["PCrSdet"]]) , 
                       PO4  = c(out$PO4flux, out$PO4deepflux),
                       FeP  = c(0, out$FePdeepflux),
                       CaP  = c(0, out$CaPdeepflux))
  if(! is.null(pF)) Fluxes <- rbind(Fluxes, 
    perturb = unlist(pF[c("FDET", "SDET", "PO4", "FeP", "CaP")])*c(Parms[["PCrFdet"]],Parms[["PCrSdet"]],1,1,1))
  else
     Fluxes <- rbind(Fluxes, perturb = c(0, 0, 0, 0))

  Fluxes           <- rbind(Fluxes, Fluxes[1,] - Fluxes[2,] + Fluxes[3,])                     
  rownames(Fluxes) <- c("surface", "bottom", "perturb", "netin")
  Fluxes$Total <- rowSums(Fluxes)
  influx <- Fluxes[1,"Total"]
  burial <- Fluxes[2,"Total"]
  Cflux   <- influx - burial

 Rates <- data.frame(
    PO4production   = out$TotPO4prod,
    FePadsorption   = out$TotFePprod,
    FePdesorption   = out$TotFePdesorp,
    MPBPO4consumption  = out$TotMPBPO4uptake,
    CaPproduction   = out$TotCaPprod,
    CaPdissolution  = out$TotCaPdiss)                   
 
    Rates$MPBPdeath <- Rates$MPBPO4consumption
    Rates$PO4prodMPBdeath <- 0
    Rates$DETPprodMPBdeath <- Rates$MPBPO4consumption
      
   rownames(Rates) <- "nmolP/cm2/d" 

     mat <- matrix (nrow = 7, ncol = 7, data = 0)
  rownames(mat) <- colnames(mat) <- c("Ext", "DET", "PO4", "FeP", "CaP", "MPB", "Burial")
  mat["Ext", "DET"] <- out$FDETflux*Parms[["PCrFdet"]]  + out$SDETflux*Parms[["PCrSdet"]] 
  mat["DET", "Burial"] <- out$FDETdeepflux*Parms[["PCrFdet"]] + out$SDETdeepflux *Parms[["PCrSdet"]]
  mat["Ext", "PO4"]    <-  out$PO4flux * (out$PO4flux > 0) 
  mat["PO4", "Ext"]    <- -out$PO4flux * (out$PO4flux < 0)
  mat["Burial", "PO4"] <- - out$PO4deepflux*(out$PO4deepflux < 0)
  mat["PO4", "Burial"] <- + out$PO4deepflux*(out$PO4deepflux > 0) 

  mat["Ext", "FeP"] <- 0
  mat["FeP", "Burial"] <- out$FePdeepflux
  mat["Ext", "CaP"] <- 0
  mat["CaP", "Burial"] <- out$CaPdeepflux
  
  mat["MPB", "DET"] <- Rates$DETPprodMPBdeath 
  mat["DET", "PO4"] <- Rates$PO4production
  mat["PO4", "MPB"] <- Rates$MPBPO4consumption
  mat["PO4", "FeP"] <- Rates$FePadsorption
  mat["FeP", "PO4"] <- Rates$FePdesorption
  mat["PO4", "CaP"] <- Rates$CaPproduction
  mat["CaP", "PO4"] <- Rates$CaPdissolution
  mat["MPB", "PO4"] <- Rates$PO4prodMPBdeath

  # derivatives
  if (is.null(dV)) 
    dC <- rowSums(mat) - colSums(mat)
  else
    dC = c(DET = dV$TotalFDET*Parms[["PCrFdet"]] + dV$TotalSDET*Parms[["PCrSdet"]], 
      PO4 = dV$TotalPO4, FeP = dV$TotalFeP, CaP = dV$TotalCaP)

  return(list(Fluxes = Fluxes, Rates = Rates, Losses = burial, dC = c(dC, sum = sum(dC)), 
     Delta = sum(dC), Fluxmat = mat))                 
}

