#####################################################################################################
######                           CNPDIA: C, N, P, O2 diagenesis                                ######
######                                     BUDGETTING                                          ######
#####################################################################################################

# Utility functions
dVal <- function(out, includeMPB)  {# takes into account unequal timing 
    toty <- c("time", "TotalFDET", "TotalSDET", "TotalO2", "TotalNO3",       
  "TotalNO2", "TotalNH3",  "TotalODU",  "TotalDIC",       
  "TotalPO4", "TotalFeP", "TotalCaP", "TotalPads") 
  if (includeMPB) 
    toty <- c(toty, c("TotalMPBC", "TotalMPBN", "TotalCHL", "TotalEPS"))
  OUT <- out[c(1, nrow(out)), toty]
  as.list((OUT[2,]-OUT[1,])/(OUT[2,1]-OUT[1,1]))
}

getMean0D <- function(out, includeMPB){
  if (includeMPB)
  OUT <- MPBDIA0D(out)
  else
  OUT <- CNPDIA0D(out)
  on <- OUT$name
  OUT <- as.list(OUT$value)
  names(OUT) <- on
  return(OUT)
}

#====================#
# budget wrapper     #
#====================#

## -----------------------------------------------------------------------------

CNPDIAbudget_all <- function(out, ..., 
   which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat"), 
   func = CNPDIAbudgetO2_one, args) {

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
  
CNPDIAbudgetO2 <- function(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  CNPDIAbudget_all(out, ..., which = which, func = CNPDIAbudgetO2_one, args = sys.call())  

CNPDIAbudgetC <- function(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  CNPDIAbudget_all(out, ..., which = which, func = CNPDIAbudgetC_one, args = sys.call())  

CNPDIAbudgetN <- function(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  CNPDIAbudget_all(out, ..., which = which, func = CNPDIAbudgetN_one, args = sys.call())  

CNPDIAbudgetP <- function(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  CNPDIAbudget_all(out, ..., which = which, func = CNPDIAbudgetP_one, args = sys.call())  

#====================#
#  O2 budgets        #
#====================#

CNPDIAbudgetO2_one <- function(out) {

  pertFlux <- attributes(out)$perturbFlux
  pF <- NULL
  if (! is.null(pertFlux)){
   pF <- colSums(pertFlux)/diff(range(out[,1]))
  }
  if ("deSolve" %in% class(out))
    includeMPB <- attributes(out)$includeMPB
  else
    includeMPB <- out$includeMPB

  dC <- NULL
  if ("deSolve" %in% class(out)){
    dC <- dVal(out, includeMPB)
    out <- getMean0D(out, includeMPB)
  }
  
  Cflux <- out$O2flux - out$O2deepflux

  Fluxes <- data.frame(O2 = c(out$O2flux, out$O2deepflux),
    ODU = c(out$ODUflux, out$ODUdeepflux))
  if(! is.null(pF)) 
     Fluxes <- rbind(Fluxes, perturb = unlist(pF[c("O2", "ODU")]))
  else
     Fluxes <- rbind(Fluxes, perturb = c(0, 0))
  Fluxes           <- rbind(Fluxes, Fluxes[1,] - Fluxes[2,] + Fluxes[3,])                     
  rownames(Fluxes) <- c("surface", "bottom", "perturb", "netin")

  Rates <- data.frame(
    Nitrification      = out$TotNitri1*1.5 + out$TotNitri2*0.5,
    ODUoxidation       = out$TotODUoxid,
    ODUoxid.dist       = out$TotODUoxsurf,
    OxicMineralisation = out$TotOxic,
    MPBO2production    = out$TotMPBO2prod,
    MPBNO3reduction    = 2*out$TotMPBNO3uptake
    
    )

  if (includeMPB){
    Rates$MPBO2respiration   = out$TotMPBO2uptake
  }
  else
    Rates$MPBO2respiration <- 0
  
  Rfac <- c(-1,-1,-1,-1,+1, 0,-1)

  Rates$Total = sum(Rates * Rfac)   
  rownames(Rates) <- "nmolO2/cm2/d" 
  
  mat <- matrix (nrow = 7, ncol = 7, data = 0)
  rownames(mat) <- colnames(mat) <- c("Ext", "O2", "NO2", "NO3", "DIC", "Oxidant", "Burial")
  mat["Ext", "O2"] <- out$O2flux * (out$O2flux > 0)      
  mat["O2", "Ext"] <- -out$O2flux * (out$O2flux < 0)      
  mat["O2", "Burial"] <- out$O2deepflux
  mat["O2", "DIC"] <- Rates$OxicMineralisation +  Rates$MPBO2respiration
  mat["O2", "NO2"] <- out$TotNitri1*1.5
  mat["O2", "NO3"] <- out$TotNitri2*0.5
  mat["O2", "Oxidant"] <- Rates$ODUoxidation + Rates$ODUoxid.dist
  mat["DIC", "O2"] <- Rates$MPBO2production

  if (is.null(dC))
  dC <- (out$O2flux - out$O2deepflux + Rates$MPBO2production 
              - Rates$Nitrification - Rates$ODUoxidation 
              - Rates$ODUoxid.dist - Rates$OxicMineralisation - Rates$MPBO2respiration)
  else 
    dC <- c(O2 = dC$TotalO2)

#    Perturb = pF
#  if (! is.null(Perturb)){
#    Perturb <- pF[c("O2", "ODU")]
#    Perturb <- c(Perturb, sum = sum(Perturb))
#    } 

  return(list(
     Fluxes = Fluxes, Rates = Rates, Losses = out$O2deepflux, dC = dC, #Perturb = Perturb,
     Delta = Cflux + Rates$Total, Fluxmat = mat))                 
 
}

## -----------------------------------------------------------------------------

#====================#
#   C budgets        #
#====================#

CNPDIAbudgetC_one <- function(out) {

  Parms <- CNPDIAparms(out, TRUE)  

  pertFlux <- attributes(out)$perturbFlux
  pF <- NULL
  if (! is.null(pertFlux)){
   pF <- colSums(pertFlux)/diff(range(out[,1]))
  }
    
  if ("deSolve" %in% class(out))
    includeMPB <- attributes(out)$includeMPB
  else
    includeMPB <- out$includeMPB

  dV <- NULL
  if ("deSolve" %in% class(out)){
    dV <- dVal(out, includeMPB)
    out <- getMean0D(out, includeMPB)
  }

  influx <- out$DICflux + out$OrgCflux 
  Fluxes <- data.frame(FDET    = c(out$FDETflux, out$FDETdeepflux), 
                       SDET    = c(out$SDETflux, out$SDETdeepflux), 
                       DIC     = c(out$DICflux,  out$DICdeepflux),
                       CinCaP  = c(out$CaPflux,  out$CaPdeepflux)*Parms[["CPrCaP"]])
  if (includeMPB)
    Fluxes <- cbind(Fluxes, MPBC = c(out$MPBCflux, out$MBCdeepflux),
                     EPS = c(out$EPSflux, out$EPSdeepflux))                       
  
  if(! is.null(pF)) Fluxes <- rbind(Fluxes, 
    perturb = unlist(pF[c("FDET", "SDET", "DIC", "CaP")])*c(1,1,1,Parms[["CPrCaP"]]))
  else
     Fluxes <- rbind(Fluxes, perturb = c(0, 0, 0, 0))

  Fluxes           <- rbind(Fluxes, Fluxes[1,] - Fluxes[2,] + Fluxes[3,])                     
  rownames(Fluxes) <- c("surface", "bottom", "perturb", "netin")

    Fluxes$Total <- rowSums(Fluxes)

  influx <- Fluxes[1,"Total"]
  burial <- Fluxes[2,"Total"]
  Cflux   <- influx - burial

  Rates <- data.frame(
    OxicMineralisation    = out$TotOxic,
    Denitrification       = out$TotDenit,
    AnoxicMineralisation  = out$TotAnoxic,
    TotalMineralisation   = out$TotMin,   
    CaPprecipitation      = out$TotCaPprod*Parms[["CPrCaP"]],
    CaPdissolution        = out$TotCaPdiss*Parms[["CPrCaP"]],
    MPBDICuptake          = out$TotMPBDICuptake
  )
    fac <- c(0,0,0,1,-1,-1,-1)
  if (includeMPB){
    Rates$MPBResp = out$TotMPBresp
    Rates$MPBCdeath = out$TotMPBCdeath
    Rates$EPSproduction = out$TotProdEPS
    Rates$EPSmineralisation = out$TotMinEps
    Rates$FDETprodMPBdeath = out$TotFDETprodMPBdeath
    Rates$DICprodMPBdeath = out$TotDICprodMPBdeath
    fac <- c(fac,1,0,0,0,-1,1)  
  } else {
    Rates$MPBResp = 0
    Rates$MPBCdeath = out$TotMPBDICuptake
    Rates$EPSproduction = 0
    Rates$EPSmineralisation = 0
    Rates$FDETprodMPBdeath = out$TotMPBDICuptake
    Rates$DICprodMPBdeath = 0
      
    }
  rownames(Rates) <- "nmolC/cm2/d"  

  # derivatives
  mat <- matrix (nrow = 7, ncol = 7, data = 0)
  rownames(mat) <- colnames(mat) <- c("Ext", "DET", "DIC", "CaP", "EPS", "MPB", "Burial")
  mat["Ext", "DET"] <- out$FDETflux      + out$SDETflux
  mat["DET", "Burial"] <- out$FDETdeepflux  + out$SDETdeepflux
  mat["Ext", "DIC"]    <-  out$DICflux * (out$DICflux > 0) 
  mat["DIC", "Ext"]    <- -out$DICflux * (out$DICflux < 0)
  mat["Burial", "DIC"] <- - out$DICdeepflux*(out$DICdeepflux < 0)
  mat["DIC", "Burial"] <- + out$DICdeepflux*(out$DICdeepflux > 0) 
  mat["Ext", "CaP"] <- 0
  mat["CaP", "Ext"] <- out$CaPdeepflux*Parms[["CPrCaP"]] 
  if (includeMPB){
   mat["Ext", "EPS"] <- 0
   mat["EPS", "Burial"] <- out$EPSdeepflux*(out$EPSdeepflux > 0) 
  }  
  mat["MPB", "DET"] <- Rates$FDETprodMPBdeath 
  mat["DET", "DIC"] <- Rates$TotalMineralisation
  mat["DIC", "MPB"] <- Rates$MPBDICuptake
  mat["DIC", "CaP"] <- Rates$CaPprecipitation
  mat["CaP", "DIC"] <- Rates$CaPdissolution
  mat["MPB", "DIC"] <- Rates$MPBResp + Rates$DICprodMPBdeath
  mat["DIC", "EPS"] <- Rates$EPSproduction
  mat["EPS", "DIC"] <- Rates$EPSmineralisation

  # derivatives
  if (is.null(dV)){
   dDETRITUSC <- (out$FDETflux + out$SDETflux 
              - out$FDETdeepflux - out$SDETdeepflux 
              + Rates$FDETprodMPBdeath - Rates$TotalMineralisation)
  dCaP <- (out$CaPflux - out$CaPdeepflux 
              + Rates$CaPprecipitation - Rates$CaPdissolution)
   if (includeMPB)
   dMPBC <- (out$MPBCflux - out$MPBCdeepflux 
           + Rates$MPBDICuptake - Rates$MPBResp 
           - Rates$MPBCdeath - Rates$EPSproduction)

   else dMPBC <- 0   
   dDIC <- (out$DICflux - out$DICdeepflux + Rates$DICprodMPBdeath
           - Rates$MPBDICuptake + Rates$MPBResp + Rates$TotalMineralisation 
       + Rates$EPSmineralisation)
   dEPS <- (out$EPSflux - out$EPSdeepflux 
          + Rates$EPSproduction - Rates$EPSmineralisation)

   dC <- c(DET = dDETRITUSC, DIC = dDIC, CaP = dCaP, MPBC = dMPBC, EPS = dEPS)
  } else{
   dC <- c(DetritusC = dV$TotalFDET+dV$TotalSDET, DIC = dV$TotalDIC, CaP = dV$TotalCaP)
   if (includeMPB)
     dC <- c(dC, MPBC = dV$TotalMPBC, EPS = dV$TotalEPS)
   }
#  Perturb = pF
#  if (! is.null(Perturb)){
#    Perturb <- pF[c("FDET", "SDET", "DIC", "CaP")]*c(1,1,1,Parms[["CPrCaP"]])
#    Perturb <- c(Perturb, sum = sum(Perturb))
#    } 


  return(list(Fluxes = Fluxes, Rates = Rates, dC = c(dC, sum = sum(dC)), #Perturb = Perturb,
       Losses = burial, Delta = influx- burial, Fluxmat = mat))                 
 
}

## -----------------------------------------------------------------------------

#====================#
#   N budgets        #
#====================#

CNPDIAbudgetN_one <- function(out) {
  Parms <- CNPDIAparms(out, TRUE)  
  pertFlux <- attributes(out)$perturbFlux
  pF <- NULL
  if (! is.null(pertFlux))
   pF <- colSums(pertFlux)/diff(range(out[,1]))
  

  if ("deSolve" %in% class(out))
    includeMPB <- attributes(out)$includeMPB
  else
    includeMPB <- out$includeMPB

  dV <- NULL
  if ("deSolve" %in% class(out)){
    dV <- dVal(out, includeMPB)
    out <- getMean0D(out, includeMPB)
  }

  Fluxes <- data.frame(FDET_N = c(out$FDETflux    *Parms[["NCrFdet"]] , 
                                  out$FDETdeepflux*Parms[["NCrFdet"]] ), 
                       SDET_N = c(out$SDETflux    *Parms[["NCrSdet"]] , 
                                  out$SDETdeepflux*Parms[["NCrSdet"]]) , 
                       NO3  = c(out$NO3flux, out$NO3deepflux),
                       NO2  = c(out$NO2flux, out$NO2deepflux),
                       NH3  = c(out$NH3flux, out$NH3deepflux),
                       N2   = c(-out$TotDenit*0.8 -2*out$TotAnammox, 0)  )

  if (includeMPB)
    Fluxes <- cbind(Fluxes, MPBN = c(out$MPBNflux, out$MBNdeepflux))

  if(! is.null(pF)) Fluxes <- rbind(Fluxes, perturb = unlist(pF[c("FDET", "SDET", "NO3", "NO2", "NH3", "NO2")])*
                                       c(Parms[["NCrFdet"]], Parms[["NCrSdet"]],1, 1, 1,0))
  else 
       Fluxes <- rbind(Fluxes, perturb = c(0, 0, 0, 0, 0, 0))

  Fluxes           <- rbind(Fluxes, Fluxes[1,] - Fluxes[2,] + Fluxes[3,])                     
  rownames(Fluxes) <- c("surface", "bottom", "perturb", "netin")
  Fluxes$Total <- rowSums(Fluxes)

  influx <- Fluxes[1,"Total"]
  burial <- Fluxes[2,"Total"]
  Nflux   <- influx - burial

  Rates <- data.frame(
    NH3production      = out$TotNH3prod,
    Denitrification    = out$TotDenit*0.8,
    MPBNO3consumption  = out$TotMPBNO3uptake,
    MPBNH3consumption  = out$TotMPBNH3uptake,
    Nitrification1     = out$TotNitri1,
    Nitrification2     = out$TotNitri2,
    Anammox            = out$TotAnammox,
    NH3adsorption      = out$TotNH3ads,
    N2production       = out$TotDenit*0.8 + 2*out$TotAnammox)
  if (includeMPB){
    Rates$MPBNdeath        = out$TotMPBNdeath
    Rates$NH3prodMPBdeath  = out$TotNH3prodMPBdeath
    Rates$DETNprodMPBdeath = out$TotFDETprodMPBdeath*Parms[["NCrFdet"]]

  } else {
    Rates$MPBNdeath = out$TotMPBNO3uptake + out$TotMPBNH3uptake
    Rates$NH3prodMPBdeath = 0
    Rates$DETNprodMPBdeath = out$TotMPBNO3uptake + out$TotMPBNH3uptake
      
    }

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

  if (is.null(dV)){
  dDETRITUSN <- (out$FDETflux*Parms[["NCrFdet"]] + out$SDETflux*Parms[["NCrSdet"]] 
              - out$FDETdeepflux*Parms[["NCrFdet"]] - out$SDETdeepflux *Parms[["NCrSdet"]]
              + Rates$DETNprodMPBdeath - Rates$NH3production)
  dNH3 <- (out$NH3flux - out$NH3deepflux 
              + Rates$NH3production - Rates$MPBNH3consumption - Rates$Nitrification1 - 
                Rates$Anammox + Rates$NH3prodMPBdeath)
  dNO3 <- (out$NO3flux - out$NO3deepflux 
              + Rates$Nitrification2 - Rates$MPBNO3consumption - Rates$Denitrification)
  dNO2 <- (out$NO2flux - out$NO2deepflux 
              + Rates$Nitrification1 - Rates$Nitrification2 - Rates$Anammox)
  if (includeMPB)
  dMPBN <- (out$MPBNflux - out$MPBNdeepflux 
           + Rates$MPBNH3consumption + Rates$MPBNO3consumption 
           - Rates$NH3prodMPBdeath - Rates$DETNprodMPBdeath)
  else dMPBN <- 0   

  dC <- c(DET = dDETRITUSN, NH3 = dNH3, NO3 = dNO3, NO2 = dNO2, MPBN = dMPBN)
  } else{
    dC = c(DET = dV$TotalFDET*Parms[["NCrFdet"]] + dV$TotalSDET*Parms[["NCrSdet"]], 
      NH3 = dV$TotalNH3*(1+Parms[["NH3Ads"]]), NO3 = dV$TotalNO3, NO2 = dV$TotalNO2)
    if (includeMPB) dC <- c(dC, MPBN = dV$TotalMPBN)
  }
#  Perturb = pF
#  if (! is.null(Perturb)){
#    Perturb <- Perturb[c("FDET", "SDET", "NH3", "NO2", "NO3", "NH3ads")]*c(Parms[["NCrFdet"]],Parms[["NCrSdet"]],1,1,1,1)
#    Perturb <- c(Perturb, sum = sum(Perturb))
#    }
  return(list(Fluxes = Fluxes, Rates = Rates, 
     Losses = burial + out$TotDenit*0.8 + 2*out$TotAnammox, #Perturb = Perturb, 
     dC = c(dC, sum = sum(dC)), Delta = Nflux, Fluxmat = mat))                 

}

## --------------------------------------------------------------------------------------------------

#====================#
#   P budgets        #
#====================#

CNPDIAbudgetP_one <- function(out) {

  pertFlux <- attributes(out)$perturbFlux
  pF <- NULL
  if (! is.null(pertFlux)){
   pF <- colSums(pertFlux)/diff(range(out[,1]))
  }
  if ("deSolve" %in% class(out))
    includeMPB <- attributes(out)$includeMPB
  else
    includeMPB <- out$includeMPB

  if (includeMPB)
    Parms <- MPBDIAparms(out, TRUE)  
  else  
    Parms <- CNPDIAparms(out, TRUE)  
  
  dV <- NULL
  if ("deSolve" %in% class(out)){
    dV <- dVal(out, includeMPB)
    out <- getMean0D(out, includeMPB)
  }

    Fluxes <- data.frame(FDET_P = c(out$FDETflux      *Parms[["PCrFdet"]] , 
                                  out$FDETdeepflux  *Parms[["PCrFdet"]] ), 
                       SDET_P = c(out$SDETflux      *Parms[["PCrSdet"]] , 
                                  out$SDETdeepflux *Parms[["PCrSdet"]]) , 
                       PO4  = c(out$PO4flux, out$PO4deepflux),
                       FeP  = c(0, out$FePdeepflux),
                       CaP  = c(0, out$CaPdeepflux))
  if (includeMPB)
    Fluxes <- cbind(Fluxes, MPBP = c(out$MPBCflux*Parms[["PCrFdet"]], out$MBCdeepflux*Parms[["PCrFdet"]]))

  if(! is.null(pF)) Fluxes <- rbind(Fluxes, perturb = unlist(pF[c("FDET", "SDET", "PO4", "FeP", "CaP")])*
                                       c(Parms[["PCrFdet"]], Parms[["PCrSdet"]],1, 1, 1))
  else      Fluxes <- rbind(Fluxes, perturb = c(0, 0, 0, 0, 0))
  Fluxes           <- rbind(Fluxes, Fluxes[1,] - Fluxes[2,] + Fluxes[3,])                     
  rownames(Fluxes) <- c("surface", "bottom", "perturb", "netin")

  Fluxes$Total <- rowSums(Fluxes)

  influx <- Fluxes[1,"Total"]
  burial <- Fluxes[2,"Total"]
  Pflux   <- influx - burial

  Rates <- data.frame(
    PO4production   = out$TotPO4prod,
    FePadsorption   = out$TotFePprod,
    FePdesorption   = out$TotFePdesorp,
    MPBPO4consumption  = out$TotMPBPO4uptake,
    CaPproduction   = out$TotCaPprod,
    CaPdissolution  = out$TotCaPdiss)                   
  if (includeMPB){
    Rates$MPBPdeath = out$TotMPBCdeath*Parms[["PCrFdet"]]
    Rates$PO4prodMPBdeath = out$TotPO4prodMPBdeath
    Rates$DETPprodMPBdeath = out$TotFDETprodMPBdeath*Parms[["PCrFdet"]]
  } else {
    Rates$MPBPdeath = Rates$MPBPO4consumption
    Rates$PO4prodMPBdeath = 0
    Rates$DETPprodMPBdeath = Rates$MPBPO4consumption
      
  }

  rownames(Rates) <- "nmolP/cm2/d" 

  # derivatives
  # derivatives
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

  if (is.null(dV)){
  dDETRITUSP <- (out$FDETflux*Parms[["PCrFdet"]] + out$SDETflux*Parms[["PCrSdet"]] 
              - out$FDETdeepflux*Parms[["PCrFdet"]] - out$SDETdeepflux *Parms[["PCrSdet"]]
              + Rates$DETPprodMPBdeath - Rates$PO4production)
  dPO4 <- (out$PO4flux - out$PO4deepflux 
              + Rates$PO4production - Rates$FePadsorption + Rates$FePdesorption 
              - Rates$CaPproduction + Rates$CaPdissolution - Rates$MPBPO4consumption 
              + Rates$PO4prodMPBdeath)
  dFeP <- (out$FePflux - out$FePdeepflux 
              + Rates$FePadsorption - Rates$FePdesorption)
  dCaP <- (out$CaPflux - out$CaPdeepflux 
              + Rates$CaPproduction - Rates$CaPdissolution)
  if (includeMPB)
  dMPBP <- (out$MPBCflux*Parms[["PCrFdet"]] - out$MPBCdeepflux*Parms[["PCrFdet"]] 
           - Rates$PO4prodMPBdeath + Rates$MPBPO4consumption
           - Rates$MPBNdeath*Parms[["PCrFdet"]]/Parms[["NCrFdet"]])
  else dMPBP <- 0   

  dC = c(DET = dDETRITUSP, PO4 = dPO4, FeP = dFeP, CaP = dCaP, MPBP = dMPBP)
  } else{
    dC = c(DET = dV$TotalFDET*Parms[["PCrFdet"]] + dV$TotalSDET*Parms[["PCrSdet"]], 
      PO4 = dV$TotalPO4, FeP = dV$TotalFeP, CaP = dV$TotalCaP)
    if (includeMPB) dC <- c(dC, MPBP = dV$TotalMPBN*Parms[["PCrFdet"]]/Parms[["NCrFdet"]])
  }

#  Perturb = pF
#  if (! is.null(pF)){
#    Perturb <- pF[c("FDET", "SDET", "PO4", "FeP", "CaP")]*c(Parms[["PCrFdet"]],Parms[["PCrSdet"]],1,1,1)
#    Perturb <- c(Perturb, sum = sum(Perturb))

#    }
  return(list(Fluxes = Fluxes, Rates = Rates, Losses = burial, dC = c(dC, sum = sum(dC)), #Perturb = Perturb, 
     Delta = Pflux , Fluxmat = mat))                 

}

# in progress

subset.CNPDIAstd <- function(x, subset = NULL, select = NULL, which = select, ...) 
{
   # deparse state variables 
    xold <- x
    x <- list()
    for (i in 1:ncol(xold$y))
      x[[i]] <- xold$y[,i]
    names(x) <- colnames(xold$y)
    x <- c(x, xold[-1])
    x$Grid <- x$porGrid <- x$model <- x$numberTries <- x$warnings <- NULL
    x$other <- x$initpar <- x$Parms <- NULL
    len <- unlist(lapply(x, FUN = length))
    
    x <- as.data.frame(x)
    
    Which <- which
    if (missing(subset)) 
        r <- TRUE
    else {
        e <- substitute(subset)
        r <- eval(e, as.data.frame(x), parent.frame())
        if (is.numeric(r)) {
            isub <- r
        }
        else {
            if (!is.logical(r)) 
                stop("'subset' must evaluate to logical or be a vector with integers")
            r <- r & !is.na(r)
        }
    }
    if (is.null(Which)) 
        return(xold$y[r, ])
    return(x[r, Which])
}

