
## --------------------------------------------------------------------------------------------------
## Perturb functions
## --------------------------------------------------------------------------------------------------

#### PERTURBATIONS ####

## --------------------------------------------------------------------------------------------------
mixSolid <- function(C, depth = 5, grid = setup.grid.1D(N = 100, L  = 30), por = 0.5) {
   N.Pert         <- length(grid$x.int[grid$x.int < depth])
   if (N.Pert == 0) return(list(C = C, Flux = 0))
   if (is.vector(por))
     por <- list(mid = rep(por, length.out = length(C)))
   Cnew <- MixDET(C, N.Pert, por, grid)
   Flux <- integrateSolid(Cnew, grid, por) - integrateSolid(C, grid, por)
   list(C = Cnew, Flux = Flux)
}

## --------------------------------------------------------------------------------------------------
# water mixed in the sediment has the overlying water concentration
mixLiquid <- function(C, depth = 5, grid = setup.grid.1D(N = 100, 
                      L  = 30), por = 0.5, BW) {
   N.Pert         <- length(grid$x.int[grid$x.int < depth])
   if (N.Pert == 0) return(list(C = C, Flux = 0))
   Cnew <- C
   Cnew[1:N.Pert] <- BW
   Flux <- integrateLiquid(Cnew, grid, por) - integrateLiquid(C, grid, por)
   list(C = Cnew, Flux = Flux)
}

## --------------------------------------------------------------------------------------------------
# water mixed in the sediment has a mix of the the overlying and interstitial concentration
mixLiquid_BW <- function(C, depth = 5, grid = setup.grid.1D(N = 100, L  = 30), por = 0.5, BW, Hwater) {
  N.Pert   <- length(grid$x.int[grid$x.int < depth])
  if (N.Pert == 0) return(list(C = C, Flux = 0))
  if (is.list(por))
    por <- por$mid
  if (is.list(grid))
    grid <- grid$dx  
  
  por <- rep(por, length.out = 100)               # expand to vector
  Cnew <- C
  ii <- 1:N.Pert
  pertvol <- sum(grid[ii]*por[ii])  # volume in porewater to exchange
  pVol   <- pertvol /(pertvol + Hwater)             # relative to total
  Conc <- (sum(Cnew[ii] + grid[ii] * por[ii])*pVol + 
          (BW * Hwater)*(1-pVol))
  Cnew[ii] <- Conc
  Flux <- integrateLiquid(Cnew, grid, por) - integrateLiquid(C, grid, por)
  list(C = Cnew, Flux = Flux)
}

## --------------------------------------------------------------------------------------------------
erodeSolid <- function(C, depth = 5, grid = setup.grid.1D(N = 100, L  = 30), por = 0.5) {
   N.Pert         <- length(grid$x.int[grid$x.int < depth])
   if (N.Pert == 0) return(list(C = C, Flux = 0))
   if (is.vector(por))
     por <- list(mid = rep(por, length.out = length(C)))
   Cnew <- Erode(C, N.Pert, por, grid) 
   Flux <- integrateSolid(Cnew, grid, por) - integrateSolid(C, grid, por)
   list(C = Cnew, Flux = Flux)
}

## --------------------------------------------------------------------------------------------------
erodeLiquid <- function(C, depth = 5, grid = setup.grid.1D(N = 100, L  = 30), por = 0.5) {
   Cnew <- erodeSolid (C, depth = depth, grid = grid, por = por)$C
   Flux <- integrateLiquid(Cnew, grid, por) - integrateLiquid(C, grid, por)
   list(C = Cnew, Flux = Flux)
  
}  

## --------------------------------------------------------------------------------------------------
depositSolid <- function(C, depth = 5, grid = setup.grid.1D(N = 100, L  = 30), por = 0.5,
  concfac = 1) {
   N.Pert         <- length(grid$x.int[grid$x.int < depth])
   if (N.Pert == 0) return(list(C = C, Flux = 0))
   if (is.vector(por))
     por <- list(mid = rep(por, length.out = length(C)))
   Cnew <- DepositDET(C, N.Pert, por, grid, concfac)
   Flux <- integrateSolid(Cnew, grid, por) - integrateSolid(C, grid, por)
   list(C = Cnew, Flux = Flux)
   
}

## --------------------------------------------------------------------------------------------------
depositLiquid <- function(C, depth = 5, grid = setup.grid.1D(N = 100, L  = 30), por = 0.5, BW) {
   N.Pert         <- length(grid$x.int[grid$x.int < depth])
   if (N.Pert == 0) return(list(C = C, Flux = 0))
   Cnew <- DepositBW(C, N.Pert, BW, grid)
   Flux <- integrateLiquid(Cnew, grid, por) - integrateLiquid(C, grid, por)
   list(C = Cnew, Flux = Flux)
}

#===============================================================================
# to estimate integrated concentrations
#===============================================================================

## --------------------------------------------------------------------------------------------------
integrateSolid <- function(C, grid, por) {
    if (is.list(por))
      por <- por$mid
    if (is.list(grid))
      grid <- grid$dx  
    sum(C * (1.-por) * grid)
}

## --------------------------------------------------------------------------------------------------
integrateLiquid <- function(C, grid, por) {
    if (is.list(por))
      por <- por$mid
    if (is.list(grid))
      grid <- grid$dx  
    sum(C * por * grid)
}    

# internal functions
IntegrateSol <- function(FDET, N.Pert, porGrid, Grid)
    sum(FDET[1:N.Pert] * (1.-porGrid$mid[1:N.Pert]) * Grid$dx[1:N.Pert])

IntegrateLiq <- function(CONC, N.Pert,porGrid, Grid)
    sum(CONC[1:N.Pert] * porGrid$mid[1:N.Pert] * Grid$dx[1:N.Pert])

#===============================================================================
# Function to mix detritus
#===============================================================================

MixDET <- function (FDET, N.Pert, porGrid, Grid) {
  # Integrated concentration
   TotalFdet <- sum(FDET[1:N.Pert] *
                   (1.-porGrid$mid[1:N.Pert]) * Grid$dx[1:N.Pert])
  # approximate mean concentration
   MeanFdet <- TotalFdet / sum(Grid$dx[1:N.Pert])

   TotalFdet2 <- sum(MeanFdet *
                  (1.-porGrid$mid[1:N.Pert]) * Grid$dx[1:N.Pert])

   if (TotalFdet2 > 0 )
     Fac <- TotalFdet / TotalFdet2
   else
     Fac <- 1
   
   # Mixed concentration
   FDET[1:N.Pert] <- MeanFdet*Fac
   return(FDET)
}

#===============================================================================
# Deposition of detritus on top of sediment
#===============================================================================

DepositDET <- function (FDET, N.Pert, porGrid, Grid, fac) {
  # Integrated concentration
   TotalFdet <- sum(FDET[1:N.Pert] *
                   (1.-porGrid$mid[1:N.Pert]) * Grid$dx[1:N.Pert])
  # approximate mean concentration
   MeanFdet <- TotalFdet / sum(Grid$dx[1:N.Pert])
   TotalFdet2 <- sum(MeanFdet *
                  (1.-porGrid$mid[1:N.Pert]) * Grid$dx[1:N.Pert])

   if (TotalFdet2 > 0 )
     Fac <- TotalFdet / TotalFdet2* fac
   else
     Fac <- 1

    NewGrid <- c(Grid$x.mid[1:N.Pert], Grid$x.mid+Grid$x.mid[N.Pert+1])
    FDETerode <- c(rep(MeanFdet*Fac, N.Pert), FDET)
    FDETn <- approx(x = NewGrid, y = FDETerode, xout = Grid$x.mid)$y

   return(FDETn)
}

#===============================================================================
# Deposition for dissolved substances -  interstitial conc = bottom water conc
#===============================================================================

DepositBW <- function (O2, N.Pert, bwO2, Grid) {
    NewGrid <- c(Grid$x.mid[1:N.Pert], Grid$x.mid+Grid$x.mid[N.Pert+1])
    O2erode <- c(rep(bwO2, N.Pert), O2)
    O2n <- approx(x = NewGrid, y = O2erode, xout = Grid$x.mid)$y

   return(O2n)
}

#===============================================================================
# Erosion of the top of the sediment
#===============================================================================

Erode <- function (FDET, N.Pert, porGrid, Grid) {

   NewGrid <- c(Grid$x.mid[-(1:N.Pert)]-Grid$x.mid[N.Pert+1], Grid$x.mid[Grid$N])
   FDETerode <- c(FDET[-(1:N.Pert)], FDET[Grid$N])
   FDETn <- approx(x = NewGrid, y = FDETerode, xout = Grid$x.mid)$y

   return(FDETn)
}

## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------
## Main CNPDIA perturbation function
## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------

CNPDIAperturb <- function (parms = list(), times = 0:365, spinup = NULL, 
               yini = NULL, gridtype = 1, Grid = NULL, porosity = NULL, 
               bioturbation = NULL, irrigation = NULL, surface = NULL, 
               diffusionfactor = NULL,  dynamicbottomwater = FALSE,   
               perturbType = "mix", perturbTimes =  seq(from = 0, to = max(times), by = 365), 
               perturbDepth = 5, concfac = 1, 
               CfluxForc   = NULL, FePfluxForc = NULL, CaPfluxForc = NULL, 
               O2bwForc    = NULL, NO3bwForc   = NULL, NO2bwForc   = NULL, 
               NH3bwForc   = NULL, ODUbwForc   = NULL, PO4bwForc   = NULL, 
               DICbwForc   = NULL, 
               wForc       = NULL, biotForc    = NULL, irrForc     = NULL, 
               rFastForc   = NULL, rSlowForc   = NULL, pFastForc   = NULL, 
               MPBprodForc = NULL, gasfluxForc = NULL, HwaterForc  = NULL, 
               ratefactor = NULL, verbose = FALSE, ...) {

## check parameter inputs
  model <- 1
  if ( dynamicbottomwater) model <- 2
  
  perttype <- match.arg(perturbType, c("mix", "erode","deposit"), several.ok = TRUE)
  
  depthPert <- rep(perturbDepth, length.out = length(perttype))
  names(depthPert) <- as.character(perttype)
  numPert <- length(depthPert)
  
  CfluxForc    <- checkforcs(CfluxForc,     "CfluxForc")
  O2bwForc     <- checkforcs(O2bwForc ,      "O2bwForc")
  NO3bwForc    <- checkforcs(NO3bwForc,     "NO3bwForc")
  NO2bwForc    <- checkforcs(NO2bwForc,     "NO2bwForc")
  NH3bwForc    <- checkforcs(NH3bwForc,     "NH3bwForc")
  ODUbwForc    <- checkforcs(ODUbwForc,     "ODUbwForc")
  PO4bwForc    <- checkforcs(PO4bwForc,     "PO4bwForc")
  DICbwForc    <- checkforcs(DICbwForc,     "DICbwForc")
  wForc        <- checkforcs(wForc,             "wForc")
  biotForc     <- checkforcs(biotForc,       "biotForc")
  irrForc      <- checkforcs(irrForc,         "irrForc")
  rFastForc    <- checkforcs(rFastForc,     "rFastForc")
  rSlowForc    <- checkforcs(rSlowForc,     "rSlowForc")
  pFastForc    <- checkforcs(pFastForc,     "pFastForc")
  MPBprodForc  <- checkforcs(MPBprodForc, "MPBprodForc")
  FePfluxForc  <- checkforcs(FePfluxForc, "FePfluxForc")
  CaPfluxForc  <- checkforcs(CaPfluxForc, "CaPfluxForc")
  gasfluxForc  <- checkforcs(gasfluxForc, "gasfluxForc")
  HwaterForc   <- checkforcs(HwaterForc,   "HwaterForc")
  ratefactor   <- checkforcs(ratefactor,   "ratefactor")
  
  Hlist <- list(Hwater = 1)
  
  if (is.null(spinup)) Times <- times else Times <- spinup
  
  if (is.null(yini)) {
    STD <- CNPDIAsolve_full(parms, gridtype = gridtype, CfluxForc = CfluxForc,
                       CaPfluxForc = CaPfluxForc, FePfluxForc = FePfluxForc,
                       O2bwForc  = O2bwForc, NO3bwForc = NO3bwForc, 
                       NO2bwForc = NO2bwForc, NH3bwForc = NH3bwForc,
                       ODUbwForc = ODUbwForc, PO4bwForc = PO4bwForc, 
                       DICbwForc = DICbwForc, gasfluxForc = gasfluxForc, 
                       wForc = wForc, biotForc = biotForc, irrForc = irrForc, 
                       ratefactor = ratefactor, rFastForc = rFastForc,
                       rSlowForc = rSlowForc, pFastForc = pFastForc, 
                       MPBprodForc = MPBprodForc, HwaterForc = HwaterForc,
                       times = Times, Grid = Grid, porosity = porosity, 
                       bioturbation = bioturbation, irrigation = irrigation, 
                       surface = surface, diffusionfactor = diffusionfactor, 
                       model = model, verbose = verbose)
    yini <- STD$y
  } else  
    STD <- initCNPDIA(parms, gridtype = gridtype, CfluxForc = CfluxForc, 
                      CaPfluxForc = CaPfluxForc, FePfluxForc = FePfluxForc,
                      O2bwForc  = O2bwForc, NO3bwForc = NO3bwForc, 
                      NO2bwForc = NO2bwForc, NH3bwForc = NH3bwForc,
                      ODUbwForc = ODUbwForc, PO4bwForc = PO4bwForc, 
                      DICbwForc = DICbwForc, gasfluxForc = gasfluxForc, 
                      wForc = wForc, biotForc = biotForc, irrForc = irrForc, 
                      ratefactor = ratefactor, rFastForc = rFastForc, 
                      rSlowForc = rSlowForc, pFastForc = pFastForc, 
                      MPBprodForc = MPBprodForc, HwaterForc = HwaterForc,
                      times = Times, Grid = Grid, porosity = porosity, 
                      bioturbation = bioturbation, irrigation = irrigation, 
                      surface = surface, diffusionfactor = diffusionfactor, 
                      model = model)
  
  Grid <- STD$Grid
  porGrid <- STD$porGrid
  
  pertCONC <- yini
  PertFlux <- PertBW <- PertBW2 <- NULL 
  events <- NULL

  Parms <- STD$Parms
  PP <- unlist(Parms)

# Number of layers that are perturbed 
# Note: there can be more than one type of perturbation 
  
  N.Pert <- rep(0, times = length(depthPert))
  
  for (i in 1:length(depthPert))
  N.Pert[i] <- length(Grid$x.int[Grid$x.int < depthPert[i]])
  names(N.Pert) <- as.character(perttype)

  nspec        <- 12

#===============================================================================
# Function to Perturb all states
#===============================================================================

   Perturb <- function (out, t) {
    
    ii <- which (N.Pert > 0)

    # integrated concentrations over the perturbed sediment layer
    if ("steady1D" %in% class(out))
      pertCONC <- out$y
    else
      pertCONC <- matrix(ncol = nspec, out, byrow = FALSE)

    colnames(pertCONC) <- .CNPDIA$svar
    
    Fluxes <-  c(FDET= 0,  SDET= 0,   O2  = 0,  NO3 = 0,
                 NO2 = 0,  NH3 = 0,  ODU = 0,  DIC = 0,
                 PO4 = 0,  FeP = 0,  CaP = 0,  Pads= 0)

    if (model == 1){
      bw_O2  <- bwO2(t)
      bw_NO3 <- bwNO3(t)
      bw_NO2 <- bwNO2(t)
      bw_NH3 <- bwNH3(t)
      bw_ODU <- bwODU(t)
      bw_DIC <- bwDIC(t)
      bw_PO4 <- bwPO4(t)
    } else {     # dynamic bottom water concentrations
      BW     <- pertCONC[1,]
      pertCONC <- pertCONC[-1,]
      bw_O2  <- BW["O2"]
      bw_NO3 <- BW["NO3"]
      bw_NO2 <- BW["NO2"]
      bw_NH3 <- BW["NH3"]
      bw_ODU <- BW["ODU"]
      bw_DIC <- BW["DIC"]
      bw_PO4 <- BW["PO4"]
    }

    for (i in ii){
    
    O2Conc  <- IntegrateLiq(pertCONC[ ,"O2"] ,  .CNPDIA$N, porGrid, Grid)
    NO3Conc <- IntegrateLiq(pertCONC[ ,"NO3"],  .CNPDIA$N, porGrid, Grid)
    NO2Conc <- IntegrateLiq(pertCONC[ ,"NO2"],  .CNPDIA$N, porGrid, Grid)
    NH3Conc <- IntegrateLiq(pertCONC[ ,"NH3"],  .CNPDIA$N, porGrid, Grid)
    PO4Conc <- IntegrateLiq(pertCONC[ ,"PO4"],  .CNPDIA$N, porGrid, Grid)
    DICConc <- IntegrateLiq(pertCONC[ ,"DIC"],  .CNPDIA$N, porGrid, Grid)
    ODUConc <- IntegrateLiq(pertCONC[ ,"ODU"],  .CNPDIA$N, porGrid, Grid)
    FDETConc<- IntegrateSol(pertCONC[ ,"FDET"], .CNPDIA$N, porGrid, Grid)
    SDETConc<- IntegrateSol(pertCONC[ ,"SDET"], .CNPDIA$N, porGrid, Grid)
    FePConc <- IntegrateSol(pertCONC[ ,"FeP"],  .CNPDIA$N, porGrid, Grid)
    CaPConc <- IntegrateSol(pertCONC[ ,"CaP"],  .CNPDIA$N, porGrid, Grid)
    PadsConc <- IntegrateSol(pertCONC[ ,"Pads"],.CNPDIA$N, porGrid, Grid)
    
    if (perttype[i] == "mix") {    # mixed
      pertCONC[,"FDET"] <- MixDET (pertCONC[,"FDET"], N.Pert[i], porGrid, Grid)
      pertCONC[,"SDET"] <- MixDET (pertCONC[,"SDET"], N.Pert[i], porGrid, Grid)
      pertCONC[,"FeP"]  <- MixDET (pertCONC[,"FeP"] , N.Pert[i], porGrid, Grid)
      pertCONC[,"CaP"]  <- MixDET (pertCONC[,"CaP"] , N.Pert[i], porGrid, Grid)
      pertCONC[,"Pads"] <- MixDET (pertCONC[,"Pads"], N.Pert[i], porGrid, Grid)
      pertCONC[1:N.Pert[i], "O2"]  <- bw_O2 
      pertCONC[1:N.Pert[i], "NO3"] <- bw_NO3
      pertCONC[1:N.Pert[i], "NO2"] <- bw_NO2
      pertCONC[1:N.Pert[i], "NH3"] <- bw_NH3
      pertCONC[1:N.Pert[i], "ODU"] <- bw_ODU
      pertCONC[1:N.Pert[i], "DIC"] <- bw_DIC
      pertCONC[1:N.Pert[i], "PO4"] <- bw_PO4

     } else if (perttype[i] == "erode") {  # erosion
       pertCONC[,"FDET"] <- Erode (pertCONC[,"FDET"], N.Pert[i], porGrid, Grid)
       pertCONC[,"SDET"] <- Erode (pertCONC[,"SDET"], N.Pert[i], porGrid, Grid)
       pertCONC[,"FeP"]  <- Erode (pertCONC[,"FeP"] , N.Pert[i], porGrid, Grid)
       pertCONC[,"CaP"]  <- Erode (pertCONC[,"CaP"] , N.Pert[i], porGrid, Grid)
       pertCONC[,"Pads"] <- Erode (pertCONC[,"Pads"], N.Pert[i], porGrid, Grid)
       pertCONC[,"O2"]   <- Erode (pertCONC[,"O2"]  , N.Pert[i], porGrid, Grid)
       pertCONC[,"NO3"]  <- Erode (pertCONC[,"NO3"] , N.Pert[i], porGrid, Grid)
       pertCONC[,"NO2"]  <- Erode (pertCONC[,"NO2"] , N.Pert[i], porGrid, Grid)
       pertCONC[,"NH3"]  <- Erode (pertCONC[,"NH3"] , N.Pert[i], porGrid, Grid)
       pertCONC[,"ODU"]  <- Erode (pertCONC[,"ODU"] , N.Pert[i], porGrid, Grid)
       pertCONC[,"DIC"]  <- Erode (pertCONC[,"DIC"] , N.Pert[i], porGrid, Grid)
       pertCONC[,"PO4"]  <- Erode (pertCONC[,"PO4"] , N.Pert[i], porGrid, Grid)

     } else { # deposit
       pertCONC[,"FDET"] <- DepositDET (pertCONC[,"FDET"], N.Pert[i], porGrid, Grid, concfac)
       pertCONC[,"SDET"] <- DepositDET (pertCONC[,"SDET"], N.Pert[i], porGrid, Grid, concfac)
       pertCONC[,"FeP"]  <- DepositDET (pertCONC[,"FeP"] , N.Pert[i], porGrid, Grid, 1)
       pertCONC[,"CaP"]  <- DepositDET (pertCONC[,"CaP"] , N.Pert[i], porGrid, Grid, 1)
       pertCONC[,"Pads"] <- DepositDET (pertCONC[,"Pads"], N.Pert[i], porGrid, Grid, 1)
       pertCONC[,"O2"]   <- DepositBW  (pertCONC[,"O2"]  , N.Pert[i], bw_O2, Grid)
       pertCONC[,"NO3"]  <- DepositBW  (pertCONC[,"NO3"] , N.Pert[i], bw_NO3, Grid)
       pertCONC[,"NO2"]  <- DepositBW  (pertCONC[,"NO2"] , N.Pert[i], bw_NO3, Grid)
       pertCONC[,"NH3"]  <- DepositBW  (pertCONC[,"NH3"] , N.Pert[i], bw_NH3, Grid)
       pertCONC[,"ODU"]  <- DepositBW  (pertCONC[,"ODU"] , N.Pert[i], bw_ODU, Grid)
       pertCONC[,"DIC"]  <- DepositBW  (pertCONC[,"DIC"] , N.Pert[i], bw_DIC, Grid)
       pertCONC[,"PO4"]  <- DepositBW  (pertCONC[,"PO4"] , N.Pert[i], bw_PO4, Grid)
     }
    Fluxes <-  Fluxes + 
               c(FDET= -(FDETConc- IntegrateSol(pertCONC[ ,"FDET"],.CNPDIA$N, porGrid, Grid)),
                 SDET= -(SDETConc- IntegrateSol(pertCONC[ ,"SDET"],.CNPDIA$N, porGrid, Grid)),
                 O2  = -(O2Conc  - IntegrateLiq(pertCONC[ ,"O2"]  ,.CNPDIA$N, porGrid, Grid)),
                 NO3 = -(NO3Conc - IntegrateLiq(pertCONC[ ,"NO3"] ,.CNPDIA$N, porGrid, Grid)),
                 NO2 = -(NO2Conc - IntegrateLiq(pertCONC[ ,"NO2"] ,.CNPDIA$N, porGrid, Grid)),
                 NH3 = -(NH3Conc - IntegrateLiq(pertCONC[ ,"NH3"] ,.CNPDIA$N, porGrid, Grid)),
                 ODU = -(ODUConc - IntegrateLiq(pertCONC[ ,"ODU"] ,.CNPDIA$N, porGrid, Grid)),
                 DIC = -(DICConc - IntegrateLiq(pertCONC[ ,"DIC"] ,.CNPDIA$N, porGrid, Grid)),
                 PO4 = -(PO4Conc - IntegrateLiq(pertCONC[ ,"PO4"] ,.CNPDIA$N, porGrid, Grid)),
                 FeP = -(FePConc - IntegrateSol(pertCONC[ ,"FeP"] ,.CNPDIA$N, porGrid, Grid)),
                 CaP = -(CaPConc - IntegrateSol(pertCONC[ ,"CaP"] ,.CNPDIA$N, porGrid, Grid)),
                 Pads= -(PadsConc- IntegrateSol(pertCONC[ ,"Pads"],.CNPDIA$N, porGrid, Grid)))
                                       
    }
    
    PertFlux <<- rbind(PertFlux, Fluxes)
    
    if (model == 2){
      PertBW   <<- rbind(PertBW  , BW)
      BW <- BW - Fluxes[.CNPDIA$svar]/bwHeight(t)
      BW <- pmax(BW, 0)
      pertCONC <- rbind(BW, pertCONC)
      PertBW2 <<- rbind(PertBW2  , BW)
    }
    
    return( pertCONC)
  }

#===============================================================================
# Event Function
#===============================================================================

  EventFunc <- function(t, y, parms){
     print (paste("event at time", t))
     pertCONC <- Perturb(y, t)
     return (as.vector(pertCONC))
  }
                        
  initpar = unlist(STD$initpar)

  events <- list(func = EventFunc, time = perturbTimes)
  if (STD$isDistReact | model == 2) {
     band   <- 0
     lrw <- 190000
  } else  {
     band   <- 1
     lrw <- 90000
  }
  
  if (model == 1) {
    func <- "omexdiamodp"
    N <- .CNPDIA$N
  } else {
    func <- "omexdiamodbw"
    N <- .CNPDIA$N+1
  }  
  
  if(length(yini) != nspec*N)
    stop ("'yini' not of correct length, should be ", nspec*N)
  
  ZZ <- NULL
  is.spinup <- ! is.null(spinup)
  if (is.spinup)
    is.spinup <- (length(spinup) > 1)
  if (is.spinup) {
    if (model == 1) 
      HWATERForc <- Setforcings (Hlist, "Hwater",     HwaterForc, spinup, fac = 1)   # not used
    else 
      HWATERForc <- Setforcings (STD$Parms, "Hwater", HwaterForc, spinup, fac = 1)
    
    forcings <- list()
    forcings[[ 1]] <- Setforcings (STD$Parms, "Cflux",     CfluxForc, spinup, fac = 1)
    forcings[[ 2]] <- Setforcings (STD$Parms, "w",             wForc, spinup, fac = 1)
    forcings[[ 3]] <- Setforcings (STD$other, "biot",       biotForc, spinup, fac = 1/STD$other[["biot"]])
    forcings[[ 4]] <- Setforcings (STD$other, "irr",         irrForc, spinup, fac = 1/STD$other[["irr"]])
    forcings[[ 5]] <- Setforcings (STD$Parms, "rFast",     rFastForc, spinup, fac = 1)
    forcings[[ 6]] <- Setforcings (STD$Parms, "rSlow",     rSlowForc, spinup, fac = 1)
    forcings[[ 7]] <- Setforcings (STD$Parms, "pFast",     pFastForc, spinup, fac = 1)
    forcings[[ 8]] <- Setforcings (STD$Parms, "MPBprod", MPBprodForc, spinup, fac = 1)
    forcings[[ 9]] <- Setforcings (STD$Parms, "FePflux", FePfluxForc, spinup, fac = 1)
    forcings[[10]] <- Setforcings (STD$Parms, "CaPflux", CaPfluxForc, spinup, fac = 1)
    forcings[[11]] <- Setforcings (STD$Parms, "gasflux", gasfluxForc, spinup, fac = 1)
    forcings[[12]] <- Setforcings (STD$Parms, "O2bw",       O2bwForc, spinup, fac = 1)
    forcings[[13]] <- Setforcings (STD$Parms, "NO3bw",     NO3bwForc, spinup, fac = 1)
    forcings[[14]] <- Setforcings (STD$Parms, "NO2bw",     NO2bwForc, spinup, fac = 1)
    forcings[[15]] <- Setforcings (STD$Parms, "NH3bw",     NH3bwForc, spinup, fac = 1)
    forcings[[16]] <- Setforcings (STD$Parms, "ODUbw",     ODUbwForc, spinup, fac = 1)
    forcings[[17]] <- Setforcings (STD$Parms, "PO4bw",     PO4bwForc, spinup, fac = 1)
    forcings[[18]] <- Setforcings (STD$Parms, "DICbw",     DICbwForc, spinup, fac = 1)
    forcings[[19]] <- HWATERForc
    forcings[[20]] <- Setforcings (list(ratefactor = 1), "ratefactor", ratefactor, spinup, fac = 1)
    
    if (model == 1){
     bwO2  <- approxfun(x = forcings[[12]], rule = 2)
     bwNO3 <- approxfun(x = forcings[[13]], rule = 2)
     bwNO2 <- approxfun(x = forcings[[14]], rule = 2)
     bwNH3 <- approxfun(x = forcings[[15]], rule = 2)
     bwODU <- approxfun(x = forcings[[16]], rule = 2)
     bwPO4 <- approxfun(x = forcings[[17]], rule = 2)
     bwDIC <- approxfun(x = forcings[[18]], rule = 2)
    } else {
      BWapprox <- approxfun(x = range(times), y = c(0,0), rule = 2)
      bwO2 <- bwNO3 <- bwNO2 <- bwNH3 <- bwODU <- bwPO4 <- bwDIC <- BWapprox
    }
    
    if (model == 2) 
      bwHeight <- approxfun(x = forcings[[19]], rule = 2)

    ZZ <- c(ZZ, capture.output(suppressWarnings(   
      DYN <- ode.1D(y = yini, names = .CNPDIA$svar, initforc = "initforcp",
                   forcings = forcings,  times = spinup, method = "lsodes",
                   func = func, initfunc = "initomexdiap", lrw = lrw,
                   initpar = initpar, dllname = "CNPDIA", bandwidth = band, 
                   nout = .CNPDIA$nout, outnames = .CNPDIA$outnames, 
                   events = events, nspec = nspec, ...)
    )))
    yini <- DYN[nrow(DYN),2:(nspec*N+1)]
  }
  if (model == 1) 
    HWATERForc <- Setforcings (Hlist, "Hwater",     HwaterForc, times, fac = 1)   # not used
  else 
    HWATERForc <- Setforcings (STD$Parms, "Hwater", HwaterForc, times, fac = 1)
  
  forcings <- list()
  forcings[[ 1]] <- Setforcings (STD$Parms, "Cflux",     CfluxForc, times, fac = 1)
  forcings[[ 2]] <- Setforcings (STD$Parms, "w",             wForc, times, fac = 1)
  forcings[[ 3]] <- Setforcings (STD$other, "biot",       biotForc, times, fac = 1/STD$other[["biot"]])
  forcings[[ 4]] <- Setforcings (STD$other, "irr",         irrForc, times, fac = 1/STD$other[["irr"]])
  forcings[[ 5]] <- Setforcings (STD$Parms, "rFast",     rFastForc, times, fac = 1)
  forcings[[ 6]] <- Setforcings (STD$Parms, "rSlow",     rSlowForc, times, fac = 1)
  forcings[[ 7]] <- Setforcings (STD$Parms, "pFast",     pFastForc, times, fac = 1)
  forcings[[ 8]] <- Setforcings (STD$Parms, "MPBprod", MPBprodForc, times, fac = 1)
  forcings[[ 9]] <- Setforcings (STD$Parms, "FePflux", FePfluxForc, times, fac = 1)
  forcings[[10]] <- Setforcings (STD$Parms, "CaPflux", CaPfluxForc, times, fac = 1)
  forcings[[11]] <- Setforcings (STD$Parms, "gasflux", gasfluxForc, times, fac = 1)
  forcings[[12]] <- Setforcings (STD$Parms,  "O2bw",      O2bwForc, times, fac = 1)
  forcings[[13]] <- Setforcings (STD$Parms, "NO3bw",     NO3bwForc, times, fac = 1)
  forcings[[14]] <- Setforcings (STD$Parms, "NO2bw",     NO2bwForc, times, fac = 1)
  forcings[[15]] <- Setforcings (STD$Parms, "NH3bw",     NH3bwForc, times, fac = 1)
  forcings[[16]] <- Setforcings (STD$Parms, "ODUbw",     ODUbwForc, times, fac = 1)
  forcings[[17]] <- Setforcings (STD$Parms, "PO4bw",     PO4bwForc, times, fac = 1)
  forcings[[18]] <- Setforcings (STD$Parms, "DICbw",     DICbwForc, times, fac = 1)
  forcings[[19]] <- HWATERForc
  forcings[[20]] <- Setforcings (list(ratefactor = 1), "ratefactor", ratefactor, times, fac = 1)
  
  bwO2  <- approxfun(x = forcings[[12]], rule = 2)
  bwNO3 <- approxfun(x = forcings[[13]], rule = 2)
  bwNO2 <- approxfun(x = forcings[[14]], rule = 2)
  bwNH3 <- approxfun(x = forcings[[15]], rule = 2)
  bwODU <- approxfun(x = forcings[[16]], rule = 2)
  bwPO4 <- approxfun(x = forcings[[17]], rule = 2)
  bwDIC <- approxfun(x = forcings[[18]], rule = 2)

  if (model == 2) 
    bwHeight <- approxfun(x = forcings[[19]], rule = 2)
  
  PertFlux <- PertBW <- PertBW2 <- NULL 
  ZZ <- c(ZZ, capture.output(suppressWarnings(   
    DYN <- ode.1D(y = yini, names = .CNPDIA$svar, initforc = "initforcp",
                   forcings = forcings,  times = times, method = "lsodes",
                   func = func, initfunc = "initomexdiap", lrw = lrw,
                   initpar = initpar, dllname = "CNPDIA", bandwidth = band, 
                   nout = .CNPDIA$nout, outnames = .CNPDIA$outnames, 
                   events = events, nspec = nspec, ...)
  )))
  colnames(DYN)[2:(12*N+1)] <- as.vector(sapply(.CNPDIA$svar, FUN = function(x) rep(x, times = N)))
  
  if (model == 2) {
    cn <- colnames(DYN)
    cn[seq(2, by = N, length.out = 12)] <- paste(.CNPDIA$svar, "bw", sep ="")
    colnames(DYN) <- cn
  }
  if (verbose) print(ZZ)
  attr(DYN, "Depth") <- STD$Depth
  attr(DYN, "dx") <- STD$dx
  attr(DYN, "porosity")   <- STD$porosity
  attr(DYN, "bioturbation") <- STD$bioturbation
  attr(DYN, "irrigation")   <- STD$irrigation
  attr(DYN, "DistReact")    <- STD$DistReact
  attr(DYN, "includeMPB") <- FALSE 
  attr(DYN, "Parms") <- STD$Parms 
  attr(DYN, "perturbFluxes")  <- data.frame(time = perturbTimes[1:nrow(PertFlux)], PertFlux)
  if (model == 2) {
    attr(DYN, "BWbeforePert") <- data.frame(time = PertBW[1:nrow(PertBW)], PertBW)
    attr(DYN, "BWafterPert")  <- data.frame(time = PertBW2[1:nrow(PertBW2)], PertBW2)
  }
  attr(DYN, "perturbSettings") <- list(perturbType = perturbType, 
               perturbTimes =  perturbTimes, perturbDepth = depthPert, concfac = concfac)
  attr(DYN, "warnings")     <- ZZ
  attr(DYN, "model") <- paste("CNPDIA_model_",model,sep="")
  class(DYN) <- c("CNPDIAdyn", class(DYN))
  
  return(DYN)
}


CNPDIAperturbFluxes <- function(out, which = NULL) {
  
  if (missing(out)) 
    stop("object 'out' should be given")
  
  if (! "CNPDIAdyn" %in% (class(out))) 
     stop("perturbation fluxes can only be obtained from a run performed with 'CNPDIAdyn' or 'CNPDIAperturb'" )
    
  W <- attributes(out)$perturbFluxes
    
  if (! is.null(W)) { 
    if (! is.null(which))
      W <- cbind(W$time, W[,which])
  }
  return(W)
}  

CNPDIAperturbSettings <- function(out) {
  
  if (missing(out)) 
    stop("object 'out' should be given")
  
  if (! "CNPDIAdyn" %in% (class(out))) 
     stop("perturbation settings can only be obtained from a run performed with 'CNPDIAdyn' or 'CNPDIAperturb'" )
    
  W <- attributes(out)$perturbSettings
  return(W)
}  
