#####################################################################################################
######                MPBDIA: C, N, P, O2 diagenesis with MicroPhytoBenthos                    ######
#####################################################################################################


## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------
## solve the steady-state condition - main function
## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------

MPBDIAsolve <- function (parms = list(), yini = NULL, Grid = NULL, porosity = NULL, 
                         bioturbation = NULL, irrigation = NULL, diffusionfactor = NULL, 
                         dynamicbottomwater = FALSE, ratefactor = NULL, verbose = FALSE, 
                         times = c(0, 1e6), method = "runsteady", ...) { 
  model <- 1
  if (dynamicbottomwater) model <- 2
  if (is.null (yini)){
  # first solve the diagenetic part - without microphytobenthos
    ParmsCNP <-  parms[names(parms) %in% names(.CNPDIA$Parms)]
    YDIA <- CNPDIAsolve_full (ParmsCNP, gridtype = 1, Grid = Grid, 
               porosity = porosity, bioturbation = bioturbation, irrigation = irrigation, 
               surface = NULL, diffusionfactor = diffusionfactor,  
               model =  model, yini = yini, ratefactor = ratefactor, verbose = verbose, 
               includeMPB = FALSE, ...)$y
  
  } else {
    YDIA <- yini
  }
  
  if (length (YDIA) == 12*.CNPDIA$N) YDIA <- c(YDIA, .MPBDIA$yini)
  
  if (! is.null(method))
    if ( method %in% c("runsteady", "mixed")) {

    DIA <- CNPDIAdyna_Full (parms = parms, times = times, gridtype = 1, yini = YDIA, 
                   Grid = Grid, porosity = porosity, bioturbation = bioturbation, 
                   irrigation = irrigation, surface = NULL, 
                   diffusionfactor = diffusionfactor, 
                   model = model, ratefactor = ratefactor, verbose = verbose, includeMPB = TRUE, ...) 
     Att <- attributes(DIA)[-1]
     DIA <- DIA[2, -1]
     SVAR <- DIA[1:(.CNPDIA$N*17)]
     if (method == "runsteady"){

       MPBdll <- CNPDIAsolve_full (parms, gridtype = 1, Grid = Grid, 
                      porosity = porosity, bioturbation = bioturbation, 
                      irrigation = irrigation, surface = NULL, 
                      diffusionfactor = diffusionfactor,  
                      model =  model, yini = SVAR, ratefactor = ratefactor, 
                      verbose = verbose, times = NA, includeMPB = TRUE, ...)
       precis <- mean(abs(MPBdll$dy))
       ewt    <- abs(SVAR*1e-8) + 1e-8
       steady <- all(precis < ewt)
     
       SVNAMES  <- unique(names(SVAR))
       STD  <- list(y = matrix(nrow = .CNPDIA$N, ncol = 17, data  = SVAR))
       colnames(STD$y) <- SVNAMES
       VAR <- DIA[-(1:(.CNPDIA$N*17))]

       OUT <- unique(names(VAR))
       Z <- lapply(OUT, FUN = function(x) {
        VV <- subset(VAR, names(VAR) == x)
        names(VV) <- NULL
        VV })
       names(Z) <- OUT   
     
       Att <- Att[-1]
       nn <- c("Depth", "dx",  "Parms", "Grid", "porosity", "bioturbation", 
         "irrigation", "isDistReact", "DistReact", "includeMPB")
       STD <- c(STD, Z, Att[nn])  
       STD$model <- paste("MPBDIA_model", model, sep = "_")
       attributes(STD) <- c(attributes(STD)[1], Att[!names(Att) %in% nn])
       attr(STD, "precis") <- precis
       attr(STD, "steady") <- steady
       class(STD) <- c("MPBDIAstd", "steady1D",  "rootSolve", "list")
       return(STD)
     } else {
       YDIA <- SVAR
       
     }   
   }                   
    STD <- CNPDIAsolve_full (parms, gridtype = 1, Grid = Grid, porosity = porosity, bioturbation = bioturbation, 
                      irrigation = irrigation, surface = NULL, diffusionfactor = diffusionfactor,  
                      model =  model, yini = YDIA,  ratefactor = ratefactor, verbose = verbose, 
                      includeMPB = TRUE, ...)
  STD                    
}


## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------
## solve the dynamic condition - main function
## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------

MPBDIAdyna <- function (parms = list(), times = 0:365, spinup = NULL, yini = NULL, 
   Grid = NULL, porosity = NULL, bioturbation = NULL, irrigation = NULL, 
   diffusionfactor = NULL, dynamicbottomwater = FALSE, 
   CfluxForc   = NULL, FePfluxForc = NULL, CaPfluxForc = NULL, O2bwForc  = NULL, 
   NO3bwForc   = NULL, NO2bwForc   = NULL, NH3bwForc   = NULL, 
   ODUbwForc   = NULL, PO4bwForc   = NULL, DICbwForc   = NULL, wForc     = NULL,  # cm/d  - advection rate
   biotForc    = NULL, irrForc     = NULL, rFastForc   = NULL, rSlowForc  = NULL, pFastForc   = NULL, 
   lightForc = NULL, gasfluxForc = NULL, HwaterForc  = NULL, ratefactor = NULL, 
   verbose = FALSE, 
    ...)  {
  
  model <- 1
  if (dynamicbottomwater) 
    model <- 2

  if (is.null(yini)) 
    yini <- MPBDIAsolve (parms = parms, yini = NULL, Grid = Grid, porosity = porosity, 
           bioturbation = bioturbation, irrigation = irrigation, 
           diffusionfactor = diffusionfactor, dynamicbottomwater = dynamicbottomwater, 
           ratefactor = ratefactor, verbose = verbose, ...)$y
  
  dyn <- CNPDIAdyna_Full (parms = parms, times = times, spinup = spinup, gridtype = 1, yini = yini, 
                   CfluxForc = CfluxForc, FePfluxForc = FePfluxForc, CaPfluxForc = CaPfluxForc, 
                   O2bwForc = O2bwForc,  NO3bwForc = NO3bwForc, NO2bwForc = NO2bwForc, NH3bwForc = NH3bwForc, 
                   ODUbwForc   = ODUbwForc, PO4bwForc   = PO4bwForc, DICbwForc   = DICbwForc, 
                   wForc       = wForc, biotForc    = biotForc, irrForc     = irrForc,
                   rFastForc   = rFastForc, rSlowForc   = rSlowForc,  pFastForc   = pFastForc, 
                   MPBprodForc = lightForc, HwaterForc = HwaterForc, gasfluxForc = gasfluxForc, 
                   Grid = Grid, porosity = porosity, bioturbation = bioturbation, 
                   irrigation = irrigation, surface = NULL, 
                   diffusionfactor = diffusionfactor, 
                   model = model, ratefactor = ratefactor, verbose = verbose, includeMPB = TRUE, ...) 
  class(dyn) <- c("MPBDIAdyn", class(dyn))
  dyn

  }
