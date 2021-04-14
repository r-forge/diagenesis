#####################################################################################################
######                           CNPDIA: C, N, P, O2 diagenesis                                ######
#####################################################################################################

## --------------------------------------------------------------------------------------------------
## Initialisation: grid generations (layers, bioturbation, irrigation), parameters
## --------------------------------------------------------------------------------------------------

initCNPDIA <- function (parms = list(), gridtype = 1, CfluxForc = NULL, 
                        FePfluxForc = NULL, CaPfluxForc = NULL, O2bwForc  = NULL, 
                        NO3bwForc = NULL, NO2bwForc = NULL, NH3bwForc = NULL,
                        ODUbwForc = NULL, PO4bwForc = NULL, DICbwForc = NULL, 
                        gasfluxForc = NULL,wForc = NULL, biotForc = NULL, irrForc = NULL, 
                        ratefactor = NULL, rFastForc = NULL, rSlowForc = NULL, 
                        pFastForc = NULL, MPBprodForc = NULL, HwaterForc = NULL,
                        Grid = NULL, porosity = NULL, bioturbation = NULL, 
                        irrigation = NULL, surface = NULL, diffusionfactor = NULL,
                        times = NULL, model = 1, includeMPB = FALSE)  {

  if (is.null(Grid))
    Grid  <- setup.grid.1D(x.up = 0, dx.1 = 0.01, N = .CNPDIA$N, L = 100)
  else {
    if (is.list (Grid)) {
      nms <- names(Grid)
      if (! "dx" %in%  nms)
          stop ("'Grid' should be a list containing 'dx' and 'dx.aux'")
      if (! "dx.aux" %in%  nms)
          stop ("'Grid' should be a list containing 'dx' and 'dx.aux'")
      if (length(Grid$dx.aux) != .CNPDIA$N +1)
            stop ("Checking 'Grid': 'dx.aux' should be a vector of length ", .CNPDIA$N+1)
      if (length(Grid$dx) != .CNPDIA$N)
          stop ("Checking 'Grid': 'dx' should be a vector of length ", .CNPDIA$N)
    } else {
      if (length(Grid) != .CNPDIA$N +1)
        stop ("Checking 'Grid': should be a vector of length ", .CNPDIA$N+1)
      Grid <- list(dx.aux = Grid, mid = 0.5*(Grid[-1] + Grid[-(.CNPDIA$N+1)]))
    }
  }   
      
## check parameter inputs
  Parms <- c(.CNPDIA$Parms, .MPBDIA$Parms)

  nms <- names(Parms)
  Parms[(namc <- names(parms))] <- parms
  if (length(noNms <- namc[!namc %in% nms]) > 0)
    warning("unknown names in parms: ", paste(noNms, collapse = ", "))
  PP <- unlist(Parms)

  # parameters need to be set based on forcing functions (if present)
  if (includeMPB) MPBprodPar <- "light" else MPBprodPar <- "MPBprod"
  
  Parms <- CreateMeanPars (Parms, "Cflux"   , CfluxForc   , times) 
  Parms <- CreateMeanPars (Parms, "FePflux" , FePfluxForc , times) 
  Parms <- CreateMeanPars (Parms, "CaPflux" , CaPfluxForc , times) 
  Parms <- CreateMeanPars (Parms, "w"      ,      wForc   , times) 
  Parms <- CreateMeanPars (Parms, "biot"    ,  biotForc   , times) 
  Parms <- CreateMeanPars (Parms, "irr"     ,   irrForc   , times) 
  Parms <- CreateMeanPars (Parms, "rFast"   , rFastForc   , times) 
  Parms <- CreateMeanPars (Parms, "rSlow"   , rSlowForc   , times) 
  Parms <- CreateMeanPars (Parms, "pFast"   , pFastForc   , times) 
  Parms <- CreateMeanPars (Parms, MPBprodPar , MPBprodForc , times) 
  Parms <- CreateMeanPars (Parms, "gasflux" , gasfluxForc , times) 
  Parms <- CreateMeanPars (Parms, "O2bw"    , O2bwForc    , times) 
  Parms <- CreateMeanPars (Parms, "NO3bw"   , NO3bwForc   , times) 
  Parms <- CreateMeanPars (Parms, "NO2bw"   , NO2bwForc   , times) 
  Parms <- CreateMeanPars (Parms, "NH3bw"   , NH3bwForc   , times) 
  Parms <- CreateMeanPars (Parms, "ODUbw"   , ODUbwForc   , times) 
  Parms <- CreateMeanPars (Parms, "PO4bw"   , PO4bwForc   , times) 
  Parms <- CreateMeanPars (Parms, "DICbw"   , DICbwForc   , times) 
  Parms <- CreateMeanPars (Parms, "Hwater"  , HwaterForc , times) 

  if (is.na(Parms[["pdepo"]]))
    Parms[["pdepo"]] <- max(0,min(1.,  0.233*(Parms[["w"]]*365)^0.336 ))

  if (Parms[["BCupLiq"]] != 3) {
    if (any(is.na(Parms[c("O2bw","NO3bw","NO2bw","NH3bw","ODUbw","PO4bw","DICbw")])))
      stop("bottom water concentrations cannot be NA if type flux or concentration")
  }
  
  if (Parms[["BCdownLiq"]] != 3) {
    if (any(is.na(Parms[[c("O2dw","NO3dw","NO2dw","NH3dw","ODUdw","PO4dw","DICdw")]])))
      stop("deep water concentrations cannot be NA if type flux or concentration")
  }

  if (! is.null(surface)) {
    Aint <- CheckProfile (surface, "surface", interface = TRUE)$int
  }  else {
    if (gridtype == 1)                        # cartesian
      Aint <- rep(1, .CNPDIA$N+1)
    else if (gridtype == 2)                   # cylindrical
      Aint <- rev(2*pi*Grid$x.int)
    else if (gridtype == 3)                   # spherical
     Aint <- rev(pi*(Grid$x.int)^2)
  }
  DF <- diffcoeff(S = Parms[["salinity"]], t = Parms[["temperature"]])
  Parms <- c(Parms, DF[c("O2","NO3","NO2","NH4","HS","H2PO4","HCO3")]*86400e4)

# porosity gradients
  exp.profile <- function(x, y.0, y.inf, x.att = 1, L = 0)
           return(y.inf + (y.0 - y.inf) * exp(-pmax(0, x-L)/x.att))

# Bioturbation profile
  if (is.null(bioturbation)) 
    bioturbation <- setup.prop.1D(func = exp.profile,
                           grid = Grid,
                           y.0 = Parms[["biot"]], y.inf = 0.,
                           L = Parms[["biotdepth"]], x.att = Parms[["biotatt"]])
  else{
    bioturbation <- CheckProfile (bioturbation, "bioturbation", interface = TRUE)
#    Parms[["biot"]] <- 1     ##### KS CHECK
  }
  Db <- bioturbation$int

# Irrigation profile  
  if (is.null(irrigation))
    irrigation <- setup.prop.1D(func = exp.profile,
                         grid = Grid,
                         y.0 = Parms[["irr"]], y.inf = 0.,
                         L = Parms[["irrdepth"]], x.att = Parms[["irratt"]])
  else{
    irrigation <- CheckProfile(irrigation,"irrigation", interface = FALSE)
#    Parms[["irr"]] <- 1
  }
  Dirr <- irrigation$mid

# porosity profile
  if (is.null(porosity))
    porGrid <- setup.prop.1D(func = exp.profile,
                         grid = Grid,
                         y.0 = Parms[["por0"]], y.inf = Parms[["pordeep"]], L = 0,
                         x.att = Parms[["porcoeff"]])

  else 
    porGrid <- CheckProfile (porosity, "porosity", interface = TRUE)

# Long-distance reactions - with oxygen in surface layer
  if (Parms[["rSurfODUox"]] > 0 )
    distreact <- setup.prop.1D(func = exp.profile,
                               grid = Grid,
                               y.0 = 1, y.inf = 0.,
                               L = Parms[["ODUoxdepth"]], x.att = Parms[["ODUoxatt"]])$mid 
  else 
    distreact <- rep(0, .CNPDIA$N)
  
# factor to multiply with the diffusion to estimate effective diffusion  
  if (is.null(diffusionfactor)){
    formationtype <- Parms[["formationtype"]]
    if (formationtype == 1) #sand
      diffusionfactor <- porGrid$int
    else if (formationtype == 2) #mud/fine sand
      diffusionfactor <- porGrid$int^2
    else   #general
      diffusionfactor <- 1/(1-log(porGrid$int^2))
  }else
    diffusionfactor <- CheckProfile (diffusionfactor, "diffusionfactor", interface = TRUE)$int
  
  toremove <- c("biot", "biotdepth", "biotatt", "irr", "irrdepth", "irratt", 
                 "ODUoxatt", "por0", "pordeep", "ODUoxdepth", "formationtype",
                "porcoeff",  "temperature", "salinity", 
                "Cflux", "FePflux", "CaPflux", "w", "pFast", "rFast", "rSlow", "pFast", 
                "MPBprod","gasflux","O2bw","NO3bw","NO2bw","NH3bw","ODUbw",
                "PO4bw","DICbw","Hwater", "light")
  parms <- unlist(Parms[-which(nms %in% toremove)])


# parameters to pass to model DLL
  initpar <- c(parms, Grid$dx, Grid$dx.aux, Grid$x.mid, Aint, porGrid$mid,  
               porGrid$int, diffusionfactor, Db, Dirr, distreact)

  list(initpar = initpar, other = Parms[c("biot", "irr")], Parms = PP, 
       Grid = Grid, porGrid = porGrid, 
       Isirr = sum(Dirr) > 0, bioturbationGrid = bioturbation,  
       Depth = Grid$x.mid, 
       dx = Grid$dx, porosity = porGrid$mid, 
       bioturbation = bioturbation$mid, irrigation = Dirr, 
       DistReact = distreact, isDistReact = (max(distreact) > 0))
}

## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------
## solve the steady-state condition - main function
## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------

CNPDIAsolve <- function (parms = list(), yini = NULL, gridtype = 1, Grid = NULL, porosity = NULL, 
                         bioturbation = NULL, irrigation = NULL, surface = NULL, diffusionfactor = NULL, 
                         dynamicbottomwater = FALSE, ratefactor = NULL, verbose = FALSE, ...) { 
  model <- 1
  if (dynamicbottomwater) model <- 2
  std <- CNPDIAsolve_full (parms = parms, yini = yini, gridtype = gridtype, Grid = Grid, 
                    porosity = porosity, bioturbation = bioturbation, irrigation = irrigation, 
                    surface = surface, diffusionfactor = diffusionfactor,  
                    model =  model, ratefactor = ratefactor, verbose = verbose, 
                    includeMPB = FALSE, ...)
  class(std) <- c("CNPDIAstd", class(std))
  std

}

## --------------------------------------------------------------------------------------------------

CNPDIAsolve_full <- function (parms = list(), gridtype = 1, CfluxForc = NULL, 
                         FePfluxForc = NULL, CaPfluxForc = NULL,
                         O2bwForc  = NULL, NO3bwForc = NULL, NO2bwForc = NULL, NH3bwForc = NULL,
                         ODUbwForc = NULL, PO4bwForc = NULL, DICbwForc = NULL, gasfluxForc = NULL,
                         wForc = NULL, biotForc= NULL, irrForc = NULL, 
                         rFastForc = NULL, rSlowForc = NULL, pFastForc = NULL, 
                         MPBprodForc = NULL, HwaterForc = NULL,
                         times = NULL, Grid = NULL, porosity = NULL, 
                         bioturbation = NULL, irrigation = NULL, surface = NULL, 
                         ratefactor = NULL, diffusionfactor = NULL, 
                         yini = NULL, model = 1, verbose = FALSE, includeMPB = FALSE, 
                         method = NULL, ...)  {
  
  P    <- initCNPDIA(parms = parms, gridtype = gridtype, CfluxForc = CfluxForc, 
                     FePfluxForc = FePfluxForc, CaPfluxForc = CaPfluxForc, 
                     O2bwForc  = O2bwForc, NO3bwForc = NO3bwForc, NO2bwForc = NO2bwForc, 
                     NH3bwForc = NH3bwForc, ODUbwForc = ODUbwForc, 
                     PO4bwForc = PO4bwForc, DICbwForc = DICbwForc, gasfluxForc = gasfluxForc, wForc = wForc, 
                     biotForc = biotForc, irrForc = irrForc,  ratefactor = ratefactor, 
                     rFastForc = rFastForc, rSlowForc = rSlowForc, pFastForc = pFastForc, 
                     MPBprodForc = MPBprodForc, HwaterForc = HwaterForc,
                     times = times, Grid = Grid, porosity = porosity, bioturbation = bioturbation, 
                     irrigation = irrigation, surface = surface, diffusionfactor = diffusionfactor, 
                     model = model, includeMPB = includeMPB)
  initfunc <- "initomexdiap"
  nspec    <- 12
  svar     <- .CNPDIA$svar
  nout     <- .CNPDIA$nout
  outnames <- .CNPDIA$outnames

  if (includeMPB) {
    
    if (model == 2) {     # DynamicBottomwater +MPB
      HWATERForc <- as.double(P$Parms["Hwater"])
      func <- "mpbdiamodbw"
      N <- .CNPDIA$N + 1
    } else {
      HWATERForc <- as.double(0)
      func <- "mpbdiamodp"
      N <- .CNPDIA$N
    }
    initfunc <- "initmpbdiap"
    nspec <- 17
    svar <- c(svar , .MPBDIA$svar)
    outnames <- c(outnames, .MPBDIA$outnames)
    outnames[outnames == "MPBprod"] <- "LightSurf"
    nout <- nout+.MPBDIA$nout

  } else if (model == 2) {  # DynamicBottomwater
    HWATERForc <- as.double(P$Parms["Hwater"])
    func <- "omexdiamodbw"
    N <- .CNPDIA$N + 1
  } else {
    HWATERForc <- as.double(0)
    func="omexdiamodp"
    N <- .CNPDIA$N
  }  
  if (includeMPB) MPBprodPar <- "light" else MPBprodPar <- "MPBprod"

  forcings <- c(as.double(P$Parms["Cflux"]),   as.double(P$Parms["w"]),
                as.double(1.),                 as.double(1.),     # for Db and irr: relative to parameter-profile
                as.double(P$Parms["rFast"]),   as.double(P$Parms["rSlow"]),   
                as.double(P$Parms["pFast"]),   as.double(P$Parms[MPBprodPar]), 
                as.double(P$Parms["FePflux"]), as.double(P$Parms["CaPflux"]),
                as.double(P$Parms["gasflux"]), as.double(P$Parms["O2bw"]),
                as.double(P$Parms["NO3bw"]),   as.double(P$Parms["NO2bw"]),   
                as.double(P$Parms["NH3bw"]), 
                as.double(P$Parms["ODUbw"]),   as.double(P$Parms["PO4bw"]), 
                as.double(P$Parms["DICbw"]),   HWATERForc, as.double(1.0))   

  Random <- is.null(yini)
  Yini <- yini  
  if (is.null (Yini))
   Yini   <- rep(10, nspec*N)
  else if(length(Yini) != nspec*N)
    stop ("'yini' not of correct length, should be ", nspec*N)

  # solution method
  sparse   <- "1D"
  jactype <- NULL
  if (is.null(method)){
#    if (includeMPB) {
#      method <- "runsteady" 
#      sparse <- "sparse"
#    } else 
      if (P$isDistReact)
      method <- "stodes"
    else if (model == 2 & P$Isirr) 
      method <- "stodes"
   else method <- "stode" 
  }

  if (method == "stodes") sparse <- "sparseint"
  
  inz <- NULL
  
  if (includeMPB) {
    if (is.null (yini))
      stop("'yini' should be properly initialised if 'includeMPB' = TRUE")

  } else if (P$isDistReact){
    sparse <- "sparseint"
    
  }  else if (model == 2 & P$Isirr) { 
    sparse <- "sparseint"
  }  
  
  ZZ <- NULL
  if (any(is.na(times))){

# expand forcings to be a "time series"
  lf <- as.list(forcings)
  forcings <- lapply(lf, FUN = function(x) cbind(t = c(0, Inf), v = x))

  ZZ <- c(ZZ, capture.output(suppressWarnings(   
     DIA  <- DLLfunc(y = as.double(Yini), func = func, initfunc = initfunc,
                    initforc = "initforcp", 
                    forcings = forcings, parms = P$initpar,
                    dllname = "CNPDIA", times = 0,
                    nout = nout, outnames = outnames, ...)
  )))
   attr(DIA,"message") <- ZZ 
   return(DIA)
  } 
  ZZ <- c(ZZ, capture.output(suppressWarnings(   
    DIA  <- steady.1D(y=as.double(Yini), func=func, initfunc=initfunc,
                    names = svar, initforc = "initforcp", 
                    jactype = sparse, method = method, inz = inz,
                    forcings = forcings, initpar=P$initpar, nspec=nspec,
                    dllname="CNPDIA", maxiter = 100,
                    nout = nout, outnames = outnames,
                    positive=TRUE, times = times, ...)
  )))
  niter <- 1
  while (! attributes(DIA)$steady & niter <= 50 & method != "runsteady")  {
    if (Random) 
      Yini <- runif(N*nspec)

    else
      Yini <- runif(N*nspec)*yini
    
    #   Yini <- runif(14 * .FESDIA$N)*niter 
    ZZ <- c(ZZ, capture.output(suppressWarnings(   
      DIA  <- steady.1D(y=as.double(Yini),func=func,initfunc=initfunc,
                        names = svar, initforc = "initforcp",
                        forcings = forcings, initpar=P$initpar,nspec=nspec,
                        dllname="CNPDIA", maxiter = 100,
                        jactype = sparse, method = method, inz = inz,
                        nout=nout, outnames = outnames,
                        positive=TRUE, ...)
    )))
    niter <- niter + 1
  }
  if (verbose) {
    if (!attributes(DIA)$steady) warning("steady-state not reached")
    print (ZZ)
  }
  
  if (method == "runsteady"){   # all output vars are in one long vector
     DD <- DIA
     names(DIA$var) <- outnames
     OUT <- unique(outnames)
     Z <- lapply(OUT, FUN = function(x) {
        VV <- subset(DIA$var, names(DIA$var) == x)
        names(VV) <- NULL
        VV })
     names(Z) <- OUT   
     DD$var <- NULL
     Att <- attributes(DIA)[-1]
     DIA <- c(DD, Z)   
     attributes(DIA) <- c(attributes(DIA)[1], Att[-1])
  }
  
  DIA$Depth        <- P$Depth
  DIA$dx           <- P$Grid$dx
  DIA$initpar      <- unlist(P$initpar)
  DIA$other        <- P$other
  DIA$Parms        <- P$Parms
  DIA$Grid         <- P$Grid
  DIA$porGrid      <- P$porGrid
  DIA$porosity     <- P$porosity
  DIA$bioturbation <- P$bioturbation
  DIA$irrigation   <- P$irrigation
  DIA$numberTries  <- niter
  DIA$warnings     <- ZZ
  DIA$model        <- paste("CNPDIA_model_",model,sep="")
  DIA$isDistReact  <- P$isDistReact
  DIA$DistReact    <- P$DistReact
  DIA$includeMPB   <- includeMPB
  return(DIA)
}


## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------
## solve the dynamic condition - main function
## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------


CNPDIAdyna <- function (parms = list(), times = 0:365, spinup = NULL, yini = NULL, 
   gridtype = 1, Grid = NULL, porosity = NULL, bioturbation = NULL, irrigation = NULL, 
   surface = NULL, diffusionfactor = NULL, dynamicbottomwater = FALSE, 
   CfluxForc   = NULL, FePfluxForc = NULL, CaPfluxForc = NULL, O2bwForc  = NULL, 
   NO3bwForc   = NULL, NO2bwForc   = NULL, NH3bwForc   = NULL, 
   ODUbwForc   = NULL, PO4bwForc   = NULL, DICbwForc   = NULL, wForc     = NULL,  # cm/d  - advection rate
   biotForc    = NULL, irrForc     = NULL, rFastForc   = NULL, rSlowForc  = NULL, pFastForc   = NULL, 
   MPBprodForc = NULL, gasfluxForc = NULL, HwaterForc  = NULL, ratefactor = NULL, verbose = FALSE, 
    ...)  {
  
  model <- 1
  if (dynamicbottomwater) 
    model <- 2
  
  dyn <- CNPDIAdyna_Full (parms = parms, times = times, spinup = spinup, gridtype = gridtype, yini = yini, 
                   CfluxForc = CfluxForc, FePfluxForc = FePfluxForc, CaPfluxForc = CaPfluxForc, 
                   O2bwForc = O2bwForc,  NO3bwForc = NO3bwForc, NO2bwForc = NO2bwForc, NH3bwForc = NH3bwForc, 
                   ODUbwForc   = ODUbwForc, PO4bwForc   = PO4bwForc, DICbwForc   = DICbwForc, 
                   wForc       = wForc, biotForc    = biotForc, irrForc     = irrForc,
                   rFastForc   = rFastForc, rSlowForc   = rSlowForc,  pFastForc   = pFastForc, 
                   MPBprodForc = MPBprodForc, HwaterForc = HwaterForc, gasfluxForc = gasfluxForc, 
                   Grid = Grid, porosity = porosity, bioturbation = bioturbation, 
                   irrigation = irrigation, surface = surface, 
                   diffusionfactor = diffusionfactor, 
                   model = model, ratefactor = ratefactor, verbose = verbose, includeMPB = FALSE, ...) 
  class(dyn) <- c("CNPDIAdyn", class(dyn))
  dyn
  }

# general dynamic function

CNPDIAdyna_Full <- function (parms = list(), times = 0:365, spinup = NULL, gridtype = 1, yini = NULL, 
                        CfluxForc   = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        FePfluxForc = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        CaPfluxForc = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        O2bwForc    = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        NO3bwForc   = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        NO2bwForc   = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        NH3bwForc   = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        ODUbwForc   = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        PO4bwForc   = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        DICbwForc   = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        wForc       = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0),  # cm/d  - advection rate
                        biotForc    = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0),  # cm2/d - bioturbation coefficient
                        irrForc     = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0),  # /d    - bio-irrigation rate
                        rFastForc   = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        rSlowForc   = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0),  
                        pFastForc   = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        MPBprodForc = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        gasfluxForc = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        HwaterForc  = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        ratefactor = list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0), 
                        
                        Grid = NULL, porosity = NULL, bioturbation = NULL, irrigation = NULL, 
                        surface = NULL, diffusionfactor = NULL,
                        model = 1, verbose = FALSE, includeMPB = FALSE, ...) 
{

  CfluxForc    <- checkforcs(CfluxForc,     "CfluxForc")
  FePfluxForc  <- checkforcs(FePfluxForc, "FePfluxForc")
  CaPfluxForc  <- checkforcs(CaPfluxForc, "CaPfluxForc")
  wForc        <- checkforcs(wForc,             "wForc")
  biotForc     <- checkforcs(biotForc,       "biotForc")
  irrForc      <- checkforcs(irrForc,         "irrForc")
  rFastForc    <- checkforcs(rFastForc,     "rFastForc")
  rSlowForc    <- checkforcs(rSlowForc,     "rSlowForc")
  pFastForc    <- checkforcs(pFastForc,     "pFastForc")
  MPBprodForc  <- checkforcs(MPBprodForc, "MPBprodForc")
  gasfluxForc  <- checkforcs(gasfluxForc, "gasfluxForc")
  O2bwForc     <- checkforcs(O2bwForc ,      "O2bwForc")
  NO3bwForc    <- checkforcs(NO3bwForc,     "NO3bwForc")
  NO2bwForc    <- checkforcs(NO2bwForc,     "NO2bwForc")
  NH3bwForc    <- checkforcs(NH3bwForc,     "NH3bwForc")
  ODUbwForc    <- checkforcs(ODUbwForc,     "ODUbwForc")
  PO4bwForc    <- checkforcs(PO4bwForc,     "PO4bwForc")
  DICbwForc    <- checkforcs(DICbwForc,     "DICbwForc")
  HwaterForc   <- checkforcs(HwaterForc,   "HwaterForc")
  ratefactor   <- checkforcs(ratefactor,   "ratefactor")
  
  if (is.null(spinup)) Times <- times else Times <- spinup
  nspec <- 12
  func <- "omexdiamodp"
  initfunc <- "initomexdiap"
  svar <- .CNPDIA$svar
  nout <- .CNPDIA$nout
  outnames <- .CNPDIA$outnames
  
  if (includeMPB) {
    func <- "mpbdiamodp"
    initfunc <- "initmpbdiap"
    nspec <- 17
    svar <- c(svar , .MPBDIA$svar)
    outnames <- c(outnames, .MPBDIA$outnames)
    outnames[outnames == "MPBprod"] <- "LightSurf"
    nout <- nout+.MPBDIA$nout
  }  
  
  if (is.null(yini)) {
    STD <- CNPDIAsolve_full(parms, gridtype = gridtype, CfluxForc = CfluxForc, 
                       FePfluxForc = FePfluxForc, CaPfluxForc = CaPfluxForc, 
                       O2bwForc  = O2bwForc, NO3bwForc = NO3bwForc, NO2bwForc = NO2bwForc, NH3bwForc = NH3bwForc,
                       ODUbwForc = ODUbwForc, PO4bwForc = PO4bwForc, DICbwForc = DICbwForc, 
                       gasfluxForc = gasfluxForc, wForc = wForc, biotForc = biotForc, 
                       irrForc = irrForc, ratefactor = ratefactor, HwaterForc = HwaterForc,
                       rFastForc = rFastForc, rSlowForc = rSlowForc, pFastForc = pFastForc,
                       MPBprodForc = MPBprodForc, 
                       times = Times, Grid = Grid, porosity = porosity,
                       bioturbation = bioturbation, irrigation = irrigation,
                       surface = surface, 
                       diffusionfactor = diffusionfactor, model = model, 
                       verbose = verbose, includeMPB = includeMPB)
    yini <- STD$y
  } else  
    STD <- initCNPDIA(parms, gridtype = gridtype, CfluxForc = CfluxForc, 
                      FePfluxForc = FePfluxForc, CaPfluxForc = CaPfluxForc, 
                      O2bwForc  = O2bwForc, NO3bwForc = NO3bwForc, NO2bwForc = NO2bwForc, NH3bwForc = NH3bwForc,
                      ODUbwForc = ODUbwForc, PO4bwForc = PO4bwForc, DICbwForc = DICbwForc, 
                      gasfluxForc = gasfluxForc, wForc = wForc, biotForc = biotForc, 
                      irrForc = irrForc, ratefactor = ratefactor, HwaterForc = HwaterForc,
                      rFastForc = rFastForc, rSlowForc = rSlowForc, pFastForc = pFastForc, 
                      MPBprodForc = MPBprodForc,
                      times = Times, Grid = Grid, porosity = porosity, 
                      bioturbation = bioturbation, irrigation = irrigation, surface = surface, 
                      diffusionfactor = diffusionfactor, model = model, includeMPB = includeMPB)

  initpar = unlist(STD$initpar)
  if (STD$isDistReact | model == 2) {
     band   <- 0
     lrw <- 190000
  } else if (! includeMPB){
    band   <- 1
    lrw <- 90000
  } else {
    band   <- 1
    lrw <- 170000
  
  }
  
  if (model == 2) {
    if (! includeMPB) func <- "omexdiamodbw" else func <- "mpbdiamodbw"
    N <- .CNPDIA$N+1  
  } else {
    N <- .CNPDIA$N
  } 
  
  if ("deSolve" %in% class(yini))   # last values of a run = initial conditions
    yini <- yini[nrow(yini), 2:(1+nspec*N)]
  
  if(length(yini) != nspec*N)
    stop ("'yini' not of correct length, should be ", nspec*N)
  
  
  is.spinup <- ! is.null(spinup)
  if (is.spinup)
    is.spinup <- (length(spinup) > 1)
  ZZ <- NULL
  if (includeMPB) MPBprodPar <- "light" else MPBprodPar <- "MPBprod"
  if (is.spinup) {

    forcings <- list()
    forcings[[ 1]] <- Setforcings (STD$Parms, "Cflux",     CfluxForc, spinup, fac = 1)
    forcings[[ 2]] <- Setforcings (STD$Parms, "w",             wForc, spinup, fac = 1)
    forcings[[ 3]] <- Setforcings (STD$other, "biot",       biotForc, spinup, fac = 1/STD$other[["biot"]])
    forcings[[ 4]] <- Setforcings (STD$other, "irr",         irrForc, spinup, fac = 1/STD$other[["irr"]])
    forcings[[ 5]] <- Setforcings (STD$Parms, "rFast",     rFastForc, spinup, fac = 1)
    forcings[[ 6]] <- Setforcings (STD$Parms, "rSlow",     rSlowForc, spinup, fac = 1)
    forcings[[ 7]] <- Setforcings (STD$Parms, "pFast",     pFastForc, spinup, fac = 1)
    forcings[[ 8]] <- Setforcings (STD$Parms, MPBprodPar, MPBprodForc, spinup, fac = 1)
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
    forcings[[19]] <- Setforcings (STD$Parms,"Hwater",    HwaterForc, spinup, fac = 1)
    forcings[[20]] <- Setforcings (list(ratefactor = 1), "ratefactor", ratefactor, spinup, fac = 1)
    
  # Note: forcings for bioturbation and irrigation are relative to the parameter   
   ZZ <- c(ZZ, capture.output(suppressWarnings(   
      DYN <- ode.1D(y = yini, names = svar, initforc = "initforcp",
                    forcings = forcings,  times = spinup, method = "lsodes",
                    func = func, initfunc = initfunc, lrw = lrw,
                    initpar = initpar, dllname = "CNPDIA", bandwidth = band, 
                    nout = nout, outnames = outnames,
                    nspec = nspec, ...)     
                     )))
     yini <- DYN[nrow(DYN),2:(nspec*N+1)]
  }
  
  forcings <- list()
  forcings[[ 1]] <- Setforcings (STD$Parms, "Cflux",     CfluxForc, times, fac = 1)
  forcings[[ 2]] <- Setforcings (STD$Parms, "w",             wForc, times, fac = 1)
  forcings[[ 3]] <- Setforcings (STD$other, "biot",       biotForc, times, fac = 1/STD$other[["biot"]])
  forcings[[ 4]] <- Setforcings (STD$other, "irr",         irrForc, times, fac = 1/STD$other[["irr"]])
  forcings[[ 5]] <- Setforcings (STD$Parms, "rFast",     rFastForc, times, fac = 1)
  forcings[[ 6]] <- Setforcings (STD$Parms, "rSlow",     rSlowForc, times, fac = 1)
  forcings[[ 7]] <- Setforcings (STD$Parms, "pFast",     pFastForc, times, fac = 1)
  forcings[[ 8]] <- Setforcings (STD$Parms, MPBprodPar, MPBprodForc, times, fac = 1)
  forcings[[ 9]] <- Setforcings (STD$Parms, "FePflux", FePfluxForc, times, fac = 1)
  forcings[[10]] <- Setforcings (STD$Parms, "CaPflux", CaPfluxForc, times, fac = 1)
  forcings[[11]] <- Setforcings (STD$Parms, "gasflux", gasfluxForc, times, fac = 1)
  forcings[[12]] <- Setforcings (STD$Parms, "O2bw",       O2bwForc, times, fac = 1)
  forcings[[13]] <- Setforcings (STD$Parms, "NO3bw",     NO3bwForc, times, fac = 1)
  forcings[[14]] <- Setforcings (STD$Parms, "NO2bw",     NO2bwForc, times, fac = 1)
  forcings[[15]] <- Setforcings (STD$Parms, "NH3bw",     NH3bwForc, times, fac = 1)
  forcings[[16]] <- Setforcings (STD$Parms, "ODUbw",     ODUbwForc, times, fac = 1)
  forcings[[17]] <- Setforcings (STD$Parms, "PO4bw",     PO4bwForc, times, fac = 1)
  forcings[[18]] <- Setforcings (STD$Parms, "DICbw",     DICbwForc, times, fac = 1)
  forcings[[19]] <- Setforcings (STD$Parms, "Hwater",   HwaterForc, times, fac = 1)
  forcings[[20]] <- Setforcings (list(ratefactor = 1), "ratefactor", ratefactor, times, fac = 1)
  
  if (model == 3)
    return(DLLfunc(y = yini, initforc = "initforcp", forcings = forcings,  
                   func = func, initfunc = initfunc, parms = initpar, 
                   times = 0, dllname = "CNPDIA", 
                   nout = nout, outnames = outnames))

  ZZ <- c(ZZ, capture.output(suppressWarnings(   
    DYN <- ode.1D(y = yini, names = svar, initforc = "initforcp",
                forcings = forcings,  times = times, method = "lsodes",
                func = func, initfunc = initfunc, lrw = lrw,
                initpar = unlist(STD$initpar), dllname = "CNPDIA", bandwidth = band, 
                nout = nout, outnames = outnames,
                nspec = nspec, ...)
  )))
  colnames(DYN)[2:(nspec*N+1)] <- as.vector(sapply(svar, FUN = function(x) rep(x, times = N)))
  if (model == 2) {
    cn <- colnames(DYN)
    cn[seq(2, by = N, length.out = nspec)] <- paste(svar, "bw", sep ="")
    colnames(DYN) <- cn
  }
  if (verbose) print(ZZ)
  attr(DYN, "Depth") <- STD$Grid$x.mid
  attr(DYN, "dx") <- STD$dx
  attr(DYN, "porosity")   <- STD$porosity
  attr(DYN, "bioturbation") <- STD$bioturbation
  attr(DYN, "irrigation")   <- STD$irrigation
  attr(DYN, "DistReact") <- STD$DistReact
  attr(DYN,"includeMPB") <- includeMPB 
  attr(DYN, "Parms") <- STD$Parms
  attr(DYN, "warnings") <- ZZ
  attr(DYN, "model") <- paste("CNPDIA_model_",model,sep="")

  return(DYN)
}

