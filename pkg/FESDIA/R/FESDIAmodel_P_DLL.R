######                       FESDIA: C, N, P, Fe, S, O2 diagenesis                             ######
#####################################################################################################

## --------------------------------------------------------------------------------------------------
## Initialisation: grid generations (layers, bioturbation, irrigation), parameters
## --------------------------------------------------------------------------------------------------

initFESDIA <- function (parms = list(), gridtype = 1, CfluxForc = NULL, 
                        FeOH3fluxForc = NULL, CaPfluxForc = NULL,
                        O2bwForc  = NULL, NO3bwForc = NULL, NO2bwForc = NULL, 
                        NH3bwForc = NULL,
                        CH4bwForc = NULL, FebwForc = NULL, H2SbwForc = NULL,
                        SO4bwForc = NULL, PO4bwForc = NULL, DICbwForc = NULL, 
                        ALKbwForc = NULL,
                        gasfluxForc = NULL, wForc = NULL, 
                        biotForc = NULL, irrForc = NULL, ratefactor = NULL, 
                        rFastForc = NULL, rSlowForc = NULL, pFastForc = NULL, 
                        MPBprodForc = NULL, HwaterForc = NULL,
                        CaCO3fluxForc = NULL, ARAGfluxForc = NULL, CabwForc = NULL,
                        Grid = NULL, porosity = NULL, bioturbation = NULL, 
                        irrigation = NULL, surface = NULL, 
                        diffusionfactor = NULL,
                        times = NULL, model = 1, dynamicpH = FALSE)  {
  

  if (is.null(Grid))
    Grid  <- setup.grid.1D(x.up = 0, dx.1 = 0.01, N = .FESDIA$N, L = 100)
  else {
    if (is.list (Grid)) {
      nms <- names(Grid)
      if (! "dx" %in%  nms)
          stop ("'Grid' should be a list containing 'dx' and 'dx.aux'")
      if (! "dx.aux" %in%  nms)
          stop ("'Grid' should be a list containing 'dx' and 'dx.aux'")
      if (length(Grid$dx.aux) != .FESDIA$N +1)
            stop ("Checking 'Grid': 'dx.aux' should be a vector of length ", .FESDIA$N+1)
      if (length(Grid$dx) != .FESDIA$N +1)
          stop ("Checking 'Grid': 'dx' should be a vector of length ", .FESDIA$N)
    } else {
      if (length(Grid) != .FESDIA$N +1)
        stop ("Checking 'Grid': should be a vector of length ", .FESDIA$N+1)
      Grid <- list(dx.aux = Grid, mid = 0.5*(Grid[-1] + Grid[-(.FESDIA$N+1)]))
    }
  }   

## check parameter inputs
  Parms <- c(.FESDIA$Parms, .PHDIA$Parms)
  nms <- names(Parms)
  Parms[(namc <- names(parms))] <- parms
  if (length(noNms <- namc[!namc %in% nms]) > 0)
    warning("unknown names in parms: ", paste(noNms, collapse = ", "))
  PP <- unlist(Parms)
  
  # parameters need to be set based on forcing functions (if present)  
  Parms <- CreateMeanPars (Parms, "Cflux"    , CfluxForc    , times) 
  Parms <- CreateMeanPars (Parms, "FeOH3flux", FeOH3fluxForc , times) 
  Parms <- CreateMeanPars (Parms, "CaPflux"  , CaPfluxForc , times) 
  Parms <- CreateMeanPars (Parms, "w"        ,  wForc   , times) 
  Parms <- CreateMeanPars (Parms, "biot"     ,  biotForc   , times) 
  Parms <- CreateMeanPars (Parms, "irr"      ,   irrForc   , times) 
  Parms <- CreateMeanPars (Parms, "rFast"    , rFastForc   , times) 
  Parms <- CreateMeanPars (Parms, "rSlow"    , rSlowForc   , times) 
  Parms <- CreateMeanPars (Parms, "pFast"    , pFastForc   , times) 
  Parms <- CreateMeanPars (Parms, "MPBprod"  , MPBprodForc , times)   
  Parms <- CreateMeanPars (Parms, "gasflux"  , gasfluxForc , times) 
  Parms <- CreateMeanPars (Parms, "O2bw"     , O2bwForc    , times) 
  Parms <- CreateMeanPars (Parms, "NO3bw"    , NO3bwForc   , times) 
  Parms <- CreateMeanPars (Parms, "NO2bw"    , NO2bwForc   , times) 
  Parms <- CreateMeanPars (Parms, "NH3bw"    , NH3bwForc   , times) 
  Parms <- CreateMeanPars (Parms, "CH4bw"    , CH4bwForc   , times) 
  Parms <- CreateMeanPars (Parms, "Febw"     , FebwForc    , times) 
  Parms <- CreateMeanPars (Parms, "H2Sbw"    , H2SbwForc   , times) 
  Parms <- CreateMeanPars (Parms, "SO4bw"    , SO4bwForc   , times) 
  Parms <- CreateMeanPars (Parms, "PO4bw"    , PO4bwForc   , times) 
  Parms <- CreateMeanPars (Parms, "DICbw"    , DICbwForc   , times) 
  Parms <- CreateMeanPars (Parms, "ALKbw"    , ALKbwForc   , times) 
  Parms <- CreateMeanPars (Parms, "Hwater"   , HwaterForc , times) 

#  if (dynamicpH) {
   Parms <- CreateMeanPars (Parms, "CaCO3flux", CaCO3fluxForc , times) 
   Parms <- CreateMeanPars (Parms, "ARAGflux" , ARAGfluxForc  , times) 
   Parms <- CreateMeanPars (Parms, "Cabw"     , CabwForc      , times) 
#  }
    
  if (Parms[["BCupLiq"]] != 3) {
    if (any(is.na(Parms[c("O2bw","NO3bw","NO2bw","NH3bw","CH4bw","Febw",
                          "SO4bw","PO4bw","DICbw","ALKbw")])))
      stop("bottom water concentrations cannot be NA if type flux or concentration")

    if (dynamicpH & is.na(Parms["Cabw"]))
      stop("bottom water concentration of Ca cannot be NA if type flux or concentration")
  }
  
  if (Parms[["BCdownLiq"]] != 3) {
    if (any(is.na(Parms[[c("O2dw","NO3dw","NO2dw","NH3dw","CH4dw","Fedw",
                           "SO4dw","PO4dw","DICdw","ALKdw")]])))
      stop("deep water concentrations cannot be NA if type flux or concentration")

    if (dynamicpH & is.na(Parms[["Cadw"]]))
      stop("deep water concentration of Ca cannot be NA if type flux or concentration")
  }
  
  Parms[which(is.na(Parms))] <- 0
  
  if (! is.null(surface)) {
    Aint <- CheckProfile (surface, "surface", interface = TRUE)$int
  }  else {
    if (gridtype == 1)                        # cartesian
      Aint <- rep(1, .FESDIA$N+1)
    else if (gridtype == 2)                   # cylindrical
      Aint <- rev(2*pi*Grid$x.int)
    else if (gridtype == 3)                   # spherical
     Aint <- rev(pi*(Grid$x.int)^2)
  }

  DF <- diffcoeff(S = Parms[["salinity"]], t = Parms[["temperature"]])           #Alkalinity7 ~ HCO3
  Parms <- c(Parms, DF[c("O2","NO3","NO2", "NH4","H2PO4","CH4","HCO3","Fe","HS","SO4","HCO3")]*86400e4) # from m2/s -> cm2/d
   

# porosity gradients
  exp.profile <- function(x, y.0, y.inf, x.att = 1, L = 0)
           return(y.inf + (y.0 - y.inf) * exp(-pmax(0, x-L)/x.att))

# Bioturbation profile
  if (is.null(bioturbation)) 
    bioturbation <- setup.prop.1D(func = exp.profile,
                           grid = Grid,
                           y.0 = Parms[["biot"]], y.inf = 0.,
                           L = Parms[["biotdepth"]], x.att = Parms[["biotatt"]])
  else
    bioturbation <- CheckProfile (bioturbation, "bioturbation", interface = TRUE)

  Db <- bioturbation$int

# Irrigation profile  
  if (is.null(irrigation))
    irrigation <- setup.prop.1D(func = exp.profile,
                         grid = Grid,
                         y.0 = Parms[["irr"]], y.inf = 0.,
                         L = Parms[["irrdepth"]], x.att = Parms[["irratt"]])
  else
    irrigation <- CheckProfile(irrigation,"irrigation", interface = FALSE)
  
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
  if (Parms[["rSurfH2Sox"]] > 0 | Parms[["rSurfCH4ox"]] > 0 )
    distreact <- setup.prop.1D(func = exp.profile,
                               grid = Grid,
                               y.0 = 1, y.inf = 0.,
                               L = Parms[["ODUoxdepth"]], x.att = Parms[["ODUoxatt"]])$mid 
  else 
    distreact <- rep(0, .FESDIA$N)

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
                "ODUoxdepth","ODUoxatt", "por0", "pordeep",   
                "porcoeff", "formationtype", "temperature", "salinity", 
                "Cflux", "FeOH3flux", "CaPflux", "w", "pFast", 
                "rFast", "rSlow", "pFast", "MPBprod","gasflux",
                "O2bw","Febw","H2Sbw","SO4bw","NO3bw","NO2bw","NH3bw",
                "PO4bw","CH4bw","DICbw","ALKbw","Hwater")
#  if (dynamicpH) {  # other diffusion coefficients, and density
    Parms <- c(Parms, DF[c("HCO3", "CO3", "Ca")]*86400e4) # from m2/s -> cm2/d
    Parms <- c(Parms, dens = sw_dens(S = Parms[["salinity"]], t = Parms[["temperature"]]))
    toremove <- c(toremove, "CaCO3flux", "ARAGflux", "Cabw")   # these are forcings
#  }

  parms <- unlist(Parms[-which(nms %in%toremove)])

  #  if (dynamicpH) 
  parms <- c(parms, temperature = Parms[["temperature"]], salinity = Parms[["salinity"]])


# parameters to pass to model DLL
  initpar <- c(parms, Grid$dx, Grid$dx.aux, Grid$x.mid, Aint, porGrid$mid,  
              porGrid$int, diffusionfactor, Db, Dirr, distreact)

  list(initpar = initpar, other = Parms[c("biot", "irr")], Parms = PP, 
       Grid = Grid, porGrid = porGrid, Isirr = sum(Dirr) > 0, 
       Dirr = Dirr, bioturbationGrid = bioturbation, 
       diffcoeff = initpar[c("O2", "NO3", "NO2", "NH4", "H2PO4", "CH4", 
         "HCO3", "Fe", "HS", "SO4", "CO3", "Ca" )], 
       Depth = Grid$x.mid, dx = Grid$dx, porosity = porGrid$mid, porint = porGrid$int, 
       bioturbation = bioturbation$mid, bioturbint = bioturbation$int, 
       irrigation = Dirr, diffusionfactor = diffusionfactor, Aint = Aint,
       DistReact = distreact, isDistReact = (max(distreact) > 0))

}

## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------
## solve the steady-state condition - main function
## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------

FESDIAsolve <- function (parms = list(), yini = NULL, gridtype = 1, Grid = NULL,  porosity = NULL, 
                         bioturbation = NULL, irrigation = NULL, surface = NULL, diffusionfactor = NULL, 
                         dynamicbottomwater = FALSE, ratefactor = NULL, calcpH = FALSE, verbose = FALSE, 
                         method = NULL, times = c(0, 1e6), ...)   {
  model <- 1
  if (dynamicbottomwater) model <- 2
  std <- FESDIAsolve_Full (parms = parms, yini = yini, gridtype = gridtype, Grid = Grid, 
                    porosity = porosity, bioturbation = bioturbation, irrigation = irrigation, 
                    surface = surface, diffusionfactor = diffusionfactor,  
                    model =  model, ratefactor = ratefactor, calcpH = calcpH, verbose = verbose, 
                    method = method, times = times, dynamicpH = FALSE, ...)
  class(std) <- unique(c("FESDIAstd", class(std)))
  std
}                         

FESDIAsolve_Full <- function (parms = list(), yini = NULL, gridtype = 1, Grid = NULL,  porosity = NULL, 
                         bioturbation = NULL, irrigation = NULL, surface = NULL, diffusionfactor = NULL, 
                         model = 1, ratefactor = NULL, calcpH = FALSE, verbose = FALSE, 
                         method = NULL, times = c(0, 1e6), dynamicpH = FALSE, ...)   {
  N <- .FESDIA$N
  if (model == 2) N <- .FESDIA$N+1
  
  if (! is.null(method))
   if ( method %in% c("runsteady", "mixed")) {
    DIA <- FESDIAdyna_Full (parms = parms, times = times, gridtype = gridtype, yini = yini, 
                   Grid = Grid, porosity = porosity, bioturbation = bioturbation, 
                   irrigation = irrigation, surface = surface, 
                   diffusionfactor = diffusionfactor, model = model, 
                   ratefactor = ratefactor, verbose = verbose, calcpH = calcpH, 
                   dynamicpH = dynamicpH, ...) 
     Att <- attributes(DIA)[-1]
     DIA <- DIA[2, -1]
     SVAR <- DIA[1:(N*17)]
     
     if (method == "runsteady"){  # finished
       SVNAMES  <- unique(names(SVAR))
       STD  <- list(y = matrix(nrow = N, ncol = 17, data  = SVAR))
       colnames(STD$y) <- SVNAMES
       VAR <- DIA[-(1:(N*17))]
       OUT <- unique(names(VAR))
       Z <- lapply(OUT, FUN = function(x) {
        VV <- subset(VAR, names(VAR) == x)
        names(VV) <- NULL
        VV })
       names(Z) <- OUT   
     
       Att <- Att[-1]
       nn <- c("Depth", "dx", "Aint", "Parms", "Grid", "porosity", "porint", 
         "diffusionfactor", "bioturbation", "bioturbint",
         "irrigation", "gridtype", "Grid",
         "diffcoeff", "isDistReact", "DistReact")
       STD <- c(STD, Z, Att[nn])  
       STD$model <- paste("FESDIA_model", model, sep = "_")
       attributes(STD) <- c(attributes(STD)[1], Att[!names(Att) %in% nn])
       class(STD) <- unique(c("FESDIAstd", "steady1D",  "rootSolve", "list"))
       return (STD)
     } else {
       yini <- SVAR
     }
  }   
  
  std <- FESDIAsolve_full (parms, gridtype, Grid = Grid, porosity = porosity, bioturbation = bioturbation, 
                    irrigation = irrigation, surface = surface, 
                    diffusionfactor = diffusionfactor, yini = yini, ratefactor = ratefactor,
                    model = model, verbose = verbose, calcpH = calcpH,
                    times = times, method = NULL, dynamicpH = dynamicpH, ...)
  class(std) <- unique(c("FESDIAstd", class(std)))
  std

}

## --------------------------------------------------------------------------------------------------

FESDIAsolve_full <- function (parms = list(), gridtype = 1, CfluxForc = NULL, 
                         FeOH3fluxForc = NULL, CaPfluxForc = NULL,
                         O2bwForc  = NULL, NO3bwForc = NULL, NO2bwForc = NULL, NH3bwForc = NULL,
                         CH4bwForc = NULL, FebwForc = NULL, H2SbwForc = NULL,
                         SO4bwForc = NULL, PO4bwForc = NULL, DICbwForc = NULL, ALKbwForc = NULL,
                         gasfluxForc = NULL, wForc = NULL, biotForc= NULL, irrForc = NULL, ratefactor = NULL,
                         rFastForc = NULL, rSlowForc = NULL, pFastForc = NULL, MPBprodForc = NULL,
                         HwaterForc = NULL, CaCO3fluxForc = NULL, ARAGfluxForc = NULL, CabwForc = NULL,
                         times = NULL, Grid = NULL, porosity = NULL, 
                         bioturbation = NULL, irrigation = NULL, surface = NULL, 
                         diffusionfactor = NULL, yini = NULL, model = 1, 
                         verbose = FALSE, calcpH = FALSE,  method = NULL, dynamicpH = FALSE, ...)  {

  P <- initFESDIA(parms = parms, gridtype = gridtype, CfluxForc = CfluxForc,  
                  FeOH3fluxForc = FeOH3fluxForc, CaPfluxForc = CaPfluxForc,
                  O2bwForc  = O2bwForc, NO3bwForc = NO3bwForc, 
                  NO2bwForc = NO2bwForc, NH3bwForc = NH3bwForc,
                  CH4bwForc = CH4bwForc, FebwForc = FebwForc, H2SbwForc = H2SbwForc,
                  SO4bwForc = SO4bwForc, PO4bwForc = PO4bwForc, DICbwForc = DICbwForc, ALKbwForc = ALKbwForc,
                  gasfluxForc = gasfluxForc, wForc = wForc, biotForc = biotForc, irrForc= irrForc,  ratefactor = ratefactor,
                  rFastForc = rFastForc, rSlowForc = rSlowForc, pFastForc = pFastForc,
                  MPBprodForc = MPBprodForc, HwaterForc = HwaterForc,
                  CaCO3fluxForc = CaCO3fluxForc, ARAGfluxForc = ARAGfluxForc, CabwForc = CabwForc,
                  times = times, Grid = Grid,
                  porosity = porosity, bioturbation = bioturbation, 
                  irrigation = irrigation, surface = surface, 
                  diffusionfactor = diffusionfactor, model = model, dynamicpH = dynamicpH)
                     
  initfunc <- "initfesdia"
  initforc <- "initfesforc"
  nspec    <- 17
  ynames   <- .FESDIA$ynames
  nout     <- .FESDIA$nout
  outnames <- .FESDIA$outnames
  
  if (dynamicpH){
    if (model == 1) {  # No dynamic bottom water
      HWATERForc <- as.double(0)
      func <- "phdiamod"
      N <- .FESDIA$N
    } else {
      HWATERForc <- as.double(P$Parms["Hwater"])
      func <- "phdiamodbw"
      N <- .FESDIA$N + 1
    } 
    initfunc <- "initph"
#    initforc <- "initphforc"
    nspec <- nspec + 3
    ynames <- c(ynames , .PHDIA$ynames)
    outnames <- c(outnames, .PHDIA$outnames)
    nout <- nout+.PHDIA$nout
 
 } else { # NO dynamic pH
    if (model == 1) {
      HWATERForc <- as.double(0)
      func <- "fesdiamod"
      N <- .FESDIA$N
    } else {
      HWATERForc <- as.double(P$Parms["Hwater"])
      func <- "fesdiamodbw"
      N <- .FESDIA$N + 1
    }  
  }
  forcings <- c(as.double(P$Parms["Cflux"]), as.double(P$Parms["FeOH3flux"]), 
                as.double(P$Parms["CaPflux"]),as.double(P$Parms["w"]),
                as.double(1.), as.double(1.),as.double(P$Parms["rFast"]),   # for Db irr ratefactor: relative
                as.double(P$Parms["rSlow"]), as.double(P$Parms["pFast"]),
                as.double(P$Parms["gasflux"]), as.double(P$Parms["MPBprod"]), 
                as.double(P$Parms["O2bw"]),
                as.double(P$Parms["NO3bw"]), as.double(P$Parms["NO2bw"]), 
                as.double(P$Parms["NH3bw"]), 
                as.double(P$Parms["CH4bw"]), as.double(P$Parms["Febw"]), 
                as.double(P$Parms["H2Sbw"]), as.double(P$Parms["SO4bw"]), 
                as.double(P$Parms["PO4bw"]), as.double(P$Parms["DICbw"]),
                as.double(P$Parms["ALKbw"]), HWATERForc, as.double(1.0))   # last = ratefactor
#  if (dynamicpH){    # 3 extra forcing functions
    forcings <- c(forcings,  
                as.double(P$Parms["CaCO3flux"]), as.double(P$Parms["ARAGflux"]), 
                as.double(P$Parms["Cabw"])) 
#  }

  Random <- is.null(yini)
  Yini <- yini  
  if (is.null(Yini)) 
     Yini   <- rep(10, nspec*N)

  else if(length(Yini) != nspec*N)
    stop ("'yini' not of correct length, should be ", nspec*N)

  sparse   <- "1D"
  jactype <- NULL
  if (is.null(method)) {
    if (P$isDistReact) 
     method <- "stodes"
    else if ( model == 2  & P$Isirr) 
      method <- "stodes"
    else  
      method <- "stode"
  }     

  if (method == "stodes") sparse <- "sparseint"

  inz <- NULL
  
  if (method == "stodes") 
    sparse <- "sparseint"
  if (P$isDistReact)      
    sparse <- "sparseint"
  if (model == 2 & P$Isirr) 
    sparse <- "sparseint"
    
  
  ZZ <- NULL
  if ((method == "runsteady" | method == "mixed") & (any(is.na(times)) | any(is.infinite(times))))
    stop ("cannot combine method = 'runsteady' or 'mixed' with 'times' that is either NA or Inf")
  
  # suppress the warnings and the messages
  if (any(is.na(times))){

# expand forcings to be a "time series"
  lf <- as.list(forcings)
  forcings <- lapply(lf, FUN = function(x) cbind(t = c(0, Inf), v = x))
    
  ZZ <- c(ZZ, capture.output(suppressWarnings(   
     DIA  <- DLLfunc(y = as.double(Yini), func = func, initfunc = initfunc,
                    initforc = initforc, forcings = forcings, parms = P$initpar,
                    dllname="FESDIA", times = 0,
                    nout = nout, outnames = outnames, ...)
  )))
   attr(DIA,"message") <- ZZ 
   return(DIA)
  } 

  ZZ <- c(ZZ, capture.output(suppressWarnings(   
    DIA  <- steady.1D(y = Yini, func = func, initfunc = initfunc,
                     names = ynames, initforc = initforc,
                     jactype = sparse, method = method, inz = inz,
                     forcings = forcings, initpar = unlist(P$initpar), nspec = nspec,
                     dllname = "FESDIA", maxiter = 100, 
                     nout = nout, outnames = outnames,
                     positive = TRUE, times = times, ...)
  )))
  
  niter <- 1
  y1 <- c(1e2, 1e4, 10, 5, 1, 1000, 10, 1e4,1e4,1e3,1,1e4,1e2,1e4,1e2,0,1e4)
  y2 <- c(1e1, 1e3, 100, 10, 10, 10, 1, 1e2,1e2,1e3,1e2,1e4,1,3e4,1,0,1e3)
#          ("FDET","SDET","O2","NO3","NO2","NH3","DIC","Fe","FeOH3","H2S","SO4","CH4","PO4","FeP","CaP","Pads","ALK")                 

  while (! attributes(DIA)$steady & niter <= 50 & method != "runsteady")  {
   if (Random) {
     aa <- runif(1)
     if(aa > 0.6)
       Yini <- runif(N*nspec)
     else if (aa > 0.3)
       Yini <- as.vector(sapply(y1*runif(nspec), FUN = function(x) rep(x, times = N)))
     else
       Yini <- as.vector(sapply(y2*runif(nspec), FUN = function(x) rep(x, times = N)))
     
   } else
     Yini <- runif( N*nspec)*yini
  
#   Yini <- runif(nspec * .FESDIA$N)*niter 
   ZZ <- c(ZZ, capture.output(suppressWarnings(   
   DIA  <- steady.1D(y=Yini, func = func, initfunc = initfunc,
                     names = ynames, initforc = initforc,
                     jactype = sparse, method = method,  inz = inz,
                     forcings = forcings, initpar = unlist(P$initpar), nspec = nspec,
                     dllname = "FESDIA", maxiter = 100, 
                     nout = nout, outnames = outnames,
                     positive = TRUE, ...)
    )))
       niter <- niter + 1
  }
  if (verbose) {
    if (!attributes(DIA)$steady) warning("steady-state not reached")
    print(ZZ)
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
  
  DIA$Depth         <- P$Depth
  DIA$dx            <- P$Grid$dx
  DIA$Aint          <- P$Aint
  DIA$initpar       <- unlist(P$initpar)
  DIA$other         <- P$other
  DIA$Parms         <- P$Parms
  DIA$porosity      <- P$porosity
  DIA$porint        <- P$porint
  DIA$diffusionfactor <- P$diffusionfactor
  DIA$bioturbation  <- P$bioturbation
  DIA$bioturbint    <- P$bioturbint
  DIA$irrigation    <- P$irrigation
  DIA$gridtype      <- gridtype
  DIA$Grid          <- P$Grid
  DIA$porGrid       <- P$porGrid
  DIA$numberTries   <- niter
  DIA$warnings      <- ZZ
  DIA$model         <- paste("FESDIA_model_",model,sep="")
  DIA$diffcoeff     <- P$diffcoeff
  DIA$diffusionfactor <- DIA$diffusionfactor
  DIA$DistReact     <- P$DistReact 
  DIA$isDistReact   <- max(P$DistReact) > 0
  DIA$DistProfile   <- P$DistProfile  
  DIA$includepH     <- dynamicpH  

  if (calcpH) {
    pH <- FESDIApH(DIA)
    
    att <- attributes(DIA)  
    cl  <- class(DIA)  
    DIA$pH <- pH
    att$names <- c(att$names, "pH")
    attributes(DIA) <- c(att, attributes(pH)) 
    class(DIA) <- cl
  }

  return(DIA)
}


## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------
## solve the dynamic condition - main function
## --------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------


FESDIAdyna <- function (parms = list(), times = 0:365, spinup = NULL, yini = NULL, 
   gridtype = 1, Grid = NULL, porosity = NULL, bioturbation = NULL, 
   irrigation = NULL, surface = NULL, diffusionfactor = NULL,  
   dynamicbottomwater = FALSE, 
   CfluxForc = NULL,FeOH3fluxForc = NULL, CaPfluxForc = NULL,  O2bwForc   = NULL,
   NO3bwForc  = NULL,  NO2bwForc  = NULL, NH3bwForc  = NULL,  FebwForc   = NULL,  H2SbwForc  = NULL,
   SO4bwForc  = NULL,  CH4bwForc  = NULL,  PO4bwForc  = NULL,  DICbwForc  = NULL, 
   ALKbwForc  = NULL,  wForc      = NULL,  biotForc   = NULL,  
   irrForc    = NULL,  rFastForc  = NULL,
   rSlowForc  = NULL,  pFastForc  = NULL,  MPBprodForc= NULL,  gasfluxForc = NULL, 
   HwaterForc = NULL,  ratefactor = NULL, calcpH = FALSE, verbose = FALSE, ...)  {

  model <- 1
  if (dynamicbottomwater) 
    model <- 2
  
  dyn <- FESDIAdyna_Full (parms = parms, times = times, spinup = spinup, yini = yini, 
     gridtype = gridtype, Grid = Grid, porosity = porosity, 
     bioturbation = bioturbation, irrigation = irrigation, surface = surface, 
     diffusionfactor = diffusionfactor, model = model, CfluxForc = CfluxForc, 
     FeOH3fluxForc = FeOH3fluxForc, CaPfluxForc = CaPfluxForc, O2bwForc = O2bwForc,
     NO3bwForc = NO3bwForc, NO2bwForc = NO2bwForc, NH3bwForc = NH3bwForc,  
     FebwForc = FebwForc, H2SbwForc = H2SbwForc, SO4bwForc = SO4bwForc, 
     CH4bwForc = CH4bwForc, PO4bwForc = PO4bwForc, DICbwForc = DICbwForc, 
     ALKbwForc = ALKbwForc, wForc = wForc, biotForc = biotForc,  
     irrForc = irrForc, rFastForc = rFastForc, rSlowForc = rSlowForc,  
     pFastForc = pFastForc, MPBprodForc = MPBprodForc, gasfluxForc = gasfluxForc, 
     HwaterForc = HwaterForc, ratefactor = ratefactor, calcpH = calcpH, 
     verbose = verbose, dynamicpH = FALSE, ...) 
     
  class(dyn) <- unique(c("FESDIAdyn", class(dyn)))
  dyn
}  

# general dynamic function

FESDIAdyna_Full <- function (parms = list(), times = 0:365, spinup = NULL, yini = NULL, 
   gridtype = 1, Grid = NULL, porosity = NULL, bioturbation = NULL, 
   irrigation = NULL, surface = NULL, diffusionfactor = NULL, model = 1, 
   CfluxForc = NULL,FeOH3fluxForc = NULL, CaPfluxForc = NULL,  O2bwForc   = NULL,
   NO3bwForc  = NULL,  NO2bwForc  = NULL, NH3bwForc  = NULL,  FebwForc   = NULL,  H2SbwForc  = NULL,
   SO4bwForc  = NULL,  CH4bwForc  = NULL,  PO4bwForc  = NULL,  DICbwForc  = NULL, 
   ALKbwForc  = NULL,  wForc      = NULL,  biotForc   = NULL,  
   irrForc    = NULL,  rFastForc  = NULL,
   rSlowForc  = NULL,  pFastForc  = NULL,  MPBprodForc= NULL,  gasfluxForc = NULL, 
   HwaterForc = NULL, CaCO3fluxForc = NULL, ARAGfluxForc = NULL, CabwForc = NULL,  
   ratefactor = NULL, calcpH = FALSE, verbose = FALSE, 
   dynamicpH = FALSE, ...)  {


  CfluxForc     <- checkforcs(CfluxForc,         "CfluxForc")
  FeOH3fluxForc <- checkforcs(FeOH3fluxForc, "FeOH3fluxForc")
  CaPfluxForc   <- checkforcs(CaPfluxForc,     "CaPfluxForc")
  O2bwForc      <- checkforcs(O2bwForc,           "O2bwForc")
  NO3bwForc     <- checkforcs(NO3bwForc,         "NO3bwForc")
  NO2bwForc     <- checkforcs(NO2bwForc,         "NO2bwForc")
  NH3bwForc     <- checkforcs(NH3bwForc,         "NH3bwForc")
  FebwForc      <- checkforcs(FebwForc,           "FebwForc")
  H2SbwForc     <- checkforcs(H2SbwForc,         "H2SbwForc")
  SO4bwForc     <- checkforcs(SO4bwForc,         "SO4bwForc")
  CH4bwForc     <- checkforcs(CH4bwForc,         "CH4bwForc")
  PO4bwForc     <- checkforcs(PO4bwForc,         "PO4bwForc")
  DICbwForc     <- checkforcs(DICbwForc,         "DICbwForc")
  ALKbwForc     <- checkforcs(ALKbwForc,         "ALKbwForc")
  wForc         <- checkforcs(wForc,                 "wForc")
  biotForc      <- checkforcs(biotForc,           "biotForc")
  irrForc       <- checkforcs(irrForc,             "irrForc")
  rFastForc     <- checkforcs(rFastForc,         "rFastForc")
  rSlowForc     <- checkforcs(rSlowForc,         "rSlowForc")
  pFastForc     <- checkforcs(pFastForc,         "pFastForc")
  MPBprodForc   <- checkforcs(MPBprodForc,     "MPBprodForc")
  gasfluxForc   <- checkforcs(gasfluxForc,     "gasfluxForc")
  HwaterForc    <- checkforcs(HwaterForc,       "HwaterForc")
  ratefactor    <- checkforcs(ratefactor,       "ratefactor")

#  if (dynamicpH)
  CaCO3fluxForc <- checkforcs(CaCO3fluxForc, "CaCO3fluxForc")
  ARAGfluxForc  <- checkforcs(ARAGfluxForc ,  "ARAGfluxForc")
  CabwForc      <- checkforcs(CabwForc,           "CabwForc")
  
  if (is.null(spinup)) Times <- times else Times <- spinup
  
  nspec <- 17 
  ynames <- .FESDIA$ynames
  initfunc <- "initfesdia"
  initforc <- "initfesforc"
  outnames <- .FESDIA$outnames
  nout <- .FESDIA$nout

  if (dynamicpH) {
     nspec <- nspec + 3
     ynames <- c(ynames , .PHDIA$ynames)
     initfunc <- "initph"
#     initforc <- "initphforc"
     outnames <- c(outnames, .PHDIA$outnames)
     nout <- nout+.PHDIA$nout
  }
  
  if (is.null(yini)) {
    STD <- FESDIAsolve_full(parms, gridtype = gridtype, CfluxForc = CfluxForc, 
                       O2bwForc  = O2bwForc, NO3bwForc = NO3bwForc, 
                       NO2bwForc = NO2bwForc, NH3bwForc = NH3bwForc,
                       CH4bwForc = CH4bwForc, FebwForc = FebwForc, H2SbwForc = H2SbwForc,
                       SO4bwForc = SO4bwForc, PO4bwForc = PO4bwForc, DICbwForc = DICbwForc, 
                       ALKbwForc = ALKbwForc, gasfluxForc = gasfluxForc, wForc = wForc, 
                       biotForc = biotForc, irrForc = irrForc, ratefactor = ratefactor, 
                       rFastForc = rFastForc, rSlowForc = rSlowForc, pFastForc = pFastForc, 
                       MPBprodForc = MPBprodForc, HwaterForc = HwaterForc, 
                       CaCO3fluxForc = CaCO3fluxForc, ARAGfluxForc = ARAGfluxForc, 
                       CabwForc = CabwForc, times = Times, 
                       Grid = Grid, porosity = porosity,  bioturbation = bioturbation, 
                       irrigation = irrigation, surface = surface, 
                       diffusionfactor = diffusionfactor, model = model, 
                       dynamicpH = dynamicpH, verbose = verbose)
    yini <- STD$y
    
  } else  
    STD <- initFESDIA(parms, gridtype = gridtype, CfluxForc = CfluxForc, 
                       O2bwForc  = O2bwForc, NO3bwForc = NO3bwForc, 
                       NO2bwForc = NO2bwForc, NH3bwForc = NH3bwForc,
                       CH4bwForc = CH4bwForc, FebwForc = FebwForc, H2SbwForc = H2SbwForc,
                       SO4bwForc = SO4bwForc, PO4bwForc = PO4bwForc, DICbwForc = DICbwForc, 
                       ALKbwForc = ALKbwForc, gasfluxForc = gasfluxForc, wForc = wForc, 
                       biotForc = biotForc, irrForc = irrForc, ratefactor = ratefactor, 
                       rFastForc = rFastForc, rSlowForc = rSlowForc, pFastForc = pFastForc, 
                       MPBprodForc = MPBprodForc, HwaterForc = HwaterForc, 
                       CaCO3fluxForc = CaCO3fluxForc, ARAGfluxForc = ARAGfluxForc,
                       CabwForc = CabwForc, times = Times, 
                       Grid = Grid, porosity = porosity, bioturbation = bioturbation, 
                       irrigation = irrigation, surface = surface, 
                       diffusionfactor = diffusionfactor, model = model, dynamicpH = dynamicpH)

  initpar = unlist(STD$initpar)
  
  band   <- ifelse(STD$isDistReact | model == 2, 0, 1)
  
  if (model == 1) {
    func <-"fesdiamod"
    N <- .FESDIA$N
  } else {
    func="fesdiamodbw"
    N <- .FESDIA$N+1  
  }
  if (dynamicpH) func <- "phdiamod"
  
  lrw <- 170000
  if(STD$isDistReact | dynamicpH) lrw <- 500000
  
  if(length(yini) != nspec*N)
    stop ("'yini' not of correct length, should be ", nspec*N)
  
  ZZ <- NULL
  is.spinup <- ! is.null(spinup)
  if (is.spinup)
    is.spinup <- (length(spinup) > 1)
  if (is.spinup) {
    forcings <- list()
    forcings[[1]]  <- Setforcings (STD$Parms, "Cflux", CfluxForc, spinup, fac = 1)
    forcings[[2]]  <- Setforcings (STD$Parms, "FeOH3flux", FeOH3fluxForc, spinup, fac = 1)
    forcings[[3]]  <- Setforcings (STD$Parms, "CaPflux",     CaPfluxForc, spinup, fac = 1)
    forcings[[4]]  <- Setforcings (STD$Parms, "w",         wForc, spinup, fac = 1)
    forcings[[5]]  <- Setforcings (STD$other, "biot",   biotForc, spinup, fac = 1/STD$other[["biot"]])
    forcings[[6]]  <- Setforcings (STD$other, "irr",     irrForc, spinup, fac = 1/STD$other[["irr"]])
    forcings[[7]]  <- Setforcings (STD$Parms, "rFast", rFastForc, spinup, fac = 1)
    forcings[[8]]  <- Setforcings (STD$Parms, "rSlow", rSlowForc, spinup, fac = 1)
    forcings[[9]]  <- Setforcings (STD$Parms, "pFast", pFastForc, spinup, fac = 1)
    forcings[[10]] <- Setforcings (STD$Parms, "gasflux", gasfluxForc, spinup, fac = 1)
    forcings[[11]] <- Setforcings (STD$Parms, "MPBprod", MPBprodForc, spinup, fac = 1)
    forcings[[12]] <- Setforcings (STD$Parms, "O2bw",   O2bwForc, spinup, fac = 1)
    forcings[[13]] <- Setforcings (STD$Parms, "NO3bw", NO3bwForc, spinup, fac = 1)
    forcings[[14]] <- Setforcings (STD$Parms, "NO2bw", NO2bwForc, spinup, fac = 1)
    forcings[[15]] <- Setforcings (STD$Parms, "NH3bw", NH3bwForc, spinup, fac = 1)
    forcings[[16]] <- Setforcings (STD$Parms, "CH4bw", CH4bwForc, spinup, fac = 1)
    forcings[[17]] <- Setforcings (STD$Parms, "Febw",  FebwForc,  spinup, fac = 1)
    forcings[[18]] <- Setforcings (STD$Parms, "H2Sbw", H2SbwForc, spinup, fac = 1)
    forcings[[19]] <- Setforcings (STD$Parms, "SO4bw", SO4bwForc, spinup, fac = 1)
    forcings[[20]] <- Setforcings (STD$Parms, "PO4bw", PO4bwForc, spinup, fac = 1)
    forcings[[21]] <- Setforcings (STD$Parms, "DICbw", DICbwForc, spinup, fac = 1)
    forcings[[22]] <- Setforcings (STD$Parms, "ALKbw", ALKbwForc, spinup, fac = 1)
    forcings[[23]] <- Setforcings (STD$Parms,"Hwater",HwaterForc, spinup, fac = 1)
    forcings[[24]] <- Setforcings (list(ratefactor = 1), "ratefactor", ratefactor, spinup, fac = 1)
    
#    if (dynamicpH){
     forcings[[25]] <- Setforcings (STD$Parms, "CaCO3flux", CaCO3fluxForc, spinup, fac = 1)
     forcings[[26]] <- Setforcings (STD$Parms,  "ARAGflux",  ARAGfluxForc, spinup, fac = 1)
     forcings[[27]] <- Setforcings (STD$Parms,      "Cabw",      CabwForc, spinup, fac = 1)
#    }

    ZZ <- c(ZZ, capture.output(suppressWarnings(   
      DYN <- ode.1D(y = yini, names = ynames, initforc = initforc,
                   forcings = forcings,  times = spinup, method = "lsodes",
                   func = func, initfunc = initfunc, lrw = lrw,
                   initpar = initpar, dllname = "FESDIA", bandwidth = band,
                   nout = nout, outnames = outnames,
                   nspec = nspec, ...)
    )))
    yini <- DYN[nrow(DYN),2:(nspec*N+1)]
  }

  forcings <- list()
    forcings[[1]]  <- Setforcings (STD$Parms, "Cflux", CfluxForc, times, fac = 1)
    forcings[[2]]  <- Setforcings (STD$Parms, "FeOH3flux", FeOH3fluxForc, times, fac = 1)
    forcings[[3]]  <- Setforcings (STD$Parms, "CaPflux",     CaPfluxForc, times, fac = 1)
    forcings[[4]]  <- Setforcings (STD$Parms, "w",         wForc, times, fac = 1)
    forcings[[5]]  <- Setforcings (STD$other, "biot",   biotForc, times, fac = 1/STD$other[["biot"]])
    forcings[[6]]  <- Setforcings (STD$other, "irr",     irrForc, times, fac = 1/STD$other[["irr"]])
    forcings[[7]]  <- Setforcings (STD$Parms, "rFast", rFastForc, times, fac = 1)
    forcings[[8]]  <- Setforcings (STD$Parms, "rSlow", rSlowForc, times, fac = 1)
    forcings[[9]]  <- Setforcings (STD$Parms, "pFast", pFastForc, times, fac = 1)
    forcings[[10]] <- Setforcings (STD$Parms, "gasflux",gasfluxForc, times, fac = 1)
    forcings[[11]] <- Setforcings (STD$Parms, "MPBprod",MPBprodForc, times, fac = 1)
    forcings[[12]] <- Setforcings (STD$Parms, "O2bw",   O2bwForc, times, fac = 1)
    forcings[[13]] <- Setforcings (STD$Parms, "NO3bw", NO3bwForc, times, fac = 1)
    forcings[[14]] <- Setforcings (STD$Parms, "NO2bw", NO2bwForc, times, fac = 1)
    forcings[[15]] <- Setforcings (STD$Parms, "NH3bw", NH3bwForc, times, fac = 1)
    forcings[[16]] <- Setforcings (STD$Parms, "CH4bw", CH4bwForc, times, fac = 1)
    forcings[[17]] <- Setforcings (STD$Parms, "Febw",  FebwForc,  times, fac = 1)
    forcings[[18]] <- Setforcings (STD$Parms, "H2Sbw", H2SbwForc, times, fac = 1)
    forcings[[19]] <- Setforcings (STD$Parms, "SO4bw", SO4bwForc, times, fac = 1)
    forcings[[20]] <- Setforcings (STD$Parms, "PO4bw", PO4bwForc, times, fac = 1)
    forcings[[21]] <- Setforcings (STD$Parms, "DICbw", DICbwForc, times, fac = 1)
    forcings[[22]] <- Setforcings (STD$Parms, "ALKbw", ALKbwForc, times, fac = 1)
    forcings[[23]] <- Setforcings (STD$Parms,"Hwater",HwaterForc, times, fac = 1)
    forcings[[24]] <- Setforcings (list(ratefactor = 1), "ratefactor", ratefactor, times, fac = 1)
    
#    if (dynamicpH){
     forcings[[25]] <- Setforcings (STD$Parms, "CaCO3flux", CaCO3fluxForc, times, fac = 1)
     forcings[[26]] <- Setforcings (STD$Parms,  "ARAGflux",  ARAGfluxForc, times, fac = 1)
     forcings[[27]] <- Setforcings (STD$Parms,      "Cabw",      CabwForc, times, fac = 1)
#    }

  ZZ <- c(ZZ, capture.output(suppressWarnings(   
  DYN <- ode.1D(y = yini, names = ynames, initforc = initforc,
                   forcings = forcings,  times = times, method = "lsodes",
                   func = func, initfunc = initfunc, lrw = lrw,
                   initpar = initpar, dllname = "FESDIA", bandwidth = band,
                   nout = nout, outnames = outnames,
                   nspec = nspec, ...)
  )))
  colnames(DYN)[2:(nspec*N+1)] <- as.vector(sapply(ynames, FUN = function(x) rep(x, times = N)))
  if (model == 2) {
    cn <- colnames(DYN)
    cn[seq(2, by = N, length.out = nspec)] <- paste(ynames, "bw", sep ="")
    colnames(DYN) <- cn
  }
  
  attr(DYN, "Depth") <- STD$Depth
  attr(DYN, "dx")    <- STD$dx
  attr(DYN, "Aint")    <- STD$Aint
  attr(DYN, "Grid")    <- STD$Grid
  attr(DYN, "porosity")   <- STD$porosity
  attr(DYN, "porint")     <- STD$porint
  attr(DYN, "bioturbation") <- STD$bioturbation
  attr(DYN, "bioturbint") <- STD$bioturbint
  attr(DYN, "irrigation")   <- STD$irrigation
  attr(DYN, "warnings") <- ZZ
  attr(DYN, "Parms") <- STD$Parms
  attr(DYN, "diffcoeff") <- STD$diffcoeff
  attr(DYN, "diffusionfactor") <- STD$diffusionfactor
  attr(DYN, "model") <- "FESDIA"

  class(DYN) <- unique(c("FESDIAdyn", class(DYN)))

  if (calcpH) {
    pH <- FESDIApH(DYN)
    
    att <- attributes(DYN)  
    cl  <- class(DYN)
    lv1 <- att$lengthvar[1]  
    i1  <- 1:(lv1+1)
    DYN <- cbind(DYN[ ,i1], pH = pH, DYN[ ,-i1])
    att$dim <- dim(DYN)
    att$nspec  <- att$nspec+1
    att$ynames <- c(att$ynames, "pH")
    att$lengthvar[1] <- lv1 + N
    attributes(DYN) <- c(attributes(DYN)[1:2], att[-(1:2)], attributes(pH)[-(1:2)]) 
    class(DYN) <- cl
  }
  
  return(DYN)
}

