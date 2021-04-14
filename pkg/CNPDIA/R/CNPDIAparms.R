##===========================================
## Interrogation Functions for CNPDIA models
##===========================================

##------------------------------------
## Get parameters and values
##------------------------------------

CNPDIAparms <- function(out = NULL, as.vector = FALSE, which = NULL) {
  if (is.null(out))
    Parms <- .CNPDIA$Parms
  else if ("MPBDIAdyn" %in% class(out) | "MPBDIAstd" %in% class(out))
    return(MPBDIAparms(out = out, as.vector= as.vector, which = which)) 
  else if ("steady1D" %in% class(out))
    Parms <- out$Parms[1:length(.CNPDIA$Parms)]
  else if ("deSolve" %in% class(out))
    Parms <- attr(out, "Parms")[1:length(.CNPDIA$Parms)]
  else stop("object 'out' not supported")
  
  if (as.vector) {
    if (! is.null(which))
      Parms <- Parms[which]
    return(Parms)
  } else {
    Units <- .CNPDIA$Parunit
    if  (Parms["BCupLiq"] == 1) {
      i1 <- which(names(Parms) == "O2bw")
      Units[i1:(i1+6)] <- "nmol/cm2/d"
    }  
    if  (Parms["BCdownLiq"] == 1) {
      i <- which(names(Parms) == "O2dw")
      Units[i1:(i1+6)] <- "nmol/cm2/d"
    }  
    Parms <- data.frame(parms = Parms, units = Units, description = .CNPDIA$Pardesc)
    if (! is.null(which))
      Parms <- Parms[which, ]
    return(Parms)
  }
}

##------------------------------------
## Grid
##------------------------------------

CNPDIAgrid <- function(out = NULL) {
  if (is.null(out))
    D <- .CNPDIA$Grid 
  else  if ("steady1D" %in% class(out))
    D <- out$Grid
  else stop("object 'out' not supported for grid calculation - try CNPDIAdepth intstead")
  D
} 

CNPDIAdepth <- function(out = NULL) {
 if (is.null(out))
   D <- .CNPDIA$Grid$x.mid
 else  if ("steady1D" %in% class(out))
   D <- out$Depth
 else if ("deSolve" %in% class(out))
   D <- attr(out, "Depth")
 else stop("object 'out' not supported")
 D
} 

CNPDIAdx <- function(out = NULL) {
  if (is.null(out))
    D <- .CNPDIA$Grid$dx
  else  if ("steady1D" %in% class(out))
    D <- out$dx
  else if ("deSolve" %in% class(out))
    D <- attr(out, "dx")
  else stop("object 'out' not supported")
  D
} 

CNPDIApor <- function(out) {
  if (missing(out))
    stop("out' needs to be given for the porosity")
  if ("steady1D" %in% class(out))
   D <- out$porosity
 else if ("deSolve" %in% class(out))
   D <- attr(out, "porosity")
 else stop("object 'out' not supported")
 D
} 

CNPDIAbiot <- function(out) {
  if (missing(out))
    stop("out' needs to be given for the bioturbation")
  
  if ("steady1D" %in% class(out))
    D <- out$bioturbation
  else if ("deSolve" %in% class(out))
    D <- attr(out, "bioturbation")
  else stop("object 'out' not supported")
  D
} 

CNPDIAirr <- function(out) {
  if (missing(out))
    stop("out' needs to be given for the irrigation")
  
  if ("steady1D" %in% class(out))
    D <- out$irrigation
  else if ("deSolve" %in% class(out))
    D <- attr(out, "irrigation")
  else stop("object 'out' not supported")
  D
} 

##------------------------------------
## Get variables 
##------------------------------------
MeanVal <- function(out)  # takes into account unequal timing 
  (colSums(diff(out[,1])*(out[-1,]+out[-nrow(out),])*0.5)/(out[nrow(out),1]-out[1,1]))[-1]

CNPDIA0D <- function(out, as.vector = FALSE, which = NULL) {
  if (missing(out)) {
    Dnames <- c(.CNPDIA$var0D,.CNPDIA$varforc)
    D <- rep(NA, times = length(Dnames))
    names(D) <- Dnames
  }  else if ("MPBDIAdyn" %in% class(out) | "MPBDIAstd" %in% class(out))
    return(MPBDIA0D(out = out, as.vector= as.vector, which = which)) 

  else if ("steady1D" %in% class(out))
   D <- unlist(out[c(.CNPDIA$var0D,.CNPDIA$varforc)])
 else if ("deSolve" %in% class(out))
   D <- MeanVal(out[, c("time",.CNPDIA$var0D,.CNPDIA$varforc)])
 else stop("object 'out' not supported")
 
 if (! as.vector)
   D <- data.frame(names = names(D), values = D, 
      units = c(.CNPDIA$unit0D,.CNPDIA$unitforc),
      description = c(.CNPDIA$descrip0D,.CNPDIA$descripforc))
 
 if (! is.null(which)){
   if (is.vector(D))
     D <- D[which]
   else D <- D[which,]  
 } 
 row.names(D) <- NULL
 D
} 

CNPDIA1D <- function(out, which = NULL) {
  if (missing(out)) {
     D <- data.frame(names = c(.CNPDIA$svar,.CNPDIA$var1D), 
            units = c(.CNPDIA$yunits, .CNPDIA$unit1D), 
            description = c(.CNPDIA$ydescrip, .CNPDIA$descrip1D))
     if (! is.null(which))
       D <- D[(D$names %in% which), ]
     return(D)
  }
  if ("MPBDIAdyn" %in% class(out) | "MPBDIAstd" %in% class(out))
    return(MPBDIA1D(out = out, which = which)) 
  if ("steady1D" %in% class(out))
   D <- cbind(out$y, as.data.frame(out[.CNPDIA$var1D]))
 else if ("deSolve" %in% class(out))  {
   D <- NULL
   for (cc in c(.CNPDIA$svar,.CNPDIA$var1D))
     D <- cbind(D,MeanVal(cbind(out[,1],subset(out, which = cc))))
   rownames(D) <- NULL
   colnames(D) <- c(.CNPDIA$svar,.CNPDIA$var1D)
   D <- as.data.frame(D)  
 }
 else stop("object 'out' not supported")
 
 D <- cbind(x = CNPDIAdepth(out), por = CNPDIApor(out), D)
 if (! is.null(which))
   D <- D[ ,c( "x", "por", which)]
 D  
} 
  
CNPDIAsvar <- function(out, which = NULL) {
  if (missing(out)) 
    return(data.frame(names = .CNPDIA$svar, units = .CNPDIA$yunits, description = .CNPDIA$ydescrip))
  
  else if ("MPBDIAdyn" %in% class(out) | "MPBDIAstd" %in% class(out))
    return(MPBDIAsvar(out = out, which = which)) 
  if ("steady1D" %in% class(out))
    D <- out$y
  else if ("deSolve" %in% class(out))  {
    D <- NULL
    for (cc in .CNPDIA$svar)
      D <- cbind(D,MeanVal(cbind(out[,1],subset(out, which = cc))))
    rownames(D) <- NULL
    colnames(D) <- .CNPDIA$svar
    D <- as.data.frame(D)  
  }
  else stop("object 'out' not supported")
  
  D <- cbind(x = CNPDIAdepth(out), por = CNPDIApor(out), D)
  if (! is.null(which))
    D <- D[ ,c( "x", "por", which)]
  D  
} 


## ============================================================================
## ============================================================================
##   Functions to extract parameters and variables from MPBDIA models
## ============================================================================
## ============================================================================

MPBDIAparms <- function(out = NULL, as.vector = FALSE, which = NULL) {
  if (is.null(out))
    Parms <- c(.CNPDIA$Parms, .MPBDIA$Parms)
  else if ("steady1D" %in% class(out))
    Parms <- out$Parms
  else if ("deSolve" %in% class(out))
    Parms <- attr(out, "Parms")
  else stop("object 'out' not supported")
  
  if (as.vector) {
    if (! is.null(which))
      Parms <- Parms[which]
    return(Parms)
  } else {
    Pn <- names(Parms)
    Units <- c(.CNPDIA$Parunit, .MPBDIA$Parunit)
    Parms <- data.frame(parms = Parms, units = Units, 
       description = c(.CNPDIA$Pardesc, .MPBDIA$Pardesc))
    row.names(Parms) <- Pn   
    if (! is.null(which))
      Parms <- Parms[which, ]
    return(Parms)
  }
}

MPBDIA0D <- function(out, as.vector = FALSE, which = NULL) {
  if (missing(out)) {
    Dnames <- c(.CNPDIA$var0D,.MPBDIA$varforc,.MPBDIA$var0D)
    D <- rep(NA, times = length(Dnames))
    names(D) <- Dnames
 } else if ("steady1D" %in% class(out))
   D <- unlist(out[c(.CNPDIA$var0D,.MPBDIA$varforc,.MPBDIA$var0D)])
 else if ("deSolve" %in% class(out))
   D <- MeanVal(out[, c("time",.CNPDIA$var0D,.MPBDIA$varforc,.MPBDIA$var0D)])
 else stop("object 'out' not supported")
 
 if (! as.vector)
   D <- data.frame(names = names(D), values = D, 
     units = c(.CNPDIA$unit0D,.MPBDIA$unitforc, .MPBDIA$unit0D), 
     description = c(.CNPDIA$descrip0D, .MPBDIA$descripforc, .MPBDIA$descrip0D))
 
 if (! is.null(which)){
   if (is.vector(D))
     D <- D[which]
#   else if ("deSolve" %in% class(out))
#     D <- D[, which]  
   else D <- D[which,]  
 } 
 D
} 

MPBDIA1D <- function(out, which = NULL) {
  if (missing(out)) 
    return(data.frame(names = c(.CNPDIA$var1D, .MPBDIA$var1D), 
      units = c(.CNPDIA$unit1D,.MPBDIA$unit1D), 
      description = c(.CNPDIA$descrip1D,.MPBDIA$descrip1D)))
  
 if ("steady1D" %in% class(out))
   D <- cbind(out$y, as.data.frame(out[c(.CNPDIA$var1D, .MPBDIA$var1D)]))
 else if ("deSolve" %in% class(out))  {
   D <- NULL
   for (cc in c(.CNPDIA$svar,.MPBDIA$svar,.CNPDIA$var1D,.MPBDIA$var1D))
     D <- cbind(D,colMeans(subset(out, which = cc)))
   rownames(D) <- NULL
   colnames(D) <- c(.CNPDIA$svar,.MPBDIA$svar,.CNPDIA$var1D,.MPBDIA$var1D)
   D <- as.data.frame(D)  
 }
 else stop("object 'out' not supported")
 
 D <- cbind(x = CNPDIAdepth(out), por = CNPDIApor(out), D)
 if (! is.null(which))
   D <- D[ ,c( "x", "por", which)]
 D  
} 
  
MPBDIAsvar <- function(out, which = NULL) {
  if (missing(out)) 
    return(data.frame(names = c(.CNPDIA$svar, .MPBDIA$svar), 
    units = c(.CNPDIA$yunits,.MPBDIA$yunits), 
    description = c(.CNPDIA$ydescrip, .MPBDIA$ydescrip)))
  
  if ("steady1D" %in% class(out))
    D <- out$y
  else if ("deSolve" %in% class(out))  {
    D <- NULL
    for (cc in c(.CNPDIA$svar,.MPBDIA$svar))
      D <- cbind(D,colMeans(subset(out, which = cc)))
    rownames(D) <- NULL
    colnames(D) <- c(.CNPDIA$svar,.MPBDIA$svar)
    D <- as.data.frame(D)  
  }
  else stop("object 'out' not supported")
  
  D <- cbind(x = CNPDIAdepth(out), por = CNPDIApor(out), D)
  if (! is.null(which))
    D <- D[ ,c( "x", "por", which)]
  D  
} 
