##===========================================
## Interrogation Functions for FESDIA models
##===========================================

##------------------------------------
## Get parameters and values
##------------------------------------

FESDIAparms <- function(out = NULL, as.vector = FALSE, which = NULL) {


 if (is.null(out))
   Parms <- .FESDIA$Parms
 else if ("PHDIAdyn" %in% class(out) | "PHDIAstd" %in% class(out))
    return(PHDIAparms(out = out, as.vector= as.vector, which = which)) 
 else if ("steady1D" %in% class(out))
   Parms <- out$Parms[1:length(.FESDIA$Parms)]
 else if ("deSolve" %in% class(out))
   Parms <- attr(out, "Parms")[1:length(.FESDIA$Parms)]
 else stop("object 'out' not supported")
 
 if (as.vector) {
   if (! is.null(which))
     Parms <- Parms[which]
   return(Parms)
 } else {
    Units <- .FESDIA$Parunit
    if  (Parms["BCupLiq"] == 1) {  # upper boundary = FLUX
      i1 <- which(names(Parms) == "O2bw")
      Units[i1:(i1+11)] <- "nmol/cm2/d"
    }  
    if  (Parms["BCdownLiq"] == 1) {
      i <- which(names(Parms) == "O2dw")
      Units[i1:(i1+6)] <- "nmol/cm2/d"
    }  
    Parms <- data.frame(parms = Parms, units = Units, description = .FESDIA$Pardesc)
    if (! is.null(which))
      Parms <- Parms[which, ]
    return(Parms)
 }  
}

FESDIAdepth <- function(out = NULL) {
 if (is.null(out))
   D <- .FESDIA$Grid$x.mid
 else if ("steady1D" %in% class(out))
   D <- out$Depth
 else if ("deSolve" %in% class(out))
   D <- attr(out, "Depth")
 else stop("object 'out' not supported")
 D
} 

FESDIAgrid <- function(out = NULL) {
  if (is.null(out))
    D <- .FESDIA$Grid 
  else  if ("steady1D" %in% class(out))
    D <- out$Grid
  else stop("object 'out' not supported for grid calculation - try FESDIAdepth intstead")
  D
} 

FESDIAdx <- function(out = NULL, mid = TRUE) {
  if (is.null(out))
    D <- .FESDIA$Grid$dx
  else if ("steady1D" %in% class(out))
    D <- out$dx
  else if ("deSolve" %in% class(out))
    D <- attr(out, "dx")
  else stop("object 'out' not supported")
  D
} 

FESDIAbiot <- function(out) {
  if (is.null(out))
    stop("out' needs to be given for the bioturbation")
  else if ("steady1D" %in% class(out))
    D <- out$bioturbation
  else if ("deSolve" %in% class(out))
    D <- attr(out, "bioturbation")
  else stop("object 'out' not supported")
  D
} 

FESDIAirr <- function(out) {
  if (missing(out))
    stop("out' needs to be given for the irrigation")
  else if ("steady1D" %in% class(out))
    D <- out$irrigation
  else if ("deSolve" %in% class(out))
    D <- attr(out, "irrigation")
  else stop("object 'out' not supported")
  D
} 

FESDIApor <- function(out) {
  if (missing(out))
    stop("out' needs to be given for the porosity")
  if ("steady1D" %in% class(out))
   D <- out$porosity
 else if ("deSolve" %in% class(out))
   D <- attr(out, "porosity")
 else stop("object 'out' not supported")
 D
} 

##------------------------------------
## Get variables 
##------------------------------------
MeanVal <- function(out)  # takes into account unequal timing 
  (colSums(diff(out[,1])*(out[-1,]+out[-nrow(out),])*0.5)/(out[nrow(out),1]-out[1,1]))[-1]


FESDIA0D <- function(out, as.vector = FALSE, which = NULL) {
  if (missing(out)) {
    Dnames <- c(.FESDIA$var0D,.FESDIA$varforc)
    D <- rep(NA, times = length(Dnames))
    names(D) <- Dnames
  } 
  
 else if ("PHDIAdyn" %in% class(out) | "PHDIAstd" %in% class(out))
   return(PHDIA0D(out = out, as.vector = as.vector, which = which))

 else if ("steady1D" %in% class(out))
   D <- unlist(out[c(.FESDIA$var0D,.FESDIA$varforc)])
 else if ("deSolve" %in% class(out))
   D <- MeanVal(out[, c("time",.FESDIA$var0D,.FESDIA$varforc)])
 else stop("object 'out' not supported")
 
 if (! as.vector)
   D <- data.frame(names = names(D), values = D, 
      units = c(.FESDIA$unit0D,.FESDIA$unitforc),
      description = c(.FESDIA$descrip0D,.FESDIA$descripforc))
 
 if (! is.null(which)){
   if (is.vector(D))
     D <- D[which]
   else D <- D[which,]  
 } 
 row.names(D) <- NULL
 D
} 

FESDIA1D <- function(out, which = NULL) {
  if (missing(out)) 
    return(data.frame(names = .FESDIA$var1D, 
            units = .FESDIA$unit1D, description = .FESDIA$descrip1D))
  
  if ("steady1D" %in% class(out))
   D <- cbind(out$y, as.data.frame(out[.FESDIA$var1D]))
 else if ("PHDIAdyn" %in% class(out) | "PHDIAstd" %in% class(out))
   return(PHDIA1D(out = out, which = which))
 else if ("deSolve" %in% class(out))  {
   D <- NULL
   for (cc in c(.FESDIA$svar,.FESDIA$var1D))
     D <- cbind(D,MeanVal(cbind(out[,1],subset(out, which = cc))))
   rownames(D) <- NULL
   colnames(D) <- c(.FESDIA$svar,.FESDIA$var1D)
   D <- as.data.frame(D)  
 }
 else stop("object 'out' not supported")
 
 D <- cbind(x = FESDIAdepth(out), por = FESDIApor(out), D)
 if (! is.null(which))
   D <- D[ ,c( "x", "por", which)]
 D  
} 
  
FESDIAsvar <- function(out, which = NULL) {
  if (missing(out)) 
    return(data.frame(names = .FESDIA$svar, units = .FESDIA$yunits, description = .FESDIA$ydescrip))
 else if ("PHDIAdyn" %in% class(out) | "PHDIAstd" %in% class(out))
   return(PHDIAsvar(out = out, which = which))

  if ("steady1D" %in% class(out))
    D <- out$y
  else if ("deSolve" %in% class(out))  {
    D <- NULL
    for (cc in .FESDIA$svar)
      D <- cbind(D,MeanVal(cbind(out[,1],subset(out, which = cc))))
    rownames(D) <- NULL
    colnames(D) <- .FESDIA$svar
    D <- as.data.frame(D)  
  }
  else stop("object 'out' not supported")
  
  D <- cbind(x = FESDIAdepth(out), por = FESDIApor(out), D)
  if (! is.null(which))
    D <- D[ ,c( "x", "por", which)]
  D  
} 



## ============================================================================
## ============================================================================
##   Functions to extract parameters and variables from PHDIA models
## ============================================================================
## ============================================================================

PHDIAparms <- function(out = NULL, as.vector = FALSE, which = NULL) {
  if (is.null(out))
    Parms <- c(.FESDIA$Parms, .PHDIA$Parms)
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
    Units <- c(.FESDIA$Parunit, .PHDIA$Parunit)
    Parms <- data.frame(parms = Parms, units = Units, 
       description = c(.FESDIA$Pardesc, .PHDIA$Pardesc))
    row.names(Parms) <- Pn   
    if (! is.null(which))
      Parms <- Parms[which, ]
    return(Parms)
  }
}

PHDIA0D <- function(out, as.vector = FALSE, which = NULL) {
  if (missing(out)) {
    Dnames <- c(.FESDIA$var0D,.PHDIA$varforc,.PHDIA$var0D)
    D <- rep(NA, times = length(Dnames))
    names(D) <- Dnames
 } else if ("steady1D" %in% class(out))
   D <- unlist(out[c(.FESDIA$var0D,.PHDIA$varforc,.PHDIA$var0D)])
 else if ("deSolve" %in% class(out))
   D <- MeanVal(out[, c("time",.FESDIA$var0D,.PHDIA$varforc,.PHDIA$var0D)])
 else stop("object 'out' not supported")
 
 if (! as.vector)
   D <- data.frame(names = names(D), values = D, 
     units = c(.FESDIA$unit0D,.PHDIA$unitforc, .PHDIA$unit0D), 
     description = c(.FESDIA$descrip0D, .PHDIA$descripforc, .PHDIA$descrip0D))
 
 if (! is.null(which)){
   if (is.vector(D))
     D <- D[which]
#   else if ("deSolve" %in% class(out))
#     D <- D[, which]  
   else D <- D[which,]  
 } 
 D
} 

PHDIA1D <- function(out, which = NULL) {
  if (missing(out)) 
    return(data.frame(names = c(.FESDIA$var1D, .PHDIA$var1D), 
      units = c(.FESDIA$unit1D,.PHDIA$unit1D), 
      description = c(.FESDIA$descrip1D,.PHDIA$descrip1D)))
  
 if ("steady1D" %in% class(out))
   D <- cbind(out$y, as.data.frame(out[c(.FESDIA$var1D, .PHDIA$var1D)]))
 else if ("deSolve" %in% class(out))  {
   D <- NULL
   for (cc in c(.FESDIA$svar,.PHDIA$svar,.FESDIA$var1D,.PHDIA$var1D))
     D <- cbind(D,colMeans(subset(out, which = cc)))
   rownames(D) <- NULL
   colnames(D) <- c(.FESDIA$svar,.PHDIA$svar,.FESDIA$var1D,.PHDIA$var1D)
   D <- as.data.frame(D)  
 }
 else stop("object 'out' not supported")
 
 D <- cbind(x = FESDIAdepth(out), por = FESDIApor(out), D)
 if (! is.null(which))
   D <- D[ ,c( "x", "por", which)]
 D  
} 
  
PHDIAsvar <- function(out, which = NULL) {
  if (missing(out)) 
    return(data.frame(names = c(.FESDIA$svar, .PHDIA$svar), 
    units = c(.FESDIA$yunits,.PHDIA$yunits), 
    description = c(.FESDIA$ydescrip, .PHDIA$ydescrip)))
  
  if ("steady1D" %in% class(out))
    D <- out$y
  else if ("deSolve" %in% class(out))  {
    D <- NULL
    for (cc in c(.FESDIA$svar,.PHDIA$svar))
      D <- cbind(D,colMeans(subset(out, which = cc)))
    rownames(D) <- NULL
    colnames(D) <- c(.FESDIA$svar,.PHDIA$svar)
    D <- as.data.frame(D)  
  }
  else stop("object 'out' not supported")
  
  D <- cbind(x = FESDIAdepth(out), por = FESDIApor(out), D)
  if (! is.null(which))
    D <- D[ ,c( "x", "por", which)]
  D  
} 
