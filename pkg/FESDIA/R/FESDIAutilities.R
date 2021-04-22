#####################################################################################################
######                       FESDIA: C, N, P, Fe, S, O2 diagenesis                             ######
#####################################################################################################

## --------------------------------------------------------------------------------------------------
## Check the forcing functions, calculate the average over "times" and update the parameter vector
## --------------------------------------------------------------------------------------------------

CreateMeanPars <- function(Parms, name, forc, times) {
  
  if (is.null(forc$data))
    return(Parms)
  else if (is.null(times))
    Parms[[name]] <- mean(forc$data[,2])
  else if (length(times) > 1){
    Data <- as.matrix(forc$data)
    ii <- which(Data[,1] >= min(times) & Data[,1] <= max(times))
    if (!length(ii)) 
      stop("cannot run model: forcing data set for", name, "is not compatible with (spinup) times")
    if (length(ii) == 1)
      res <- Data[ii,2]
    else res <- approx(x = Data[ ,1], y = Data[ ,2], xout = times)$y
    Parms[[name]] <- mean(res)
  } else 
    Parms[[name]] <- approx(x = forc$data[,1], y=forc$data[,2], xout = times, rule = 2)$y
  Parms 
}

## --------------------------------------------------------------------------------------------------
## Checks the forcing functions and creates time series in case they are not provides as such
## --------------------------------------------------------------------------------------------------

Setforcings <- function(Parms, name, forc, times, fac = 1){
  if (is.nan(fac) | is.infinite(fac)) fac <- 1
  
  if (is.null(forc$data)) {
    meanforc <- Parms[[name]]*fac
    if(forc$amp[1] == 0 | forc$pow[1] == 0)
      forcings <- cbind(range(times), meanforc)
    else {
      Val <- pmax(forc$min[1], meanforc*(1 + forc$amp[1]*sin((times - forc$phase[1])/forc$period[1]*2*pi))^forc$pow[1])
      if (forc$n >1)
          for (i in  2:forc$n) 
            Val <- Val + pmax(forc$min[i], meanforc*(1 + forc$amp[i]*sin((times - forc$phase[i])/forc$period[i]*2*pi))^forc$pow[i])
          
      Val <- meanforc / mean(Val) * Val
      forcings <- cbind(times, Val)
    }  
  } else {
    forcings  <-  as.matrix(forc$data)
    forcings[,2] <- forcings[,2]*fac
  }
  forcings  
}

## --------------------------------------------------------------------------------------------------
## Checks a forcing function declaration
## --------------------------------------------------------------------------------------------------

checkforcs <- function(forc, name) {
  formal <- list(data = NULL, amp = 0, period = 365, phase = 0, pow = 1, min = 0)
  if (is.null(forc)) 
    return(formal)

  if (is.data.frame(forc) | is.matrix(forc))
    forc <- list(data = forc)
  
  if (! is.list(forc))
    stop ("forcing ", name, " should be a list or NULL")
  
  nms <- names(formal)
  formal[(namc <- names(forc))] <- forc
  if (length(noNms <- namc[!namc %in% nms]) > 0)
    warning("unknown names in forcing ",name, ": ", paste(noNms, collapse = ", "))

  NN <- c("amp", "period", "phase", "pow", "min")
  formal$n <- max(unlist(lapply(formal[NN], FUN = length)))
  if (formal$n > 1) 
    for (i in NN) formal[[i]] <- rep(formal[[i]], length.out = formal$n)
  formal
}

## --------------------------------------------------------------------------------------------------
## Checks a profile and returns a list with int and mid.
## --------------------------------------------------------------------------------------------------

CheckProfile <- function(prof, Name, interface = TRUE) {  
  if (is.null(prof)) prof <- 0 
  nms <- names(prof)
  if (is.list (prof)) {
    if (! "int" %in%  nms)
      stop (Name, " should be a list containing 'int' and 'mid'")
    if (! "mid" %in%  nms)
      stop (Name, " should be a list containing 'int' and 'mid'")
    if (length(prof$int) != .FESDIA$N +1)
      stop ("Checking", Name, ": 'int' should be a vector of length ", .FESDIA$N +1)
    if (length(prof$mid) != .FESDIA$N )
      stop ("Checking", Name, ": 'mid' should be a vector of length ", .FESDIA$N )
    return(prof)
  }
  if (interface) {
    if (length(prof) == 1) prof <- rep(prof, .FESDIA$N+1)
    if (length(prof) != .FESDIA$N+1)
      stop ("Checking", Name, ": 'should be a vector of length ", .FESDIA$N +1)
    new <- list(int = prof, mid = 0.5*(prof[-1] + prof[-(.FESDIA$N+1)]))
  } else {
    if (length(prof) == 1) prof <- rep(prof, .FESDIA$N)
    if (length(prof) != .FESDIA$N)
      stop ("Checking", Name, ": 'should be a vector of length ", .FESDIA$N)
    new <- list(mid = prof, ind = c(prof[1], 0.5*(prof[-1] + prof[-(.FESDIA$N)]),prof[.FESDIA$N]))
  }  
  return (new)
}  

