# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------
# Periodic function to generate e.g. variation in light (daily) or tidal height
# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

StepFunction <- function(value = 1, avg = NULL, 
                         period = 1, 
                         fraction = 0.5, 
                         phase = 0,
                         times = seq(0, 5, length.out = 100)){
  
  if (! is.null(avg) & ! is.null(value))
     stop("'avg' and 'value' cannot both have a value - select one")
  if (!is.null(avg)) value <- avg/fraction
  Z <- ((times-phase)%%period  < fraction*period) * value
  cbind(times, Z)
}

PeriodicFunction <- function(avg = 1, max = NULL,   # mean and max y-value, only one should be specified
                             period = 1,            # period in days
                             fraction = 0.5,        # fraction of period in which y is not 0
                             phase = 0,             # time offset in days
                             min = 0,               # minimum y-value
                             pow = 1,               # power
                             times = seq(0, 5, length.out = 100)){  # times, days for which function needs to be applied
  
  if (! is.null(avg) & ! is.null(max))
     stop("'avg' and 'max' cannot both have a value - select one")

  if (is.null(avg) & is.null(max))
     stop("'avg' and 'max' cannot both be NULL - one should have a value")
  
  N <- 1000
  tseq <- seq(from = 0, to = period, length.out = N)
  
  Lfun <- function(x, M = 1, min = 0) 
    M*pmax(0,(0+ (sin(2*pi*tseq/period)+x)))^pow + min

  SW <- function(x = 0){
    sum (Lfun(x) > 0)/N - fraction
  }
  Shift <- uniroot(f= SW, interval = c(-1,1)) 
  x <- Shift$root
  
  if (! is.null(avg))
    M <- avg/mean(Lfun(x = x, min = 0)) - min 
  else
    M <- (max+min)/max(Lfun(x = x, min = 0)) - min
  
  cbind(times,
        M*pmax(0,(0+ (sin(2*pi*(times-phase)/period)+x)))^pow  + min)
}
