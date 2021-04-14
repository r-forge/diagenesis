\name{Periodic functions}
\alias{PeriodicFunction}
\alias{StepFunction}

\title{
  Utility functions for use with the CNPDIA and MPBDIA model.
}
\description{

  \code{PeriodicFunction} generates a periodic function that is a truncated sine wave.
  
  \code{StepFunction} generates a periodic step function
}
\usage{
PeriodicFunction (avg = 1, max = NULL, period = 1, fraction = 0.5, 
                  phase = 0, min = 0, pow = 1,
                  times = seq(0, 5, length.out = 100))

StepFunction (value = 1, avg = NULL, period = 1, 
              fraction = 0.5, phase = 0,
              times = seq(0, 5, length.out = 100))
  
}
\arguments{

  \item{avg, max }{mean and maximum function-value, only one of those should be specified. 
  }
  \item{value }{function-value when not 0, only one of \code{avg, value} should be specified.
  }
  \item{period }{periodicity, in [days].
  }
  \item{fraction }{fraction of period in which the function value is not 0, dimensionless.
  }
  \item{phase }{time offset, in [days].
  }
  \item{min }{minimum y-value.
  }
  \item{pow }{power to which the function value is raised
  }
  \item{times }{time values for which the function needs to be applied, in [days]
  }
}

\details{
 \code{StepFunction} creates a periodic on-off function, where the non-0 (on)-value as given by \code{value} occurs \code{fraction}*\code{period} of the time, and the 0-value (1-\code{fraction})*\code{period} of the time.
For example, this can be used to mimic tides, swiching between inundated and dry.
 
 
 \code{PeriodicFunction} creates a periodic sinusoidal function with a give periodicity, and where the fraction of time that the function is nonzero can be specified. For instance to generate light within a day or within a year.
}

\author{
  Karline Soetaert
}


\examples{

#=====================================
# Tides
#=====================================

# tides with a periodicity of 12.4 hours, inundated half of the time.
# inundated 1/4 of the time

 times     <- seq(0, 5, length.out = 1000)
 Tides     <- StepFunction(period = 12.4/24, fraction = 0.5, 
                           times = times)
 Tides2    <- StepFunction(period = 12.4/24, fraction = 0.25, 
                           times = times, value = 2)
 matplot(x = times, y = cbind(Tides[,2], Tides2[,2]), type = "l", lty = 1)
 c(mean(Tides[,2]), mean(Tides2[,2]))

#=====================================
# Light within a day - 
# different daylength
#=====================================
 
 times     <- seq(0, 5, length.out = 1000)

# daylength = 6, 12, 18 hours (fraction = 6/24, 0.5, 18/24) 
 LightDay  <- PeriodicFunction(period = 1, fraction = 6/24, 
                       times = times, avg = NULL, max = 200)
 LightDay2 <- PeriodicFunction(period = 1, fraction = 0.5, 
                       times = times, avg = NULL, max = 200)
 LightDay3 <- PeriodicFunction(period = 1, fraction = 18/24, 
                       times = times, avg = NULL, max = 200)

 matplot(x = times, y = cbind(LightDay[,2], LightDay2[,2], LightDay3[,2]), 
         type = "l", lty = 1)
 c(mean(LightDay[,2]), mean(LightDay2[,2]), mean(LightDay3[,2]))


#=====================================
# Light within a day 
# sharpness of the peak
#=====================================
 
 times     <- seq(0, 5, length.out = 1000)

# daylength = 12 hours (fraction = 0.5) 
 LightDay  <- PeriodicFunction(period = 1, fraction = 0.5, 
                       times = times, avg = 200)
 LightDay2 <- PeriodicFunction(period = 1, fraction = 0.5, 
                       times = times, avg = 200, pow = 2)

 LightDay3 <- PeriodicFunction(period = 1, fraction = 0.5, 
                       times = times, avg = 200, pow = 0.5)

 matplot(x = times, y = cbind(LightDay[,2], LightDay2[,2], LightDay3[,2]), 
         type = "l", lty = 1)
 c(mean(LightDay[,2]), mean(LightDay2[,2]), mean(LightDay3[,2]))


#=====================================
# Yearly light variations
#=====================================

 LightYear <- PeriodicFunction(period = 365, fraction = 1, times = 1:365, 
                          max = 600, avg = NULL, min = 10, phase = 100)
 plot(LightYear, type = "l")
 
 mean(LightYear[,2])

#=====================================
# Year + daily light variations
# constant daylength
#=====================================

 times <- seq(0, 365, length.out = 10000)
 LightYear <- PeriodicFunction(period = 365, fraction = 1, times = times, 
                          max = 600, avg = NULL, min = 10, phase = 100)
 Dayvar    <- PeriodicFunction(period = 1, fraction = 0.5, times = times, 
                          max = 1, avg = NULL, phase = 0.5)
 LightDay  <- cbind(times, LightYear[,2]*Dayvar[,2])                          

 plot(LightDay, type = "l")
 plot(LightDay[1:1000,], type = "l")
 

#=====================================
# Year + daily light variations
# variable daylength
#=====================================

 times <- seq(0, 365, length.out = 10000)
 LightYear <- PeriodicFunction(period = 365, fraction = 1, times = times, 
                          max = 600, avg = NULL, min = 10, phase = 100)
                          
 Daylength <- PeriodicFunction(period = 365, fraction = 1, times = 0:365, 
                          max = 18, avg = NULL, min = 10, phase = 100)


 Dayvar    <- NULL
 for (i in 1:nrow(Daylength)){
   Dayvar <- rbind( Dayvar, 
      PeriodicFunction(period = 1,  
                       times = times[trunc(times) == Daylength[i,1]], 
                       max = 1, avg = NULL, fraction = Daylength[i,2]/24)
                  )
 } 
 
 LightDay  <- cbind(times, LightYear[,2]*Dayvar[,2])                          

 mf <- par (mfrow = c(2,2))
 plot(Daylength, type = "l", main = "Daylength")
 plot(LightDay, type = "l", main = "Light")
 plot(LightDay[1:100,], type = "l", main = "light")
 plot(LightDay[5000:5100,], type = "l", main = "light")
 
 par(mfrow = mf)
}
\value{
 A matrix with two columns, the \code{times} and the calculated function value.
}

\keyword{ utilities }
