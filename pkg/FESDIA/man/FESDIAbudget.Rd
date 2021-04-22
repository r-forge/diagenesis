\name{FESDIAbudget}
\alias{FESDIAbudgetO2}
\alias{FESDIAbudgetC}
\alias{FESDIAbudgetN}
\alias{FESDIAbudgetP}
\alias{FESDIAbudgetS}
\alias{FESDIAbudgetFe}
\title{
  Budgets for the CNPDIA or MPBDIA model.
}
\description{
 \code{FESDIAbudgetO2, FESDIAbudgetC, FESDIAbudgetN, FESDIAbudgetP} estimate mass 
 budgets from FESDIA model solutions.
}
\usage{
  FESDIAbudgetO2(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudgetC(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudgetN(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudgetP(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudgetS(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  FESDIAbudgetFe(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
}
\arguments{

  \item{out }{an output object returned by \link{FESDIAsolve}
   or \link{FESDIAdyna}. 
  }

  \item{which }{if not \code{NULL}, a vector with names of the items to return.
  }
  \item{... }{unused.
  }
}
\value{
 \code{FESDIAbudgetx} returns the element budget (C, N, P, S, Fe, O2) of the solution as a \code{list}, with the following items. 
 
 \itemize{
  \item \code{Fluxes}, the boundary fluxes at the surface and bottom of the sediment, the perturbation fluxes (only when the model was solved with \link{FESDIAperturb}) and the net input. Positive fluxes are directed into the sediment for the surface, and out of the sediment at the bottom. Negative perturbation fluxes are directed out of the sediment. For dynamic runs, fluxes are averaged over the simulation period.
  \item \code{Rates}, the integrated process rates, in nmol/cm2/d
  \item \code{dC}, the rate of change of integrated values of state vairables, in nmol/cm2/d, defined as ([concentration at the end] - [concentration at beginning])/[length of simulation]
  \item \code{Losses}, the total amount lost from the system (burial, and removal e.g. N2 production for N-budget)
  \item \code{Delta}, the difference between total fluxes in and fluxes out, i.e. the deviation from steady-state.   
  \item \code{Fluxmat}, the flux matrix (rows: from, columns:to). For N, C, Fe, S and P budget the column sums - row sums is equal to the rate of change. For the O2 budget, this only applies to O2.
 }

If more than one FESDIA object is passed to these functions, a matrix is returned, one column for each object. 

For dynamic runs, the budget is taken over the mean of the run; Delta is then the integrated mean rate of change; for steady-state runs, Delta should be very small.

}

\author{
  Karline Soetaert
}

\examples{

# some runs to work with  
  defsteady <- FESDIAsolve()
  defdyn    <- FESDIAdyna()

# altered steady-state run
  out <- FESDIAsolve(parms = list(Cflux = 1000))
  cbind(default = FESDIAparms(), altered = FESDIAparms(out))

# budgets
  B1 <- FESDIAbudgetO2(out)  
  B1$Fluxmat
  colSums(B1$Fluxmat) - rowSums(B1$Fluxmat)  # Small only for O2
  FESDIAbudgetO2(out, defsteady)$Rates  
  
  B2 <- FESDIAbudgetC(out)
  colSums(B2$Fluxmat) - rowSums(B2$Fluxmat)  # Small for all states
  
  FESDIAbudgetC(out,defsteady)
  FESDIAbudgetC(out,defsteady,defdyn, which = "Rates")
}

\references{
  Soetaert K, PMJ Herman and JJ Middelburg, 1996a.
  A model of early diagenetic processes from the shelf to abyssal depths.
  Geochimica Cosmochimica Acta, 60(6):1019-1040.

  Soetaert K, PMJ Herman and JJ Middelburg, 1996b.
  Dynamic response of deep-sea sediments to seasonal variation: a model.
  Limnol. Oceanogr. 41(8): 1651-1668.
}

\keyword{ utilities }
