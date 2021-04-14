\name{CNPDIAbudget}
\alias{CNPDIAbudgetO2}
\alias{CNPDIAbudgetC}
\alias{CNPDIAbudgetN}
\alias{CNPDIAbudgetP}
\title{
  Budgets for the CNPDIA or MPBDIA model.
}
\description{
 \code{CNPDIAbudgetO2, CNPDIAbudgetC, CNPDIAbudgetN, CNPDIAbudgetP} estimate mass 
 budgets from CNPDIA model solutions.
}
\usage{
  CNPDIAbudgetO2(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  CNPDIAbudgetC(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  CNPDIAbudgetN(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  CNPDIAbudgetP(out, ..., which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
}
\arguments{

  \item{out }{an output object returned by \link{CNPDIAsolve}, \link{MPBDIAsolve},
   or \link{CNPDIAdyna}, \link{MPBDIAdyna}. 
  }

  \item{which }{if not \code{NULL}, a vector with names of the items to return.
  }
  \item{... }{unused.
  }
}
\value{
 \code{CNPDIAbudgetx} returns the element budget (C, N, P, O2) of the solution as a \code{list}, with the following items. 
 
 \itemize{
  \item \code{Fluxes}, the boundary fluxes at the surface and bottom of the sediment, the perturbation fluxes (only when the model was solved with \link{CNPDIAperturb}) and the net input. Positive fluxes are directed into the sediment for the surface, and out of the sediment at the bottom. Negative perturbation fluxes are directed out of the sediment. For dynamic runs, fluxes are averaged over the simulation period.
  \item \code{Rates}, the integrated process rates, in nmol/cm2/d
  \item \code{dC}, the rate of change of integrated values of state vairables, in nmol/cm2/d, defined as ([concentration at the end] - [concentration at beginning])/[length of simulation]
  \item \code{Losses}, the total amount lost from the system (burial, and removal e.g. N2 production for N-budget)
  \item \code{Delta}, the difference between total fluxes in and fluxes out, i.e. the deviation from steady-state.   
  \item \code{Fluxmat}, the flux matrix (rows: from, columns:to). For N, C and P budget the column sums - row sums is equal to the rate of change. For the O2 budget, this only applies to O2.
 }

If more than one CNPDIA object is passed to these functions, a matrix is returned, one column for each object. 

For dynamic runs, the budget is taken over the mean of the run; Delta is then the integrated mean rate of change; for steady-state runs, Delta should be very small.

}

\author{
  Karline Soetaert
}

\details{
  
}

\examples{

# some runs to work with  
  defsteady <- CNPDIAsolve()
  defdyn    <- CNPDIAdyna()

# altered steady-state run
  out <- CNPDIAsolve(parms = list(Cflux = 1000))
  cbind(default = CNPDIAparms(), altered = CNPDIAparms(out))

# budgets
  B1 <- CNPDIAbudgetO2(out)  
  B1$Fluxmat
  colSums(B1$Fluxmat) - rowSums(B1$Fluxmat)  # Small only for O2
  CNPDIAbudgetO2(out, defsteady)$Rates  
  
  B2 <- CNPDIAbudgetC(out)
  colSums(B2$Fluxmat) - rowSums(B2$Fluxmat)  # Small for all states
  
  CNPDIAbudgetC(out,defsteady)
  CNPDIAbudgetC(out,defsteady,defdyn, which = "Rates")
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
