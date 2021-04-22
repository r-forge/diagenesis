\name{FESDIAdefault}
\alias{FESDIAdefault}
\alias{FESDIAdefault2}
\title{
  Model solutions to aid in solving new applications
}
\description{

  \code{FESDIAdefault, FESDIAdefaults2} are two different FESDIA model solutions.
}
\usage{
  data(FESDIAdefault)
}
\examples{

# defaults
  head(FESDIAparms(FESDIAdefault, as.vector = TRUE)) 
#        FESDIAparms(FESDIAdefault2))
  plot(FESDIAdefault, xyswap = TRUE, ylim = c(30, 0),
    mfrow = c(4,4), grid = FESDIAdepth(FESDIAdefault))
  
}

\references{
  Soetaert K, PMJ Herman and JJ Middelburg, 1996a.
  A model of early diagenetic processes from the shelf to abyssal depths.
  Geochimica Cosmochimica Acta, 60(6):1019-1040.

  Soetaert K, PMJ Herman and JJ Middelburg, 1996b.
  Dynamic response of deep-sea sediments to seasonal variation: a model.
  Limnol. Oceanogr. 41(8): 1651-1668.
}

\keyword{ data }
