## ====================================================================
## A local environment for non user-visible data,
## for pH part of the model
## ====================================================================

.PHDIA <- new.env()

##------------------------------------
## state variables
##------------------------------------

.PHDIA$ynames   <- c("CaCO3", "ARAG", "Ca")
.PHDIA$svar   <- .PHDIA$ynames

.PHDIA$yunits   <- c("mmolC/m3 solid", "mmolC/m3 solid", "mmolCa/m3 liquid")
.PHDIA$ydescrip <- c("Calcium carbonate",  "Aragonite", "Calcium (Ca2+)")

.PHDIA$ynamesall <- as.vector(sapply(.PHDIA$ynames, FUN = function(x) 
                                      rep(x, times = .FESDIA$N)))
                                      
##------------------------------------
## Parameters
##------------------------------------
                                      
.PHDIA$Parms <- c(

  CaCO3flux = 0.39*1e5/365  , # nmolC/cm2/d   CaCO3     deposition: 0.39 M/m2/yr
  ARAGflux  = 0.21*1e5/365  , # nmolC/cm2/d   Aragonite deposition: 0.21 M/m2/yr
  rCadiss   = 40/365        , # /d            Carbonate dissolution rate
  pCa       = 2.74          , # -             Carbonate dissolution/precipitation power
  rArdiss   = 110/365       , # /d            Aragonite dissolution
  pAr       = 2.43          , # -             Aragonite dissolution/precipitation power
  rCaprod   = 0.            , # /nmol/cm3/d   Carbonate precipitation rate
  
  Cabw      = 10000         , # mmol/m3       Ca++ conc in bottom water
  Cadw      = NA            , # mmol/m3       deep water concentration

  CaCO3fall       = 100     , # cm/day        fall speed of suspended CaCO3
  ARAGfall        = 100       # cm/day        fall speed of suspended aragonite
    )

.PHDIA$Parunit <- c("nmolC/cm2/d", "nmol/cm2/d", "/nmol/cm3/d", "-", "/nmol/cm3/d", "-", "/nmol/cm3/d", "mmol/m3", "mmol/m3", "cm/d", "cm/d")

.PHDIA$Pardesc <- c(
  "CaCO3 deposition rate", "Aragonite deposition rate", 
  "Carbonate dissolution rate", "Carbonate dissolution/precipitation power", 
  "Aragonite dissolution/precipitation rate", "Aragonite dissolution/precipitation power", 
  "Carbonate precipitation rate", 
  "Ca++ conc in bottom water", "deep water Ca2+ concentration", 
  "fall speed of suspended CaCO3", "fall speed of suspended aragonite")


##------------------------------------
## forcing functions 
##------------------------------------

.PHDIA$varforc <- c("CaCO3flux", "ARAGflux", "bwCa")
     
.PHDIA$unitforc <- c("nmolC/cm2/d", "nmolC/cm2/d", "mmol/m3") 

.PHDIA$descripforc <- c("CaCO3 sediment-water deposition",
   "Aragonite sediment-water deposition", "Ca2+ bottom water concentration")

##------------------------------------
## Variables
##------------------------------------

.PHDIA$var0D <- c("Caflux", "Cadeepflux", "CaCO3deepflux", "ARAGdeepflux",
                  "TotCaCO3diss", "TotARAGdiss", "TotCaCO3prod", 
                  "TotalCaCO3", "TotalARAG", "TotalCa",
                  "KH2O", "KHF", "KH2CO3", "KHCO3", "KBOH3", 
                  "KNH4", "KH2S", "KHS", "KH3PO4",
                  "KH2PO4", "KHPO4", "KHSO4", "KH2SO4", "K0CO2", 
                  "KspCalcite", "KspAragonite", "TotToFree",
                  "SWSToFree", "Borate", "Fluoride")

.PHDIA$unit0D <- c("nmolCa/cm2/d", "nmolCa/cm2/d", "nmolC/cm2/d", "nmolC/cm2/d", 
      "nmolC/cm2/d", "nmolC/cm2/d", "nmolC/cm2/d", "nmolC/cm2", "nmolC/cm2", "nmolCa/cm2",
      rep("mol/kg-soln", 16), "-", "-", "micromol/kg", "micromol/kg")

.PHDIA$descrip0D <- c("Ca2+ influx sediment-water",  "Ca2+ efflux lower boundary",
  "CaCO3 efflux lower boundary", "Aragonite efflux lower boundary",
  "Integrated Carbonate dissolution",  "Integrated Aragonate dissolution",  
  "Integrated Carbonate production",  "Integrated Ca Carbonate concentration",  
  "Integrated Aragonite concentration",  "Integrated Ca2+ (ion) concentration",  
  "ion product of water/dissociation ct of H2O",
  "dissociation constant of HF", "dissociation constant of H2CO3",
  "dissociation constant of HCO3", "dissociation constant of BOH3",
  "dissociation constant of NH4", "dissociation constant of H2S", 
  "dissociation constant of HS", "dissociation constant of H3PO4",
  "dissociation constant of H2PO4", "dissociation constant of HPO4",
  "dissociation constant of HSO4", "dissociation constant of H2SO4", 
  "CO2 solubility, Henrys ct", "solubility product for calcite",
  "solubility product for aragonite", "conversion from total to free pH scale",
  "conversion from SWS to free pH scale", "Borate concentration", 
  "Fluoride concentration")
 
.PHDIA$var1D <- c("CaCO3diss", "ARAGdiss", "CaCO3prod", "pH", 
                  "omegaCa", "omegaAr", "CO3", "CO2") 
.PHDIA$unit1D <- c("nmolC/cm3 solid/d", "nmolC/cm3 solid/d", "nmolC/cm3 liquid/d", "-",  "-", "-",
                    "mmolC/m3","mmolC/m3")
.PHDIA$descrip1D <- c("CaCO3 dissolution profile",  "Aragonite dissolution profile", 
              "CaCO3 production profile", "pH profile", 
              "CaCO3 saturation state profile", "Aragonite saturation state profile", 
              "CO3 concentration profile", "CO2 concentration profile")

.PHDIA$var1Dall <- as.vector(sapply(.PHDIA$var1D, FUN = function(x) rep(x, times = .FESDIA$N)))

.PHDIA$outnames <- c(.PHDIA$var0D, .PHDIA$var1Dall, .PHDIA$varforc)

.PHDIA$nout <- length(.PHDIA$outnames)
