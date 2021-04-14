## ============================================================================
## ====================================================================
## A local environment for non user-visible data, 
## for microphytobenthos (MPB) part of the model
## ====================================================================
## ============================================================================

.MPBDIA <- new.env()

##------------------------------------
## extra state variables
##------------------------------------

.MPBDIA$ynames <- c("MPBC", "MPBN", "CHL", "EPS", "PSIIin")
.MPBDIA$svar   <- .MPBDIA$ynames

.MPBDIA$yunits <- c("mmolC/m3 solid", "mmolN/m3 solid",
                    "mgChl/m3 solid", "mmolC/m3 solid", "-")
.MPBDIA$ydescrip <- c("Microphytobenthos C","Microphytobenthos N",
   "Microphytobenthos Chlorophyll", "Extracellular Polymeric Substances", 
   "Fraction of PSII in inhibited state")
.MPBDIA$ynamesall <- as.vector(sapply(.MPBDIA$ynames, FUN = function(x) 
                                      rep(x, times = .CNPDIA$N)))

##------------------------------------
## Parameters
##------------------------------------

.MPBDIA$Parms <- c(
  ## Light and deposition
    MPBflux   = 0       , # nmolC/cm2/day  - input of MPB from the water
    light     = 100     , # uEinst/m2/s    - light at surface 
    kdChl     = 0.0002  , # /(mgChl/m3)/cm - Chl-specific extinction
      
  ## Light and deposition
    QmaxCHLN  = 3.8     , # mgChl:mmolN    - maximum Chl:N ratio
    QmaxNC    = 0.2     , # mmolC:mmolN,     maximum N:C ratio
    QminNC    = 0.05    , # mmolC:mmolN,     minimum N:C ratio

    sigmaPSII = 3.31    , # m2/umolQuanta    absorption cross-section of PSII 
    kdH       = 4.5e-8  , # -                Dimensionless damage rate constant
    tau       = 6.48e-8 , # d                Turnover time of the electron transport chain              
    kr        = 22.464  , # /d               Rate of PSII damage repair                                        
    phiM      = 5.2e-5  , # mmolC/umolQuanta Maximum quantum yield of photosynthesis            
    a         = 0.042   , # m2/mgChl         Spectrally-integrated Chla-specific absorption coeff  
  
# Nutrient dynamics
    maxUpNO3  = 0.2     , # mmolN/mmolC/d    maximal uptake rate of NO3
    maxUpNH3  = 0.2     , # mmolN/mmolC/d    maximal uptake rate of NH3
  
# Respiration
    respMPB   = 0.1     , # /d               basal respiration rate of MPB
    rG        = 0.25    , # mmolC/mmolC      Cost of growth
    rVNO3     = 0.387   , # mmolC/mmolN      Cost of nitrate uptake
    rVN       = 0.198   , # mmolC/mmolN      Cost of ammonium/don uptake
    rRNO3     = 1.98    , # mmolC/mmolN      Cost of nitrate reduction

    kEps      = 0.2     , # -                fraction of Carbohydrate exudated
    maxGrazing= 0.1     , # /d               mortality rate induced by grazing
    deepMort  = 0.01    , # /d               mortality rate in anoxic zone
    rEPS      = 0.26    , # /d               rate of EPS remineralization
    MPBfall   = 0.01    , # cm/d             sinking rate of suspended MPB
    ksWater   = 6e-5      # /cm              extinction coeff of light in water
)

.MPBDIA$Parunit <- c("nmolC/cm2/day", "uEinst/m2/s", "/(mgChl/m3)/cm", 
 "mgChl:mmolN", "mmolC:mmolN", "mmolC:mmolN",              
 "m2/umolQuanta", "-", "d", "/d", "mmolC/umolQuant", "m2/mgChl",   
 "mmolN/mmolC/d","mmolN/mmolC/d", "/d", "mmolC/mmolC", "mmolC/mmolN", 
 "mmolC/mmolN", "mmolC/mmolN", "-", "/d", "/d", "/d", "cm/d", "/cm")

.MPBDIA$Pardesc <- c("deposition of MPB from the water",
    "light intensity at surface",
    "Chl-specific extinction coefficient",
    "maximum Chl:N ratio of MPB",
    "maximum N:C ratio of MPB",
    "minimum N:C ratio of MPB",

    "absorption cross-section of Photosystem II of MPB",
    "dimensionless damage rate constant of PSII",
    "turnover time of the electron transport chain of PSII",
    "Rate of PSII damage repair",
    "Maximum quantum yield of photosynthesis",
    "Spectrally-integrated Chla-specific absorption coeff",
    "maximal uptake rate of NO3 by MPB",
    "maximal uptake rate of NH3 by MPB",

    "basal respiration rate of MPB",
    "Cost of growth",
    "Cost of nitrate uptake",
    "Cost of ammonium/don uptake",
    "Cost of nitrate reduction",

    "fraction of photosynthesis exudated as Carbohydrate",
    "mortality rate of MPB induced by grazing",
    "mortality rate of MPB in anoxic zone",
    "rate of EPS remineralization",
    
    "sinking rate of suspended MPB",
    "extinction coeff of light in water"
    )

##------------------------------------
## forcing functions - overrules MPDIAforc
##------------------------------------

.MPBDIA$varforc  <- .CNPDIA$varforc
iL <- which (.MPBDIA$varforc == "MPBprod")
.MPBDIA$varforc[iL]  <- "LightSurf"

.MPBDIA$unitforc <- .CNPDIA$unitforc
.MPBDIA$unitforc[iL] <- "uEinst/m2/s"

.MPBDIA$descripforc <- .CNPDIA$descripforc
.MPBDIA$descripforc[iL] <- "Light intensity on sediment-water interface"

##------------------------------------
## Variables
##------------------------------------

.MPBDIA$var0D <- c("MPBCflux", "MPBCdeepflux", "MPBNflux", "MPBNdeepflux", 
    "Chlflux", "Chldeepflux",  "EPSflux", "EPSdeepflux", 
    "TotMPBprod", "TotMPBresp", "TotMPBCdeath", "TotMPBNdeath", 
    "TotProdEPS", "TotChlProd", 
    "TotPrimProd", "TotMinEps", "TotFDETprodMPBdeath", "TotDICprodMPBdeath", 
    "TotNH3prodMPBdeath", "TotPO4prodMPBdeath", "TotMPBO2uptake", "TotalMPBC", 
    "TotalMPBN", "TotalCHL", "TotalEPS")

.MPBDIA$unit0D <- c("nmolC/cm2/d", "nmolC/cm2/d", "nmolN/cm2/d", 
     "nmolN/cm2/d", "ngChl/cm2/d", "ngChl/cm2/d", "nmolC/cm2/d",
     "nmolC/cm2/d", "nmolC/cm2/d", "nmolC/cm2/d", "nmolC/cm2/d", 
     "nmolN/cm2/d","nmolC/cm2/d","ngChl/cm2/d", "molC/gChl/d",
     "nmolC/cm2/d", "nmolC/cm2/d", "nmolC/cm2/d", "nmolN/cm2/d",
     "nmolP/cm2/d", "nmolO2/cm2/d","nmolC/cm2", "nmolN/cm2",
                    "ngChl/cm2", "nmolC/cm2")

.MPBDIA$descrip0D <- c("MPBC influx sediment-water", 
      "MPBC efflux lower boundary", 
      "MPBN influx sediment-water", 
      "MPBN efflux lower boundary", 
      "MPB Chl influx sediment-water", 
      "MPB Chl efflux lower boundary",   
      "EPS influx sediment-water", 
      "EPS efflux lower boundary", 
  
      "Integrated MPB C production (DIC uptake)", 
      "Integrated MPB respiration", 
      "Integrated MPB C-mortality", 
      "Integrated MPB N-mortality", 
      "Integrated MPB EPS production", 
      "Integrated MPB Chl production",
      "Integrated Chl-specific primary production", 
      "Integrated EPS mineralisation", 
      "Integrated flux from MPBdeath to FDET", 
      "Integrated stoichiometric flux from MPBdeath to DIC", 
      "Integrated stoichiometric flux from MPBdeath to NH3", 
      "Integrated stoichiometric flux from MPBdeath to PO4",
      "Integrated net oxygen uptake by MPB (respiration)",
      "Integrated Microphytobenthos C",
      "Integrated Microphytobenthos N",
      "Integrated Microphytobenthos Chlorophyll", 
      "Integrated Extracellular Polymeric Substances")

.MPBDIA$var1D <- c("PrimProd",   
      "Chlproduction", "MPBResp", "MPBCdeath", "MPBNdeath", 
      "Chldeath", "ProdEPS", "MinEPS", "ChlCrMPB", "NCrMPB", 
      "ChlNrMPB", "CNrMPB", "Light", "FDETprodMPBdeath", "DICprodMPBdeath", 
      "NH3prodMPBdeath", "PO4prodMPBdeath")
.MPBDIA$unit1D <- c("molC/gChl/d",
       "ngChl/cm3 solid/d", "nmolC/cm3 solid/d", "nmolC/cm3 solid/d", 
       "nmolN/cm3 solid/d", "ngChl/cm3 solid/d", "nmolC/cm3 solid/d",  
       "nmolC/cm3 solid/d", "gChl/molC", "molN/molC", "gChl/molN", 
       "molC/molN", "uEinst/m2/d", "nmolC/cm3 solid/d", "nmolC/cm3 liquid/d", 
       "nmolN/cm3 liquid/d", "nmolP/cm3 liquid/d")

.MPBDIA$descrip1D <- c("Chl-specific MPB C production profile", 
       "MPB Chl production profile", 
       "MPB C-respiration profile", 
       "MPB C-mortality profile", 
       "MPB N-mortality profile", 
       "MPB Chl-mortality profile", 
       "MPB EPS production profile", 
       "EPS mineralisation profile", 
       "Chlorophyll:carbon ratio MPB", 
       "Nitrogen:carbon ratio MPB", 
       "Chlorophyll:nitrogen ratio MPB", 
       "Carbon:nitrogen ratio MPB", 
       "Light profile (daily dose)", 
       "Stoichiometric flux from MPB death to FDET", 
       "Stoichiometric flux from MPB death to DIC", 
       "Stoichiometric flux from MPB death to NH3", 
       "Stoichiometric flux from MPB death to PO4")
      
.MPBDIA$var1Dall <- as.vector(sapply(.MPBDIA$var1D, FUN = function(x) 
                                    rep(x, times = .CNPDIA$N)))
.MPBDIA$outnames <- c(.MPBDIA$var0D, .MPBDIA$var1Dall)
.MPBDIA$nout <- length(.MPBDIA$outnames)



# initial conditions: MPB present in the first upper cms, below nothing
.MPBDIA$MPBCini   <- c(rep(100, 25), rep(0, .CNPDIA$N-25))       
# Internal N in stiochiometric proportion to MPBC
.MPBDIA$MPBNini   <- .MPBDIA$MPBCini * 0.5*(.MPBDIA$Parms["QminNC"] + .MPBDIA$Parms["QmaxNC"])
# Internal Chla in stoichiometric proportion to MPBN
.MPBDIA$CHLini    <- .MPBDIA$MPBNini*2       
.MPBDIA$yini    <- c(.MPBDIA$MPBCini, .MPBDIA$MPBNini, .MPBDIA$CHLini, 
                     rep(0., .CNPDIA$N), rep(0, .CNPDIA$N))

