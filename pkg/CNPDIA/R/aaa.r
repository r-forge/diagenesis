## ====================================================================
## A local environment for non user-visible data,
## ====================================================================
.CNPDIA <- new.env()

.CNPDIA$N     <- 100
.CNPDIA$Grid  <- setup.grid.1D(x.up=0, dx.1=0.01, N = .CNPDIA$N, L = 100)

##------------------------------------
## Parameters
##------------------------------------

.CNPDIA$Parms <- c(
  ## organic matter dynamics  #
    Cflux       = 500            ,  # nmolC/cm2/d - Carbon deposition: around 20gC/m2/yr
    pFast       = 0.9            ,  # -           - fraction fast detritus in flux
    FePflux     = 0              ,  # nmolP/cm2/d - deposition rate of FeP
    CaPflux     = 0              ,  # nmolP/cm2/d - deposition rate of CaP
    rFast       = 25/365         ,  # /d          - decay rate fast decay detritus
    rSlow       = 0.05/365       ,  # /d          - decay rate slow decay detritus
    NCrFdet     = 16/106         ,  # molN/molC   - NC ratio fast decay detritus
    NCrSdet     = 16/106         ,  # molN/molC   - NC ratio slow decay detritus
    PCrFdet     = 1/106          ,  # molP/molC   - PC ratio fast decay det.
    PCrSdet     = 1/106          ,  # molP/molC   - PC ratio slow decay det.
    
  ## Nutrient bottom water conditions
    BCupLiq         = 2        ,    # upper boundary condition for liquid; 1= flux, 2=conc, 3 = 0-grad
    BCdownLiq       = 3        ,    # lower boundary condition; default = 0-gradient
    O2bw            = 300      ,    # mmol/m3     if BC = 2: Oxygen conc in bottom water; else: flux 
    NO3bw           = 10       ,    # mmol/m3
    NO2bw           = 0        ,    # mmol/m3
    NH3bw           = 1        ,    # mmol/m3
    ODUbw           = 0        ,    # mmol/m3
    PO4bw           = 0.5      ,    # mmol/m3
    DICbw           = 2200     ,    # mmol/m3
    O2dw            = NA       ,    # mmol/m3  deep-water boundary concentration
    NO3dw           = NA       ,
    NO2dw           = NA        ,    # mmol/m3
    NH3dw           = NA       ,
    ODUdw           = NA       ,
    PO4dw           = NA       ,
    DICdw           = NA       ,

  ## Bioturbation, advection, bio-irrigation
    w               = 0.1/365000  , # cm/d       - advection rate
    biot            = 1/365       , # cm2/d      - bioturbation coefficient
    biotdepth       = 5        ,    # cm         - depth of mixed layer
    biotatt         = 1        ,    # cm         - depth attenuation coefficient below biotdepth
    irr             = 0        ,    # /d         - bio-irrigation rate
    irrdepth        = 5        ,    # cm         - depth of irrigates layer
    irratt          = 1        ,    # cm         - depth attenuation coefficient below irrdepth
    gasflux         = 0        ,    # cm/d       - piston velocity for dry flats - exchange of O2 and DIC only+deposition 
  
  ## Nutrient parameters
    NH3Ads          = 1.3      ,    #-           Adsorption coeff ammonium
    rnitri1         = 20       ,    #/d          Max nitrification rate (step1) - oxidation of ammonium  
    rnitri2         = 20.      ,    #/d          Max nitrification rate (step2) - oxidatino of nitrite
    ksO2nitri       = 1.       ,    #mmolO2/m3   half-sat O2 in nitrification
    ranammox        = 0.1      ,    #/d          Anammox rate
    ksNO2anammox    = 0.1      ,    #mmolNO2/m3  half-sat NO2 in annamox
    rODUox          = 20.      ,    #/d          Max rate oxidation of ODU in one layer
    rSurfODUox      = 0.       ,    #/d          Max rate oxidation of deep ODU with surface O2 - 0 to toggle it off 
    ODUoxdepth      = 5.       ,    #cm          Depth where oxidation of ODU with surface water O2 is possible
    ODUoxatt        = 1.       ,    #cm          the depth attenuation coefficient of ODU oxidation below ODUoxdepth.
    ksO2oduox       = 1.       ,    #mmolO2/m3   half-sat O2 in oxidation of ODU
    ksO2oxic        = 3.       ,    #mmolO2/m3   half-sat O2 in oxic mineralisation
    ksNO3denit      = 30.      ,    #mmolNO3/m3  half-sat NO3 in denitrification
    kinO2denit      = 1.       ,    #mmolO2/m3   half-sat O2 inhib denitrification
    kinNO3anox      = 1.       ,    #mmolNO3/m3  half-sat NO3 inhib anoxic degr
    kinO2anox       = 1.       ,    #mmolO2/m3   half-sat O2 inhib anoxic min
    pdepo           = NA       ,    #-           part ODU produced that is lost - if NA, calculated from w
    rdepo           = 0        ,    #/d          ODU removal rate
   
  ## To estimate diffusion coefficients
    temperature     = 10       ,    # temperature
    salinity        = 35       ,    # 

  ## Other
    TOC0            = 0.5      ,
    rFePdesorp      = 0.01     ,    # /d
    rFePadsorp      = 0.3      ,    # /day    rFePadsorp      = 0.3      
    rCaPprod        = 0        ,
    rCaPdiss        = 0        ,
    CPrCaP          = 1.32/4.6 ,    # Ca10(PO4)4.6(CO3)1.32F1.87(OH)1.45   Jilbert and Slomp
    rPads           = 0.0      ,    # /day - adsorption rate of phosphate
    rPdes           = 0.0      ,    # /day - desorption rate of adsorbed P
    maxPads         = 1e3      ,    # Max adsorbed P concentration - mmolP/m3 solid
    por0            = 0.9      ,    #-  surface porosity
    pordeep         = 0.5      ,    #-  deep porosity
    porcoeff        = 0.3      ,    #cm porosity coefficient
    formationtype   = 1        ,    # to estimate effective diffusion, 1=sand, 2=mud, 3=general
  
  ## parameters if bottom water is dynamically described
    dilution        = 0        ,    # /day   - relaxation towards background conc dissolved
    Hwater          = 10       ,    # cm     - height of water over core
    Cfall           = 100      ,    # cm/day - fall speed of organic carbon (FDET, SDET)
    FePfall         = 100      ,    # cm/day - fall speed of FeP
    CaPfall         = 100      ,    # cm/day - fall speed of CaP
    
  ## MPB production
    MPBprod         =   0      ,    # mmol/m3/d  maximal production rate - range: 5000-5e4
    kdSed           =  20      ,    # /cm,       exponential decay - kd of sediment ~ 2 to 4 /mm
    kNO3upt         =   3      ,    # mmolN/m3   half-saturation conc NO3 uptake
    kNH3upt         =   3      ,    # mmol/m3,   half-saturation concentration for N uptake MPB
    kPO4upt         = 0.1      ,    # mmolP/m3   half-saturation conc PO4 uptake
    kDICupt         =   1           # mmolC/m3   half-saturation conc DIC uptake
  
)

.CNPDIA$Parunit <- c("nmolC/cm2/d","-","nmolP/cm2/d","nmolP/cm2/d",
  "/d","/d","molN/molC","molN/molC","molP/molC","molP/molC",
  "-", "-",rep("mmol/m3",times = 7), rep("mmol/m3",times = 7),
  "cm/d","cm2/d","cm","cm","/d","cm","cm",
  "cm/d", "-","/d","/d","mmolO2/m3","/d","mmolN/m3", 
  "/d","/d","cm","cm","mmolO2/m3","mmolO2/m3","mmolNO3/m3","mmolO2/m3","mmolNO3/m3",
  "mmolO2/m3","-","/d","dgC","psu","%","/d",
  "/d","/d","/d","molC/molP","/d","/d","mmolP/m3 solid", "-","-","cm","-", "/d","cm","cm/d","cm/d","cm/d","mmol/m3/d","/cm", "mmolN/m3", "mmolN/m3", "mmolP/m3", "mmolC/m3")

.CNPDIA$Pardesc <- c("total organic C deposition", "part FDET in carbon flux", 
  "deposition rate of FeP", "deposition rate of CaP",
  "decay rate FDET", "decay rate SDET", "NC ratio FDET", "NC ratio SDET",
  "PC ratio FDET", "PC ratio SDET",
  "upper boundary liq. 1:flux, 2:conc, 3:0-grad",
  "lower boundary liq. 1:flux, 2:conc, 3:0-grad",
  "upper boundary O2  -if BC=1: flux, 2:conc",
  "upper boundary NO3 -if BC=1: flux, 2:conc",
  "upper boundary NO2 -if BC=1: flux, 2:conc",
  "upper boundary NH3 -if BC=1: flux, 2:conc",
  "upper boundary ODU -if BC=1: flux, 2:conc",
  "upper boundary PO4 -if BC=1: flux, 2:conc",
  "upper boundary DIC -if BC=1: flux, 2:conc",
  "lower boundary O2  -if BC=1: flux, 2:conc",
  "lower boundary NO3 -if BC=1: flux, 2:conc",
  "lower boundary NO2 -if BC=1: flux, 2:conc",
  "lower boundary NH3 -if BC=1: flux, 2:conc",
  "lower boundary ODU -if BC=1: flux, 2:conc",
  "lower boundary PO4 -if BC=1: flux, 2:conc",
  "lower boundary DIC -if BC=1: flux, 2:conc",
  "advection rate", "bioturbation coefficient", "depth of mixed layer", 
  "attenuation coeff below biotdepth", "bio-irrigation rate", 
  "depth of irrigated layer", "attenuation coeff below irrdepth", 
  "piston velocity for dry flats", "Adsorption coeff ammonium",
  "Max nitrification rate step1 (NH3ox)", "Max nitrification rate step2 (NO2ox)", "half-sat O2 in nitrification",
  "Anammox rate", "half-sat NO2 in anammox", "Max rate ODU oxidation in one layer", "Max rate ODU oxidation with BW O2",
  "Max depth ODU oxidation with BW O2", "depth attenuation ODU oxidation",
  "half-sat O2 in oxidation of ODU", "half-sat O2 in oxic mineralisation",
  "half-sat NO3 in denitrification", "half-sat O2 inhib denitrification",
  "half-sat NO3 inhib anoxic degr", "half-sat O2 inhib anoxic min",
  "part ODU prod lost (NA:estimated from w)", "ODU removal rate",
  "temperature", "salinity", "refractory Carbon conc",
  "rate FeP desorption", "rate FeP adsorption",
  "rate CaP production", "rate CaP dissolution","C:Pratio in CaP",
  "adsorption rate PO4", "desorption rate of adsorbed P", "Max adsorbed P concentration", 
  "surface porosity", "deep porosity", "porosity decay coefficient",
  "formationfactor, 1=sand,2=fine sand,3=general",
  "relaxation towards background conc ", "height of overlying water",
  "fall speed of organic C (FDET, SDET)", "fall speed of FeP", "fall speed of CaP",
  "maximal MPB production rate", "sedimentary light extinction coefficient",
  "NO3 limitation constant MPB", "NH3 limitation constant MPB", 
  "P limitation constant MPB", 
  "C limitation constant MPB")

##------------------------------------
## State variables
##------------------------------------

.CNPDIA$ynames <- c("FDET", "SDET", "O2", 
                    "NO3", "NO2", "NH3", 
                    "ODU", "DIC",
                    "PO4", "FeP", "CaP", "Pads")
.CNPDIA$svar <- .CNPDIA$ynames

.CNPDIA$yunits <- c("mmolC/m3 solid", "mmolC/m3 solid", "mmolO/m3 liquid",
                   "mmolN/m3 liquid", "mmolN/m3 liquid", "mmolN/m3 liquid",
                   "mmolO/m3 liquid", "mmolC/m3 liquid",
                   "mmolP/m3 liquid", "mmolP/m3 solid", "mmolP/m3 solid", 
                   "mmolP/m3 solid")

.CNPDIA$ydescrip <- c("Fast decaying Detritus (solid)", 
                      "Slow decaying Detritus (solid)",
                      "Oxygen (liquid)", 
                      "Nitrate (liquid)", 
                      "Nitrite (liquid)",
                      "Ammonium/ammonia (liquid)", 
                      "Oxygen Demand Units (liquid)",
                      "Dissolved Inorganic Carbon (liquid)",
                      "Phosphate (liquid)",
                      "Iron-bound P (solid)", 
                      "Ca-bound P (solid)", 
                      "Adsorbed P (solid)")
.CNPDIA$ynamesall <- as.vector(sapply(.CNPDIA$ynames, FUN = function(x) rep(x, times = .CNPDIA$N)))

##------------------------------------
## 0-D Variables
##------------------------------------

.CNPDIA$var0D <- c("O2flux", "O2deepflux", 
      "NO3flux", "NO3deepflux", 
      "NO2flux", "NO2deepflux",
      "NH3flux", "NH3deepflux", 
      "ODUflux", "ODUdeepflux",
      "PO4flux", "PO4deepflux", 
      "DICflux", "DICdeepflux", 
      "FDETflux", "FDETdeepflux", 
      "SDETflux", "SDETdeepflux",
      "FePdeepflux", "CaPdeepflux", 
      "OrgCflux", "OrgNflux", "OrgPflux",
      "O_Cflux", "ODUO_Cflux", 
      "DINDIPflux", "DINDIPmean", "DINDIPdeep",
      "TotMin", "TotOxic", "TotDenit", "TotAnoxic", 
      "PartOxic", "PartDenit", "PartAnoxic", 
      "TotNitri1", "TotNitri2", "TotAnammox", "TotODUoxid",
      "TotFePprod", "TotCaPprod", "TotFePdesorp", "TotCaPdiss", "TotPads",
      "TotODUoxsurf", "TotNH3ads", "PartPremoved", "PartNremoved",
      "TotNH3prod", "TotPO4prod",
      "TotMPBNO3uptake", "TotMPBNH3uptake", "TotMPBPO4uptake", "TotMPBDICuptake",
      "TotMPBO2prod", "TotalFDET", "TotalSDET", "TotalO2", 
      "TotalNO3", "TotalNO2", "TotalNH3", 
      "TotalODU", "TotalDIC", "TotalPO4", "TotalFeP", 
      "TotalCaP", "TotalPads")

.CNPDIA$unit0D <- c("nmolO2/cm2/d", "nmolO2/cm2/d", 
      "nmolN/cm2/d", "nmolN/cm2/d", 
      "nmolN/cm2/d", "nmolN/cm2/d", 
      "nmolN/cm2/d", "nmolN/cm2/d", 
      "nmolO2/cm2/d", "nmolO2/cm2/d",
      "nmolP/cm2/d", "nmolP/cm2/d", 
      "nmolC/cm2/d", "nmolC/cm2/d", 
      "nmolC/cm2/d", "nmolC/cm2/d", 
      "nmolC/cm2/d", "nmolC/cm2/d",
      "nmolP/cm2/d", "nmolP/cm2/d", 
      "nmolC/cm2/d", "nmolN/cm2/d", "nmolP/cm2/d",
      "molO/molC", "molO/molC", 
      "molN/molP", "molN/molP", "molN/molP",
      "nmolC/cm2/d", "nmolC/cm2/d", "nmolC/cm2/d", "nmolC/cm2/d",
      "-", "-", "-",
      "nmolN/cm2/d", "nmolN/cm2/d", "nmolN/cm2/d", "nmolO2/cm2/d",
      "nmolP/cm2/d", "nmolP/cm2/d", "nmolP/cm2/d", "nmolP/cm2/d", "nmolP/cm2/d",
      "nmolO2/cm2/d", "nmolN/cm2/d", "-", "-",
      "nmolN/cm2/d", "nmolP/cm2/d",
      "nmolN/cm2/d", "nmolN/cm2/d", "nmolP/cm2/d", "nmolC/cm2/d",
      "nmolO2/cm2/d", "nmolC/cm2", "nmolC/cm2", "nmolO/cm2", 
      "nmolN/cm2", "nmolN/cm2", "nmolN/cm2", 
      "nmolO/cm2", "nmolC/cm2", "nmolP/cm2", "nmolP/cm2",
      "nmolP/cm2", "nmolP/cm2")
      
.CNPDIA$descrip0D <- c("O2 influx sediment-water", "O2 efflux lower boundary", 
      "NO3 influx sediment-water", "NO3 efflux lower boundary", 
      "NO2 influx sediment-water", "NO2 efflux lower boundary", 
      "NH3 influx sediment-water", "NH3 efflux lower boundary", 
      "ODU influx sediment-water", "ODU efflux lower boundary",
      "PO4 influx sediment-water", "PO4 efflux lower boundary", 
      "DIC influx sediment-water", "DIC efflux lower boundary", 
      "FDET flux to sediment", "FDET efflux lower boundary",
      "SDET flux to sediment", "SDET efflux lower boundary", 
      "FeP efflux lower boundary", "CaP efflux lower boundary",
      "OrgC influx to sediment", "OrgN influx to sediment",
      "OrgP influx to sediment", 
      "O2:DIC ratio flux sediment-water", 
      "(O2+ODUflux): (DIC flux) sediment-water",
      "DIN:DIP ratio flux sediment-water", 
      "DIN:DIP mean concentration",
      "DIN:DIP deep concentration", 
      "Vertically integrated Mineralisation", 
      "Vertically integrated oxic Mineralisation", 
      "Vertically integrated Denitrification",
      "Vertically integrated anoxic Mineralisation",
      "Part of mineralisation by oxic min", 
      "Part of mineralisation by denitrification", 
      "Part of mineralisation by anoxic min", 
      "Vertically integrated nitrification step 1 (NH3 oxidation)",
      "Vertically integrated nitrification step 2 (NO2 oxidation)",
      "Vertically integrated anammox",
      "Vertically integrated ODU oxidation",
      "Vertically integrated FeP production", 
      "Vertically integrated CaP production", 
      "Vertically integrated FeP desorption",
      "Vertically integrated CaP dissolution", 
      "Vertically integrated P adsorption",
      "Vertically integrated ODU oxid by BW O2", 
      "Vertically integrated NH3 adsorption",
      "Part P removed", "Part N removed",
      "Vertically integrated NH3 production",
      "Vertically integrated PO4 production",
      "Vertically integrated MPB NO3 uptake",
      "Vertically integrated MPB NH3 uptake",
      "Vertically integrated MPB PO4 uptake",
      "Vertically integrated MPB DIC uptake",
      "Vertically integrated MPB O2 production", 
      "Vertically integrated Fast decaying Detritus",
      "Vertically integrated Slow decaying Detritus",
      "Vertically integrated Oxygen",
      "Vertically integrated Nitrate",
      "Vertically integrated Nitrite",
      "Vertically integrated Ammonium/ammonia",
      "Vertically integrated Oxygen Demand Units",
      "Vertically integrated Dissolved Inorganic Carbon",
      "Vertically integrated Phosphate",
      "Vertically integrated Iron-bound P", 
      "Vertically integrated Ca-bound P", 
      "Vertically integrated Adsorbed P")
      
##------------------------------------
## forcing functions 
##------------------------------------

.CNPDIA$varforc <- c("Cflux", "w", "biotfac", "irrfac",
       "rFast", "rSlow", "pFast", "MPBprod", 
       "CaPflux", "FePflux", "gasflux", 
       "bwO2", "bwNO3", "bwNO2", "bwNH3", "bwODU", "bwPO4", "bwDIC", 
       "Hwater", "ratefac")

.CNPDIA$unitforc <- c("nmolC/cm2/d", "cm/d", "-", "-",
       "/d", "/d", "-", "mmol/m3/d",
       "nmolP/cm2/d", "nmolP/cm2/d", "cm/d",
       rep("mmol/m3", times=7),
       "cm", "-")

.CNPDIA$descripforc <- c( "Carbon flux to sediment", 
      "Sedimentation rate", 
      "Bioturbation multiplication factor",
      "Irrigation multiplication factor", 
      "Decay rate FDET",
      "Decay rate SDET",
      "Part FDET in flux", 
      "MicroPhytoBenthos forcing",
      "CaP deposition flux", 
      "FeP deposition flux", 
      "Gas exchange flux (piston velocity)",
      "Bottom water O2 concentration", 
      "Bottom water NO3 concentration", 
      "Bottom water NO2 concentration", 
      "Bottom water NH3 concentration", 
      "Bottom water ODU concentration", 
      "Bottom water PO4 concentration", 
      "Bottom water DIC concentration", 
      "Height of water above the sediment", 
      "Rate multiplication factor"
      )
      
##------------------------------------
## 1D variables
##------------------------------------
      
.CNPDIA$var1D <- c("TOC", "DICprod", "DINprod", 
      "DIPprod", "O2prod",
      "Oxicmin", "Denitrific", "Anoxicmin",
      "Nitri1", "Nitri2", "Anammox", 
      "ODUox", "ODUoxsurf","ODUdepo",
      "FePadsorp", "FePdesorp",  
      "CaPprod", "CaPdiss", "Padsorp", 
      "MPBproduction", "MPBNO3uptake", "MPBNH3uptake", 
      "MPBPO4uptake", "MPBDICuptake")

.CNPDIA$unit1D <- c("%", "nmolC/cm3 liquid/d", "nmolN/cm3 liquid/d", 
      "nmolP/cm3 liquid/d", "nmolO/cm3 liquid/d",
      "nmolC/cm3 liquid/d", "nmolC/cm3 liquid/d", "nmolC/cm3 liquid/d",
      "nmolN/cm3 liquid/d", "nmolN/cm3 liquid/d", "nmolN/cm3 liquid/d", 
      "nmolO/cm3 liquid/d", "nmolO/cm3 liquid/d", "nmolO/cm3 liquid/d",
      "nmolP/cm3 liquid/d", "nmolP/cm3 solid/d",  
      "nmolP/cm3 liquid/d", "nmolP/cm3 solid/d", "nmolP/cm3 liquid/d", 
      "nmolC/cm3 solid/d", "nmolN/cm3 liquid/d", "nmolN/cm3 liquid/d", 
      "nmolP/cm3 liquid/d", "nmolC/cm3 liquid/d")
      
.CNPDIA$descrip1D <- c("Total Organic Carbon % profile", 
      "DIC production profile (C-mineralisation)", 
      "DIN production profile (N-mineralisation)",  
      "DIP production profile (P-mineralisation)", 
      "O2 production profile (microphytobenthos)", 
      "Oxic mineralisation profile", 
      "Denitrification profile", 
      "Anoxic mineralisation profile",
      "Nitrification step 1 profile (NH3 oxidation)", 
      "Nitrification step 2 profile (NO2 oxidation)", 
      "Anammox profile", 
      "ODU oxidation profile",  "ODU oxidation with surface O2 profile", 
      "ODU deposition profile",
      "FeP adsorption profile", "FeP desorption profile",  
      "CaP production profile", "CaP dissolution profile", 
      "P adsorption profile", 
      "MPB production profile",
      "MPB NO3 uptake profile", "MPB NH3 uptake profile", 
      "MPB PO4 uptake profile", "MPB DIC uptake profile")


      
.CNPDIA$var1Dall <- as.vector(sapply(.CNPDIA$var1D, FUN = function(x) rep(x, times = .CNPDIA$N)))

.CNPDIA$outnames <- c(.CNPDIA$var0D, .CNPDIA$var1Dall,.CNPDIA$varforc)

.CNPDIA$nout     <- length(.CNPDIA$outnames)

# text used for labeling plots
.CNPDIA$labels <- data.frame(
    Units = c("mmolC/m3 solid", "mmolO/m3 liquid", 
        "mmolN/m3 liquid", "mmolC/m3 liquid",
        "mmolP/m3 liquid", "mmolP/m3 solid", 
        "mmolN/m3 solid", "mgChl/m3 solid", "-", "%",
        "nmolC/cm3 liquid/d", "nmolN/cm3 liquid/d", "nmolP/cm3 liquid/d", 
        "nmolO/cm3 liquid/d", "nmolP/cm3 solid/d",  "nmolC/cm3 solid/d", 
        "molC/gChl/d", "ngChl/cm3 solid/d", 
        "nmolN/cm3 solid/d", "gChl/molC",
        "molN/molC", "gChl/molN", "molC/molN", 
        "uEinst/m2/d", "nmolC/cm2/d", 
        "nmolP/cm2/d", "/d", "molP/molC", "mmol/m3", 
        "cm/d", "cm2/d", "cm",                
        "mmolO2/m3", "mmolN/m3", "mmolNO3/m3", 
        "dgC", "psu", "molC/molP",          
        "mmol/m3/d", "/cm", "mmolP/m3", "mmolC/m3", 
        "nmolO2/cm2/d", "nmolN/cm2/d",
        "molO/molC", "molN/molP", "ngChl/cm2/d"),
    Labels = c("mmol/m3.s", "mmol/m3.l", 
        "mmol/m3.l", "mmol/m3.l", 
        "mmol/m3.l", "mmol/m3.s", 
        "mmol/m3.s", "mg/m3.s", "-", "%", 
        "mmol/m3.l/d", "mmol/m3.l/d", "mmol/m3.l/d", 
        "mmol/m3.l/d", "mmol/m3.s/d", "mmol/m3.s/d", 
        "mol/g/d", "mg/m3.s/d", 
        "mmol/m3.s/d",  "g/mol", 
        "mol/mol", "g/mol", "mol/mol", 
        "uEinst/m2/d", "nmol/cm2/d", 
        "nmol/cm2/d", "/d", "mol/mol", "mmol/m3", 
        "cm/d", "cm2/d", "cm", 
        "mmol/m3", "mmol/m3", "mmol/m3", 
        "dgC", "", "mol/mol", 
        "mmol/m3/d", "/cm", "mmol/m3", "mmol/m3", 
        "nmol/cm2/d", "nmol/cm2/d", 
        "mol/mol", "mol/mol", "ng/cm2/d")    )

row.names(.CNPDIA$labels) <- .CNPDIA$labels$Units

.CNPDIA$getplot1D <- 
    rbind(data.frame(names = .CNPDIA$svar, units = .CNPDIA$labels[.CNPDIA$yunits,2]),
          data.frame(names = .CNPDIA$var1D, units = .CNPDIA$labels[.CNPDIA$unit1D,2]))
