## ====================================================================
## A local environment for non user-visible data,
## ====================================================================

.FESDIA <- new.env()

.FESDIA$N     <- 100
.FESDIA$Grid  <- setup.grid.1D(x.up=0, dx.1 = 0.01, N = .FESDIA$N, L = 100)

##------------------------------------
## Parameters
##------------------------------------

.FESDIA$Parms <- c(
  ## organic matter dynamics  #
    Cflux       = 20*1e5/12/365 ,  # nmolC/cm2/d  - Carbon deposition: 20gC/m2/yr
    pFast       = 0.9         , # -            fraction fast detritus in flux
    FeOH3flux   = 1           , # nmol/cm2/d   FeOH3 deposition rate
    CaPflux     = 0           , # nmolP/cm2/d deposition rate of CaP
    rFast       = 25/365      , # /d          decay rate fast decay detritus
    rSlow       = 0.05/365    , # /d          decay rate slow decay detritus
    NCrFdet  = 16/106         , # molN/molC   NC ratio fast decay detritus
    NCrSdet  = 16/106         , # molN/molC   NC ratio slow decay detritus
    PCrFdet  = 1/106          , # molP/molC   PC ratio fast decay det.
    PCrSdet  = 1/106          , # molP/molC   PC ratio slow decay det.

  ## Nutrient bottom water conditions
    BCupLiq   = 2             , # upper boundary condition for liquid; 1= flux, 2=conc, 3 = 0-grad
    BCdownLiq = 3             , # lower boundary condition; default = 0-gradient
    O2bw      = 300           , # mmol/m3    Oxygen conc in bottom water
    NO3bw     = 10            , # mmol/m3
    NO2bw     = 0             , # mmol/m3
    NH3bw     = 1             , # mmol/m3
    CH4bw     = 0             , # mmol/m3
    PO4bw     = 0.5           , # mmol/m3
    DICbw     = 2100          ,
    Febw      = 0             ,
    H2Sbw     = 0             ,
    SO4bw     = 31000         ,
    ALKbw     = 2400          ,

    O2dw      = NA            , # mmol/m3    deep boundary concentrations
    NO3dw     = NA            , # mmol/m3
    NO2dw     = NA            , # mmol/m3
    NH3dw     = NA            , # mmol/m3
    CH4dw     = NA            , # mmol/m3
    PO4dw     = NA            , # mmol/m3
    DICdw     = NA            ,
    Fedw      = NA            ,
    H2Sdw     = NA            ,
    SO4dw     = NA            ,
    ALKdw     = NA            ,

  ## Bioturbation, advection, bio-irrigation
    w         = 1/365000      , # cm/d        advection rate
    biot      = 1/365         , # cm2/d       bioturbation coefficient
    biotdepth = 5             , # cm          depth of mixed layer
    biotatt   = 1             , # cm          depth attenuation coefficient below biotdepth
    irr       = 0             , # /d          bio-irrigation rate
    irrdepth  = 5             , # cm          depth of irrigates layer
    irratt    = 1             , # cm          depth attenuation coefficient below irrdepth
    gasflux   = 0             , # cm/d        piston velocity for dry flats - exchange of O2 and DIC only+deposition 
  
  ## Nutrient parameters
    NH3Ads         = 1.3      , # -           Adsorption coeff ammonium
    rnitri1        = 20       , #/d           Max nitrification rate (step1) - oxidation of ammonium  
    rnitri2        = 20.      ,  #/d          Max nitrification rate (step2) - oxidation of nitrite
    ranammox       = 0.1      , # /(mmol/m3)/d   Anammox rate
    ksO2nitri      = 1.       , # mmolO2/m3   half-sat O2 in nitrification
    ksO2oxic       = 3.       , # mmolO2/m3   half-sat O2 in oxic mineralisation
    ksNO3denit     = 30.      , # mmolNO3/m3  half-sat NO3 in denitrification
    kinO2denit     = 1.       , # mmolO2/m3   half-sat O2 inhib denitrification
    kinNO3anox     = 1.       , # mmolNO3/m3  half-sat NO3 inhib anoxic degr
    kinO2anox      = 0.001    , # mmolO2/m3   half-sat O2 inhib anoxic min

  ## To estimate diffusion coefficients
    temperature    = 10       , # temperature
    salinity       = 35       , # 
    TOC0           = 0.5      , # %           TOC concentration at depth 
    rFePadsorp     = 1e-6     , # /(nmolLiq/cm3)/day
    rCaPprod       = 0.0      ,
    rCaPdiss       = 0.0      ,
    CPrCaP         = 1.32/4.6 , # Ca10(PO4)4.6(CO3)1.32F1.87(OH)1.45   Jilbert and Slomp
    rPads          = 0.0      , # /day -  adsorption rate of phosphate
    rPdes          = 0.0      , # /day - desorption rate of adsorbed P
    maxPads        = 1e3      , # Max adsorbed P concentration - mmolP/m3 solid
    ksFeOH3        = 12500.   , # mmolFeOH3/m3 half-sat FeOH3 in iron red  
    kinFeOH3       = 12500.   , # mmolFeOH3/m3 half-sat FeOH3 inhib BSR
    ksSO4BSR       = 1600.    , # mmolSO4/m3   half-sat SO4 in sulfate reduction
    kinSO4Met      = 1000     , # mmolSO4/m3   half-sat SO4 inhibition for methanogenesis
    rFeox          = 0.3      , # /d/nmol/cm3  oxidation constant for iron by O2 (bimolecular rate law)
    rH2Sox         = 5e-4     , # /d/nmol/cm3  oxidation constant for diss Sulfide by O2 (bimolecular rate law)
    rFeS           = 1e-3     , # /d/nmol/cm3  oxidation constant for diss Sulfide by O2 (bimolecular rate law)
    rCH4ox         = 27       , # /d/nmol/cm3  oxidation constant for CH4 by O2 (bimolecular rate law)
    rAOM           = 3e-5     , # /d/nmol/cm3  oxidation constant for AOM CH4 by SO4 (bimolecular rate law)
    rSurfH2Sox     = 0.       , # /d           Max rate oxidation of deep H2S with surface O2 - 0 to toggle it off
    rSurfCH4ox     = 0.       , # /d           Max rate oxidation of deep CH4 with surface O2 - 0 to toggle it off
    ksSurfALK      = 3000.    , # mmol/m3      half-sat alkalinity in reoxidation of CH4 or H2S with O2
    ksO2reox       = 1.       , # mmolO2/m3    half-sat O2 in reoxidation of CH4 or H2S with O2
    ODUoxdepth     = 5.      , # cm           Depth where oxidation of H2S/ODU with surface water O2 is possible
    ODUoxatt       = 1.      , # /cm          the depth attenuation coefficient of H2S/ODU oxidation below ODUoxdepth.
# NOT USED  rFeS2      = 8.9^10-6,       # cm3 liquid/nmol/d      -  Rickard (1997a)
  
    por0           = 0.9      , # -            surface porosity
    pordeep        = 0.5      , # -            deep porosity
    porcoeff       = 0.3      , # cm           porosity coefficient
    formationtype   = 1       ,    # to estimate effective diffusion, 1=sand, 2=mud, 3=general
  
# parameters if bottom water is dynamically described
    dilution        = 0       , # /day         relaxation towards background conc dissolved
    Hwater          = 10      , # cm           height of water over core
    Cfall           = 100     , # cm/day       fall speed of organic carbon (FDET, SDET)
    FePfall         = 100     , # cm/day       fall speed of FeP
    FeOH3fall       = 100     , # cm/day       fall speed of FeP
    CaPfall         = 100     , # cm/day       fall speed of FeP
    addalk          = 1.      , # if 1: alkalinity dyanmically modeled 

  ## MPB production
    MPBprod   =     0         , # mmol/m3/d   maximal rate - range: 5000-5e4
    kMPB      =     4         , # /cm,        exponential decay 
    kDINupt   =  0.01         , # mmol/m3,    DIN limitation constant
    kPO4upt   = 0.001         , # mmol/m3,    P limitation constant
    kDICupt   =     1           # mmol/m3,    C limitation constant
  
  
    )

# parameter units
.FESDIA$Parunit <- c("nmolC/cm2/d", "-", "nmol/cm2/d", "nmol/cm2/d", "/d", "/d",      
    "molN/molC", "molN/molC", "molP/molC", "molP/molC", "-", "-", 
    rep("mmol/m3", times=11), rep("mmol/m3",times=11), "cm/d",
    "cm2/d", "cm", "/cm", "/d", "cm", "cm", "cm/d", 
    "-", "/d","/d", "/(mmol/m3)/d", "mmolO2/m3", "mmolO2/m3",
    "mmolNO3/m3", "mmolO2/m3", "mmolNO3/m3", "mmolO2/m3",
    "dgC", "psu", "%",
    "/d", "/d", "/d", "mol/mol", "/d", "/d", "mmolP/m3solid", "mmolFeOH3/m3","mmolFeOH3/m3",
    "mmolS/m3", "mmolS/m3", "/(mmol/m3)/d", "/(mmol/m3)/d", "/(mmol/m3)/d",
    "/(mmol/m3)/d", "/(mmol/m3)/d", "/d", "/d", "mmol/m3", "mmolO2/m3", "cm", "/cm", "-", "-", "cm", 
    "-", "/d", "cm", "cm/d", "cm/d", "cm/d", "cm/d", "-", "mmol/m3/d", "/cm", "mmol/m3",
     "mmol/m3", "mmol/m3")


.FESDIA$Pardesc <- c("total organic C deposition","part FDET in carbon flux", 
  "deposition rate of FeOH3", "deposition rate of CaP",
  "decay rate FDET", "decay rate SDET", "NC ratio FDET", "NC ratio SDET",
  "PC ratio FDET", "PC ratio SDET",
  "upper boundary liq. 1:flux, 2:conc, 3:0-grad",
  "lower boundary liq. 1:flux, 2:conc, 3:0-grad",
  "upper boundary O2  -if BC=1: flux, 2:conc",
  "upper boundary NO3 -if BC=1: flux, 2:conc",
  "upper boundary NO2 -if BC=1: flux, 2:conc",
  "upper boundary NH3 -if BC=1: flux, 2:conc",
  "upper boundary CH4 -if BC=1: flux, 2:conc",
  "upper boundary PO4 -if BC=1: flux, 2:conc",
  "upper boundary DIC -if BC=1: flux, 2:conc",
  "upper boundary Fe2+ -if BC=1: flux, 2:conc",
  "upper boundary H2S -if BC=1: flux, 2:conc",
  "upper boundary SO4 -if BC=1: flux, 2:conc",
  "upper boundary alkalinity -if BC=1: flux, 2:conc",
  "lower boundary O2  -if BC=1: flux, 2:conc",
  "lower boundary NO3 -if BC=1: flux, 2:conc",
  "lower boundary NO2 -if BC=1: flux, 2:conc",
  "lower boundary NH3 -if BC=1: flux, 2:conc",
  "lower boundary CH3 -if BC=1: flux, 2:conc",
  "lower boundary PO4 -if BC=1: flux, 2:conc",
  "lower boundary DIC -if BC=1: flux, 2:conc",
  "lower boundary Fe2+ -if BC=1: flux, 2:conc",
  "lower boundary H2S -if BC=1: flux, 2:conc",
  "lower boundary SO4 -if BC=1: flux, 2:conc",
  "lower boundary alkalinity -if BC=1: flux, 2:conc",
  
  "advection rate", "bioturbation coefficient", "depth of mixed layer", 
  "attenuation coeff below biotdepth", "bio-irrigation rate", 
  "depth of irrigated layer", "attenuation coeff below irrdepth", 
  "piston velocity for dry flats", "Adsorption coeff ammonium",
  "Max nitrification rate step1 (NH3ox)", "Max nitrification rate step2 (NO2ox)",
  "Anammox rate", "half-sat O2 in nitrification", "half-sat O2 in oxic mineralisation",
  "half-sat NO3 in denitrification", "half-sat O2 inhib denitrification",
  "half-sat NO3 inhib anoxic degr", "half-sat O2 inhib anoxic min",
  "temperature", "salinity", "refractory Carbon conc",
  "rate FeP adsorption",
  "rate CaP production", "rate CaP dissolution","C:Pratio in CaP",
  "adsorption rate PO4", "desorption rate of adsorbed P", "Max adsorbed P concentration", 
  
  "half-sat FeOH3 conc in iron reduction", "half-sat FeOH3 inhibition S reduction",
  "half-sat SO4 conc in sulphate reduction", "half-sat SO4 inhibition methanogenesis",
  "Max rate Fe oxidation", "Max rate H2S oxidation", "maximum rate FeS production",
  "Max rate CH4 oxidation with O2", "Max rate anaerobic oxidation Methane", 
  "Max rate H2S oxidation with BW O2", 
  "Max rate CH4 oxidation with BW O2", 
  "half-sat Alkalinity in oxidation of H2S/CH4 with bwO2", 
  "half-sat Oxygen in oxidation of H2S/CH4 with bwO2", 
  "Max depth H2S/CH4 oxidation with BW O2", "depth attenuation ODU oxidation",
  
  "surface porosity", "deep porosity", "porosity decay coefficient",
  "formationfactor, 1=sand,2=fine sand,3=general",
  "relaxation towards background conc ", "height of water over core",
  "fall speed of organic C (FDET, SDET)", "fall speed of FeP", "fall speed of FeOH3", 
  "fall speed of CaP", "solve for alkalinity", 
  "maximal MPB production rate", "sedimentary light extinction coefficient",
  "DIN limitation constant MPB", 
  "P limitation constant MPB", 
  "C limitation constant MPB")

##------------------------------------
## State variables
##------------------------------------

.FESDIA$ynames <- c("FDET", "SDET", "O2",
                    "NO3", "NO2", "NH3", "DIC",
                     "Fe", "FeOH3", "H2S",
                     "SO4", "CH4", "PO4",
                     "FeP", "CaP", "Pads", "ALK")
.FESDIA$svar <- .FESDIA$ynames 
.FESDIA$yunits <- c("mmolC/m3 solid", "mmolC/m3 solid", "mmolO/m3 liquid",
                    "mmolN/m3 liquid", "mmolN/m3 liquid", "mmolN/m3 liquid", "mmolC/m3 liquid", 
                    "mmolFe/m3 liquid", "mmolFe/m3 solid", "mmolS/m3 liquid",
                    "mmolS/m3 liquid", "mmolC/m3 liquid", "mmolP/m3 liquid",
                    "mmolP/m3 solid", "mmolP/m3 solid", "mmolP/m3 solid", "mmol/m3 liquid")
.FESDIA$ydescrip <- c("Fast decaying Detritus (solid)", 
                      "Slow decaying Detritus (solid)",
                      "Oxygen (liquid)", 
                      "Nitrate (liquid)", "Nitrite (liquid)", 
                      "Ammonium/ammonia (liquid)", 
                      "Dissolved Inorganic Carbon (liquid)",
                      "Fe2+ (liquid)",
                      "Fe-oxide (solid)",
                      "Sulphide (liquid)",
                      "Sulphate (liquid)",
                      "Methane (liquid)",
                      "Phosphate (liquid)",
                      "Iron-bound P (solid)", 
                      "Ca-bound P (solid)", 
                      "Adsorbed P (solid)",
                      "Alpkalinity (liquid)")

.FESDIA$ynamesall <- as.vector(sapply(.FESDIA$ynames, FUN = function(x) rep(x, times = .FESDIA$N)))


##------------------------------------
## 0-D Variables
##------------------------------------

.FESDIA$var0D <- c("O2flux", "O2deepflux",
    "NO3flux", "NO3deepflux",
    "NO2flux", "NO2deepflux",
    "NH3flux", "NH3deepflux",
    "PO4flux", "PO4deepflux",
    "DICflux", "DICdeepflux",
    "Feflux", "Fedeepflux",
    "H2Sflux", "H2Sdeepflux",
    "SO4flux", "SO4deepflux",
    "CH4flux", "CH4deepflux",
    "ALKflux", "ALKdeepflux",
    "FDETflux", "FDETdeepflux",
    "SDETflux", "SDETdeepflux",
    "FePsurfflux", "FePdeepflux",
    "CaPsurfflux", "CaPdeepflux",
    "FeOH3surfflux", "FeOH3deepflux",
    "OrgCflux", "OrgNflux", "OrgPflux",
    "DINDIPflux", "DINDIPmean", "DINDIPdeep",
    "TotMin", "TotOxic", "TotDenit",
    "TotFered", "TotBSR", "TotMeth",
    "PartOxic", "PartDenit", 
    "PartFered", "PartBSR", "PartMethano",
    "TotNitri1", "TotNitri2", "TotAnammox", "TotFeoxid",
    "TotH2Soxid", "TotCH4oxid", "TotAOM",
    "TotFeSprod", "TotFePprod", "TotCaPprod", 
    "TotFePdesorp", "TotCaPdiss", "TotPadsorb",
    "TotNH3prod", "TotPO4prod", "TotNH3ads", "TotO2prod", "TotH2Soxsurf",
    "TotCH4oxsurf", "TotALkprod", "PartPremoved", "PartNremoved",
    "TotMPBNO3uptake", "TotMPBNH3uptake", "TotMPBPO4uptake", "TotMPBDICuptake",
    "TotMPBO2prod","TotalFDET", "TotalSDET", "TotalO2", 
    "TotalNO3", "TotalNO2", "TotalNH3", 
    "TotalDIC", "TotalFe", "TotalFeOH3", "TotalH2S", "TotalSO4", 
    "TotalCH4", "TotalPO4", "TotalFeP", "TotalCaP", "TotalPads")

.FESDIA$unit0D <- c("nmolO2/cm2/d", "nmolO2/cm2/d", "nmolN/cm2/d", "nmolN/cm2/d", 
      "nmolN/cm2/d", "nmolN/cm2/d", "nmolN/cm2/d", "nmolN/cm2/d", 
      "nmolP/cm2/d", "nmolP/cm2/d", 
      "nmolC/cm2/d", "nmolC/cm2/d", "nmolFe/cm2/d", "nmolFe/cm2/d",
      "nmolS/cm2/d", "nmolS/cm2/d", "nmolS/cm2/d", "nmolS/cm2/d",
      "nmolC/cm2/d", "nmolC/cm2/d", "nmol/cm2/d",  "nmol/cm2/d",
      "nmolC/cm2/d", "nmolC/cm2/d", "nmolC/cm2/d", "nmolC/cm2/d",
      "nmolP/cm2/d", "nmolP/cm2/d", "nmolP/cm2/d", "nmolP/cm2/d",
      "nmolFe/cm2/d", "nmolFe/cm2/d", "nmolC/cm2/d", "nmolN/cm2/d", 
      "nmolP/cm2/d", "molN/molP", "molN/molP", "molN/molP",
      "nmolC/cm2/d", "nmolC/cm2/d", "nmolC/cm2/d", 
      "nmolC/cm2/d", "nmolC/cm2/d", "nmolC/cm2/d",
      "-", "-", "-", "-", "-", 
      "nmolN/cm2/d", "nmolN/cm2/d", "nmolN/cm2/d", "nmolFe/cm2/d", 
      "nmolS/cm2/d", "nmolC/cm2/d", "nmolS/cm2/d",
      "nmolFe/cm2/d", "nmolFe/cm2/d", "nmolP/cm2/d", 
      "nmolP/cm2/d", "nmolP/cm2/d", "nmolP/cm2/d",
      "nmolN/cm2/d", "nmolP/cm2/d", "nmolN/cm2/d", "nmolO/cm2/d",
      "nmolS/cm2/d", "nmolC/cm2/d", "nmol/cm2/d", "-", "-",
      "nmolN/cm2/d", "nmolN/cm2/d", "nmolP/cm2/d", "nmolC/cm2/d",
      "nmolO2/cm2/d", "nmolC/cm2", "nmolC/cm2",   
       "nmolO/cm2",  "nmolN/cm2",  "nmolN/cm2",  "nmolN/cm2",  "nmolC/cm2", 
      "nmolFe/cm2",  "nmolFe/cm2",  "nmolS/cm2",  "nmolS/cm2",  
       "nmolC/cm2",  "nmolP/cm2",  "nmolP/cm2",  
       "nmolP/cm2",  "nmolP/cm2") 
 
.FESDIA$descrip0D <- c("O2 influx sediment-water", "O2 efflux lower boundary", 
      "NO3 influx sediment-water", "NO3 efflux lower boundary", 
      "NO2 influx sediment-water", "NO2 efflux lower boundary", 
      "NH3 influx sediment-water", "NH3 efflux lower boundary", 
      "PO4 influx sediment-water", "PO4 efflux lower boundary", 
      "DIC influx sediment-water", "DIC efflux lower boundary", 
      "Fe2+ influx sediment-water", "Fe2+ efflux lower boundary", 
      "H2S influx sediment-water", "H2S efflux lower boundary", 
      "SO4 influx sediment-water", "SO4 efflux lower boundary", 
      "CH4 influx sediment-water", "CH4 efflux lower boundary", 
      "Alkalinity influx sediment-water", "Alkalinity efflux lower boundary", 
      "FDET flux to sediment", "FDET efflux lower boundary",
      "SDET flux to sediment", "SDET efflux lower boundary", 
      "FeP flux upper boundary", "FeP efflux lower boundary",
      "CaP flux upper boundary", "CaP efflux lower boundary",
      "FeOH3 flux upper boundary", "FeOH3 efflux lower boundary",
      "OrgC influx to sediment", "OrgN influx to sediment",
      "OrgP influx to sediment", 
#      "O2:DIC ratio flux sediment-water", 
#      "(O2+ODUflux): (DIC flux) sediment-water",
      "DIN:DIP ratio flux sediment-water", 
      "DIN:DIP mean concentration",
      "DIN:DIP deep concentration", 
      "Vertically integrated Mineralisation", 
      "Vertically integrated oxic Mineralisation", 
      "Vertically integrated Denitrification",
      "Vertically integrated Iron reduction",
      "Vertically integrated Sulphate reduction",
      "Vertically integrated Methanogenesis",
      "Part of mineralisation by oxic min", 
      "Part of mineralisation by denitrification", 
      "Part of mineralisation by iron reduction", 
      "Part of mineralisation by sulphate reduction", 
      "Part of mineralisation by methanogenisis", 
      "Vertically integrated nitrification step 1 (NH3 ox)",
      "Vertically integrated nitrification step 2 (NO2 ox)",
      "Vertically integrated anammox",
      "Vertically integrated Fe2+ oxidation",
      "Vertically integrated H2S oxidation",
      "Vertically integrated CH4 oxidation",
      "Vertically integrated Anaerobic oxidation methane", 
      "Vertically integrated FeS production", 
      "Vertically integrated FeP production", 
      "Vertically integrated CaP production", 
      "Vertically integrated FeP desorption",
      "Vertically integrated CaP dissolution", 
      "Vertically integrated P adsorption",
      "Vertically integrated NH3 production",
      "Vertically integrated PO4 production",
      "Vertically integrated NH3 adsorption",
      "Vertically integrated O2 production (?)",
      "Vertically integrated H2S oxidation by surface O2",
      "Vertically integrated CH4 oxidation by surface O2",
      "Total alkalinity production",
      "Part P removed", "Part N removed",
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
      "Vertically integrated Dissolved Inorganic Carbon",
      "Vertically integrated Fe",
      "Vertically integrated FeOH3",
      "Vertically integrated H2S",
      "Vertically integrated SO4",
      "Vertically integrated CH4",
      "Vertically integrated Phosphate",
      "Vertically integrated Iron-bound P", 
      "Vertically integrated Ca-bound P", 
      "Vertically integrated Adsorbed P")

##------------------------------------
## forcing functions 
##------------------------------------

.FESDIA$varforc <- c("Cflux", "FeOH3flux", "CaPflux", "w", "biotfac", "irrfac",
                     "rFast", "rSlow", "pFast", "MPBprod", 
                      "gasflux", "bwO2", "bwNO3", "bwNO2", "bwNH3", 
                     "bwCH4", "bwFe", "bwH2S", "bwSO4",
                     "bwPO4", "bwDIC", "bwALK", "Hwater", "Ratefactor")

.FESDIA$unitforc <- c("nmolC/cm2/d", "nmolP/cm2/d", "nmolP/cm2/d",
                     "cm/d", "-", "-", "/d", "/d",
                     "-", "mmol/m3/d", "cm/d", rep("mmol/m3", times=11), "cm", "-")

.FESDIA$descripforc <- c( "Carbon flux to sediment", 
     "FeOH3 flux to sediment", "CaP flux to sediment", 
      "Sedimentation rate", 
      "Bioturbation multiplication factor",
      "Irrigation multiplication factor", 
      "Decay rate FDET",
      "Decay rate SDET",
      "Part FDET in flux", 
      "MicroPhytoBenthos forcing",
      "Gas exchange flux (piston velocity)",
      "Bottom water O2 concentration", 
      "Bottom water NO3 concentration", 
      "Bottom water NO2 concentration", 
      "Bottom water NH3 concentration", 
      "Bottom water CH4 concentration", 
      "Bottom water Fe concentration", 
      "Bottom water H2S concentration", 
      "Bottom water SO4 concentration", 
      "Bottom water PO4 concentration", 
      "Bottom water DIC concentration", 
      "Bottom water alkalinity concentration", 
      "Height of water above the sediment", 
      "Rate multiplication factor"
      )
      
##------------------------------------
## 1D variables
##------------------------------------
      
.FESDIA$var1D <- c("TOC", "DICprodMin", "DINprodMin", "DIPprodMin", "O2prod", 
                   "Oxicmin", "Denitrific", "Feredmin", 
                   "BSRmin", "Methmin", "nitri1", "nitri2", "Anammox", 
                   "Feoxid", "H2Soxid", "CH4oxid",
                   "AOM", "FeSprod", "FePadsorp", "FePdesorp", 
                   "CaPprod", "CaPdiss","Padsorb",
                   "H2Soxsurf", "CH4oxsurf",
                   "O2distConsump", "ALKprod", "DICprodCH4",
                   "MPBCprod", "MPBuptakeNO3", "MPBuptakeNH3", 
                   "MPBuptakePO4", "MPBuptakeDIC")                                          

.FESDIA$unit1D <- c("%", "nmolC/cm3 liquid/d", "nmolN/cm3 liquid/d", 
      "nmolP/cm3 liquid/d", "nmolO/cm3 liquid/d", 
      "nmolC/cm3 liquid/d", "nmolC/cm3 liquid/d", "nmolC/cm3 liquid/d", 
      "nmolC/cm3 liquid/d", "nmolC/cm3 liquid/d",
      "nmolN/cm3 liquid/d", "nmolN/cm3 liquid/d", "nmolN/cm3 liquid/d", 
      "nmolFe/cm3 liquid/d", "nmolS/cm3 liquid/d", 
      "nmolC/cm3 liquid/d", "nmolS/cm3 liquid/d", "nmolFe/cm3 liquid/d",      
      "nmolFe/cm3 liquid/d", "nmolP/cm3 solid/d",  "nmolP/cm3 liquid/d", 
      "nmolP/cm3 solid/d", "nmolP/cm3 solid/d",
      "nmolS/cm3 liquid/d",  "nmolC/cm3 liquid/d",  "nmolO/cm3 liquid/d",  
      "nmol/cm3 liquid/d", "nmolC/cm3 liquid/d", "nmolC/cm3 solid/d", 
      "nmolN/cm3 liquid/d", "nmolN/cm3 liquid/d", 
      "nmolP/cm3 liquid/d", "nmolC/cm3 liquid/d")
     
.FESDIA$descrip1D <- c("Total Organic Carbon % profile", 
      "DIC production profile (mineralisation)", 
      "DIN production profile (mineralisation)",  
      "DIP production profile (mineralisation)", 
      "O2 production profile (microphytobenthos)", 
      "Oxic mineralisation profile", "Denitrification profile", 
      "Fe reduction mineralisation profile", "Sulphate reduction mineralisation profile",
      "Methanogensis mineralisation profile",
      "Nitrification step 1 profile (NH3 oxidation)", 
      "Nitrification step 2 profile (NO2 oxidation)", 
      "Anammox profile", "Fe2+ oxidation profile", "H2S oxidaton profile",
      "CH4 oxidation profile",  
      "Anaerobic oxidation of methane profile", "FeS production profile",
      "FeP adsorption profile", "FeP desorption profile",  
      "CaP production profile", "CaP dissolution profile", 
      "P adsorption profile", 
      "H2S oxidation with surface O2 profile", 
      "CH4 oxidation with surface O2 profile", 
      "O2 uptake oxidation with surface O2 profile", 
      "Alkalinity production profile", "DIC production via Methane profile",
      "MPB production profile",
      "MPB NO3 uptake profile", "MPB NH3 uptake profile", 
      "MPB PO4 uptake profile", "MPB DIC uptake profile")

  
.FESDIA$var1Dall <- as.vector(sapply(.FESDIA$var1D, FUN = function(x) rep(x, times = .FESDIA$N)))

.FESDIA$outnames <- c(.FESDIA$var0D, .FESDIA$var1Dall,.FESDIA$varforc)

.FESDIA$nout <- length(.FESDIA$outnames)

# how to plot the units
.FESDIA$labels <- data.frame(
    Units = c("%", "nmolC/cm3 liquid/d", "nmolN/cm3 liquid/d", 
       "nmolP/cm3 liquid/d", "nmolO/cm3 liquid/d", "nmolFe/cm3 liquid/d",
       "nmolS/cm3 liquid/d", "nmolP/cm3 solid/d", "nmolC/cm3 solid/d",
       "nmol/cm3 liquid/d",
       "nmolO2/cm2/d", "nmolN/cm2/d", "nmolP/cm2/d", "nmolC/cm2/d",
       "nmolFe/cm2/d", "nmolS/cm2/d", "nmol/cm2/d", "molN/molP", "-",                  
       "nmolO/cm2/d" , "cm/d", "/d", "mmol/m3/d", "mmol/m3", 
       "cm", "mmolC/m3 solid",                 
       "mmolO/m3 liquid", "mmolN/m3 liquid", "mmolC/m3 liquid", 
       "mmolFe/m3 liquid", "mmolFe/m3 solid",  
       "mmolS/m3 liquid", "mmolP/m3 liquid", "mmolP/m3 solid", "mmol/m3 liquid",
       "nmolC/cm2", "nmolO/cm2", "nmolN/cm2", "nmolFe/cm2", "nmolS/cm2"),
    Labels = c("%", "mmol/m3.l/d", "mmol/m3.l/d", 
       "mmol/m3.l/d", "mmol/m3.l/d", "mmol/m3.l/d",
       "mmol/m3.l/d", "mmol/m3.s/d", "mmol/m3.s/d", "mmol/m3.l/d",
       "nmol/cm2/d", "nmol/cm2/d", "nmol/cm2/d", "nmol/cm2/d",
       "nmol/cm2/d", "nmol/cm2/d", "nmol/cm2/d", "mol/mol", "-",                  
       "nmol/cm2/d" , "cm/d", "/d", "mmol/m3/d", "mmol/m3", "cm", "mmol/m3.s",                 
       "mmol/m3.l", "mmol/m3.l", "mmol/m3.l", "mmol/m3.l", "mmol/m3.s",  
       "mmol/m3.l", "mmol/m3.l", "mmol/m3.s", "mmol/m3.l",
       "nmol/cm2", "nmol/cm2", "nmol/cm2", "nmol/cm2", "nmol/cm2")  )

row.names(.FESDIA$labels) <- .FESDIA$labels$Units

.FESDIA$getplot1D <- 
    rbind(data.frame(names = .FESDIA$ynames, units = .FESDIA$labels[.FESDIA$yunits,2]),
          data.frame(names = .FESDIA$var1D, units = .FESDIA$labels[.FESDIA$unit1D,2]))
