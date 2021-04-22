!==========================================================================
! THE OMEXDIA model with P, Fe and S, implemented in FORTRAN
!
! This file contains the module declarations and common subroutines
! Karline Soetaert, nioz-yerseke
!==========================================================================

!==========================================================================
! MODULES
!========================================================================== 

      MODULE dimFESDIA
      IMPLICIT NONE
      INTEGER, PARAMETER :: N = 100, Np1 = N+1
      INTEGER, PARAMETER :: nparmsdia = 5*N + 5*(N+1) + 69
      INTEGER, PARAMETER :: nforcsdia = 24
      INTEGER, PARAMETER :: nparmspH = 14
      INTEGER, PARAMETER :: nforcspH = 3 

      INTEGER, PARAMETER :: nparms = nparmsdia + nparmspH
      INTEGER, PARAMETER :: nforc  = nforcsdia + nforcspH
      INTEGER, PARAMETER :: Noutdia = 2758+634 

      LOGICAL ::  DynamicpH 

      END MODULE dimFESDIA

!==========================================================================

      MODULE commonFESDIA
      use dimFESDIA
      IMPLICIT NONE

! -----------------------------------
! State variables and derivatives
! -----------------------------------

      DOUBLE PRECISION  :: Fdet(N),Sdet(N),O2(N),NO3(N),NO2(N),NH3(N)
      DOUBLE PRECISION  :: DIC(N),Fe(N),FeOH3(N),H2S(N),SO4(N),CH4(N)
      DOUBLE PRECISION  :: PO4(N),FeP(N),CaP(N),Pads(N),ALK(N)

      DOUBLE PRECISION  :: dFdet(N),dSdet(N),dO2(N),dNO3(N),dNO2(N),            &
     &              dNH3(N),dDIC(N),dFe(N),dFeOH3(N),dH2S(N),dSO4(N),           &     
     &              dCH4(N),dPO4(N),dFeP(N),dCaP(N),dPads(N),dALK(N)

! -----------------------------------
! parameters
! -----------------------------------

! Parameter profiles
! 
      DOUBLE PRECISION :: intpor(N+1), por(N) ! - porosity (at interface, middle),
      DOUBLE PRECISION :: porfac(N+1)         !        porosity multiplication factor (for diffusion)
      DOUBLE PRECISION :: Db0(N+1)            ! cm2/d  bioturbation profile (interface)
      DOUBLE PRECISION :: dx(N), dxInt(N+1)   ! cm     thickness and distance from mid to mid
      DOUBLE PRECISION :: x(N)                ! cm     depth below surface
      DOUBLE PRECISION :: Aint(N+1)           ! cm2    surface area at interface
      DOUBLE PRECISION :: Dirr0(N)            ! /d     irrigation profile
      DOUBLE PRECISION :: distreact(N)        ! -      profile for distance reoxidation with O2

! N:C and P:C ratios of organic matter 
      DOUBLE PRECISION  :: NCrFdet, NCrSdet, PCrFDET, PCrSDET

! boundary conditions for solutes
      DOUBLE PRECISION  :: BCupLiq, BCdownLiq

! deep water concentrations for solutes
      DOUBLE PRECISION  :: dwO2, dwNO3, dwNO2, dwNH3, dwCH4, dwPO4,             &
     &  dwDIC, dwFe, dwH2S, dwSO4, dwALK

! rate coefficients
      DOUBLE PRECISION ::  NH3Ads, rnit1, rnit2,                                &
     &  ranammox, ksO2nitri, ksO2oxic, ksNO3denit, kinO2denit, rPads,           &
     &  rPdes, maxPads, kinNO3anox, kinO2anox, TOC0, rFePadsorp,                &
     &  rCaPprod, rCaPdiss, CPrCaP, ksFeOH3, kinFeOH3, ksSO4,                   &
     &  kinSO4, rFeOx, rH2Sox, rFeS, rCH4ox, rAOM, rH2Soxsurf,                  &
     &  rCH4oxsurf, ksALK, ksO2reox, relax, Cfall, FePfall, FeOH3fall,          &
     &  CaPfall

! diffusion coefficients     
      DOUBLE PRECISION :: DispO2, DispNO3, DispNO2, DispNH3, DispPO4,           &
     &  DispCH4, DispDIC, DispFe, DispH2S, dispSO4, DispAlk, addalk

! (simple) MPB dynamics     
      DOUBLE PRECISION :: kMPB, ksDIN, ksPO4, ksDIC

! Parameters relating to pH dynamics
      DOUBLE PRECISION :: rCaDiss    ! /d   Carbonate dissolution rate
      DOUBLE PRECISION :: pCa        ! -             Carbonate dissolution/precipitation power
      DOUBLE PRECISION :: rArDiss    ! /nmol/cm3/d   Aragonite dissolution/precipitation rate
      DOUBLE PRECISION :: pAr        ! -             Aragonite dissolution/precipitation power
      DOUBLE PRECISION :: rCalc      ! d   Carbonate precipitation rate
      DOUBLE PRECISION :: dwCa       ! mmol/m3       deep Ca2+ concentration
      DOUBLE PRECISION :: CaCO3fall  ! cm/day        fall speed of suspended CaCO3
      DOUBLE PRECISION :: ARAGfall   ! cm/day        fall speed of suspended aragonite
      DOUBLE PRECISION :: DispHCO3, DispCO3, DispCa  ! cm2/d  diffusion coefficients
      DOUBLE PRECISION :: dens       ! kg/m3         seawater density
      DOUBLE PRECISION :: temp, sal  ! dgC, -        temperature, salinity

! Parameter common block: a combination of FESDIA and pHDIA parameters

      COMMON /myparmsFES  /NCrFdet,NCrSdet,PCrFDET,PCrSDET,BCupLiq,             &
     &  BCdownLiq,dwO2,dwNO3,dwNO2,dwNH3,dwCH4,dwPO4,dwDIC,dwFe,                &
     &  dwH2S,dwSO4,dwALK,NH3Ads,rnit1, rnit2,ranammox,                         &
     &  ksO2nitri,ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,                    &
     &  kinO2anox,TOC0,rFePadsorp,rCaPprod,rCaPdiss,CPrCaP,rPads,               &
     &  rPdes, maxPads, ksFeOH3, kinFeOH3, ksSO4, kinSO4, rFeOx,                &
     &  rH2Sox, rFeS, rCH4ox, rAOM, rH2Soxsurf, rCH4oxsurf, ksALK,              &
     &  ksO2reox,relax, Cfall, FePfall, FeOH3fall, CaPfall, addalk,             &
     &  kMPB, ksDIN, ksPO4, ksDIC, rCaDiss, pCa,                                &
     &  rArDiss, pAr, rCalc, dwCa, CaCO3fall, Aragfall,                         &
     &  DispO2, DispNO3, DispNO2, DispNH3, DispPO4, DispCH4, DispDIC,           &
     &  DispFe, DispH2S, dispSO4, DispAlk, DispHCO3, DispCO3, DispCa,           &
     &  dens, temp, Sal, dx, dxint, x, Aint, por, intpor, porfac,               &
     &  Db0, Dirr0, distreact

! -----------------------------------
! forcings
! -----------------------------------

      DOUBLE PRECISION CarbonFlux,BWO2,bwNO3,bwNO2,bwNH3,bwCH4,bwFe,            &
     &   bwH2S,bwSO4,bwPO4,bwDIC,bwALK,w,biotfac,irrfac,rFast,rSlow,            &
     &   pFast,FeOH3flux,CaPflux, gasflux,MPBforc,Hwater, ratefac 

! Forcings relating to pH dynamics
      DOUBLE PRECISION :: CaCO3flux  ! nmolC/cm2/d   CaCO3     deposition: 0.39 M/m2/yr
      DOUBLE PRECISION :: ARflux     ! nmolC/cm2/d   Aragonite deposition: 0.21 M/m2/yr
      DOUBLE PRECISION :: bwCa       ! mmol/m3       Ca++ conc in bottom water

      COMMON /myforcsFES/ CarbonFlux, FeOH3flux, CaPflux, w, biotfac,           &
     &   irrfac, rFast, rSlow, pFast, gasflux, MPBforc, bwO2, bwNO3,            &
     &   BWno2, bwNH3, bwCH4, bwFe, bwH2S, bwSO4, bwPO4, bwDIC,                 &
     &   bwALK, Hwater, ratefac, CaCO3flux, ARflux, bwCa

! -----------------------------------
! output variables
! -----------------------------------

      DOUBLE PRECISION :: O2flux,O2deepflux,NO3flux,NO3deepflux,                &
     &  NO2flux, NO2deepflux, NH3flux,  NH3deepflux,                            &
     &  PO4flux,PO4deepflux,DICflux,DICdeepflux,Feflux,                         &
     &  Fedeepflux,H2Sflux,H2Sdeepflux,SO4flux,SO4deepflux,CH4flux,             &
     &  CH4deepflux,partOxic,partDenit,partFered,partBSR,partMeth,              &
     &  TotMin,TotOxic,TotDenit,TotFered,TotBSR,TotMeth,NPmean,NPflux,          &
     &  NPdeep,TotFePprod,TotCaPprod,Premoved,totO2prod,Nremoved,               &
     &  Cflux,Nflux,Pflux,TotH2Soxsurf,TotCH4oxsurf,TotALkprod,                 &
     &  ALKflux,ALKdeepflux

      DOUBLE PRECISION :: sumPO4,sumDIN,FDETdeepflux,SDETdeepflux,              & 
     &  FePsurfflux,CaPsurfflux,FeOH3surfflux,FePdeepflux,CaPdeepflux,          &
     &  FeOH3deepflux,totNitri1,totNitri2,totAnammox,totFeoxid,                 &
     &  totH2Soxid,totCH4oxid,totAOM,totNH3ads,totFeSprod,totFePdesorp,         &
     &  totCaPdiss,totPadsorp,TotNprod,TotPprod,FDETflux,SDETflux,              &
     &  TotMPBNO3uptake, TotMPBNH3uptake, TotMPBPO4uptake,                      &
     &  TotMPBDICuptake, TotMPBO2prod, TotFDET, TotSDET, TotO2,                 &
     &  TotNO3, TotNO2, TotNH3, TotDIC, TotFe, TotFeOH3, TotH2S,                &
     &  TotSO4, TotCH4, TotPO4, TotFeP, TotCaP, TotPads 
      
      DOUBLE PRECISION  :: Cprod(N), Nprod(N), Pprod(N), O2prod(N),             &
     &  TOC(N), Oxicminlim(N), Denitrilim(N), Feredlim(N), BSRlim(N),           &
     &  Methlim(N), Oxicmin(N), Denitrific(N), nitri1(N), nitri2(N),            &
     &  Anammox(N), FeSprod(N), FePadsorp(N), NetadsorpP(N),                    &
     &  FePdesorp(N), CaPprod(N), CaPdiss(N), Alkprod(N), DICprodCH4(N),        &
     &  FeredMin(N), BSRMin(N), MethMin(N), FeOxid(N), H2SOxid(N),              &
     &  CH4Oxid(N), AOM(N), H2Soxsurf(N), CH4oxsurf(N),O2distConsump(N),        &
     &  MPBproduction(N), NO3uptake(N), NH3uptake(N), PO4uptake(N),             &
     &  DICuptake(N)

      COMMON /myoutFES/O2flux, O2deepflux, NO3flux, NO3deepflux,                &
     &  NO2flux, NO2deepflux, NH3flux, NH3deepflux,                             &
     &  PO4flux, PO4deepflux, DICflux, DICdeepflux, Feflux,                     &
     &  Fedeepflux, H2Sflux, H2Sdeepflux, SO4flux, SO4deepflux,                 &
     &  CH4flux, CH4deepflux, ALKflux, ALKdeepflux, FDETflux,                   &
     &  FDETdeepflux, SDETflux, SDETdeepflux, FePsurfflux, FePdeepflux,         &
     &  CaPsurfflux, CaPdeepflux, FeOH3surfflux, FeOH3deepflux, Cflux,          &
     &  Nflux, Pflux, NPflux, NPmean, NPdeep, TotMin, TotOxic,                  & 
     &  TotDenit, TotFered, TotBSR, TotMeth, partOxic, partDenit,               &
     &  partFered, partBSR, partMeth, totNitri1, totNitri2, totAnammox,         &
     &  totFeoxid, totH2Soxid, totCH4oxid, totAOM, totFeSprod,                  &
     &  TotFePprod, TotCaPprod, totFePdesorp, totCaPdiss, totPadsorp,           &
     &  TotNprod, TotPprod, totNH3ads, totO2prod, TotH2Soxsurf,                 &
     &  TotCH4oxsurf, TotALkprod,Premoved,Nremoved,TotMPBNO3uptake,             &
     &  TotMPBNH3uptake, TotMPBPO4uptake, TotMPBDICuptake, TotMPBO2prod,        &
     &  TotFDET, TotSDET, TotO2, TotNO3, TotNO2, TotNH3, TotDIC, TotFe,         &
     &  TotFeOH3, TotH2S, TotSO4, TotCH4, TotPO4, TotFeP,TotCaP,TotPads,        & 
     &  TOC, Cprod, Nprod, Pprod, O2prod, Oxicmin, Denitrific, Feredmin,        &
     &  BSRmin, MethMin, nitri1, nitri2, Anammox, Feoxid, H2Soxid,              &
     &  CH4oxid, AOM, FeSprod, FePadsorp, FePdesorp, CaPprod, CaPdiss,          &
     &  NetadsorpP, H2Soxsurf, CH4oxsurf, O2distConsump, AlkProd,               &
     &  DICprodCH4, MPBproduction, NO3uptake, NH3uptake, PO4uptake,             &
     &  DICuptake

! -----------------------------------

      DOUBLE PRECISION  :: Flux(N+1)

      END MODULE commonFESDIA
     
!==========================================================================
!==========================================================================
! subroutine for calculating the biogeochemical rates of
! the fesdia model 
!==========================================================================
!==========================================================================
 
      SUBROUTINE FESDIAbiochem
      USE commonFESDIA
      IMPLICIT NONE

      DOUBLE PRECISION :: Rescale(N), mPads
      DOUBLE PRECISION :: Sum, TotalO2, pO2(N)
      INTEGER :: I,J

! --------------------------------------------------------------------------
! Rate of change due to biogeochemistry 
! --------------------------------------------------------------------------

! Production of DIC and DIN, expressed per cm3 LIQUID/day

       Cprod= (rFast*FDET        +rSlow*SDET        ) * (1.d0-por)/por
       Nprod= (rFast*FDET*NCrFdet+rSlow*SDET*NCrSdet) * (1.d0-por)/por
       Pprod= (rFast*FDET*PCrFdet+rSlow*SDET*PCrSdet) * (1.d0-por)/por

! Oxic mineralisation, denitrification, anoxic mineralisation
! limitation terms
       Oxicminlim= O2/(O2+ksO2oxic)
       Denitrilim= (1.d0-O2/(O2+kinO2denit))*NO3/(NO3+ksNO3denit)
       Feredlim  = (1.d0-O2/(O2+kinO2denit))*                                   &
     &         (1.d0-NO3/(NO3+ksNO3denit))*FeOH3/(FEOH3+ksFEOH3)
       BSRlim    = (1.d0-O2/(O2+kinO2anox))*(1.d0-NO3/(NO3+kinNO3anox))*        &
     &         (1.D0-FeOH3/(FeOH3+kinFeOH3)) * SO4/ (SO4 + ksSO4)
       Methlim= (1.d0-O2/(O2+kinO2anox))*(1.d0-NO3/(NO3+kinNO3anox))*           &
     &         (1.D0-FeOH3/(FeOH3+kinFeOH3))*(1.D0-SO4/(SO4+kinSO4))

       Rescale   = 1.d0/(Oxicminlim+Denitrilim+Feredlim+BSRlim+Methlim)

! mineralisation rates - per liquid
       OxicMin    = Cprod * Oxicminlim * Rescale    ! oxic mineralisation
       Denitrific = Cprod * Denitrilim * Rescale    ! Denitrification
       FeredMin   = Cprod * FeRedlim   * Rescale    ! FeOH3 reduction, /LIQUID
       BSRMin     = Cprod * BSRlim     * Rescale    ! bacterial sulfate reduction
       MethMin    = Cprod * Methlim    * Rescale    ! Methanogenesis

! reoxidation rates
       Nitri1     = rnit1   * NH3 * O2/(O2+ksO2nitri)
       Nitri2     = rnit2   * NO2 * O2/(O2+ksO2nitri)
       Anammox    = ranammox* NH3 * NO2
       Feoxid     = rFeox  * Fe  * O2
       H2Soxid    = rH2Sox * H2S * O2
       CH4oxid    = rCH4ox * CH4 * O2
       AOM        = rAOM   * CH4 * SO4

       pO2(:) = 0.D0
       IF (rH2Soxsurf + rCH4oxsurf .GT. 0) THEN ! long-distance reoxidation
        TotalO2 = 0.D0
        Sum = 0.D0
        DO I = 1, N
            PO2(I) = O2(I)*dx(I)*por(I)
            Sum = Sum + dx(I)*por(I)
            TotalO2 = TotalO2 + O2(I)*dx(I)*por(I)
        ENDDO
        pO2 = pO2/TotalO2
       ENDIF

       TotH2Soxsurf = 0.D0
       H2Soxsurf(:) = 0.d0
       IF (rH2Soxsurf .GT. 0) THEN ! long-distance reoxidation
        H2Soxsurf = rH2Soxsurf*H2S * ALK**2/(ALK**2 + ksALK**2) *                       &
     &    (1.d0-O2/(O2+0.01)) * TotalO2/(TotalO2+ksO2reox) * distreact
        DO I = 1, N
          TotH2Soxsurf = TotH2Soxsurf+H2Soxsurf(I) * por(I)*dx(I)
        ENDDO
       ENDIF

       TotCH4oxsurf = 0.D0
      
       IF (rCH4oxsurf .GT. 0) THEN ! long-distance reoxidation
        CH4oxsurf = rCH4oxsurf*CH4 * ALK**2/(ALK**2 + ksALK**2) *                       &
     &    (1.d0-O2/(O2+0.01)) * TotalO2/(TotalO2+ksO2reox) * distreact
        DO I = 1, N
          TotCH4oxsurf = TotCH4oxsurf  + CH4oxsurf(I) * por(I)*dx(I)
        ENDDO
       ENDIF


! precipitation of FeS
       FeSprod    = rFeS * Fe * H2S

! P adsorption to Fe-oxides if sufficient O2, P desorption other way around
       FePadsorp  = 0.d0

       DO I = 1, N
        if (FeOH3(I) > 0.1) THEN
         FePadsorp(I)  = rFePadsorp * FeOH3(I) * PO4(I)            ! nmol liquid/cm3/d
        endif
        FePdesorp(I)=4.D0*FeredMin(I)*FeP(I)/max(1.d-8, FeOH3(I))  ! nmol liquid/cm3/d
       ENDDO

! P adsorption
       mPads = max(1D-8, maxPads)
         DO I = 1, N
           IF (Pads(I) <  maxPads) THEN
             NetadsorpP(I) = rPads *PO4(I) *(1.d0 - Pads(I)/ mPads)     ! nmol liquid/cm3/d
           ELSE
             NetadsorpP(I) = -rPdes*Pads(I)*(1.d0-por(I))/por(I)        ! nmol liquid/cm3/d
           ENDIF   
         ENDDO   

! P binding to Ca and P-release
       CaPprod   = rCaPprod * PO4* DIC/(DIC+1.D0)     ! nmol liquid/cm3/d
       CaPdiss   = rCaPdiss * CaP                     ! nmol solid/cm3/d

! --------------------------------------------------------------------------
! Update the rate of change with rates due to biogeochemical processes
! Here the increase as in rate factor (e.g. due to temperature) is used!

       IF (ratefac .NE. 1) THEN
         Cprod        = Cprod       *ratefac
         Nprod        = Nprod       *ratefac
         Pprod        = Pprod       *ratefac
         Nitri1       = Nitri1      *ratefac
         Nitri2       = Nitri2      *ratefac
         Anammox      = Anammox     *ratefac
         OxicMin      = OxicMin     *ratefac
         FeredMin     = FeredMin    *ratefac
         MethMin      = MethMin     *ratefac
         BSRMin       = BSRMin      *ratefac
         Feoxid       = Feoxid      *ratefac
         H2Soxid      = H2Soxid     *ratefac
         CH4oxid      = CH4oxid     *ratefac
         FePdesorp    = FePdesorp   *ratefac
         AOM          = AOM         *ratefac
         TotH2Soxsurf = TotH2Soxsurf*ratefac
         TotCH4oxsurf = TotCH4oxsurf*ratefac
         H2Soxsurf    = H2Soxsurf   *ratefac
         CH4oxsurf    = CH4oxsurf   *ratefac
       ENDIF
       
       dFDET = dFDET  - rFast*FDET*ratefac 

       dSDET = dSDET  - rSlow*SDET*ratefac
      
       dO2   = dO2    - OxicMin -1.5d0* Nitri1 - 0.5D0* Nitri2                  &
     &                - 0.25*Feoxid - 2*H2Soxid - 2*CH4oxid 
     
       dNH3  = dNH3  +  (Nprod - Nitri1 - Anammox)  / (1.d0+NH3Ads)
     
       dNO3  = dNO3  - 0.8d0*Denitrific  + Nitri2 
     
       dNO2  = dNO2  - Anammox + Nitri1 - Nitri2 

       DICprodCH4 = - 0.5D0*MethMin + CH4oxid + AOM  

       dCH4  = dCH4  - DICprodCH4

       dDIC  = dDIC  + Cprod + DICprodCH4 +                                     &
     &                 CPrCaP*(CaPdiss*(1.d0-por)/por-CaPprod)   
     
       dPO4  = dPO4  + Pprod + FePdesorp - NetadsorpP                           &
     &                - FePadsorp - CaPprod + CaPdiss*(1.d0-por)/por 
     
       dFeP  = dFeP  + (FePadsorp - FePdesorp)*por/(1.d0-por)
       
       dCaP  = dCaP  + CaPprod*por/(1.d0-por) - CaPdiss
       
       dFe   = dFe   + 4.d0*FeredMin  - Feoxid - FeSprod
       
       dFeOH3= dFeOH3+ (FeOxid -4.d0*FeredMin) * por/(1.d0-por)
       
       dH2S  = dH2S  + 0.5D0 *BSRMin - H2Soxid - FeSprod + AOM
       
       dSO4  = dSO4   -0.5D0 *BSRMin + H2Soxid           - AOM
       
       dPads = dPads  + NetadsorpP*por/(1.d0-por)


! longdistance HSoxidation: effect on alkalinity is split in two parts
! O2 + 4e + 4H+ -> 2H2O                  dALK = +4 (per mmol O)   OR +8 per mmol S
! H2S + 4H2O    -> SO42- + 2e + 10H+     dALK = -10 (per mmol S)

       O2distConsump = 2.D0*(TotH2Soxsurf+TotCH4oxsurf)*pO2
       
       IF (rH2Soxsurf .GT. 0) THEN         ! long-distance reoxidation
        dO2  = dO2  - 2.D0*TotH2Soxsurf*pO2/dx/por
        dH2S = dH2S - H2Soxsurf
        dSO4 = dSO4 + H2Soxsurf
       ENDIF
       
       IF (rCH4oxsurf .GT. 0) THEN         ! long-distance reoxidation
        dO2  = dO2 - 2.D0*TotCH4oxsurf*pO2/dx/por
        dCH4 = dCH4 - CH4oxsurf
        dDIC = dDIC + CH4oxsurf
       ENDIF
       
       IF (AddAlk > 0.1) THEN

! note: dTA = dNH3 - dNO3 - dNO2 - dPO4 - 2*dSO4 + dFe - dF

       AlkProd  = Nprod - 2* Nitri1 - Pprod +                                  &
     &    0.8d0*Denitrific + 8.d0*FeredMin  + BSRMin + 2.d0*AOM                &
     &  - 2.d0 * H2Soxid - 2.d0*Feoxid - 2.d0*FeSprod                          &
     &  - FePdesorp + NetadsorpP + FePadsorp                                   & 
     &  + (CaPprod - CaPdiss*(1.d0-por)/por)*(1.d0+1.8/6.4) 
       !  4.6*H3PO4+1.3H2CO3+1.8*HF+1.45H2O -> CaP     =>> 1pO4 ~ 1.8/4.6 F

        IF (rH2Soxsurf .GT. 0) THEN         ! long-distance reoxidation
          Alkprod = Alkprod +8.D0*TotH2Soxsurf*pO2/dx/por                      &
     &                      -8.D0*H2Soxsurf                 ! or -10???
        ENDIF
       
        dAlk = dAlk + Alkprod

       ELSE  
        Alkprod = 0.d0
        dAlk = 0.d0
       ENDIF       

      RETURN 
      END SUBROUTINE FESDIAbiochem



!==============================================================================
!==============================================================================
! Transport of solid substances
!==============================================================================
!==============================================================================

      SUBROUTINE FESDIAtransolid(FDETdepo, SDETdepo, FePdepo, Padsdepo,         &
     &                           FeOH3depo, CaPdepo)
      USE commonFESDIA
      IMPLICIT NONE

      DOUBLE PRECISION :: FDETdepo,SDETdepo,FePdepo,Padsdepo,                   &
     &                    FeOH3depo,CaPdepo
      DOUBLE PRECISION :: Db(N+1), irrf      
      
      DOUBLE PRECISION zero(N)
      COMMON /myzero/zero
! --------------------------------------------------------------------------
! Rate of change due to transport 
! --------------------------------------------------------------------------

        Db   = Db0   * biotfac

        CALL diff1D (N, FDET, FDETdepo, 0.d0, 0.d0, 0.d0,                       &
     &   1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,           &
     &   Flux, dFDET, irrf)
        FDETdeepflux = Flux(N+1)

        CALL diff1D (N, SDET, SDETdepo, 0.d0, 0.d0, 0.d0,                       &
     &   1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,           &
     &   Flux, dSDET, irrf)
        SDETdeepflux = Flux(N+1)

        FePsurfflux = FePdepo
        CALL diff1D (N, FeP, FePdepo, 0.d0, 0.d0, 0.d0,                         &
     &   1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,           &
     &   Flux, dFeP, irrf)
        FePdeepflux = Flux(N+1)

        CALL diff1D (N, Pads, Padsdepo, 0.d0, 0.d0, 0.d0,                       &
     &   1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,           &
     &   Flux, dPads, irrf)
!        Padsdeepflux = Flux(N+1)

        CaPsurfflux = CaPdepo
        CALL diff1D (N, CaP, CaPdepo, 0.d0, 0.d0, 0.d0,                         &
     &   1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,           &
     &   Flux, dCaP, irrf)
        CaPdeepflux = Flux(N+1)

        CALL diff1D (N, Pads, 0.d0, 0.d0, 0.d0, 0.d0,                         &
     &   1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,           &
     &   Flux, dPads, irrf)
!        Padsdeepflux = Flux(N+1)

        FeOH3surfflux = FeOH3depo
        CALL diff1D (N, FeOH3, FeOH3depo, 0.d0, 0.d0, 0.d0,                     &
     &   1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,           &
     &   Flux, dFeOH3, irrf)
        FeOH3deepflux = Flux(N+1)

      END SUBROUTINE FESDIAtransolid


 !==============================================================================
! Transport of liquid substances
!==============================================================================

      SUBROUTINE FESDIAtranliquid(O2BW, NO3bw, NO2bw, NH3bw, CH4bw,              &
     &    PO4bw,Febw, H2Sbw, SO4bw, DICbw, ALKbw)
     
      USE commonFESDIA
      IMPLICIT NONE

      DOUBLE PRECISION :: O2BW, NO3bw, NO2bw, NH3bw, CH4bw, PO4bw,                     &
     &    Febw, H2Sbw, SO4bw, DICbw, ALKbw
      DOUBLE PRECISION :: Dirr(N)
      DOUBLE PRECISION :: Sum, DS(N+1), irrf
      INTEGER :: BCup, BCDwn
! ----------------------------------------------------------------------------
! transport of O2 and DIC depends on gasflux; 
! if gasflux > 0: dry flat and piston-like exchange

        Dirr = Dirr0 * irrfac

        BCup  = INT(BCupLiq   + 0.1)   ! 2
        BCdwn = INT(BCdownLiq + 0.1)

        if (gasflux > 0) BCup = 4

        Ds = dispO2 * porfac  ! effective diffusion coefficient
        CALL diff1D (N, O2, O2bw , dwO2, gasflux, 0.d0,                         &
     &   BCup,  3, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dO2,irrf)
        O2flux      = Flux(1) + irrf
        O2deepflux  = Flux(N+1)
                        
        Ds  = dispDIC*porfac
        CALL diff1D (N, DIC, DICbw , dwDIC, gasflux, 0.d0,                      &
     &   BCup,  3, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dDIC, irrf)
        DICflux     = Flux(1) + irrf
        DICdeepflux = Flux(N+1)

! Other substances: if dry flat: zero-flux
        BCup  = INT(BCupLiq + 0.1)   ! 2
        if (gasflux > 0) BCup = 5

        Ds  = dispNO3*porfac 
        CALL diff1D (N, NO3, NO3bw ,dwNO3, 0.d0, 0.d0,                          &
     &   BCup,  3, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dNO3,irrf)
        NO3flux     = Flux(1) + irrf
        NO3deepflux = Flux(N+1)

        Ds  = dispNO2*porfac 
        CALL diff1D (N, NO2, NO2bw ,dwNO2, 0.d0, 0.d0,                          &
     &   BCup,  3, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dNO2,irrf)
        NO2flux     = Flux(1) + irrf
        NO2deepflux = Flux(N+1)

        Ds = dispNH3/(1.D0+NH3Ads)*porfac
        CALL diff1D (N, NH3, NH3bw ,dwNH3, 0.d0, 0.d0,                          &
     &   BCup,  3, w, Ds, Dirr/ (1.d0+NH3Ads), Aint, intpor, por,               &
     &   dx, dxint, Flux, dNH3, irrf)
        NH3flux     = (Flux(1)+ irrf)*(1.D0+NH3Ads) 
        NH3deepflux = Flux(N+1)*(1.D0+NH3Ads)

        Ds  = dispPO4*porfac
        CALL diff1D (N, PO4, PO4bw ,dwPO4, 0.d0, 0.d0,                          &
     &   BCup,  3, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dPO4, irrf)
        PO4flux     = Flux(1) + irrf
        PO4deepflux = Flux(N+1)

        Ds = dispCH4*porfac
        CALL diff1D (N, CH4, CH4bw , dwCH4, 0.d0, 0.d0,                         &
     &   BCup,  3, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dCH4, irrf)
        CH4flux     = Flux(1) + irrf
        CH4deepflux = Flux(N+1)

        Ds  = dispFe*porfac
        CALL diff1D (N, Fe, Febw , dwFe, 0.d0, 0.d0,                            &
     &   BCup,  3, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dFe, irrf)
        Feflux     = Flux(1) + irrf
        Fedeepflux = Flux(N+1)

        Ds  = dispH2S*porfac
        CALL diff1D (N, H2S, H2Sbw , dwH2S, 0.d0, 0.d0,                         &
     &   BCup,  3, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dH2S, irrf)
        H2Sflux     = Flux(1) + irrf
        H2Sdeepflux = Flux(N+1)

        Ds  = dispSO4*porfac
        CALL diff1D (N, SO4, SO4bw , dwSO4, 0.d0, 0.d0,                         &
     &   BCup,  3, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dSO4, irrf)
        SO4flux     = Flux(1) + irrf
        SO4deepflux = Flux(N+1)

        Ds  = dispAlk*porfac
        CALL diff1D (N, Alk, ALKbw ,dwAlk, 0.d0, 0.d0,                          &
     &   BCup,  3, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dAlk, irrf)
        Alkflux     = Flux(1) + irrf
        Alkdeepflux = Flux(N+1)

      END SUBROUTINE FESDIAtranliquid

!==========================================================================
!==========================================================================
! subroutine for calculating integrated rates and writing the output
!==========================================================================
!==========================================================================
      
      SUBROUTINE FESDIAout(yout)
      USE commonFESDIA
      IMPLICIT NONE
      
      DOUBLE PRECISION  ::  yout(*), liqfac, solfac
      INTEGER :: I

! ------------------------------------------------------------------------

       TOC = (FDET + SDET)*1200d0*1e-9/2.5 + TOC0
       Cflux = CarbonFlux
       Nflux = cFlux*pFast*NCrFDET + (1.d0-pFast)*cflux*NCrSDET
       Pflux = cFlux*pFast*PCrFDET + (1.d0-pFast)*cflux*PCrSDET

! calculate integrated rates
       totNitri1  = 0.D0
       totNitri2  = 0.D0
       totAnammox = 0.D0
       totFeoxid  = 0.D0
       totH2Soxid = 0.D0
       totCH4oxid = 0.D0
       totNH3ads  = 0.D0
       totDenit   = 0.D0
       totOxic    = 0.D0
       totFered   = 0.D0
       totBSR     = 0.D0
       totMeth    = 0.D0
       totAOM     = 0.D0
       totFeSprod = 0.D0
       totO2prod  = 0.D0
       sumDIN     = 0.D0
       sumPO4     = 0.D0
       totFePprod = 0.D0
       totPadsorp = 0.D0
       totCaPprod = 0.D0
       Premoved   = 0.D0
       totPprod   = 0.D0
       totNProd   = 0.D0
       totFePdesorp = 0.D0
       totCaPdiss   = 0.D0  
       totALkprod   = 0.D0
       totMPBO2prod    = 0.D0
       totMPBPO4uptake = 0.D0
       totMPBNH3uptake = 0.D0
       totMPBNO3uptake = 0.D0
       TotFDET   = 0.D0
       TotSDET   = 0.D0 
       TotFeOH3  = 0.D0
       TotFeP    = 0.D0 
       TotCaP    = 0.D0 
       TotPads   = 0.D0      
       TotO2     = 0.D0 
       TotNO3    = 0.D0  
       TotNO2    = 0.D0 
       TotNH3    = 0.D0 
       TotDIC    = 0.D0 
       TotFe     = 0.D0 
       TotH2S    = 0.D0 
       TotSO4    = 0.D0 
       TotCH4    = 0.D0 
       TotPO4    = 0.D0 
       DO I = 1, N
         liqfac =  por(I)      *dx(I)    ! from /cm3 liquid -> cm2 bulk
         solfac = (1.d0-por(I))*dx(I)    ! from /cm3 solid  -> cm2 bulk
       
         totPprod     = totPprod   + Pprod(I)   * liqfac
         totNprod     = totNprod   + Nprod(I)   * liqfac
         totO2prod    = totO2prod  + O2prod(I)  * liqfac
         totNitri1    = totNitri1  + Nitri1(I)  * liqfac
         totNitri2    = totNitri2  + Nitri2(I)  * liqfac
         totAnammox   = totAnammox + Anammox(I) * liqfac
         totFeoxid    = totFeoxid  + Feoxid(I)  * liqfac
         totH2Soxid   = totH2Soxid + H2Soxid(I) * liqfac
         totCH4oxid   = totCH4oxid + CH4oxid(I) * liqfac
         totALkprod   = TotALkprod + AlkProd(I) * liqfac
         totAOM       = totAOM     + AOM(I)     * liqfac
         totFeSprod   = totFeSprod + FeSprod(I) * liqfac
         totDenit     = totDenit   + Denitrific(I) * liqfac
         totFeRed     = totFeRed   + FeredMin(I)   * liqfac
         totBSR       = totBSR     + BSRMin(I)     * liqfac
         totMeth      = totMeth    + MethMin(I)    * liqfac
         totOxic      = totOxic    + OxicMin(I)    * liqfac
         sumDIN       = sumDIN + (NH3(I)+NO3(I)+NO2(I)) * liqfac
         sumPO4       = sumPO4     + PO4(I)        * liqfac
         TotFePprod   = TotFePprod + FePadsorp(I)  * liqfac
         TotCaPprod   = TotCaPprod + CaPprod(I)    * liqfac
         totFePdesorp = totFePdesorp + FePdesorp(I) * solfac
         TotPadsorp   = TotPadsorp + NetadsorpP(I)  * liqfac
         totCaPdiss   = totCaPdiss + CaPdiss(I)     * solfac
         Premoved = Premoved + (FePdesorp(I) + CaPdiss(I))* solfac
         
         TotMPBO2prod    = TotMPBO2prod    + O2prod(I)    * liqfac
         TotMPBNO3uptake = TotMPBNO3uptake + NO3uptake(I) * liqfac
         TotMPBNH3uptake = TotMPBNH3uptake + NH3uptake(I) * liqfac
         TotMPBPO4uptake = TotMPBPO4uptake + PO4uptake(I) * liqfac

         TotFDET   = TotFDET  + FDET(I)  * solfac
         TotSDET   = TotSDET  + SDET(I)  * solfac
         TotFeOH3  = TotFeOH3 + FeOH3(I) * solfac
         TotFeP    = TotFeP   + FeP(I)   * solfac
         TotCaP    = TotCaP   + CaP(I)   * solfac
         TotPads   = TotPads  + Pads(I)  * solfac
         TotO2     = TotO2    + O2(I)  * liqfac
         TotNO3    = TotNO3   + NO3(I) * liqfac
         TotNO2    = TotNO2   + NO2(I) * liqfac
         TotNH3    = TotNH3   + NH3(I) * liqfac
         TotDIC    = TotDIC   + DIC(I) * liqfac
         TotFe     = TotFe    + Fe(I)  * liqfac
         TotH2S    = TotH2S   + H2S(I) * liqfac
         TotSO4    = TotSO4   + SO4(I) * liqfac
         TotCH4    = TotCH4   + CH4(I) * liqfac
         TotPO4    = TotPO4   + PO4(I) * liqfac

     
       ENDDO
       Premoved = (TotCaPprod + TotFePprod+totPads- Premoved)/TotPprod
       TotMin = totDenit+ totFered + totMeth + totBSR + totOxic
       TotMPBDICuptake = TotMPBO2prod
       
       Nremoved = (totDenit*0.8+TotAnammox*2)/TotNProd
       NPdeep = (NH3(100) + NO3(100)) / max(1e-8,PO4(100))
       partDenit = totDenit / TotMin
       partOxic = totOxic / TotMin
       partFered = totFered / TotMin
       partBSR = totBSR / TotMin
       partMeth = totMeth / TotMin
       NPmean = sumDIN/max(1D-8,sumPO4)
       NPflux = (NO3flux + NH3flux) / PO4flux
       totNH3ads = (totNprod-totNitri1-totAnammox-TotMPBNH3uptake) *          &
     &               (1.d0-1.d0/(1.D0+NH3Ads))

       CALL getoutFES(yout)
      
      RETURN
      END SUBROUTINE FESDIAout

 

!==========================================================================
! put output variables in one vector
!==========================================================================

      SUBROUTINE getoutFES(yout)
      USE dimFESDIA
      DOUBLE PRECISION :: yout(*), out(noutdia), forc(nforc)
      INTEGER :: i

      COMMON /myoutFES  /out
      COMMON /myforcsFES/forc
      
      DO i = 1, Noutdia
       yout(i) = out (i)
      ENDDO       
      DO i = 1, Nforcsdia
       yout(noutdia+i) = forc (i)
      ENDDO       
 
      END SUBROUTINE getoutFES
       

!==============================================================================
! Diffusion in a 1-dimensional finite difference grid
! all inputs are vectors
! subroutine from ReacTran in isnt\doc\fortran directory
!==============================================================================

      SUBROUTINE diff1d (N, C, Cup, Cdown, aup, adown,BcUp, BcDown,         &
     &            v, D, Dirr, Aint, VF, VFmid, dx, dxaux,                   &
     &            Flux, dC, irrf)
      IMPLICIT NONE
      INTEGER N                  ! length of C
C input
      DOUBLE PRECISION C(N)

C Boundary concentrations (used if Bc..=2,4), fluxes (used if Bc= 1)
C and convection coeff (used if Bc=4) Cup and Cdown: concentration or flux
      DOUBLE PRECISION Cup, Cdown, aup, adown

C Diffusion, volume fraction, advection
      DOUBLE PRECISION D(N+1), Dirr(N), Aint(N+1), VF(N+1), VFmid(N), v

C grid size, distance from mid to mid
      DOUBLE PRECISION dx(N), dxaux(N+1)

C boundary concitions (1= flux, 2=conc, 3 = 0-grad, 4=convect)
      INTEGER BcUp, BcDown

C output: fluxes and rate of change
      DOUBLE PRECISION Flux(N+1), dC(N), irrf

C locals
      INTEGER I
      DOUBLE PRECISION AVF, Amid, irrigation, Cbnd

C -------------------------------------------------------------------------------

C Flux - first internal cells

      IF (v >= 0) THEN
       DO I = 2,N
        Flux(I) = -VF(I)*D(I) * (C(I)-C(I-1)) /dxaux(I)                         &
     &           + VF(I)*v*C(I-1)
       ENDDO
      ELSE
       DO I = 2,N
        Flux(I) = -VF(I)*D(I) * (C(I)-C(I-1)) /dxaux(I)                         &
     &           + VF(I)*v*C(I)
       ENDDO
      ENDIF
      
C Then the outer cells
C upstream boundary
      IF (v >= 0) THEN
        Cbnd = Cup
      ELSE
        Cbnd = C(1)
      ENDIF
      
      IF (BcUp .EQ. 1) THEN
        Flux(1) = Cup

      ELSE IF (BcUp .EQ. 2) THEN
        Flux(1) = -VF(1)*D(1) * (C(1)-Cup) /dxaux(1)                            &
     &           + VF(1)*v*Cbnd

      ELSE IF (BcUp .EQ. 3) THEN
        Flux(1) = VF(1)*v*Cbnd

      ELSE IF (BcUp .EQ. 4) THEN
        Flux(1) = aup * (Cup - C(1))
      ELSE
      
        Flux(1) = 0.D0
      ENDIF


C downstream boundary
      IF (v >= 0 .OR. BcDown .eq. 3) THEN
        Cbnd = C(N)
      ELSE
        Cbnd = Cdown
      ENDIF

      IF (BcDown .EQ. 1) THEN
        Flux(N+1) = Cdown

      ELSE IF (BcDown .EQ. 2) THEN
        Flux(N+1) = -VF(N+1)*D(N+1) * (Cdown-C(N)) /dxaux(N+1)                  &
     &              + VF(N+1) * v * Cbnd

      ELSE IF (BcDown .EQ. 3) THEN
        Flux(N+1) = VF(N+1) * v * Cbnd

      ELSE IF (BcDown .EQ. 3) THEN
        Flux(N+1) = -adown * (Cdown-C(N))

      ELSE
        Flux(N+1) = 0.D0
      ENDIF


C Rate of change = negative flux gradient
      DO I = 1,N
        Amid  = 0.5 * (Aint(I)+Aint(I+1))
        dC(I) = -(Aint(I+1)*Flux(I+1) - Aint(I)*Flux(I))/                       &
     &   Amid/VFmid(I)/dx(I)
      ENDDO

! bioirrigation
      irrf = 0.d0
      IF (sum(Dirr) > 0) THEN
        DO I = 1, N
          irrigation = Dirr(I)*(Cup - C(I)) 
          irrf = irrf + Irrigation*dx(I)*VFmid(I)
          dC(I) = dC(I) + irrigation
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE diff1D

