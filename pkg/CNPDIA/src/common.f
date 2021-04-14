!==========================================================================
! THE OMEXDIA model with P, implemented in FORTRAN
!
! This file contains the module declarations and common subroutines
! Karline Soetaert, nioz-yerseke
!==========================================================================

!==========================================================================
! Modules with declarations
!==========================================================================

      MODULE dimCNPDIA
      IMPLICIT NONE

      INTEGER, PARAMETER :: N=100, Np1 = 101 
      INTEGER, PARAMETER :: nparmsMPB = 24
      INTEGER, PARAMETER :: nparms = 5*N + 5*(N+1) + 54 + nparmsMPB
      INTEGER, PARAMETER :: nout    = 2467         ! number of output variables without forcings
      INTEGER, PARAMETER :: noutMPB = 1725         ! output for dynamic MPB

! Flag for MPB type 
      LOGICAL MPBdynamic

      END MODULE dimCNPDIA

!==========================================================================

      MODULE commonCNPDIA
      USE dimCNPDIA

      IMPLICIT NONE

      
! Common state variables      
      DOUBLE PRECISION  :: Fdet(N),Sdet(N),O2(N),NO3(N),NH3(N),ODU(N)
      DOUBLE PRECISION  :: DIC(N),PO4(N),FeP(N),CaP(N),Pads(N),NO2(N)
      DOUBLE PRECISION  :: dFdet(N),dSdet(N),dO2(N),dNO3(N),dNH3(N),            &
     &      dODU(N),dDIC(N),dPO4(N),dFeP(N),dCaP(N),dPads(N),dNO2(N)

! parameters
      DOUBLE PRECISION  :: intpor(N+1),por(N),porfac(N+1),Db0(N+1),             &   
     &  dx(N),dxInt(N+1),x(N),Aint(N+1),Dirr0(N),distreact(N)

      DOUBLE PRECISION  :: NCrFdet,NCrSdet,PCrFDET,PCrSDET,                     &
     &  BCupLiq,BCdownLiq,dwO2,dwNO3,dwNO2,dwNH3,dwODU,dwDIC,dwPO4,             &
     &  NH3Ads,rnitrif1,rnitrif2,ksO2nitri,ranammox,ksNO2anammox,               &
     &  rODUox,rODUoxsurf,ksO2oduox,ksO2oxic,ksNO3denit,kinO2denit,             &
     &  kinNO3anox,kinO2anox,pdepo,rdepo,TOC0,rFePdesorp,rFePadsorp,            &
     &  rCaPprod,rCaPdiss,CPrCaP,rPads,rPdes,maxPads,relax,Cfall,               &
     &  FePfall,CaPfall,kdLight,ksNO3,ksNH3,ksPO4,ksDIC,                        &
     &  DispO2,DispNO3,DispNO2,DispNH3,DispODU,DispPO4,DispDIC

! parameters relating to MPB dynamics
      DOUBLE PRECISION :: MPBflux, Kdchl, QmaxCHLN, QmaxNC, QminNC,             &
     &  sigmaPSII, kdH, tau, kr, phiM, a, maxUpNO3, maxUpNH3,                   &
     &  respMPB, rG, rVNO3, rVN, rRNO3, kEps, maxGrazing, deepMort,             &
     &  rEPS, MPBfall, kwLight                      

      COMMON /myparmsP/NCrFdet,NCrSdet,PCrFDET,PCrSDET,                         &
     &  BCupLiq,BCdownLiq,dwO2,dwNO3,dwNO2,dwNH3,dwODU,dwPO4,dwDIC,             &
     &  NH3Ads,rnitrif1,rnitrif2,ksO2nitri,rAnammox,ksNO2anammox,               &
     &  rODUox,rODUoxsurf,ksO2oduox,ksO2oxic,ksNO3denit,kinO2denit,             &
     &  kinNO3anox,kinO2anox,pdepo,rdepo,TOC0,rFePdesorp,rFePadsorp,            &
     &  rCaPprod,rCaPdiss,CPrCaP,rPads,rPdes,maxPads,relax,Cfall,               &
     &  FePfall,CaPfall,kdLight,ksNO3,ksNH3,ksPO4,ksDIC, MPBflux,               & 
     &  Kdchl, QmaxCHLN, QmaxNC, QminNC, sigmaPSII, kdH, tau, kr,               &
     &  phiM, a, maxUpNO3, maxUpNH3, respMPB, rG, rVNO3, rVN, rRNO3,            &
     &  kEps, maxGrazing, deepMort, rEPS, MPBfall, kwLight,                     &
     &  DispO2,DispNO3,DispNO2,DispNH3,DispODU,DispPO4,                         &
     &  DispDIC,dx,dxint,x,Aint,por,intpor,porfac,Db0,Dirr0,distreact

! Forcing functions: includes water column properties and deposition.
      DOUBLE PRECISION CarbonFlux,bwO2,bwNO3,bwNO2,bwNH3,bwODU,bwPO4,           &
     &   bwDIC,w, biotfac, irrfac, rFast, rSlow ,pFast, MPBforc,                &
     &   CaPflux, FePflux, gasflux, Hwater, ratefac, forc(20)  
      COMMON /myforcsP/ CarbonFlux, w, biotfac, irrfac, rFast,  rSlow,          &
     &   pFast, MPBforc, CaPflux, FePflux, gasflux,                             &
     &   BWO2, bwNO3,bwNO2,bwNH3, bwODU, bwPO4, bwDIC, Hwater, ratefac   

! Other variables
      DOUBLE PRECISION  :: Oxicminlim(N),Denitrilim(N),Anoxiclim(N)
      DOUBLE PRECISION  :: Flux(N+1)
      
! output variables
      DOUBLE PRECISION  :: Cprod(N),Nprod(N),Pprod(N),O2prod(N),TOC(N),         &
     & Oxicmin(N),anoxicmin(N),Denitrific(N),nitri1(N),nitri2(N),               &
     & anammox(N),oduox(N),odudepo(N),FePadsorp(N),FePdesorp(N),                &
     & CaPprod(N),CaPdiss(N),netPadsorp(N),ODUoxsurf(N),TotODUoxsurf,           &
     & TotCProd,MPBproduction(N),NO3uptake(N),NH3uptake(N),                     &
     & PO4uptake(N),DICuptake(N) 

      DOUBLE PRECISION  :: O2flux, O2deepflux, NO3flux,NO3deepflux,             &              
     & NO2flux, NO2deepflux, NH3flux,                                           &  
     & NH3deepflux, ODUflux, ODUdeepflux, PO4flux, PO4deepflux,                 &
     & DICflux, DICdeepflux, FDETflux, FDETdeepflux, SDETflux,                  &
     & SDETdeepflux, FePdeepflux, CaPdeepflux, TotMPBNO3uptake,                 & 
     & totMPBNH3uptake, TotMPBPO4uptake, TotMPBDICuptake, TotNH3ads
      DOUBLE PRECISION  :: TotMin, TotOxic, TotDenit, TotAnoxic,                &
     & partDenit, partAnoxic, partOxic, NPdeep, NPmean, NPflux,                 & 
     & Cflux, Pflux, Nflux, TotFePprod, TotCaPprod, Premoved,                   &
     & Nremoved, totFePdesorp, totCaPdiss, TotNitri1, TotNitri2,                &
     & TotODUox, TotMPBO2prod, TotAnammox, totPAdsorp, OCflux,                  &
     & ODUOCflux, totNProd, totPProd, TotFDET, TotSDET, TotO2, TotNO3,          &
     & TotNO2, TotNH3, TotODU, TotDIC, TotPO4, TotFeP, TotCaP, TotPads 
     
      COMMON /myoutP      /O2flux, O2deepflux, NO3flux, NO3deepflux,            &
     &  NO2flux, NO2deepflux, NH3flux, NH3deepflux, ODUflux,                    &
     &  ODUdeepflux, PO4flux,PO4deepflux, DICflux, DICdeepflux,                 &
     &  FDETflux, FDETdeepflux, SDETflux, SDETdeepflux, FePdeepflux,            & 
     &  CaPdeepflux, Cflux, Nflux, Pflux, OCflux, ODUOCflux, NPflux,            &
     &  NPmean, NPdeep, TotMin, TotOxic, TotDenit, TotAnoxic,                   &
     &  partOxic, partDenit, partAnoxic, TotNitri1, TotNitri2,                  &
     &  TotAnammox, TotODUox, TotFePprod, TotCaPprod, totFePdesorp,             &
     &  totCaPdiss, totPadsorp, TotODUoxsurf, TotNH3ads, Premoved,              &
     &  Nremoved, TotNprod, TotPprod, TotMPBNO3uptake, TotMPBNH3uptake,         &
     &  TotMPBPO4uptake, TotMPBDICuptake, TotMPBO2prod, TotFDET,                &
     &  TotSDET, TotO2, TotNO3, TotNO2, TotNH3, TotODU, TotDIC, TotPO4,         &
     &  TotFeP, TotCaP, TotPads, TOC, Cprod, Nprod, Pprod, O2prod,              &
     &  Oxicmin, Denitrific, anoxicmin, nitri1, nitri2, anammox, oduox,         &
     &  ODUoxsurf, odudepo, FePadsorp, FePdesorp, CaPprod, CaPdiss,             &
     &  netPadsorp, MPBproduction, NO3uptake, NH3uptake, PO4uptake,             &
     &  DICuptake, forc


      END MODULE commonCNPDIA

!==========================================================================
!==========================================================================
! subroutine for calculating the biogeochemical rates of
! the omexdia model 
!==========================================================================
!==========================================================================
      
      SUBROUTINE CNPDIAbiochem 

      USE commonCNPDIA
      IMPLICIT NONE
      
!......................... declaration section.............................
      DOUBLE PRECISION :: Rescale(N),mPads
      DOUBLE PRECISION :: Sum, pO2(N), O2lim(N), TotalO2
      
      INTEGER :: I
      
      CHARACTER(len=80) msg
!............................ statements ..................................

! --------------------------------------------------------------------------
! Rate of change due to biogeochemistry 
! --------------------------------------------------------------------------

! Production of DIC and DIN, expressed per cm3 LIQUID/day

      Cprod= (rFast*FDET        +rSlow*SDET        ) * (1.d0-por)/por
      Nprod= (rFast*FDET*NCrFdet+rSlow*SDET*NCrSdet) * (1.d0-por)/por
      Pprod= (rFast*FDET*PCrFdet+rSlow*SDET*PCrSdet) * (1.d0-por)/por

! Oxic mineralisation, denitrification, anoxic mineralisation

! first the limitation terms
      Oxicminlim = O2/(O2+ksO2oxic)                ! limitation terms
      Denitrilim = (1.d0-O2/(O2+kinO2denit)) * NO3/(NO3+ksNO3denit)
      Anoxiclim  = (1.d0-O2/(O2+kinO2anox)) * (1-NO3/(NO3+kinNO3anox))
      Rescale    = 1.d0/(Oxicminlim+Denitrilim+Anoxiclim)

! then the mineralisation rates
      OxicMin    = Cprod*Oxicminlim*Rescale        ! oxic mineralisation
      Denitrific = Cprod*Denitrilim*Rescale        ! Denitrification
      AnoxicMin  = Cprod*Anoxiclim *Rescale        ! anoxic mineralisation

! reoxidation and ODU deposition
      Nitri1     = rnitrif1* NH3    *O2/(O2+ksO2nitri)
      Nitri2     = rnitrif2* NO2    *O2/(O2+ksO2nitri)
      
      Anammox    = ranammox* NH3 * NO2/(NO2+ksNO2anammox)
      OduOx      = rODUox  * ODU *  O2/(O2 +ksO2oduox)

      TotODUoxsurf = 0.D0
      
      IF (rODUoxsurf .GT. 0) THEN ! long-distance reoxidation
        ODUoxsurf = rODUoxsurf*ODU*O2(1)/(O2(1)+ksO2oduox)*distreact
        TotalO2 = 0.D0
        DO I = 1, N
          TotODUoxsurf = TotODUoxsurf  + ODUoxsurf(I) * por(I)*dx(I)
          pO2(I) = O2(I) *por(I)*dx(I)
          TotalO2 = TotalO2 + pO2(I)
        ENDDO
        pO2 = pO2/TotalO2
      ENDIF

      OduDepo    = AnoxicMin*pDepo + ODU*rdepo

! P adsorption
       mPads = max(1D-8, maxPads)
         DO I = 1, N
           IF (Pads(I) <  maxPads) THEN
             NetPadsorp(I) = rPads *PO4(I) *(1.d0 - Pads(I)/ mPads)     ! nmol liquid/cm3/d
           ELSE
             NetPadsorp(I) = -rPdes*Pads(I)*(1.d0-por(I))/por(I)         ! nmol liquid/cm3/d
           ENDIF   
         ENDDO   

! P adsorption to Fe-oxides if sufficient O2, P desorption other way around
! Note: for very small negative O2 this creates a problem as then the Monod = negative
      DO I = 1, N
        O2lim(I) = max(0.D0, O2(I)/(O2(I)+0.1D0))
      ENDDO  
      
      FePadsorp  = rFePadsorp * O2lim* PO4                ! nmol liquid/cm3/d
      FePdesorp  = rFePdesorp * (1.d0-O2lim) *FeP         ! nmol solid/cm3/d

! P binding to Ca and P-release
      CaPprod    = rCaPprod * PO4* DIC/(DIC+1.D0)         ! nmol liquid/cm3/d
      CaPdiss    = rCaPdiss * CaP                         ! nmol solid/cm3/d

! --------------------------------------------------------------------------
! Update the rate of change with rates due to biogeochemical processes
! --------------------------------------------------------------------------

       IF (ratefac .NE. 1) THEN
         Cprod        = Cprod          *ratefac
         Nprod        = Nprod          *ratefac
         Pprod        = Pprod          *ratefac
         Nitri1       = Nitri1         *ratefac
         Nitri2       = Nitri2         *ratefac
         Anammox      = Anammox        *ratefac
         OxicMin      = OxicMin        *ratefac
         Denitrific   = Denitrific     *ratefac
         AnoxicMin    = AnoxicMin      *ratefac
         ODUox        = ODUox          *ratefac
         TotODUoxsurf = TotODUoxsurf   *ratefac
         ODUoxsurf    = ODUoxsurf      *ratefac
       ENDIF
       

      dFDET = dFDET - rFast*FDET*ratefac      
      dSDET = dSDET - rSlow*SDET*ratefac 

      dO2   = dO2   - OxicMin -1.5d0* Nitri1 - 0.5d0*Nitri2 - OduOx

      dNH3  = dNH3  + (Nprod -Nitri1- Anammox)/(1.d0+NH3Ads)
      dNO2  = dNO2  -  Anammox  + Nitri1 - Nitri2  
      dNO3  = dNO3  - 0.8d0*Denitrific   + Nitri2 
      dODU  = dODU  + AnoxicMin                   - OduOx - OduDepo

      dDIC  = dDIC  + Cprod + CPrCaP*(CaPdiss*(1.d0-por)/por-CaPprod)

      dPO4  = dPO4  + Pprod - FePadsorp - CaPprod                         &
     &              + (FePdesorp+CaPdiss)*(1.d0-por)/por - NetPadsorp
      dFeP  = dFeP  + FePadsorp*por/(1.d0-por) - FePdesorp
      dCaP  = dCaP  + CaPprod*por/(1.d0-por) - CaPdiss
      dPads = dPads  + NetPadsorp*por/(1.d0-por)

      IF (rODUoxsurf .GT. 0) THEN         ! long-distance reoxidation
        dO2 = dO2 - TotODUoxsurf*pO2/dx/por
        dODU = dODU - ODUoxsurf
      ENDIF


      END SUBROUTINE CNPDIAbiochem

      
!==========================================================================
!==========================================================================
! subroutine calculating integrated rates and writing the output
!==========================================================================
!==========================================================================
      
      SUBROUTINE CNPDIAout 
      USE commonCNPDIA
      IMPLICIT NONE
      
      INTEGER :: I
      DOUBLE PRECISION :: solfac, liqfac
! c--------------------------------------------------------------------

      Cflux = CarbonFlux
      Nflux = cFlux*pFast*NCrFDET + (1.d0-pFast)*cflux*NCrSDET
      Pflux = cFlux*pFast*PCrFDET + (1.d0-pFast)*cflux*PCrSDET
      totDenit = 0.D0
      totOxic = 0.D0
      totAnoxic = 0.D0
      TotNitri1 = 0.D0
      TotNitri2 = 0.D0
      totAnammox = 0.D0
      TotODUox = 0.D0      
      TotFePprod = 0.D0
      TotCaPprod = 0.D0
      TotFePdesorp = 0.D0 
      TotCaPdiss = 0.D0
      TotPadsorp = 0.D0
      Premoved = 0.D0
      TotNprod = 0.D0
      TotPprod = 0.D0
      TotNH3ads = 0.D0
      TotFDET  = 0.D0
      TotSDET  = 0.D0 
      TotO2    = 0.D0
      TotNO3   = 0.D0
      TotNO2   = 0.D0 
      TotNH3   = 0.D0 
      TotODU   = 0.D0 
      TotDIC   = 0.D0 
      TotPO4   = 0.D0 
      TotFeP   = 0.D0 
      TotCaP   = 0.D0 
      TotPads  = 0.D0
      DO I = 1, N
         liqfac =  por(I)      *dx(I)    ! from /cm3 liquid -> cm2 bulk
         solfac = (1.d0-por(I))*dx(I)    ! from /cm3 solid  -> cm2 bulk
         
         totPprod   = totPprod  + Pprod(I)        * liqfac
         totNprod   = totNprod  + Nprod(I)        * liqfac
         totDenit   = totDenit  + Denitrific(I)   * liqfac
         totAnoxic  = totAnoxic + AnoxicMin(I)    * liqfac
         totOxic    = totOxic   + OxicMin(I)      * liqfac
         TotNitri1  = TotNitri1 + nitri1(I)       * liqfac
         TotNitri2  = TotNitri2 + nitri2(I)       * liqfac
         TotAnammox = TotAnammox + Anammox(I)     * liqfac
         TotODUox   = TotODUox   + ODUox(I)       * liqfac
         TotFePprod = TotFePprod + FePadsorp(I)   * liqfac
         TotCaPprod = TotCaPprod + CaPprod(I)     * liqfac
         TotPadsorp = TotPadsorp + netPadsorp(I)  * liqfac

         Premoved     = Premoved + (FePdesorp(I)+CaPdiss(I)) * solfac  
         TotFePdesorp = TotFePdesorp + FePdesorp(I)          * solfac  
         TotCaPdiss   = TotCaPdiss + CaPdiss(I)              * solfac  

         TotFDET  = TotFDET  + FDET (I)  * solfac  
         TotSDET  = TotSDET  + SDET (I)  * solfac  
         TotO2    = TotO2    + O2 (I)    * liqfac
         TotNO3   = TotNO3   + NO3(I)    * liqfac
         TotNO2   = TotNO2   + NO2(I)    * liqfac
         TotNH3   = TotNH3   + NH3(I)    * liqfac
         TotODU   = TotODU   + ODU(I)    * liqfac
         TotDIC   = TotDIC   + DIC(I)    * liqfac
         TotPO4   = TotPO4   + PO4(I)    * liqfac
         TotFeP   = TotFeP   + FeP  (I)  * solfac  
         TotCaP   = TotCaP   + CaP  (I)  * solfac  
         TotPads  = TotPads  + Pads (I)  * solfac  

      ENDDO
      
      TotMin    = totDenit + totAnoxic + totOxic
      partDenit = totDenit / TotMin
      partOxic  = totOxic / TotMin
      partAnoxic = totAnoxic / TotMin

      Premoved = (TotCaPprod+TotFePprod+TotPadsorp-Premoved)/TotPprod
      Nremoved = (totDenit*0.8+TotAnammox*2)/TotNProd
      
      NPdeep = (NH3(100) + NO3(100) + NO2(100)) / max(1e-8,PO4(100))
      NPmean = (totNO3 + totNO2 + totNH3)/max(1D-8, totPO4)
      NPflux = (NO3flux + NH3flux + NO2flux) / PO4flux

      OCflux    = O2flux/DICflux
      ODUOCflux = (O2flux+ODUflux)/DICflux 
      totNH3ads = (totNprod-totNitri1-totAnammox) *                            &
     &               (1.d0-1.d0/(1.D0+NH3Ads))
      
      RETURN
      END SUBROUTINE CNPDIAout

!==========================================================================
! put output variables in one vector
!==========================================================================

      SUBROUTINE getoutP(yout)
      USE dimCNPDIA
      
      IMPLICIT NONE
      INTEGER :: i
      DOUBLE PRECISION :: yout(*), out(nout+20), forc(20)

      COMMON /myoutP   /out
      COMMON /myforcsP /forc
      
      DO i = 1, nout
       yout(i) = out (i)
      ENDDO       
      DO i = 1, 20
       yout(nout+i) = forc (i)
      ENDDO       
 
      END SUBROUTINE getoutP
       
!==============================================================================
! Transport of solid substances
!==============================================================================

      SUBROUTINE CNPDIAtransolid(FDETdepo, SDETdepo, FePdepo, CaPdepo)
      USE commonCNPDIA
      IMPLICIT NONE

      DOUBLE PRECISION :: FDETdepo,SDETdepo,FePdepo,CaPdepo
      DOUBLE PRECISION  :: Db(N+1), irrf      
      
      DOUBLE PRECISION zero(N)
      COMMON /myzero/zero
! ----------------------------------------------------------------------------
      
      TOC = (FDET + SDET)*1200d0*1e-9/2.5 + TOC0
      Db   = Db0   * biotfac

      CALL diff1D (N, FDET, FDETdepo, 0.d0, 0.d0, 0.d0,                         &
     & 1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,             &
     & Flux, dFDET, irrf)
      FDETdeepflux = Flux(N+1)

      CALL diff1D (N, SDET, SDETdepo, 0.d0, 0.d0, 0.d0,                         &
     & 1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,             &
     & Flux, dSDET, irrf)
      SDETdeepflux = Flux(N+1)

      CALL diff1D (N, FeP,FePdepo, 0.d0, 0.d0, 0.d0,                            &
     & 1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,             &
     & Flux, dFeP, irrf)
      FePdeepflux = Flux(N+1)

      CALL diff1D (N, CaP, CaPdepo, 0.d0, 0.d0, 0.d0,                           &
     & 1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,             &
     & Flux, dCaP, irrf)
      CaPdeepflux = Flux(N+1)

      CALL diff1D (N, Pads, 0.d0, 0.d0, 0.d0, 0.d0,                             &
     & 1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,             &
     & Flux, dPads, irrf)
!      Padsdeepflux = Flux(N+1)

      END SUBROUTINE 

!==============================================================================
! Transport of liquid substances
!==============================================================================

      SUBROUTINE CNPDIAtranliquid(O2BW, NO3bw, NH3bw, ODUbw, PO4bw,             &
     &    DICbw, NO2bw)
     
      USE commonCNPDIA
      IMPLICIT NONE

      DOUBLE PRECISION :: O2BW, NO3bw, NH3bw, ODUbw, PO4bw, DICbw, NO2bw
      DOUBLE PRECISION :: Dirr(N)
      DOUBLE PRECISION :: Sum, DS(N+1), irrf
      INTEGER :: BCup, BCDwn
! ----------------------------------------------------------------------------

      Dirr = Dirr0 * irrfac

! transport of O2 and DIC depends on gasflux; if > 0: dry flat and piston-like exchange
      BCup  = INT(BCupLiq + 0.1)   ! 2
      BCdwn = INT(BCdownLiq + 0.1)
      if (gasflux > 0) BCup = 4

      Ds = dispO2 * porfac  ! effective diffusion coefficient
      CALL diff1D (N, O2, O2bw , dwO2, gasflux, 0.d0,                              &
     &   BCup, BCdwn, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dO2, irrf)
      O2flux      = Flux(1) + irrf
      O2deepflux  = Flux(N+1)
                        
      Ds  = dispDIC * porfac
      CALL diff1D (N, DIC, DICbw, dwDIC, gasflux, 0.d0,                            &
     &   BCup, BCdwn, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dDIC, irrf)
      DICflux     = Flux(1) + irrf
      DICdeepflux = Flux(N+1)

! Other substances: if dry flat: imposed flux = 0
      BCup  = INT(BCupLiq + 0.1)   ! 2
      if (gasflux > 0) BCup = 5

      Ds  = dispNO3*porfac 
      CALL diff1D (N, NO3, NO3bw , dwNO3, 0.d0, 0.d0,                              &
     &   BCup, BCdwn, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dNO3, irrf)
      NO3flux     = Flux(1) + irrf
      NO3deepflux = Flux(N+1)

      Ds  = dispNO2*porfac 
      CALL diff1D (N, NO2, NO2bw , dwNO2, 0.d0, 0.d0,                              &
     &   BCup, BCdwn, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dNO2, irrf)
      NO2flux     = Flux(1) + irrf
      NO2deepflux = Flux(N+1)

      Ds = dispNH3/(1.D0+NH3Ads)*porfac
      CALL diff1D (N, NH3, NH3bw , dwNH3, 0.d0, 0.d0,                              &
     &   BCup, BCdwn, w, Ds, Dirr/ (1.d0+NH3Ads), Aint, intpor, por,               &
     &   dx, dxint, Flux, dNH3, irrf)
      NH3flux     = (Flux(1)+ irrf)*(1.D0+NH3Ads) 
      NH3deepflux = Flux(N+1)*(1.D0+NH3Ads)

      Ds  = dispPO4*porfac
      CALL diff1D (N, PO4, PO4bw , dwPO4, 0.d0, 0.d0,                              &
     &   BCup, BCdwn, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dPO4, irrf)
      PO4flux     = Flux(1) + irrf
      PO4deepflux = Flux(N+1)

      Ds = dispODU*porfac
      CALL diff1D (N, ODU, ODUbw , dwODU, 0.d0, 0.d0,                              &
     &   BCup, BCdwn, w, Ds, Dirr, Aint, intpor, por, dx, dxint,                   &
     &   Flux, dODU, irrf)
      ODUflux     = Flux(1) + irrf
      ODUdeepflux = Flux(N+1)

      END SUBROUTINE
      
!==============================================================================
! Diffusion in a 1-dimensional finite difference grid
! all inputs are vectors
! subroutine from ReacTran in isnt\doc\fortran directory
!==============================================================================

      SUBROUTINE diff1d (N, C, Cup, Cdown,  aup, adown, BcUp, BcDown,        &
     &            v, D, Dirr, Aint, VF, VFmid, dx, dxaux,                    &
     &            Flux, dC, irrf)
      IMPLICIT NONE
      INTEGER N                  ! length of C
C input
      DOUBLE PRECISION C(N)

C Boundary concentrations (used if Bc.. = 2,4), fluxes (used if Bc= 1)
C and convection coeff (used if Bc=4) Cup, Cdown: either conc or flux
      DOUBLE PRECISION Cup, Cdown, aup, adown

C Diffusion, volume fraction, advection
      DOUBLE PRECISION D(N+1), Dirr(N), Aint(N+1), VF(N+1), VFmid(N), v

C grid size, distance from mid to mid
      DOUBLE PRECISION dx(N), dxaux(N+1)

C boundary concitions (1= flux, 2=conc, 3 = 0-grad, 4 = convection)
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

      ELSE IF (BcDown .EQ. 4) THEN
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
