!==========================================================================
! Microphytobenthos dynamics, to be used with CNPDIA
!
! Karline Soetaert, nioz-yerseke
! Main modules and transport functions are in file common.f
!==========================================================================

      MODULE commonMPBDIA
      use dimCNPDIA  !N, Nparms, Nout etc...
      IMPLICIT NONE

! state variables and derivatives: 
! Microphytobenthos C, N, Chlorophyll (MPBC, MPBN, Chl), 
! Extracellular polymeric substances  (EPS),  
! Inhibited state of photosynthesis system II (PSIIin)

      DOUBLE PRECISION ::  MPBC(N), MPBN(N), CHL(N), EPS(N), PSIIin(N)
      DOUBLE PRECISION :: dMPBC(N),dMPBN(N),dCHL(N),dEPS(N),dPSIIin(N)


! Output variables
      DOUBLE PRECISION  :: PrimProd(N),Chlproduction(N),MPBResp(N),             &
     & MinEPS(N),ProdEPS(N),MPBCdeath(N),MPBNdeath(N),Chldeath(N),              &
     & ChlCrMPB(N),NCrMPB(N),ChlNrMPB(N), CNrMPB(N),Light(N),                   &
     & FDETprodMPB(N),DICprodMPB(N),NH3prodMPB(N),PO4prodMPB(N),                &
     & totMPBprod,totMPBresp,                                                   &
     & totMPBCdeath,totMPBNdeath,totProdEPS,totChlProd, TotPrimProd,            &
     & totMinEps,totFDETprodMPB,totDICprodMPB,totNH3prodMPB,                    &
     & totPO4prodMPB, TotMPBO2uptake, TotalMPBC, TotalMPBN, TotalChl,           &
     & TotalEPS

      DOUBLE PRECISION  :: MPBCflux,MPBCdeepflux,MPBNflux,MPBNdeepflux,         &              
     & Chlflux, Chldeepflux, EPSflux, EPSdeepflux 

      COMMON /myoutMPB    /MPBCflux,MPBCdeepflux,MPBNflux,MPBNdeepflux,         &              
     & Chlflux, Chldeepflux, EPSflux, EPSdeepflux, totMPBprod,                  &
     & totMPBresp, totMPBCdeath, totMPBNdeath, totProdEPS, totChlProd,          &
     & totPrimProd, totMinEps, totFDETprodMPB, totDICprodMPB,                   &
     & totNH3prodMPB, totPO4prodMPB, TotMPBO2uptake, TotalMPBC,                 &
     & TotalMPBN, TotalChl, TotalEPS, PrimProd,                                 &
     & Chlproduction, MPBResp, MPBCdeath, MPBNdeath, Chldeath,                  &
     & ProdEPS, MinEPS, ChlCrMPB, NCrMPB, ChlNrMPB, CNrMPB,                     &
     & Light, FDETprodMPB,DICprodMPB, NH3prodMPB,PO4prodMPB 

      DOUBLE PRECISION :: O2uptake(N)

      END MODULE commonMPBDIA


!==========================================================================
! initialise the common block with parameter values, 
! followed by thicknesses, surface, porosities, bioturbation values
!==========================================================================

      SUBROUTINE initmpbdiap (steadyparms)
      USE dimCNPDIA
      
      IMPLICIT NONE
      EXTERNAL steadyparms

      DOUBLE PRECISION parms(nparms)
      COMMON /myparmsP/parms

      DOUBLE PRECISION zero(N)
      COMMON /myzero/zero

       Zero(:) = 0.D0
       CALL steadyparms(nparms, parms)
       MPBdynamic = .TRUE.
      
      RETURN
      END SUBROUTINE initmpbdiap

!==========================================================================
!==========================================================================
! subroutine calculating the rate of change of
! the MPBdia model - here the quantities in the common blocks are named
!==========================================================================
!==========================================================================
      
      SUBROUTINE mpbdiamodp (neq, t, Conc, dConc, yout, ip)

      USE commonCNPDIA
      USE commonMPBDIA
      IMPLICIT NONE
      
!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i

      DOUBLE PRECISION  :: t, Conc(17*N), dConc(17*N), yout(*)

      CHARACTER(len=80) msg
!............................ statements ..................................

!     check memory allocated to output variables
      IF (ip(1) < nout)  CALL rexit("nout not large enough") 

! from Conc to fdet, sdet, o2,...
      DO I = 1, N
        Fdet(I) = Conc(I)
        Sdet(I) = Conc(N+I)
        O2(I)   = Conc(2*N+I)
        NO3(I)  = Conc(3*N+I)
        NO2(I)  = Conc(4*N+I)
        NH3(I)  = Conc(5*N+I)
        ODU(I)  = Conc(6*N+I)
        DIC(I)  = Conc(7*N+I)
        PO4(I)  = Conc(8*N+I)
        FeP(I)  = Conc(9*N+I)
        CaP(I)  = Conc(10*N+I)
        Pads(I) = Conc(11*N+I)
        MPBC(I) = Conc(12*N+I)
        MPBN(I) = Conc(13*N+I)
        CHL(I)  = Conc(14*N+I)
        EPS(I)  = Conc(15*N+I)
        PSIIin(I) = Conc(16*N+I)
        
      ENDDO
      
! --------------------------------------------------------------------------
      FDETFLux = Carbonflux*pFast
      SDETflux = Carbonflux - FDETflux
      
      CALL CNPDIAtransolid(FDETFLux,SDETFLux,FePflux,CaPflux)

      MPBCflux = MPBflux
      MPBNflux = MPBflux*NCrFDET
      Chlflux  = MPBflux*NCrFDET   ! assumes a Chl:N ratio of 1
      EPSflux  = 0.D0
      CALL MPBDIAtransolid(MPBCflux,MPBNflux,CHLflux,EPSflux)
      CALL CNPDIAtranliquid(bwO2,bwNO3,bwNH3,bwODU,bwPO4,bwDIC,bwNO2)

      CALL CNPDIAbiochem 
      CALL MPBDIAdynamic

      CALL CNPDIAout         ! integrated vars
      CALL MPBDIAout     
      CALL getoutP(yout)     ! store output vars in yout
      CALL getoutMPB(yout)   


! --------------------------------------------------------------------------

! from dfdet, dsdet, do2,... to dconc

      DO I = 1, N
         dConc(I)      = dFdet(I)
         dConc(N+I)    = dSdet(I) 
         dConc(2*N+I)  = dO2(I)  
         dConc(3*N+I)  = dNO3(I) 
         dConc(4*N+I)  = dNO2(I)
         dConc(5*N+I)  = dNH3(I) 
         dConc(6*N+I)  = dODU(I)
         dConc(7*N+I)  = dDIC(I)
         dConc(8*N+I)  = dPO4(I)
         dConc(9*N+I)  = dFeP(I)
         dConc(10*N+I) = dCaP(I)
         dConc(11*N+I) = dPads(I)
         dConc(12*N+I) = dMPBC(I)
         dConc(13*N+I) = dMPBN(I)
         dConc(14*N+I) = dCHL(I)
         dConc(15*N+I) = dEPS(I)
         dConc(16*N+I) = dPSIIin(I)
         
      ENDDO

      RETURN
      END SUBROUTINE mpbdiamodp
      
!==============================================================================
! Transport of MPB solid substances
!==============================================================================

      SUBROUTINE MPBDIAtransolid(MPBCdepo,MPBNdepo,CHLdepo,EPSdepo)
      USE commonCNPDIA
      USE commonMPBDIA
      IMPLICIT NONE

      DOUBLE PRECISION :: Db(N+1), irrf      
      DOUBLE PRECISION :: MPBCdepo,MPBNdepo,CHLdepo,EPSdepo
      
      DOUBLE PRECISION zero(N)
      COMMON /myzero/zero
! ----------------------------------------------------------------------------
      
      Db   = Db0   * biotfac

      CALL diff1D (N, MPBC, MPBCdepo, 0.d0, 0.d0, 0.d0,                        &
     & 1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,            &
     & Flux, dMPBC, irrf)
      MPBCdeepflux = Flux(N+1)

      CALL diff1D (N, MPBN, MPBNdepo, 0.d0, 0.d0, 0.d0,                        &
     & 1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,            &
     & Flux, dMPBN, irrf)
      MPBNdeepflux = Flux(N+1)

      CALL diff1D (N, CHL, Chldepo, 0.d0, 0.d0, 0.d0,                          &
     & 1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,            &
     & Flux, dCHL, irrf)
      Chldeepflux = Flux(N+1)

      CALL diff1D (N, EPS, EPSdepo, 0.d0, 0.d0, 0.d0,                          &
     & 1, 3, w, Db, Zero, Aint, 1.d0-intpor, 1.d0 - por, dx, dxint,            &
     & Flux, dEPS, irrf)
      EPSdeepflux = Flux(N+1)

      END SUBROUTINE MPBDIAtransolid

!==============================================================================
! MPB dynamics
!==============================================================================

      SUBROUTINE MPBDIAdynamic
      USE commonCNPDIA
      USE commonMPBDIA
      IMPLICIT NONE
      
! Local variables
      DOUBLE PRECISION :: limNH3(N), limNO3(N), limDIC(N), limPO4(N)
      DOUBLE PRECISION :: PSIIo(N), PSIIc(N), PSDamage(N), PSrepair(N)
      DOUBLE PRECISION :: DroopNC(N), Nut_quota(N), DroopChlN(N)
      DOUBLE PRECISION :: rhoChl(N), sol2liq(N), Nresp(N), Chlresp(N)
      
      DOUBLE PRECISION :: Small_val, Eint, Kdsum, fac
      DOUBLE PRECISION :: ChlCr, ChlNr, PNrMPB, QmaxChlC
      INTEGER :: I

      small_val = 1.D-8

      sol2liq = (1.d0-por)/por            ! from /solid to /liquid

! The quota
      ChlCrMPB = CHL  /(MPBC + small_val)
      NCrMPB   = MPBN /(MPBC + small_val)
      ChlNrMPB = CHL  /(MPBN + small_val)
      CNrMPB   = MPBC /(MPBN + small_val) 
      PNrMPB   = PCrFdet / NCrFdet        ! PN of MPB: same as FDET = 1/16

! Respiration initialised to basal respiration rate
      MPBResp = respMPB   ! /day

! =======================================================
!  Nutrient uptake
! =======================================================
      limNH3 = NH3/(NH3 + ksNH3)
      limNO3 = NO3/(NO3 + ksNO3) * exp(-1.5 * NH3) 
      limPO4 = PO4/(PO4 + ksPO4)
      
      DO I = 1, N
        limNH3(I) = min(limNH3(I),limPO4(I))
        limNO3(I) = min(limNO3(I),limPO4(I))
      ENDDO
      
! feedback from N status (from Flynn 2001 and Baklouti 2006)
      Nut_Quota = (QmaxNC - NCrMPB)/(QmaxNC - QminNC)   ! 0 for max NCrMPB, 1 for min

      NH3uptake = maxUpNH3 * limNH3 * Nut_Quota         ! mmolN/mmolC/d
      NO3uptake = maxUpNO3 * limNO3 * Nut_Quota

! Respiration cost related to N uptake and NO3 reduction: 
      DO I = 1, N
        MPBResp(I) = MPBResp(I) + (rRNO3+rVNO3)*max(0.D0, NO3uptake(I))
        MPBResp(I) = MPBResp(I) +  rVN         *max(0.D0, NH3uptake(I))  ! /day
      ENDDO
      
! Total nutrient uptake (can be negative)
      NH3uptake = NH3uptake * MPBC                      ! mmolN/m3solid/d
      NO3uptake = NO3uptake * MPBC                      ! mmolN/m3solid/d
      PO4uptake = (NH3uptake + NO3uptake)*PNrMPB        ! mmolP/m3solid/d
      
! ============================================================
!                   LIGHT DEPENDENT PHOTOSYNTHESIS 
! ============================================================
      
! ============================================================
! Light profile
! ============================================================

! Light at sediment-water interface - includes extinction by water
      Eint = MPBforc*86400D0*exp(-kwLight*Hwater)     ! [uEinst/m2/d]

      DO I = 1, N
        kdsum = kdLight + CHL(I)*kdChl                
        Light(I) = Eint * exp(-(kdsum * dx(i)/2))    ! par in grid middle
        Eint     = Eint * exp(-(kdsum * dx(i)))      ! par at interface
      ENDDO
      
! Fraction of PSII states: open, closed, inhibited  (sum = 1)
      PSIIo  = (1.d0 - PSIIin)/(1.d0+sigmaPSII*tau*Light) ! fraction PSII open                            
      PSIIc  = 1.d0 - PSIIin - PSIIo                      ! fraction PSII closed                

! Rates that change the fraction inhibited      
      PSDamage = kdH * sigmaPSII * Light * PSIIc                ! (c -> i) photoinhibition of PSII system       
      PSRepair = kr * PSIIin                                    ! (i -> c) damage repair                        
      
! Chl-a-specific primary production rate              
      PrimProd = PSIIo * phiM * a * Light                       ! [mmolC/mgChl/d]       

! Dependency of Photsynthesis and chl Synthesis on NCratio          
      DO I = 1, N
        DroopNC(I)  = max(0.D0, 1.d0 - QminNC/(NCrMPB(I) + 1D-8))                  ! [-]   
      ENDDO      
! Limitation by DIC      
      limDIC = DIC/(DIC + ksDIC)
      
! C-biomass -specific primary production rate     
      QmaxChlC = QmaxCHLN*QmaxNC
      DO I = 1, N
        ChlCr = max(0.D0, min (QmaxChlC, ChlCrMPB(I))) 
        MPBproduction(I) = PrimProd(I)*ChlCr*limDIC(I)*DroopNC(I)    ! [/d]  
      ENDDO

! Cost of photosynthesis added to respiration      
      MPBResp = MPBResp + rG*MPBproduction                           ! [/d]

! Total primary production rate              
      MPBproduction = MPBproduction * MPBC                           ! [mmolC/m3solid/d]
      
! ======================================================
! Chlorophyll dynamics
! ======================================================

! Dependency on N:C ratio and open state of PSII
      rhoChl    = QmaxCHLN * PSIIo * DroopNC                             ! [mgChl/mmolN]           
! Dependency on Chl:n ratio   
      DO I = 1, N
        ChlNr = min(ChlNrMPB(I), QmaxCHLN)
        DroopChlN(I) = (1.D0-ChlNr/QmaxCHLN)/(1.D0-ChlNr/QmaxCHLN+0.05) ! [-] 
      ENDDO

! Total Chl production      
      CHLproduction = maxUpNO3 * limNO3 + maxUpNH3 * limNH3              ! [/day]
      CHLproduction = rhoChl * CHLproduction * DroopChlN * MPBC          ! [mgChl/m3solid/d]

! ======================================================
! Carbohydrate exudation 
! ======================================================
      ProdEPS = kEps * MPBproduction
      
! ========================================================
! Algal Death and other loss term
! MPB dies at a constant rate - in the oxic zone: by grazing
! in the anoxic zone additional mortality by oxygen stress
! ===========================================================
      
      MPBCdeath = maxGrazing * MPBC  + deepMort * MPBC *(1.-O2/(O2+1))
      MPBNdeath = maxGrazing * MPBN  + deepMort * MPBN *(1.-O2/(O2+1))
      CHLdeath  = maxGrazing * CHL   + deepMort * CHL * (1.-O2/(O2+1))

      ! the respiration will consume proteins, affecting MPBN and CHL
      ! depending on the NCratio
      
      DO I = 1, N  
         fac = NCrMPB(I)**5.d0/(NCrMPB(I)**5.d0+ QmaxNC**5.d0)
         Nresp(I)    = MPBResp(I)*MPBN(I)*fac
         fac = ChlCrMPB(I)**5.d0/(ChlCrMPB(I)**5.d0+ QmaxChlC**5.d0)
         CHLresp(I)  = MPBResp(I)*CHL(I)*fac
      ENDDO
      
! =======================================================
! Total algal respiration
! =======================================================

      MPBResp = MPBResp*MPBC                    !   [mmolC/m3 solid/d]    B06_18
      
! ==========================================================
! EPS mineralisation
! ==========================================================

      MinEPS = rEPS*EPS  
      
! ==========================================================
! ACCOUNT FOR DECOUPLING OF DIFFERENT COMPARTMENTS
! Stoichiometric fluxes when MPB dies:
! from variable C:N in MPB to constant C:N in Detritus fractions 
! ==========================================================
      DO I = 1, N
       fac = sol2liq(I)

       FDETprodMPB(I)= min(MPBCdeath(I), MPBNdeath(I)/ NCrFdet)          ! mmolC/m3 solid/d
       DICprodMPB(I) = max(0.d0, MPBCdeath(I)-FDETprodMPB(I))*fac        ! mmolC/m3 liquid/d
       NH3prodMPB(I) = max(0.d0,MPBNdeath(I)-FDETprodMPB(I)*NCrFdet)*fac ! mmolN/m3 liquid/d
       NH3prodMPB(I) = NH3prodMPB(I) + Nresp(I)*fac
       PO4prodMPB(I) = NH3prodMPB(I)*PNrMPB                              ! mmolP/m3 liquid/d
      ENDDO 

      
! ==========================================================
! Mass balance equations
! ==========================================================
      
! Microphytobenthos mass balances      
      dMPBC  = dMPBC + MPBproduction - MPBresp - MPBCdeath - ProdEps
      dMPBN  = dMPBN + NH3uptake + NO3uptake - Nresp - MPBNdeath
      dCHL   = dCHL  + CHLproduction - Chldeath - Chlresp
      dEPS   = dEPS  + Prodeps - MinEPS
      
      dPSIIin  = PSDamage - PSRepair               ! NOTE: No transport of PSIIin

      NO3uptake = NO3uptake    *sol2liq            ! convert to nmol/cm3 liquid/d
      NH3uptake = NH3uptake    *sol2liq
      PO4uptake = PO4uptake    *sol2liq
      DICuptake = MPBproduction*sol2liq

      ! Note extra oxygen production related to the reduction of nitrate (2 moles per NO3 uptake)
      ! but NOT when NO3uptake is negative
      
      O2uptake  = MPBresp*O2/(O2+0.1)*sol2liq 
      O2prod    = MPBproduction*sol2liq 
      DO I = 1, N
        O2prod(I) = O2prod(I) + 2.d0*max(0.D0, NO3uptake (I))
      ENDDO

! Updating biogeochemical mass balances      
      dFDET = dFDET+ FDETprodMPB                  ! Mortality MPB
      dO2   = dO2  + O2prod - O2uptake
      dNH3  = dNH3 + (NH3prodMPB - NH3uptake)/(1.d0+NH3Ads)
      dNO3  = dNO3 - NO3uptake
      dDIC  = dDIC + DICprodMPB - DICuptake + (MinEPS+MPBresp)*sol2liq
      dPO4  = dPO4 + PO4prodMPB - PO4uptake

      END SUBROUTINE MPBDIAdynamic
      
! ==========================================================
! ==========================================================
! Integrated quantities (for output)
! ==========================================================
! ==========================================================

      SUBROUTINE MPBDIAout
      USE commonCNPDIA
      USE commonMPBDIA
      IMPLICIT NONE
      
      INTEGER :: I
      DOUBLE PRECISION :: solfac, liqfac      

      totMPBNO3uptake = 0.d0
      totMPBNH3uptake = 0.d0
      totMPBPO4uptake = 0.d0
      totMPBDICuptake = 0.d0
      totMPBprod   = 0.d0
      totMPBresp   = 0.d0
      totMPBCdeath = 0.d0
      totMPBNdeath = 0.d0
      totProdEPS   = 0.d0
      totChlProd   = 0.d0
      totPrimProd  = 0.d0
      totMinEps    = 0.d0
      totFDETprodMPB = 0.d0
      totDICprodMPB  = 0.d0
      totNH3prodMPB  = 0.d0
      totPO4prodMPB  = 0.d0
      totMPBO2prod   = 0.d0
      totMPBO2uptake = 0.d0
      TotalMPBC      = 0.d0
      TotalMPBN      = 0.d0 
      TotalChl       = 0.d0 
      TotalEPS       = 0.d0

      DO I = 1, N
        liqfac =  por(I)      *dx(I)    ! from /cm3 liquid -> cm2 bulk
        solfac = (1.d0-por(I))*dx(I)    ! from /cm3 solid  -> cm2 bulk

        totMPBO2prod    = totMPBO2prod    + O2prod(I)       * liqfac
        totMPBO2uptake  = totMPBO2uptake  + O2uptake(I)     * liqfac

        totMPBNO3uptake = totMPBNO3uptake + NO3uptake(I)    * liqfac
        totMPBNH3uptake = totMPBNH3uptake + NH3uptake(I)    * liqfac
        totMPBPO4uptake = totMPBPO4uptake + PO4uptake(I)    * liqfac
        totMPBDICuptake = totMPBDICuptake + DICuptake(I)    * liqfac
        totDICprodMPB   = totDICprodMPB   + DICprodMPB(I)   * liqfac
        totNH3prodMPB   = totNH3prodMPB   + NH3prodMPB(I)   * liqfac
        totPO4prodMPB   = totPO4prodMPB   + PO4prodMPB(I)   * liqfac
        
        totMPBprod      = totMPBprod      + MPBproduction(I)* solfac
        totMPBresp      = totMPBresp      + MPBresp(I)      * solfac
        totMPBCdeath    = totMPBCdeath    + MPBCdeath(I)    * solfac
        totMPBNdeath    = totMPBNdeath    + MPBNdeath(I)    * solfac
        totProdEPS      = totProdEPS      + ProdEps(I)      * solfac
        totChlProd      = totChlProd      + Chlproduction(I)* solfac
        totPrimProd     = totPrimProd     + PrimProd(I)     * solfac
        totMinEps       = totMinEps       + MinEps(I)       * solfac
        totFDETprodMPB  = totFDETprodMPB  + FDETprodMPB(I)  * solfac
       
        TotalMPBC      = TotalMPBC        + MPBC(I)         * solfac
        TotalMPBN      = TotalMPBN        + MPBN(I)         * solfac
        TotalChl       = TotalChl         + CHL (I)         * solfac
        TotalEPS       = TotalEPS         + EPS (I)         * solfac
      ENDDO
      
        totNH3ads  = totNH3ads  + (totNH3prodMPB - totMPBNH3uptake)*                      &
     &    (1.d0-1.d0/(1.D0+NH3Ads))
      
      END SUBROUTINE MPBDIAout
      
!==========================================================================
! put output variables in one vector
!==========================================================================

      SUBROUTINE getoutMPB(yout)
      USE dimCNPDIA
      
      IMPLICIT NONE
      INTEGER :: i
      DOUBLE PRECISION :: yout(*), MPBout(noutMPB)

      COMMON /myoutMPB   /MPBout

      DO i = 1, noutMPB
       yout(nout+20+i) = MPBout (i)      !nout FESDIA output vars, 20 forcings
      ENDDO       

      END SUBROUTINE getoutMPB
       
