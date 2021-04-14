!==========================================================================
! THE MPBDIA model, omexdia + microphytobenthos + P-dynamics
! version with dynamic water concentrations
! Karline Soetaert, nioz-yerseke
!==========================================================================

!==========================================================================
! initialisation = same as mpbdiamod - see file mpbdia_P.f
!==========================================================================

!==========================================================================
!==========================================================================
! subroutine calculating the rate of change
! - here the quantities in the common blocks are named
!==========================================================================
!==========================================================================

      SUBROUTINE mpbdiamodbw (neq, t, Conc, dConc, yout, ip)

      USE commonCNPDIA
      USE commonMPBDIA
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i

      DOUBLE PRECISION  :: t,Conc(17*Np1),dConc(17*Np1),yout(*)
      
      DOUBLE PRECISION  :: FdetBW,SdetBW,FePBW,CaPBW,O2BW,NO3bw,NH3bw,          &
     &                   ODUbw,PO4bw,DICbw,PadsBW,NO2BW,                        &
     &                   MPBCbw,MPBNbw,CHLbw,EPSbw,PSIIinbw
      DOUBLE PRECISION  :: dFdetBW,dSdetBW,dO2BW,dNO3BW,dNH3BW,dODUBW,          &
     &                   dPO4BW,dFePBW,dCaPBW,dDICBW,dPadsBW,dNO2BW,            &
     &                   dMPBCbw,dMPBNbw,dCHLbw,dEPSbw,dPSIIinbw

      DOUBLE PRECISION :: pF, FDETdepo, SDETdepo, FePdepo, CaPdepo 
      DOUBLE PRECISION :: MPBCdepo, MPBNdepo, CHLdepo, EPSdepo
      CHARACTER(len=120) msg
!............................ statements ..................................

!     check memory allocated to output variables
      IF (ip(1) < nout)  CALL rexit("nout not large enough") 

! from Conc to fdet, sdet, o2,...
      DO I = 1, N
        Fdet(I) = Conc(      1+I)
        Sdet(I) = Conc(  Np1+1+I)
        O2(I)   = Conc(2*Np1+1+I)
        NO3(I)  = Conc(3*Np1+1+I)
        NO2(I)  = Conc(4*Np1+1+I)
        NH3(I)  = Conc(5*Np1+1+I)
        ODU(I)  = Conc(6*Np1+1+I)
        DIC(I)  = Conc(7*Np1+1+I)
        PO4(I)  = Conc(8*Np1+1+I)
        FeP(I)  = Conc(9*Np1+1+I)
        CaP(I)  = Conc(10*Np1+1+I)
        Pads(I) = Conc(11*Np1+1+I)
        MPBC(I) = Conc(12*Np1+1+I)
        MPBN(I) = Conc(13*Np1+1+I)
        CHL(I)  = Conc(14*Np1+1+I)
        EPS(I)  = Conc(15*Np1+1+I)
        PSIIin(I) = Conc(16*Np1+1+I)
        
      ENDDO
      FdetBW = Conc(      1)
      SdetBW = Conc(  Np1+1)
      O2BW   = Conc(2*Np1+1)
      NO3BW  = Conc(3*Np1+1)
      NO2BW  = Conc(4*Np1+1)
      NH3BW  = Conc(5*Np1+1)
      ODUBW  = Conc(6*Np1+1)
      DICBW  = Conc(7*Np1+1)
      PO4BW  = Conc(8*Np1+1)
      FePBW  = Conc(9*Np1+1)
      CaPBW  = Conc(10*Np1+1)
      PadsBW = Conc(11*Np1+1)
      MPBCBW = Conc(12*Np1+1)
      MPBNBW = Conc(13*Np1+1)
      CHLBW  = Conc(14*Np1+1)
      EPSBW  = Conc(15*Np1+1)
      PSIIinBW = Conc(16*Np1+1)

! --------------------------------------------------------------------------
      FDETFLux = Carbonflux*pFast
      SDETflux = Carbonflux - FDETflux

      FDETdepo  = FdetBW*Cfall
      SDETdepo  = SdetBW*Cfall
      FePdepo   = FePBW*FePfall
      CaPdepo   = CaPBW*CaPfall

      CALL CNPDIAtransolid(FDETdepo,SDETdepo,FePdepo,CaPdepo)

      MPBCdepo  = MPBCBW*MPBfall
      MPBNdepo  = MPBNBW*MPBfall
      CHLdepo   = CHLBW*MPBfall
      EPSdepo   = 0.D0

      CALL MPBDIAtransolid(MPBCdepo,MPBNdepo,CHLdepo,EPSdepo)

      MPBCflux = MPBflux   ! input TO Bottom water
      
! input TO overlying water = parameter, output = estimated based on sinking
      if (Hwater > 0) THEN
        dFDETBW = (FDETflux - FDETdepo)/ Hwater
        dSDETBW = (SDETflux - SDETdepo)/ Hwater
        dFePBW  = (FePflux  - FePdepo )/ Hwater
        dCaPBW  = (CaPflux  - CaPdepo )/ Hwater
        dMPBCBW = (MPBCflux - MPBCdepo)/ Hwater
        dMPBNBW = (MPBCflux*NCrFDET - MPBNdepo)/ Hwater
        dCHLBW  = (MPBCflux*NCrFDET - CHLdepo)/ Hwater
        dEPSBW   = 0.D0
        dPSIIinBW = 0.D0
      ELSE
        dFDETBW = 0.D0
        dSDETBW = 0.D0
        dFePBW  = 0.D0
        dCaPBW  = 0.D0
        dMPBCBW = 0.D0
        dMPBNBW = 0.D0
        dCHLBW  = 0.D0
        dEPSBW    = 0.D0
        dPSIIinBW = 0.D0
      ENDIF
      
      dPadsBW = 0.D0
      CALL CNPDIAtranliquid(O2BW,NO3bw,NH3bw,ODUbw,PO4bw,DICbw,NO2bw)
      
      if (Hwater > 0) THEN
        dO2BW   = (bwO2  -  O2BW) * relax - O2flux  /Hwater
        dNO3BW  = (bwNO3 - NO3BW) * relax - NO3flux /Hwater
        dNO2BW  = (bwNO2 - NO2BW) * relax - NO2flux /Hwater
        dNH3BW  = (bwNH3 - NH3BW) * relax - NH3flux /Hwater
        dODUBW  = (bwODU - ODUBW) * relax - ODUflux /Hwater
        dPO4BW  = (bwPO4 - PO4BW) * relax - PO4flux /Hwater
        dDICBW  = (bwDIC - DICBW) * relax - DICflux /Hwater
      ELSE
        dO2BW   = (bwO2  -  O2BW) * relax
        dNO3BW  = (bwNO3 - NO3BW) * relax
        dNO2BW  = (bwNO2 - NO2BW) * relax
        dNH3BW  = (bwNH3 - NH3BW) * relax
        dODUBW  = (bwODU - ODUBW) * relax
        dPO4BW  = (bwPO4 - PO4BW) * relax
        dDICBW  = (bwDIC - DICBW) * relax
      ENDIF
      
      CALL CNPDIAbiochem 
      CALL MPBDIAdynamic

      CALL CNPDIAout         ! integrated vars
      CALL MPBDIAout     

      CALL getoutP(yout)     ! store output vars in yout
      CALL getoutMPB(yout)   
      
!      write(MSG,'(3(A4, 2(F8.3, 1X)))') "O2 ", O2BW, O2flux, "NO3 ",              &
!     &  NO3BW, NO3flux, "NH3 ", NH3BW, NH3flux
!      CALL rexit(MSG) 
      
! --------------------------------------------------------------------------

! from dfdet, dsdet, do2,... to dconc
      DO I = 1, N
         dConc(      1+I)  =  dFdet(I)
         dConc(  Np1+1+I)  =  dSdet(I) 
         dConc(2*Np1+1+I)  =  dO2(I)  
         dConc(3*Np1+1+I)  =  dNO3(I) 
         dConc(4*Np1+1+I)  =  dNO2(I)
         dConc(5*Np1+1+I)  =  dNH3(I) 
         dConc(6*Np1+1+I)  =  dODU(I)
         dConc(7*Np1+1+I)  =  dDIC(I)
         dConc(8*Np1+1+I)  =  dPO4(I)
         dConc(9*Np1+1+I)  =  dFeP(I)
         dConc(10*Np1+1+I) =  dCaP(I)
         dConc(11*Np1+1+I) =  dPads(I)
         dConc(12*Np1+1+I) =  dMPBC(I)  
         dConc(13*Np1+1+I) =  dMPBN(I)
         dConc(14*Np1+1+I) =  dCHL(I) 
         dConc(15*Np1+1+I) =  dEPS(I)   
         dConc(16*Np1+1+I) =  dPSIIin(I) 
      ENDDO
      dConc(      1)  =  dFdetBW
      dConc(  Np1+1)  =  dSdetBW 
      dConc(2*Np1+1)  =  dO2BW  
      dConc(3*Np1+1)  =  dNO3BW 
      dConc(4*Np1+1)  =  dNO2BW
      dConc(5*Np1+1)  =  dNH3BW 
      dConc(6*Np1+1)  =  dODUBW
      dConc(7*Np1+1)  =  dDICBW
      dConc(8*Np1+1)  =  dPO4BW
      dConc(9*Np1+1)  =  dFePBW
      dConc(10*Np1+1) =  dCaPBW
      dConc(11*Np1+1) =  dPadsBW
      dConc(12*Np1+1) =  dMPBCBW  
      dConc(13*Np1+1) =  dMPBNBW
      dConc(14*Np1+1) =  dCHLBW 
      dConc(15*Np1+1) =  dEPSBW   
      dConc(16*Np1+1) =  dPSIIinBW
 
      RETURN
      END SUBROUTINE
