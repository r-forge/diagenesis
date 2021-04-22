!==========================================================================
! THE OMEXDIA model with P, Fe and S, implemented in FORTRAN
!
! Karline Soetaert, nioz-yerseke
! Modules and transport functions are in file FEScommon.f
!==========================================================================


!==========================================================================
! initialise the common block with parameter values, 
! followed by thicknesses, surface, porosities, bioturbation values
!==========================================================================
      SUBROUTINE initfesdia (steadyparms)
      USE dimFESDIA
      IMPLICIT NONE
      EXTERNAL steadyparms

      DOUBLE PRECISION parms(nparms)
      COMMON /myparmsfes/parms

      DOUBLE PRECISION zero(N)
      COMMON /myzero/zero

        Zero(:) = 0.D0

        CALL steadyparms(nparms, parms)
        DynamicpH = .FALSE.

       
      RETURN
      END SUBROUTINE

!==========================================================================
! Initialise the forcing function common block (Carbon flux)
!==========================================================================

      SUBROUTINE initfesforc (steadyforcs)
      USE dimFESDIA
      IMPLICIT NONE
      EXTERNAL steadyforcs

      DOUBLE PRECISION forcs(nforc)
      COMMON /myforcsfes/forcs

       CALL steadyforcs(nforc, forcs)
       
      RETURN
      END SUBROUTINE

!==========================================================================
!==========================================================================
! subroutine calculating the rate of change of
! the fesdia model - here the quantities in the common blocks are named
!==========================================================================
!==========================================================================

      SUBROUTINE fesdiamod (neq, t, Conc, dConc, yout, ip)
      use commonFESDIA
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i

      DOUBLE PRECISION  :: t,Conc(17*N),dConc(17*N),yout(*)

      CHARACTER(len=80) msg

!............................ statements ..................................

!     check memory allocated to output variables
       IF (ip(1)< noutdia+nforcsdia) THEN
         CALL rexit("nout should be at least nout")
       ENDIF  

! from Conc to fdet, sdet, o2,...
       DO I = 1, N
        Fdet(I) = Conc(I)
        Sdet(I) = Conc(N+I)
        O2(I)   = Conc(2*N+I)
        NO3(I)  = Conc(3*N+I)
        NO2(I)  = Conc(4*N+I)
        NH3(I)  = Conc(5*N+I)
        DIC(I)  = Conc(6*N+I)
        Fe(I)   = Conc(7*N+I)
        FeOH3(I)= Conc(8*N+I)
        H2S(I)  = Conc(9*N+I)
        SO4(I)  = Conc(10*N+I)
        CH4(I)  = Conc(11*N+I)
        PO4(I)  = Conc(12*N+I)
        FeP(I)  = Conc(13*N+I)
        CaP(I)  = Conc(14*N+I)
        Pads(I) = Conc(15*N+I)
        ALK(I)  = Conc(16*N+I)
       ENDDO
      
! --------------------------------------------------------------------------
       FDETFLux = Carbonflux*pFast
       SDETflux = Carbonflux - FDETflux
      
       CALL FESDIAtransolid(FDETFLux,SDETFLux,0.D0,0.D0,                        &
     &                      FeOH3flux,CaPflux)

       CALL FESDIAtranliquid(bwO2, bwNO3, bwNO2, bwNH3, bwCH4, bwPO4,           &
     &                      bwFe, bwH2S, bwSO4, bwDIC, bwALK)

       CALL FESDIAbiochem 
       
       CALL MPBFESDIAsimple 

       CALL FESDIAout     (yout)

! from dfdet, dsdet, do2,... to dconc
       DO I = 1, N
         dConc(I)      =  dFdet(I)
         dConc(N+I)    =  dSdet(I) 
         dConc(2*N+I)  =  dO2(I)  
         dConc(3*N+I)  =  dNO3(I) 
         dConc(4*N+I)  =  dNO2(I) 
         dConc(5*N+I)  =  dNH3(I) 
         dConc(6*N+I)  =  dDIC(I)
         dConc(7*N+I)  =  dFe(I)
         dConc(8*N+I)  =  dFeOH3(I)
         dConc(9*N+I)  =  dH2S(I)
         dConc(10*N+I) =  dSO4(I)
         dConc(11*N+I) =  dCH4(I)
         dConc(12*N+I) =  dPO4(I)
         dConc(13*N+I) =  dFeP(I)
         dConc(14*N+I) =  dCaP(I)
         dConc(15*N+I) =  dPADS(I)
         dConc(16*N+I) =  dALK(I)
       ENDDO 
      END SUBROUTINE fesdiamod
