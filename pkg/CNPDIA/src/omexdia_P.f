!==========================================================================
! THE OMEXDIA model with P, implemented in FORTRAN
!
! Karline Soetaert, nioz-yerseke
! Modules and transport functions are in file common.f
!==========================================================================

!==========================================================================
! initialise the common block with parameter values, 
! followed by thicknesses, surface, porosities, bioturbation values
!==========================================================================

      SUBROUTINE initomexdiaP (steadyparms)
      USE dimCNPDIA
      
      IMPLICIT NONE
      EXTERNAL steadyparms

      DOUBLE PRECISION parms(nparms)
      COMMON /myparmsP/parms

      DOUBLE PRECISION zero(N)
      COMMON /myzero/zero

       Zero(:) = 0.D0
       CALL steadyparms(nparms, parms)
       MPBdynamic = .FALSE.
       
      RETURN
      END SUBROUTINE

!==========================================================================
! Initialise the forcing function common block (Carbon flux)
!==========================================================================

      SUBROUTINE initforcP (steadyforcs)
      IMPLICIT NONE
      EXTERNAL steadyforcs

      INTEGER,PARAMETER :: N=20

      DOUBLE PRECISION forcs(N)
      COMMON /myforcsP/forcs

       CALL steadyforcs(N, forcs)
       
      RETURN
      END SUBROUTINE


!==========================================================================
!==========================================================================
! subroutine calculating the rate of change of
! the omexdia model - here the quantities in the common blocks are named
!==========================================================================
!==========================================================================
      
      SUBROUTINE omexdiamodP (neq, t, Conc, dConc, yout, ip)

      USE commonCNPDIA
      IMPLICIT NONE
      
!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i

      DOUBLE PRECISION  :: t, Conc(12*N), dConc(12*N), yout(*)

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
      ENDDO
      
! --------------------------------------------------------------------------
      FDETFLux = Carbonflux*pFast
      SDETflux = Carbonflux - FDETflux
      
      CALL CNPDIAtransolid(FDETFLux,SDETFLux,FePflux,CaPflux)

      CALL CNPDIAtranliquid(bwO2,bwNO3,bwNH3,bwODU,bwPO4,bwDIC,bwNO2)

      CALL CNPDIAbiochem 
      CALL MPBDIAsimple 

      CALL MPBDIAsimpleout     
      CALL CNPDIAout     
      CALL getoutP(yout)

! --------------------------------------------------------------------------

! from dfdet, dsdet, do2,... to dconc

      DO I = 1, N
         dConc(I)      =  dFdet(I)
         dConc(N+I)    =  dSdet(I) 
         dConc(2*N+I)  =  dO2(I)  
         dConc(3*N+I)  =  dNO3(I) 
         dConc(4*N+I)  =  dNO2(I)
         dConc(5*N+I)  =  dNH3(I) 
         dConc(6*N+I)  =  dODU(I)
         dConc(7*N+I)  =  dDIC(I)
         dConc(8*N+I)  =  dPO4(I)
         dConc(9*N+I)  =  dFeP(I)
         dConc(10*N+I) =  dCaP(I)
         dConc(11*N+I) =  dPads(I)
      ENDDO

      RETURN
      END SUBROUTINE omexdiamodP

