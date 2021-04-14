!==========================================================================
!==========================================================================
! subroutine for calculating the biogeochemical rates of
! the MPBDIA model (omexdia + microphytobenthos + P-dynamics)
!==========================================================================
!==========================================================================
      
      SUBROUTINE MPBDIAsimple 

      USE commonCNPDIA
      IMPLICIT NONE
      
!......................... declaration section.............................
      
      DOUBLE PRECISION :: DIN,MPBlim,pNH3(N)
      
      INTEGER :: I
      
      CHARACTER(len=80) msg
!............................ statements ..................................

! --------------------------------------------------------------------------
! Rate of change due to biogeochemistry 
! --------------------------------------------------------------------------
! MPB production without explicitly describing MPB - no inhibition of NO3 uptake by NH3

      IF (.NOT. MPBdynamic .AND. MPBforc > 0.D0) THEN
        DO I = 1, N
         DIN  = NH3(I) + NO3(I)
         pNH3(I) = NH3(I)/ max(DIN, 1D-10)
       
         MPBlim  = min (DIN   / (DIN+ksNH3), PO4(I) / (PO4(I)+ksPO4),          &
     &                 DIC(I) / (DIC(I)+ksDIC))
         MPBproduction(I) =max(0.D0,MPBforc* exp(-kdLight*x(I))*MPBlim)   ! per solid
        ENDDO   
        
        IF (ratefac .NE. 1) THEN                      ! multiplication factor (eg temperature)
          MPBproduction = MPBproduction*ratefac
        ENDIF
        
        O2prod    = MPBproduction*(1.d0-por)/por      ! change units: /cm3 liquid/d	
        NO3uptake = O2prod*NCrFdet*(1.d0-pNH3)
        NH3uptake = O2prod*NCrFdet*pNH3
        PO4uptake = O2prod*PCrFdet
        DICuptake = O2prod 

        dFDET = dFDET + MPBproduction

      ! Note extra oxygen production related to the reduction of nitrate 
      ! (2 moles per NO3 uptake)
        dO2   = dO2   + O2prod + 2.D0*NO3uptake

        dNH3  = dNH3  - NH3uptake/(1.D0+NH3Ads)
        dNO3  = dNO3  - NO3uptake 

        dDIC  = dDIC - DICuptake

        dPO4  = dPO4 - PO4uptake
      ELSE
        MPBproduction = 0.D0
        pNH3          = 0.D0
        O2prod        = 0.D0
        NO3uptake     = 0.D0
        NH3uptake     = 0.D0
        PO4uptake     = 0.D0
        DICuptake     = 0.D0
      ENDIF

      END SUBROUTINE MPBDIAsimple

!==========================================================================
! Integrated rates for MPB production
!==========================================================================
      
      SUBROUTINE MPBDIAsimpleout
      USE commonCNPDIA
      IMPLICIT NONE
      
      INTEGER :: I
      
! c--------------------------------------------------------------------

      TotMPBO2prod = 0.D0
      TotMPBPO4uptake = 0.D0
      totMPBNH3uptake = 0.D0
      totMPBNO3uptake = 0.D0

      DO I = 1, N
        TotMPBO2prod    = TotMPBO2prod    + O2prod(I)    *por(I)*dx(I)
        TotMPBNO3uptake = TotMPBNO3uptake + NO3uptake(I) *por(I)*dx(I)
        TotMPBNH3uptake = TotMPBNH3uptake + NH3uptake(I) *por(I)*dx(I)
        TotMPBPO4uptake = TotMPBPO4uptake + PO4uptake(I) *por(I)*dx(I)
      ENDDO
      TotMPBDICuptake = TotMPBO2prod
      RETURN
      END SUBROUTINE MPBDIAsimpleout

