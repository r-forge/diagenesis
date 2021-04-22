!==========================================================================
! THE OMEXDIA model with P, Fe and S, implemented in FORTRAN
!
! Karline Soetaert, nioz-yerseke
!==========================================================================


!==========================================================================
!==========================================================================
! subroutine for calculating the biogeochemical rates of
! microphytobenthos (simple)
!==========================================================================
!==========================================================================
 
      SUBROUTINE MPBFESDIAsimple
      USE commonFESDIA
      IMPLICIT NONE

      DOUBLE PRECISION :: DIN, MPBlim, pNH3(N)
      INTEGER :: I

! --------------------------------------------------------------------------
! Rate of change due to Microphytobenthos
! --------------------------------------------------------------------------

      IF (MPBforc > 0.D0) THEN

       DO I = 1, N
        DIN     = NH3(I) + NO3(I)
        pNH3(I) = NH3(I)/ max(DIN, 1D-10)
        MPBlim  = min (DIN   / (DIN+ksDIN), PO4(I) / (PO4(I)+ksPO4),            &
     &                 DIC(I) / (DIC(I)+ksDIC))
        MPBproduction(I) = max(0.D0, MPBforc * exp(-kMPB*x(I)) * MPBlim)   ! per solid
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
        
! note: dTA = dNH3 - dNO3 - dNO2 - dPO4 - 2*dSO4 

        IF (AddAlk > 0.1) THEN
          dAlk  = dAlk + NO3uptake - NH3uptake + PO4uptake
          AlkProd  = AlkProd + NO3uptake - NH3uptake + PO4uptake 
        ENDIF  

       ELSE
       
        MPBproduction = 0.D0
        pNH3          = 0.D0
        O2prod        = 0.D0
        NO3uptake     = 0.D0
        NH3uptake     = 0.D0
        PO4uptake     = 0.D0
        DICuptake     = 0.D0
       
       ENDIF
      END SUBROUTINE MPBFESDIAsimple
     
