
!------------------------------------------------------------------------------------------
! based on code originally created by Karline Soetaert and 
! altered by Andreas Hofmann 
!------------------------------------------------------------------------------------------

!##########################################################################
! Module with constants
!##########################################################################

MODULE Environment
IMPLICIT NONE
  
  ! physical and conversion constants
  double precision, parameter :: R         = 83.14472d0  ! [(bar cm3) / (mol K)] gas constant
  double precision, parameter :: uMolToMol = 1d-6
    
  ! equilibrium constants assumed to be constant (scale and concentration neglected)
  double precision, parameter :: KSiOOH3  = 1.61d-13
  double precision, parameter :: KHNO2    = 1.584893d-3 ! pKHNO2 = 2.8
  double precision, parameter :: KHNO3    = 23.44d0
  double precision, parameter :: KH2SO4   = 1d2
  double precision, parameter :: KHS      = 1.1d-12

  ! physical quantities
  double precision :: Tdeg     ! [dg C]                temperature
  double precision :: S        ! [psu]                 salinity
  double precision :: T        ! [dg K]                absolute temerature
  double precision :: P        ! [bar]                 applied pressure: real pressure -1 atm
  double precision :: IS       ! [mol/(kg-H2O)]        ionic strength
  double precision :: S2       ! -                     squared salinity
  double precision :: SQRTS    ! -                     square root of salinity
  double precision :: SSQRTS   ! -                     S times SQRTS 
  double precision :: P2       ! -                     squared pressure 
  double precision :: LNT      ! -                     natural logarithm of absolute temp
  double precision :: T2       ! -                     squared absolute temperature 
  double precision :: Tdeg2    ! -                     squared temperature in degC
  double precision :: I2       ! -                     squared ionic strength
  double precision :: SQRTI    ! -                     square root of ionic strength
  double precision :: ISQRTI   ! -                     IS times SQRTI 

  double precision :: TA       ! [mol/kg-soln]         total alkalinity
  double precision :: SumB     ! [mol/kg-soln]         lump sum of borate
  double precision :: SumC     ! [mol/kg-soln]         lump sum of CO
  double precision :: SumP     ! [mol/kg-soln]         lump sum of phosphates
  double precision :: SumN     ! [mol/kg-soln]         lump sum of ammonium
  double precision :: SumHS    ! [mol/kg-soln]         lump sum of sulfides
  double precision :: SumS     ! [mol/kg-soln]         lump sum of sulfates
  double precision :: SumF     ! [mol/kg-soln]         lump sum of fluorides
  
  ! factors calculated out of physical quantities stripped of units
  ! conversion summands (without units)
  ! e^(lnK) has property X -> e^(lnK + lnXToY) has property Y

  double precision :: lnMolalToMolin  ! X = [mol/(kg-H2O)], Y = [mol/kg-soln]
  double precision :: lnSWSToFree     ! X = SWS pH        , Y = FREE pH scale
  double precision :: lnTotToFree     ! X = TOTAL pH scale, Y = FREE pH scale

  ! conversion factors (without units) 
  ! K has property X -> K * XToY has property Y
  double precision :: MolalToMolin    ! X = [mol/(kg-H2O)], Y = [mol/kg-soln]
  double precision :: SWSToFree       ! X = SWS pH        , Y = FREE pH scale
  double precision :: TotToFree       ! X = TOTAL pH scale, Y = FREE pH scale

  ! calculated equilibrium constants, converted to FREE pH scale and molinity
  double precision :: KHSO4  ! [mol/kg-soln], FREE pH scale
  double precision :: KHF    ! [mol/kg-soln], FREE pH scale
  double precision :: KH2CO3 ! [mol/kg-soln], FREE pH scale
  double precision :: KHCO3  ! [mol/kg-soln], FREE pH scale
  double precision :: KH2O   ! [mol/kg-soln], FREE pH scale
  double precision :: KBOH3  ! [mol/kg-soln], FREE pH scale
  double precision :: KNH4   ! [mol/kg-soln], FREE pH scale
  double precision :: KH2S   ! [mol/kg-soln], FREE pH scale
  double precision :: KH3PO4 ! [mol/kg-soln], FREE pH scale
  double precision :: KH2PO4 ! [mol/kg-soln], FREE pH scale
  double precision :: KHPO4  ! [mol/kg-soln], FREE pH scale
  double precision :: KSiOH4 ! [mol/kg-soln], FREE pH scale
  double precision :: K0CO2  ! Henrys constant
  double precision :: KspCalcite ! solubility ct
  double precision :: KspAragonite 

  ! contributions to alkalinity
  DOUBLE PRECISION :: Walk, Balk, PO4alk, HSalk, NH3alk, SO4alk, Falk, Sialk

END MODULE Environment

!##########################################################################
! Initialisation
!##########################################################################

SUBROUTINE initialize (Temp, Sal, Depth, Cts)
USE Environment
IMPLICIT NONE
DOUBLE PRECISION :: Temp, Sal   ! [dgC], [-]
DOUBLE PRECISION :: Depth       ! [m]     mean depth
DOUBLE PRECISION :: Cts(12)

  Tdeg   = Temp
  S      = Sal
  Tdeg2  = Tdeg*Tdeg

  T      = Tdeg + 273.15d0 
  T2     = T*T
  LNT    = dlog(T)                         !"dlog" = natural logarithm ("ln") for doubles
  
  P      = 0.1d0 * depth                   ! + 1.01325d0
  P2     = P*P
  
  S2     = S*S
  SQRTS  = dsqrt(S)                        !"dsqrt" = "sqrt" for doubles
  SSQRTS = S * SQRTS
  
  IS     = 19.924d0*S/(1000d0-1.005d0*S)  ! [mol/kg-H20], DOE [9]
  I2     = IS*IS
  SQRTI  = dsqrt(IS)              
  ISQRTI = IS * SQRTI
 
  MolalToMolin   = 1d0 - 0.001005d0*S      ! DOE [9]
  lnMolalToMolin = dlog(MolalToMolin)      ! ln
  
  Cts(1) = T
  Cts(2) = P + 1.01325d0
  Cts(3) = IS
  Cts(4) = MolalToMolin

END SUBROUTINE initialize

!##########################################################################
! Pressure dependence function
!##########################################################################

DOUBLE PRECISION FUNCTION deltaPlnK(a0, a1, a2, b0, b1, b2)
  use Environment
  implicit none
  double precision, intent(in) :: a0, a1, a2, b0, b1, b2
  double precision             :: deltaV, deltaK

  IF (P .EQ. 0) THEN
    deltaPlnK = 0.D0
  ELSE
    deltaV     =  a0 + a1 * Tdeg + a2 * Tdeg2
    deltaK     = (b0 + b1 * Tdeg + b2 * Tdeg2)/1000d0
    deltaPlnK = -(deltaV/(R*T))*P + (0.5D0*(deltaK/(R*T)))*(P2)
  ENDIF
END FUNCTION deltaPlnK

! ##########################################################################
! Dissociation constants
! ##########################################################################

SUBROUTINE Dissociation (KS)

USE environment
IMPLICIT NONE
DOUBLE PRECISION :: KS(20)
DOUBLE PRECISION :: deltaPlnk
DOUBLE PRECISION :: A, B, C, D, E, M
character (len = 200) :: MSG

! Dissociation of HSO4^-
!--------------------------

  A = 324.57d0*SQRTI - 771.54d0*IS + 141.328d0  
  B = 35474d0*IS + 1776d0*I2 - 13856d0*SQRTI - 2698d0*ISQRTI - 4276.1d0
  C = 114.723d0*IS - 47.986d0*SQRTI - 23.093d0

  M = lnMolalToMolin + deltaPlnK(-18.03d0, 0.0466d0, 0.316d-3,-4.53d0, 0.09d0, 0d0)

  KHSO4 = DEXP(A + (B/T) + C*LNT + M)

! Dissociation of HF                                      
!-------------------

  A = 1.525d0 * SQRTI - 12.641d0
  B = 1590.2d0
  M = lnMolalToMolin + deltaPlnK(-9.78d0,-0.009d0,-0.9420d-3,-3.91d0, 0.054d0,0d0)

  KHF = DEXP(A + (B/T) + M)
 
 ! Zeebe/Wolf-Gladrow [25], p.57
  TotToFree    = 1.d0/(1d0 + (SumS/KHSO4))                
  lnTotToFree  = dlog(TotToFree)                         
  SWSToFree    = 1.d0/(1d0 + (SumS/KHSO4) + (SumF/KHF))    
  lnSWSToFree  = dlog(SWSToFree)                         

! "Dissociation" of CO2(aq) via the intermediate H2CO3
!---------------------------------------------------------

  if (S .gt. 5) then    ! two different formulations for high and low salinity
     A = 2.83655d0    - 0.20760841d0*SQRTS + 0.08468345d0*S - 0.00654208d0*SSQRTS
     B = -2307.1266d0 -     4.0484d0*SQRTS
     C = -1.5529413d0
  else
     A = 290.9097d0  -  228.39774d0*SQRTS +  54.20871d0*S -    3.969101*SSQRTS - 0.00258768*S2
     B = -14554.21d0 + 9714.36839d0*SQRTS - 2310.48919d0*S + 170.22169d0*SSQRTS
     C = -45.0575d0  +  34.485796d0*SQRTS -   8.19515d0*S +   0.60367d0*SSQRTS
  end if
  M = lnTotToFree + lnMolalToMolin + deltaPlnK(-25.5d0, 0.1271d0, 0d0,-3.08d0, 0.0877d0,0d0)

  KH2CO3 = DEXP(A + (B/T) + C*LNT + M)

! Dissociation of HCO3^-
!---------------------------------------------------------

  if (S .gt. 5) then    ! two different formulations for high and low salinity
     A = -9.226508d0  - 0.106901773d0*SQRTS + 0.1130822d0*S - 0.00846934d0*SSQRTS 
     B = -3351.6106d0 -     23.9722d0*SQRTS 
     C = -0.2005743d0 
  else
     A = 207.6548d0  -  167.69908d0*SQRTS +   39.75854d0*S -   2.892532d0*SSQRTS - 0.00613142d0*S2
     B = -11843.79d0 + 6551.35253d0*SQRTS - 1566.13883d0*S + 116.270079d0*SSQRTS
     C = -33.6485d0  +  25.928788d0*SQRTS -   6.171951d0*S + 0.45788501d0*SSQRTS
  end if
  M = lnTotToFree + lnMolalToMolin + deltaPlnK(-15.82d0,-0.0219d0, 0.0d0, 1.13d0,-0.1475d0,0d0)

  KHCO3 = DEXP(A + (B/T) + C*LNT + M)

! Dissociation of H20:
!---------------------

  A = 148.9652d0  -  5.977d0*SQRTS - 0.01615d0*S
  B = -13847.26d0 + 118.67d0*SQRTS
  C = -23.6521d0  + 1.0495d0*SQRTS
  M = lnTotToFree + deltaPlnK(-25.60d0, 0.2324d0,-3.6246d-3,-5.13d0, 0.0794d0,0d0)

  KH2O = DEXP(A + (B/T) + C*LNT + M)

! Dissociation of B(OH)3:
!------------------------

  A = 148.0248d0 + 137.1942d0*SQRTS + 1.62142d0*S
  B = -8966.90d0 -  2890.53d0*SQRTS -  77.942d0*S + 1.728d0*SSQRTS - 0.0996d0*S2
  C = -24.4344d0 -     25.085*SQRTS -  0.2474d0*S
  D =              0.053105d0*SQRTS
  M = lnTotToFree + deltaPlnK(-29.48d0, 0.1622d0, 2.608d-3,-2.84d0, 0.0d0, 0d0)

  KBOH3 = DEXP(A + (B/T) + C*LNT + D*T + M)

! Dissociation of NH4^+:
!------------------------

  A = -0.25444d0 +  0.46532d0*SQRTS - 0.01992d0 * S
  B = -6285.33d0 - 123.7184d0*SQRTS + 3.17556d0 * S
  C = 0d0
  D = 0.0001635d0
  E = 0d0
  M = lnSWSToFree + deltaPlnK(-26.43d0, 0.0889d0,-0.9050d-3,-5.03d0, 0.0814d0,0d0)

  KNH4 = DEXP(A + (B/T) + C*LNT + D*T + E*T2 + M)

! Dissociation of H2S:
!---------------------

  A = 225.838d0  + 0.3449d0*SQRTS - 0.0274d0*S
  B = -13275.3d0
  C = -34.6435d0
  M = lnTotToFree + deltaPlnK(-14.80d0, 0.0020d0,-0.4000d-3, 2.89d0, 0.0540d0,0d0)

  KH2S = DEXP(A + (B/T) + C*LNT + M)

! Dissociation of H3PO4:
!-----------------------
  
  A = 115.525d0   + 0.69171d0*SQRTS - 0.01844d0*S
  B = -4576.752d0 - 106.736d0*SQRTS - 0.65643d0*S
  C = -18.453d0
  M = lnTotToFree + deltaPlnK(-14.51d0, 0.1211d0,-0.3210d-3,-2.67d0, 0.0427d0,0d0)

  KH3PO4 = DEXP(A + (B/T) + C*LNT + M)

! Dissociation of H2PO4^-:
!-------------------------
  
  A = 172.0883d0  +  1.3566d0*SQRTS - 0.05778d0*S
  B = -8814.715d0 - 160.340d0*SQRTS + 0.37335d0*S
  C = -27.927d0
  M = lnTotToFree + deltaPlnK(-23.12d0, 0.1758d0,-2.6470d-3,-5.15d0, 0.0900d0,0d0)

  KH2PO4 = DEXP(A + (B/T) + C*LNT + M)

! Dissociation of HPO4^2-:
!-------------------------
  
  A = -18.141d0  +  2.81197d0*SQRTS -  0.09984d0*S
  B = -3070.75d0 + 17.27039d0*SQRTS - 44.99486d0*S
  M = lnTotToFree + deltaPlnK(-26.57d0, 0.2020d0,-3.0420d-3,-4.08d0, 0.0714d0,0d0)

  KHPO4 = DEXP(A + (B/T) + M)

! Henrys constant
!-------------------------

  A = 0.023517d0*S - 167.81077d0
  B = 9345.17d0
  C = 23.3585d0
  D = -2.3656d-4*S
  E = 4.7036d-7*S
  K0CO2 = DEXP(A + (B/T) + C*LNT + D*T + E*T2)
  
! solubility products
!-------------------------
    
  A = -171.9065d0 - 0.77712d0*SQRTS - 0.07711d0*S + 0.0041249d0*SSQRTS
  B = 2839.319d0  +  178.34d0*SQRTS
  C = 71.595d0
  D = -0.077993d0 + 0.0028426d0*SQRTS
  E = 0.d0
  M =  deltaPlnK(-48.76d0, 0.5304d0,0.d0,-11.76d0, 0.3692d0,0d0)
    
  KspCalcite = 10.d0**(A + (B/T) + C*log10(T) + D*T) *DEXP(M)

  A = -171.945d0  - 0.068393d0*SQRTS - 0.10018d0*S + 0.0059415d0*SSQRTS
  B = 2903.293d0  +   88.135d0*SQRTS
  C = 71.595d0
  D = -0.077993d0 + 0.0017276d0*SQRTS
  M =  deltaPlnK(-45.96d0, 0.5304d0,0.d0,-11.76d0, 0.3692d0,0d0)
    
  KspAragonite = 10.d0**(A + (B/T) + C*log10(T) + D*T) *DEXP(M)  



  KS(1) = KH2O   ! [mol/kg-soln], FREE pH scale
  KS(2) = KHF    ! [mol/kg-soln], FREE pH scale
  KS(3) = KH2CO3 ! [mol/kg-soln], FREE pH scale
  KS(4) = KHCO3  ! [mol/kg-soln], FREE pH scale
  KS(5) = KBOH3  ! [mol/kg-soln], FREE pH scale
  KS(6) = KNH4   ! [mol/kg-soln], FREE pH scale
  KS(7) = KH2S   ! [mol/kg-soln], FREE pH scale
  KS(8) = KHS    ! [mol/kg-soln], FREE pH scale
  KS(9) = KH3PO4 ! [mol/kg-soln], FREE pH scale
  KS(10) = KH2PO4 ! [mol/kg-soln], FREE pH scale
  KS(11) = KHPO4  ! [mol/kg-soln], FREE pH scale
  KS(12) = KHSO4  ! [mol/kg-soln], FREE pH scale
  KS(13) = KH2SO4 !
  KS(14) = K0CO2
  KS(15) = KspCalcite
  KS(16) = KspAragonite

END SUBROUTINE Dissociation 

!##########################################################################
!
! Estimate pH from summed concentrations, total alkalinity, 
! temperature and salinity
!
! as implemented by Karline Soetaert, 
! Used definition of alkalinity: (note: assumed that H2SO4 and HNO3 = 0)
! ALK = HCO3- + 2CO3-- + OH- + BOH4- + HPO4-- + 2 PO4.3- + 
!       NH3 + HS- +2S— -H+ - H3PO4 – HSO4- - HF   (-HNO3 - 2H2SO4)

!##########################################################################

SUBROUTINE EstimatePH (Temp, Sal, depth, Alk, SumCO2, SumH3PO4, SumNH4, SumH2S,  &
                       SumH2SO4, maxIter, numvals, pH, Ks, cts)

USE Environment
IMPLICIT NONE

  INTEGER          :: maxIter, numvals
  DOUBLE PRECISION :: Temp, Sal, depth
  DOUBLE PRECISION :: pH(*), Alk(*), SumCO2(*), SumH3PO4(*), SumH2S(*) 
  DOUBLE PRECISION :: SumNH4(*), SumH2SO4(*)
  DOUBLE PRECISION :: pHI, Hguess, H
  DOUBLE PRECISION :: Ks(20), Cts(12)
  LOGICAL          :: finished
  INTEGER          :: I, II
  DOUBLE PRECISION, EXTERNAL  :: EstimateH, RootFun, Zeroin, COAlk, NonCOAlk
  DOUBLE PRECISION, PARAMETER :: tolerance = 1D-12
  CHARACTER (LEN = 120) :: msg
!--------------------------------------------------------------------

  SumB      = 1.1878788D-5 * Sal
  SumF      = 1.952167D-06 * Sal
  CALL initialize (Temp, Sal, depth, Cts)
  PHI = PH(1)
  KS(17) = TotToFree
  KS(18) = SWSToFree
  KS(19) = SumB*1e6
  KS(20) = SumF*1e6
  
  DO I = 1, numvals
    TA        = Alk(I)
    SumS      = SumH2SO4(I)
    SumC      = SumCO2(I)
    SumN      = SumNH4(I)
    SumP      = SumH3PO4(I)
    SumHS     = SumH2S(I)

    CALL dissociation (KS)

! first guess of H+ concentration

    IF (phI .LT. 1 .OR. pHI .GT. 10.5) pHI = 8.D0
    H        = 10.d0**(-phI)
    finished = .False.       ! Will be true when convergence is reached
    II = 0

    DO WHILE (.NOT. Finished)
       II = II+1
       Hguess = EstimateH(H)
       IF (ABS(HGuess - H)*1e12 .LT. Tolerance) Finished = .TRUE.
       IF (II .GE. maxiter)                     Finished = .TRUE.
       IF (Hguess .LE. 0)                       Finished = .TRUE.
       H = Hguess
     ENDDO

     IF (H .GT. 0.d0 .AND. II .LE. maxiter .AND. .NOT. Finished) THEN
       maxiter = II
     ELSE
       H = zeroin(10.D0**(-10), 10.D0**(-2), rootfun, Tolerance)
     ENDIF

       IF (H .GT. 0.d0) THEN
         pH(I)   = -LOG10(H)
       ELSE
         pH(I) = -99.
       ENDIF

!    IF (I .EQ. 1) THEN
!      Cts(5) = TotToFree
!      Cts(6) = SWSToFree
!      Cts(7) = SumB
!      Cts(8) = SumS
!      Cts(9) = SumF
!      Cts(10) = COAlk(H)
!      Cts(11) = NonCOAlk(H)       
!      Cts(12) = rootFun(H)                        
!     ENDIF
   
     pHI   = pH(I)
     
   ENDDO
  END SUBROUTINE

!##########################################################################
!  Netlib root solver
!##########################################################################

  double precision function zeroin(ax,bx,f,tol)
      double precision ax,bx,f,tol

!      a zero of the function  f(x)  is computed in the interval ax,bx .
!
!  input..
!
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result (.ge.0.)
!
!  output..
!
!  zeroin abscissa approximating a zero of  f  in the interval ax,bx
!
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  this is checked, and an error message is printed if this is not
!  satisfied.   zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
!  the  relative machine precision defined as the smallest representable
!  number such that  1.+macheps .gt. 1.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice-hall, inc. (1973).
!
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision  dabs, d1mach

   10 eps = 1.d-12
      tol1 = eps+1.0d0
!


      a=ax
      b=bx
      fa=f(a)
      fb=f(b)
!     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0d0 .or. fb .eq. 0.0d0) go to 20
      if (fa * (fb/dabs(fb)) .le. 0.0d0) go to 20
!         write(6,2500)
!2500     format(1x,'f(ax) and f(bx) do not have different signs,',
!     1             ' zeroin is aborting')
        zeroin = -99.
         return
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if (dabs(fc).ge.dabs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0d0*eps*dabs(b)+0.5d0*tol
      xm = 0.5d0*(c-b)
      if ((dabs(xm).le.tol1).or.(fb.eq.0.0d0)) go to 150
!
! see if a bisection is forced
!
      if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
      if (a.ne.c) go to 60
!
! linear interpolation
!
      p=2.0d0*xm*s
      q=1.0d0-s
      go to 70
!
! inverse quadratic interpolation
!
   60 q=fa/fc
      r=fb/fc
      p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
      q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
   70 if (p.le.0.0d0) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0d0*p).ge.(3.0d0*xm*q-dabs(tol1*q))).or.            &
                   (p.ge.dabs(0.5d0*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if (dabs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0d0) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
  140 fb=f(b)
      if ((fb*(fc/dabs(fc))).gt.0.0d0) go to 20
      go to 30
  150 zeroin=b
      return
      end

!============================================================================
! Alkalinity as a function of H
!============================================================================
      
DOUBLE PRECISION FUNCTION CalcTA (H) 

DOUBLE PRECISION, EXTERNAL :: COAlk, NonCOAlk
DOUBLE PRECISION :: H, TA1, TA2

CHARACTER(LEN = 100) :: MSG

   CalcTA = COAlk(H) + NonCOAlk(H)

END FUNCTION CalcTA

!============================================================================
! Alkalinity component due to CO species 
!============================================================================
      
SUBROUTINE CalcCO3(H, DIC, CO3)
USE Environment
IMPLICIT NONE
DOUBLE PRECISION :: H, DIC, CO3
   
   IF (H > 10) THEN
     CO3 = 0.D0
   ELSE   
     CO3 = KH2CO3*KHCO3/(H*H+KH2CO3*H+KH2CO3*KHCO3)*DIC  
   ENDIF
   
END SUBROUTINE CalcCO3

DOUBLE PRECISION FUNCTION COAlk(H) 

USE Environment
IMPLICIT none 
DOUBLE PRECISION :: H
    !Habi + 2*Abi
    COalk = (KH2CO3*H + 2.D0*KH2CO3*KHCO3)/(H*H+KH2CO3*H+KH2CO3*KHCO3)*SumC 

END FUNCTION COAlk

!============================================================================
! Alkalinity component not due to CO species 
! BOH4 + OH + HPO4 + 2*PO4 + SiOOH3 + HS + 2*S2min + NH3 - H - HSO4 - HF - H3PO4)  
!============================================================================
      
DOUBLE PRECISION FUNCTION NonCOAlk(H) 

USE Environment
IMPLICIT none 
DOUBLE PRECISION :: H, H2, H3, Denom

    H2   = H*H
    H3   = H2*H
    Walk = KH2O  / H                                      !OH
    Balk = KBOH3 / (KBOH3 + H) * SumB                     !BOH  Auni

    ! HPO4 + 2*PO4 - H3PO4                  Hatri , Atr, H3atri
    Denom  = H3 + H2*KH3PO4 + H*KH3PO4*KH2PO4 + KH3PO4*KH2PO4*KHPO4
   
    PO4Alk = (KH3PO4*KH2PO4*H + 2*KH3PO4*KH2PO4*KHPO4 - H3)/Denom *SumP

    Denom = H2 + KH2S*H + KH2S*KHS                        
    HSalk =  (KH2S*H + 2*KH2S*KHS)/Denom *SumHS           !HS- +2S2-  HAbi, Abi

    NH3alk = KNH4/(H + KNH4) *SumN                        !NH3 Auni

    SO4alk = -KH2SO4*H/(H2+KH2SO4*H+KH2SO4*KHSO4)*SumS    !-HSO4   Habi        

    Falk = - H/(H + KHF)*SumF                             !-HF     HAuni

    Sialk = 0.D0 !(No silicium)
    NonCOAlk =  Walk + Balk+ PO4alk + HSAlk + NH3alk +   &
                 SO4alk  + Falk + SiAlk - H  

END FUNCTION NonCOAlk

!============================================================================
! Function to give new value of H 
!============================================================================

DOUBLE PRECISION FUNCTION EstimateH(H) 
USE Environment
IMPLICIT none 
DOUBLE PRECISION :: NonCOalk, H, H2, H3, Denom
DOUBLE PRECISION :: AA, BB, CC, COalkalinity, NonCOalkalinity

!----------------------------------
! 1. Estimate the contribution of the less important species directy

     NONCOAlkalinity = NonCOAlk(H)
     
! We estinate H+ by solving the equation for the DIC species
! carbonate alkalinity
     COAlkalinity = TA - NonCOAlkalinity

! COalkalinity equals HCO3-  +2* CO3--, with
!HCO3 = (KH2CO3*H    /(H2+KH2CO3*H+KH2CO3*KHCO3))*SumCO2
!CO3  = (KH2CO3*KHCO3/(H2+KH2CO3*H+KH2CO3*KHCO3))*SumCO2        

!Calk = HCO3 + 2*CO3
!Calk = (KH2CO3*H + 2.D0*KH2CO3*KHCO3)/(H2+KH2CO3*H+KH2CO3*KHCO3)*SumCO2

! This gives rise to a quadratic equation
     aa = COAlkalinity
     bb = (COAlkalinity - 1.d0*SumC)*KH2CO3
     cc = (COAlkalinity - 2.d0*SumC)*KH2CO3*KHCO3

! new guess of H+
     EstimateH = (-bb + SQRT(BB*BB-4.d0*aa*cc) )/(2.d0*aa)
     RETURN

END FUNCTION EstimateH

DOUBLE PRECISION FUNCTION rootFun(H) 
USE Environment
IMPLICIT NONE
DOUBLE PRECISION :: H
DOUBLE PRECISION, EXTERNAL :: CalcTA

     rootFun = (calcTA(H) - TA)*1.D6  ! in micromol/kg rather than mol/kg

END FUNCTION rootfun



       
     