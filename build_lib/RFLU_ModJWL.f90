










!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
! ******************************************************************************
!
! Purpose: Collection of routines for JWL equation of state. 
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModJWL.F90,v 1.4 2016/02/08 17:00:26 fred Exp $
!
! Copyright: (c) 2006-2007 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModJWL

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input, &
                        t_mixt
  USE ModError
  USE ModMPI
  USE ModTools, ONLY: IsNan
  USE ModInterfaces, ONLY: MixtPerf_C_DGP, &
                           MixtPerf_Eo_DGPVm, &
                           MixtPerf_P_DEoGVm2, &
                           MixtPerf_T_DPR
  ! Josh - 1-Eqtn Model
  USE ModSpecies, ONLY: t_spec_input
  USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex 
  USE SPEC_RFLU_ModPBAUtils, ONLY: SPEC_RFLU_PBA_ComputeY


  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModJWL.F90,v $'
  REAL(RFREAL), PRIVATE :: ATNT,BTNT,CTNT,wTNT,R1TNT,R2TNT,rhoTNT,shcvTNT,ETNT
    !ATNT      = 371.213E9_RFREAL, & ! Pa
    !BTNT      = 3.2306E9_RFREAL, &  ! Pa
    !CTNT      = 1.04527E9_RFREAL, & ! Pa
    !wTNT      = 0.30_RFREAL, &
    !R1TNT     = 4.15_RFREAL, &
    !R2TNT     = 0.95_RFREAL, &
    !rhoTNT    = 1.63E3_RFREAL, &            ! Kg/m3
    !shcvTNT   = 613.4969325153374_RFREAL, & ! 1e6/rhoTNT, J/(K-kg), Cv'=rhoTNT*Cv=1.0D6 (pa/K)
    !kTNT      = 0.297_RFREAL, &
    !muTNT     = 1.74E-6_RFREAL, &         ! NOTE: take air viscosity temporary
    !PrTNT     = 0.4672491789056_RFREAL, & ! shcvTNT*(wTNT+1.0_RFREAL)*muTNT/kTNT
    !MwTNT     = 24.6789_RFREAL, &
    !mlfN2TNT  = 0.2_RFREAL, & ! 1.5_RFREAL/7.5_RFREAL
    !mlfCOTNT  = 0.466666666666667_RFREAL, & ! 3.5_RFREAL/7.5_RFREAL
    !mlfH2OTNT = 0.333333333333333, & ! 2.5_RFREAL/7.5_RFREAL
    !MwAir     = 28.8501_RFREAL, &
    !mlfN2Air  = 0.79_RFREAL, &
    !mlfO2Air  = 0.21_RFREAL
  
! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: &
      RFLU_JWL_C_ER, &
      RFLU_JWL_ComputeEnergyMixt, &
      RFLU_JWL_ComputePressureMixt, &
      RFLU_JWL_E_PR, &
      RFLU_JWL_P_ER, &
      RFLU_JWL_T_PR, &
      RFLU_JWL_Yfix

! ==============================================================================
! Private functions
! ==============================================================================
  
  PRIVATE :: RFLU_JWL_FindInverseMatrix

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  




! ******************************************************************************
!
! Purpose: Compute speed of sound from JWL equation of state.
!
! Description: None.
!
! Input:
!   e             Energy
!   r             Density
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION RFLU_JWL_C_ER(pRegion,e,r,Y)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region),POINTER :: pRegion
    REAL(RFREAL), INTENT(IN) :: e,r,Y
    REAL(RFREAL) :: RFLU_JWL_C_ER
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    
    REAL(RFREAL) :: a,CTNT,w,mp,ma,mb,ra,eJ,An,Bn,p,eamb
    
    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3
    ETNT   = pRegion%mixtInput%prepRealVal13  

! ******************************************************************************
!   Start
! ******************************************************************************

   ! ma = 1.655246553e+05_RFREAL ! (LX-10-1 LLNL 1985 Exp. Handbook)
   ! mb = 3.450701342E+03_RFREAL 
   ! ra = 1.2581927_RFREAL

    ra = pRegion%mixtInput%prepRealVal3  

    !eJ = ETNT/rhoTNT
    eJ = ETNT/1758.000_RFREAL!1860.0000_RFREAL
    eamb = 2.55750E+05_RFREAL 

    ma = (ATNT-0.00_RFREAL)/(eJ-eamb)
    mb = (BTNT-0.00_RFREAL)/(eJ-eamb)

    IF (r < ra) THEN
    mp = 0.0_RFREAL
    ELSE 
  !  mp = -1.0731E-05_RFREAL
    mp = (wTNT-0.4_RFREAL)/(rhoTNT-ra) 
  !  mp = (wTNT-0.4_RFREAL)/(1860.000_RFREAL-ra) 
    END IF
    
    w = max(mp*(r-ra)+0.4_RFREAL,wTNT)

    !IF (r > 1.865E+03_RFREAL) THEN
    !IF (r > 1860.000_RFREAL) THEN
    IF (r > rhoTNT) THEN
    mp = 0.0_RFREAL
    w = wTNT 
    END IF

    IF ( Y .le. 0.99_RFREAL ) THEN
    IF (e .GE. eJ) THEN
     An = ATNT
     Bn = BTNT
     ma = 0.0_RFREAL
     mb = 0.0_RFREAL
    !ELSE IF (e .LE. 2.55750E+05_RFREAL) THEN
    ELSE IF (e .LE. eamb) THEN
     An = 0.0_RFREAL
     Bn = 0.0_RFREAL
     ma = 0.0_RFREAL
     mb = 0.0_RFREAL
    ELSE
     An = ATNT+ma*(e-eJ)
     Bn = BTNT+mb*(e-eJ)
    END IF !Linear fit + constants
    ELSE
     An = ATNT
     Bn = BTNT
     ma = 0.0_RFREAL
     mb = 0.0_RFREAL
     w = wTNT
     mp = 0.00_RFREAL
    ENDIF ! Y .le. 0.99

    p = RFLU_JWL_P_ER(pRegion,e,r,Y)

    CTNT = exp(-R1TNT*rhoTNT/r)*(An*(R1TNT*rhoTNT/r/r - w/r - w/R1TNT/rhoTNT &
                               -r*mp/R1TNT/rhoTNT)+(ma*p/r/r)*(1.0_RFREAL &
                                - (w*r/R1TNT/rhoTNT)))  + &
           exp(-R2TNT*rhoTNT/r)*(Bn*(R2TNT*rhoTNT/r/r - w/r - w/R2TNT/rhoTNT &
                               -r*mp/R2TNT/rhoTNT)+(mb*p/r/r)*(1.0_RFREAL &
                                - (w*r/R2TNT/rhoTNT)))  + & 
           w*(e + p/r) + mp*r*e
  
    a    =  SQRT(CTNT)
 
     IF ( Y .le. 0.01_RFREAL ) THEN

        a = sqrt(1.4_RFREAL*p/r)

     end if 
                               
    RFLU_JWL_C_ER = a

! ******************************************************************************
!   End  
! ******************************************************************************

  END FUNCTION RFLU_JWL_C_ER








! ******************************************************************************
!
! Purpose: Compute speed of sound, energy, and temperature
!          from JWL equation of state
!          for mixture of JWL and perfect gas.
!
! Description: None.
!
! Input:
!   icg           Cell Number 
!   g             Gamma, ratio of specific heats for mixture
!   gc            R, gas constant for mixture
!   p             Pressure
!   r             Density
!   Y             Mass fraction of explosive (JWL) gas
!
! Output:
!   a             Speed of sound
!   e             Energy
!   T             Temperature
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_JWL_ComputeEnergyMixt(pRegion,icg,g,gc,p,r,Y,a,e,T)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region),POINTER :: pRegion
    INTEGER :: icg
    REAL(RFREAL), INTENT(IN) :: g,gc,p,r,Y
    REAL(RFREAL), INTENT(OUT) :: a,e,T
    
! ==============================================================================
!   Locals
! ==============================================================================
    TYPE(t_global),POINTER :: global
    
    global => pRegion%global

    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3
    shcvTNT  = pRegion%mixtInput%prepRealVal20


    IF (IsNaN(p) .EQV. .TRUE. .OR. &
        IsNaN(r) .EQV. .TRUE. ) THEN
      WRITE(*,*) 'Input variables to energy function are NaN'
      WRITE(*,*) 'p,r', p,r,pRegion%irkStep
      CALL ErrorStop(global,ERR_INVALID_VALUE,354,'Invalid quantity in ModJWL')
    END IF

    IF (p .LE. 0.0_RFREAL .OR. &
        r .LE. 0.0_RFREAL ) THEN
      WRITE(*,*) 'Input variables to energy function are Negative'
      WRITE(*,*) 'p,r', p,r,pRegion%irkStep
      CALL ErrorStop(global,ERR_INVALID_VALUE,361,'Invalid quantity in ModJWL')
    END IF

   ! e = RFLU_JWL_E_PR(pRegion,p,r)
   ! T = RFLU_JWL_T_PR(pRegion,p,r,e)
   ! a = RFLU_JWL_C_ER(pRegion,e,r)

    e = RFLU_JWL_E_PR(pRegion,p,r,Y)
    T = RFLU_JWL_T_PR(pRegion,p,r,e,Y)
    a = RFLU_JWL_C_ER(pRegion,e,r,Y)

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_JWL_ComputeEnergyMixt




! ******************************************************************************
!
! Purpose: Compute speed of sound, mass fractions of perfect gas and detonation
!          gas (JWL), pressure, and temperature from JWL equation of state
!          for mixture of JWL and perfect gas.
!
! Description: None.
!
! Input:
!   icg           Cell Number
!   g             Gamma, ratio of specific heats for mixture
!   gc            R, gas constant for mixture
!   e             Energy
!   r             Density
!   Y             Mass fraction of explosive (JWL) gas
!
! Output:
!   a             Speed of sound
!   eJWL          Mass fraction for detonation gas (JWL)
!   ePerf         Mass fraction for perfect gas
!   p             Pressure
!   T             Temperature
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_JWL_ComputePressureMixt(pRegion,icg,g,gc,e,r,Y,a,eJWL,ePerf,p,T)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    INTEGER :: icg
    REAL(RFREAL), INTENT(IN) :: g,gc,e,r,Y
    REAL(RFREAL), INTENT(OUT) :: a,eJWL,ePerf,p,T
    
! ==============================================================================
!   Locals
! ==============================================================================
    TYPE(t_global),POINTER :: global
    global => pRegion%global
    
    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3
    shcvTNT  = pRegion%mixtInput%prepRealVal20
    
    IF (IsNaN(e) .EQV. .TRUE. .OR. &
        IsNaN(r) .EQV. .TRUE. ) THEN
      WRITE(*,*) 'Input variables to pressure function are NaN'
      WRITE(*,*) 'e,r,rkstep,icg,iReg',e,r,pRegion%irkStep,icg,pRegion%iRegionGlobal
      CALL ErrorStop(global,ERR_INVALID_VALUE,443,'Invalid quantity in ModJWL')
    END IF

    IF (e .LE. 0.0_RFREAL .OR. &
        r .LE. 0.0_RFREAL ) THEN
      WRITE(*,*) 'Input variables to pressure function are Negative'
      WRITE(*,*) 'e,r,VF,rkstep,icg,iReg,xcoord,ycoord,time',e,r,pRegion%mixt%piclVF(icg),pRegion%irkStep,icg,pRegion%iRegionGlobal,pRegion%grid%cofg(XCOORD,icg),&
                pRegion%grid%cofg(YCOORD,icg),global%currentTime
      CALL ErrorStop(global,ERR_INVALID_VALUE,451,'Invalid quantity in ModJWL')
    END IF

   ! p = RFLU_JWL_P_ER(pRegion,e,r)
   ! T = RFLU_JWL_T_PR(pRegion,p,r,e)
   ! a = RFLU_JWL_C_ER(pRegion,e,r)
    p = RFLU_JWL_P_ER(pRegion,e,r,Y)
    T = RFLU_JWL_T_PR(pRegion,p,r,e,Y)
    a = RFLU_JWL_C_ER(pRegion,e,r,Y)
    ePerf = e
    eJWL = ePerf
    
! ******************************************************************************
!   End  
! ******************************************************************************

  END SUBROUTINE RFLU_JWL_ComputePressureMixt


! ******************************************************************************
!
! Purpose: Correction for Explosive Mass Fraction for JWL equation of state.
!
! Description: None.
!
! Input:
!   opt           Function where correction is applied - 1=Flux, 2=Dependent Var 
!   c             cell index
!   e             Energy
!   r             Density
!   p             Pressure
!   Y             Explosive Mass Fraction
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_JWL_Yfix(pRegion,c,opt,ri,ei,pi,Yi,rf,ef,pf,Yf)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region),POINTER :: pRegion
    INTEGER, INTENT(IN) :: opt,c
    REAL(RFREAL), INTENT(IN) :: ei,ri,pi,Yi
    REAL(RFREAL), INTENT(OUT) :: ef,rf,pf,Yf

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvOld
! ******************************************************************************
!   Start
! ******************************************************************************

   pCv    => pRegion%mixt%cv 
   pCvOld => pRegion%mixt%cvOld

   ef = ei
   rf = ri
   pf = pi
   Yf = Yi !Default is to do nothing to values if Y is acceptable

   IF (opt==1) THEN
    IF ((IsNan(pi) .EQV. .TRUE.) .OR. (IsNan(ri) .EQV. .TRUE.) .OR. &
             (IsNan(Yi) .EQV. .TRUE.)) THEN
            rf = pCvOld(CV_MIXT_DENS,c)
            pf = rf*0.4_RFREAL*2.5E5_RFREAL
            Yf = 0.0_RFREAL
    END IF

   ELSEIF (opt==2) THEN

 
   IF(IsNan(Yi) .EQV. .TRUE. .OR. Yi .LE. 0.0_RFREAL) THEN
     Yf = 0.0_RFREAL
   END IF

   IF (Yi .GE. 1.0_RFREAL) THEN
     Yf = 1.0_RFREAL
   END IF !Fred - fixing Y update issue where the numerical solver pushes Y
          !slightly past its limit points


  END IF 
! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_JWL_Yfix






! ******************************************************************************
!
! Purpose: Compute internal energy from JWL equation of state.
!
! Description: None.
!
! Input:
!   p             Pressure
!   r             Density
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION RFLU_JWL_E_PR(pRegion,p,r,Y)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region),POINTER :: pRegion
    REAL(RFREAL), INTENT(IN) :: p,r,Y
    REAL(RFREAL) :: RFLU_JWL_E_PR
     
! ==============================================================================
!   Locals
! ==============================================================================
        
    ! Josh
     
  CHARACTER(LEN=25) :: str2 
  LOGICAL :: I_EXIST
  REAL(RFREAL) :: Ynot
    
    TYPE(t_global),POINTER :: global
    REAL(RFREAL) :: e,mp,ma,mb,ra,eJ,w,C1,C2,An,Bn,eamb
    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3
    ETNT   = pRegion%mixtInput%prepRealVal13  
    


    global => pRegion%global
! ******************************************************************************
!   Start
! ******************************************************************************

 
  !  ma = 1.655246553e+05_RFREAL
  !  mb = 3.4507013E+03_RFREAL
  !  ra = 1.258_RFREAL
    Ynot = Y 
    ra = pRegion%mixtInput%prepRealVal3  

   ! eJ = ETNT/rhoTNT
    eJ = ETNT/1758.000_RFREAL!1860.000_RFREAL
    eamb = 2.55750E+05_RFREAL 

    ma = (ATNT-0.00_RFREAL)/(eJ-eamb)
    mb = (BTNT-0.00_RFREAL)/(eJ-eamb)


    IF (r < ra) THEN
    mp = 0.0_RFREAL
    ELSE 
  !  mp = -1.0731E-05_RFREAL
  !  mp = (wTNT-0.4_RFREAL)/(1860.0000_RFREAL-ra) 
    mp = (wTNT-0.4_RFREAL)/(rhoTNT-ra) 
    END IF
    
    w = max(mp*(r-ra)+0.4_RFREAL,wTNT)

    !IF (r > 1.865E+03_RFREAL) THEN
    !IF (r > 1860.000_RFREAL) THEN
    IF (r > rhoTNT) THEN
    mp = 0.0_RFREAL
    w = wTNT 
    END IF
    
        ! Josh - Write Slopes mA,mB,mW
   !   write (str2, "(A14)") "JWL_slopes.dat"
   !     str2 = trim(str2)
   !    INQUIRE (FILE=str2, EXIST=I_EXIST)
   !     IF ( I_EXIST ) THEN
   !     ELSE
   !         OPEN(1117,file=str2,form='formatted',status='unknown')
   !         WRITE(1117,*) 'ma,mb,mW',ma,mb,mp
   !     ENDIF

!    IF (r < 1.258_RFREAL) THEN
!    mp = 0.0_RFREAL
!    ELSE
!    mp = -1.0731E-05_RFREAL
!    END IF
!
!    eJ = 5.576E+06_RFREAL
!
!    IF (r .LE. 1865_RFREAL) THEN
!    w = max(mp*(r-ra)+0.4_RFREAL,wTNT)
!    ELSE
!    w = wTNT
!    END IF
120 CONTINUE
    ! Josh
    !IF ( Y .le. 0.99_RFREAL ) THEN
    IF ( Ynot .le. 0.99_RFREAL ) THEN
    C1 = (1.0_RFREAL - ( (w*r)/(R1TNT*rhoTNT) ))*exp(-R1TNT*rhoTNT/r)
    C2 = (1.0_RFREAL - ((w*r)/(R2TNT*rhoTNT)))*exp(-R2TNT*rhoTNT/r)

    e = ( p + (ma*eJ - ATNT)*C1 + (mb*eJ - BTNT)*C2)/(ma*C1 + mb*C2 + w*r)
    !First attempt

    !IF (e .LE. eJ .AND. e .GE. 2.5575E+05) THEN
    IF (e .LE. eJ .AND. e .GE. eamb) THEN
!     RFLU_JWL_E_PR = e
    ELSE
     An = ATNT
     Bn = BTNT
    
     e = (p/r/w) - (An*((1.0_RFREAL/(w*r))-(1.0_RFREAL/rhoTNT/R1TNT))*exp(-R1TNT*rhoTNT/r)) &
                 - (Bn*((1.0_RFREAL/(w*r))-(1.0_RFREAL/rhoTNT/R2TNT))*exp(-R2TNT*rhoTNT/r))
     IF (e .LT. eJ) THEN
      e = p/w/r
     END IF ! e .LE. eJ Correction step in case assumption that e is in linear range is wrong
    END IF ! e .LE. eJ .AND. e .GE. eAir
    ELSE
     An = ATNT
     Bn = BTNT
     w = wTNT
    
     e = (p/r/w) - (An*((1.0_RFREAL/(w*r))-(1.0_RFREAL/rhoTNT/R1TNT))*exp(-R1TNT*rhoTNT/r)) &
                 - (Bn*((1.0_RFREAL/(w*r))-(1.0_RFREAL/rhoTNT/R2TNT))*exp(-R2TNT*rhoTNT/r))
    ENDIF ! Y .le. 0.99
    IF (IsNaN(e) .EQV. .TRUE.) THEN
      WRITE(*,*) 'Energy in RFLU_JWL_E_PR function is NaN'
      WRITE(*,*) 'e,r,p', e,r,p,pRegion%irkStep
      CALL ErrorStop(global,ERR_INVALID_VALUE,705,'Invalid quantity in ModJWL')
    END IF


    IF ( e .LE. 0.0_RFREAL ) THEN
                IF ( Ynot .gt. 0.99_RFREAL ) THEN
                      !Ynot = Y
                      Ynot = 0.981230_RFREAL
                      GOTO 120
                ELSE
      WRITE(*,*) 'Energy in RFLU_JWL_E_PR function is negative'
      WRITE(*,*) 'e,r,p', e,r,p,pRegion%irkStep
      CALL ErrorStop(global,ERR_INVALID_VALUE,717,'Invalid quantity in ModJWL')
                 ENDIF
    ENDIF

     IF ( Ynot .le. 0.01_RFREAL ) THEN

        e = p/r/(0.4_RFREAL)

     end if   

    RFLU_JWL_E_PR = e

! ******************************************************************************
!   End  
! ******************************************************************************

  END FUNCTION RFLU_JWL_E_PR








! ###############################################
! TO DO: Manoj-JWL, Clean up following code    
!Subroutine to find the inverse of a square matrix
!Author : Yoshifumi Nozaki

  SUBROUTINE RFLU_JWL_FindInverseMatrix(matrix,inverse,n,nmax)
  implicit none
        
!---Declarations
        INTEGER, INTENT(IN) :: n,nmax
        double precision, INTENT(IN), DIMENSION(n,n) :: matrix  !Input A matrix
        double precision, INTENT(OUT), DIMENSION(nmax,nmax) :: inverse !Inverted matrix
        
        INTEGER :: i, j, k, l
        double precision :: m
        double precision, DIMENSION(n,2*n) :: augmatrix !augmented matrix

    
        !Augment input matrix with an identity matrix
        DO i = 1,n
          DO j = 1,2*n
            IF (j <= n ) THEN
              augmatrix(i,j) = matrix(i,j)
            ELSE IF ((i+n) == j) THEN
              augmatrix(i,j) = 1.0d0
            Else
              augmatrix(i,j) = 0.0d0
            ENDIF
          END DO
        END DO   
		
        !Ensure diagonal elements are non-zero
        DO k = 1,n-1
          DO j = k+1,n
            IF (augmatrix(k,k) == 0) THEN
!			  write(*,*) 'There exists diagonal elements'
               DO i = k+1, n
                 IF (augmatrix(i,k) /= 0) THEN
                   DO  l = 1, 2*n
                     augmatrix(k,l) = augmatrix(k,l)+augmatrix(i,l)
                   END DO
                 ENDIF
               END DO
            ENDIF
          END DO
        END DO   
		
        !Reduce augmented matrix to upper traingular form
        DO k =1, n-1
          DO j = k+1, n   
            m = augmatrix(j,k)/augmatrix(k,k)
			!write(*,*) k, j, m
            DO i = k, 2*n
              augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
            END DO
          END DO
        END DO
        
        !Test for invertibility
        DO i = 1, n
          IF (augmatrix(i,i) == 0) THEN
            write(*,*) "ERROR-Matrix is non-invertible"
            inverse = 0.d0
            WRITE(11,'(4(E23.16,1X))') matrix(1,1:4)
            WRITE(11,'(4(E23.16,1X))') Matrix(2,1:4)
            WRITE(11,'(4(E23.16,1X))') Matrix(3,1:4)
            WRITE(11,'(4(E23.16,1X))') Matrix(4,1:4)
            STOP 
            return
          ENDIF
        END DO
       		
        !Make diagonal elements as 1
        DO i = 1 , n
          m = augmatrix(i,i)
          DO j = i , (2 * n)                                
            augmatrix(i,j) = (augmatrix(i,j) / m)
          END DO
        END DO
		
        !Reduced right side half of augmented matrix to identity matrix
        DO k = n-1, 1, -1
          DO i =1, k
            m = augmatrix(i,k+1)
            DO j = k, (2*n)
              augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
            END DO
          END DO
        END DO        
		
		!store answer
        DO i =1, n
          DO j = 1, n
            inverse(i,j) = augmatrix(i,j+n)
          END DO
        END DO

        !do i=1,n
	    !write(26,99) augmatrix (i,1:2*n)
	    !end do
		
!99  format(100(1X,E30.15)) 

  END SUBROUTINE RFLU_JWL_FindInverseMatrix 
! ###############################################








! ******************************************************************************
!
! Purpose: Compute pressure from JWL equation of state.
!
! Description: None.
!
! Input:
!   e             Energy
!   r             Density
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION RFLU_JWL_P_ER(pRegion,e,r,Y)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    REAL(RFREAL), INTENT(IN) :: e,r,Y
    REAL(RFREAL) :: RFLU_JWL_P_ER
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    REAL(RFREAL) :: p,mp,ma,mb,ra,eJ,w,An,Bn,eamb
    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3
    ETNT   = pRegion%mixtInput%prepRealVal13  
    
! ******************************************************************************
!   Start
! ******************************************************************************

    ra = pRegion%mixtInput%prepRealVal3  

   ! eJ = ETNT/rhoTNT
    eJ = ETNT/1758.000_RFREAL!1860.000_RFREAL
    eamb = 2.55750E+05_RFREAL 

    ma = (ATNT-0.00_RFREAL)/(eJ-eamb)
    mb = (BTNT-0.00_RFREAL)/(eJ-eamb)

    IF (r < ra) THEN
    mp = 0.0_RFREAL
    ELSE 
  !  mp = -1.0731E-05_RFREAL
  !  mp = (wTNT-0.4_RFREAL)/(1860.000_RFREAL-ra) 
    mp = (wTNT-0.4_RFREAL)/(rhoTNT-ra) 
    END IF
    
    w = max(mp*(r-ra)+0.4_RFREAL,wTNT)

    !IF (r > 1.865E+03_RFREAL) THEN
    !IF (r > 1860.000_RFREAL) THEN
    IF (r > rhoTNT) THEN
    mp = 0.0_RFREAL
    w = wTNT 
    END IF
    
    IF ( Y .le. 0.99_RFREAL ) THEN
    IF (e .GE. eJ) THEN
     An = ATNT
     Bn = BTNT
     ma = 0.0_RFREAL
     mb = 0.0_RFREAL
    !ELSE IF (e .LE. 2.55750E+05_RFREAL) THEN
    ELSE IF (e .LE. eamb) THEN
     An = 0.0_RFREAL
     Bn = 0.0_RFREAL
     ma = 0.0_RFREAL
     mb = 0.0_RFREAL
    ELSE
     An = ATNT+ma*(e-eJ)
     Bn = BTNT+mb*(e-eJ)
    END IF !Linear fit + constants
    ELSE
     An = ATNT
     Bn = BTNT
     ma = 0.0_RFREAL
     mb = 0.0_RFREAL
     w = wTNT
    ENDIF ! Y .le. 0.99 
    
    p = An*(1.0_RFREAL - (w*r)/(R1TNT*rhoTNT))*exp(-R1TNT*rhoTNT/r) + &
        Bn*(1.0_RFREAL - (w*r)/(R2TNT*rhoTNT))*exp(-R2TNT*rhoTNT/r) + &
        w*e*r
  
     IF ( Y .le. 0.01_RFREAL ) THEN

        p = e*r*(0.4_RFREAL)

     end if   




    RFLU_JWL_P_ER = p

! ******************************************************************************
!   End  
! ******************************************************************************

  END FUNCTION RFLU_JWL_P_ER


! ******************************************************************************
!
! Purpose: Compute temperature from JWL equation of state.
!
! Description: None.
!
! Input:
!   p             Pressure
!   r             Density
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION RFLU_JWL_T_PR(pRegion,p,r,e,Y)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region),POINTER :: pRegion
    REAL(RFREAL), INTENT(IN) :: p,r,e,Y
    REAL(RFREAL) :: RFLU_JWL_T_PR
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    REAL(RFREAL) :: T,mp,mcv,ma,mb,ra,eJ,w,cv,An,Bn,eamb,shcvair
    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3
    shcvTNT  = pRegion%mixtInput%prepRealVal20
    ETNT   = pRegion%mixtInput%prepRealVal13  
       
 
! ******************************************************************************
!   Start
! ******************************************************************************

    ra = pRegion%mixtInput%prepRealVal3  

   ! eJ = ETNT/rhoTNT
    eJ = ETNT/1758.000_RFREAL!1860.000_RFREAL
    eamb = 2.55750E+05_RFREAL
    shcvair = 717.00_RFREAL 

    ma = (ATNT-0.00_RFREAL)/(eJ-eamb)
    mb = (BTNT-0.00_RFREAL)/(eJ-eamb)

    IF (r < ra) THEN
    mp = 0.0_RFREAL
    mcv = 0.0_RFREAL
    ELSE 
  !  mp = -1.0731E-05_RFREAL
  !  mcv = -0.10731802_RFREAL
   ! mp = (wTNT-0.4_RFREAL)/(1860.000_RFREAL-ra)
    mp = (wTNT-0.4_RFREAL)/(rhoTNT-ra)
   ! mcv = (shcvTNT-shcvair)/(1860.000_RFREAL-ra) 
    mcv = (shcvTNT-shcvair)/(rhoTNT-ra) 
    END IF
    
    w = max(mp*(r-ra)+0.4_RFREAL,wTNT)
    cv = max(mcv*(r-ra)+717.0_RFREAL,shcvTNT)

    !IF (r > 1.865E+03_RFREAL) THEN
    !IF (r > 1860.000_RFREAL) THEN
    IF (r > rhoTNT) THEN
    mp = 0.0_RFREAL
    mcv = 0.00_RFREAL
    w = wTNT
    cv = shcvTNT 
    END IF
    
    IF ( Y .le. 0.99_RFREAL ) THEN
    IF (e .GE. eJ) THEN
     An = ATNT
     Bn = BTNT
     ma = 0.0_RFREAL
     mb = 0.0_RFREAL
    !ELSE IF (e .LE. 2.55750E+05_RFREAL) THEN
    ELSE IF (e .LE. eamb) THEN
     An = 0.0_RFREAL
     Bn = 0.0_RFREAL
     ma = 0.0_RFREAL
     mb = 0.0_RFREAL
    ELSE
     An = ATNT+ma*(e-eJ)
     Bn = BTNT+mb*(e-eJ)
    END IF !Linear fit + constants
    ELSE
     An = ATNT
     Bn = BTNT
     ma = 0.0_RFREAL
     mb = 0.0_RFREAL
     mp = 0.0_RFREAL
     mcv = 0.00_RFREAL
     w = wTNT
     cv = shcvTNT
    ENDIF ! Y .le. 0.99
    
    T = (p - An*EXP(-R1TNT*rhoTNT/r)  &
           - Bn*EXP(-R2TNT*rhoTNT/r))*1.0_RFREAL/w/cv/r    

    IF ( Y .le. 0.01_RFREAL ) THEN

        T = p/r/(287.04_RFREAL)

     end if 


    RFLU_JWL_T_PR = T

! ******************************************************************************
!   End  
! ******************************************************************************

  END FUNCTION RFLU_JWL_T_PR

!BRAD ADDED A JWL JUNCTION

! ******************************************************************************
!
! Purpose: Compute Pressure from JWL equation of state.
!
! Description: None.
!
! Input:
!   T             Temperature
!   r             Density
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION RFLU_JWL_P_TR(pRegion,t,r,e,Y)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region),POINTER :: pRegion
    REAL(RFREAL), INTENT(IN) :: t,r,e,Y
    REAL(RFREAL) :: RFLU_JWL_P_TR
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    REAL(RFREAL) :: P,mp,mcv,ma,mb,ra,eJ,w,cv,An,Bn,eamb,shcvair
    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3
    shcvTNT  = pRegion%mixtInput%prepRealVal20
    ETNT   = pRegion%mixtInput%prepRealVal13  
       
 
! ******************************************************************************
!   Start
! ******************************************************************************

    ra = pRegion%mixtInput%prepRealVal3  

   ! eJ = ETNT/rhoTNT
    eJ = ETNT/1758.000_RFREAL!1860.000_RFREAL
    eamb = 2.55750E+05_RFREAL
    shcvair = 717.00_RFREAL 

    ma = (ATNT-0.00_RFREAL)/(eJ-eamb)
    mb = (BTNT-0.00_RFREAL)/(eJ-eamb)

    IF (r < ra) THEN
    mp = 0.0_RFREAL
    mcv = 0.0_RFREAL
    ELSE 
  !  mp = -1.0731E-05_RFREAL
  !  mcv = -0.10731802_RFREAL
   ! mp = (wTNT-0.4_RFREAL)/(1860.000_RFREAL-ra)
    mp = (wTNT-0.4_RFREAL)/(rhoTNT-ra)
   ! mcv = (shcvTNT-shcvair)/(1860.000_RFREAL-ra) 
    mcv = (shcvTNT-shcvair)/(rhoTNT-ra) 
    END IF
    
    w = max(mp*(r-ra)+0.4_RFREAL,wTNT)
    cv = max(mcv*(r-ra)+717.0_RFREAL,shcvTNT)

    !IF (r > 1.865E+03_RFREAL) THEN
    !IF (r > 1860.000_RFREAL) THEN
    IF (r > rhoTNT) THEN
    mp = 0.0_RFREAL
    mcv = 0.00_RFREAL
    w = wTNT
    cv = shcvTNT 
    END IF
    
    IF ( Y .le. 0.99_RFREAL ) THEN
    IF (e .GE. eJ) THEN
     An = ATNT
     Bn = BTNT
     ma = 0.0_RFREAL
     mb = 0.0_RFREAL
    !ELSE IF (e .LE. 2.55750E+05_RFREAL) THEN
    ELSE IF (e .LE. eamb) THEN
     An = 0.0_RFREAL
     Bn = 0.0_RFREAL
     ma = 0.0_RFREAL
     mb = 0.0_RFREAL
    ELSE
     An = ATNT+ma*(e-eJ)
     Bn = BTNT+mb*(e-eJ)
    END IF !Linear fit + constants
    ELSE
     An = ATNT
     Bn = BTNT
     ma = 0.0_RFREAL
     mb = 0.0_RFREAL
     mp = 0.0_RFREAL
     mcv = 0.00_RFREAL
     w = wTNT
     cv = shcvTNT
    ENDIF ! Y .le. 0.99
    
    !T = (p - An*EXP(-R1TNT*rhoTNT/r) - Bn*EXP(-R2TNT*rhoTNT/r))*1.0_RFREAL/w/cv/r    

    P = (T * w*cv*r) + An*EXP(-R1TNT*rhoTNT/r) + Bn*EXP(-R2TNT*rhoTNT/r)  

    RFLU_JWL_P_TR = P

! ******************************************************************************
!   End  
! ******************************************************************************

  END FUNCTION RFLU_JWL_P_TR



END MODULE RFLU_ModJWL

! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
!
! ******************************************************************************

