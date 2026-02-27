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
!******************************************************************************
!
! Purpose: 
!
! Description: none.
!
! Input: 
!
! Output:
!
! Notes: 
!
!******************************************************************************
!
! $Id: PICL_F90,v 1.0 2022/05/08 bdurant Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PICL_TEMP_Runge(pRegion)

!  USE 

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_level,t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt

  USE RFLU_ModDifferentiationCells
  USE RFLU_ModLimiters, ONLY: RFLU_CreateLimiter, &
                              RFLU_ComputeLimiterBarthJesp, &
                              RFLU_ComputeLimiterVenkat, &
                              RFLU_LimitGradCells, &
                              RFLU_LimitGradCellsSimple, &
                              RFLU_DestroyLimiter
  USE RFLU_ModWENO, ONLY: RFLU_WENOGradCellsWrapper, &
                          RFLU_WENOGradCellsXYZWrapper
#ifdef PICL
USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                             RFLU_ConvertCvPrim2Cons

 USE ModInterfaces, ONLY: RFLU_DecideWrite !BRAD added for picl
 
 use ppiclf_solve, only: ppiclf_solve_GetProFld, ppiclf_solve_InterpFieldUser, ppiclf_solve_IntegrateParticle
#endif



#ifdef PICL
!DEC$ NOFREEFORM
#include "PPICLF_STD.h"
!DEC$ FREEFORM
#endif


  IMPLICIT NONE


! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString


TYPE(t_global), POINTER :: global
TYPE(t_level), POINTER :: levels(:)
TYPE(t_region), POINTER :: pRegion
TYPE(t_grid), POINTER :: pGrid
!INTEGER :: errorFlag

#ifdef PICL
  LOGICAL :: doWrite      
  INTEGER(KIND=4) :: i,piclIO,nCells
  INTEGER :: errorFlag,icg      
  REAL(KIND=8) :: piclDtMin,piclCurrentTime, &
          temp_dudtMixt,temp_dvdtMixt,temp_dwdtMixt,energydotg
  REAL(KIND=8) :: dudx,dudy,dudz
  REAL(KIND=8) :: dvdx,dvdy,dvdz
  REAL(KIND=8) :: dwdx,dwdy,dwdz
  REAL(KIND=8) :: vFrac

  REAL(KIND=8), DIMENSION(3) :: ug      
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhoF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: uxF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: uyF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: uzF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: csF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tpF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ppF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: vfP
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dpxF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dpyF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dpzF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: SDRX
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: SDRY
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: SDRZ
  REAL(KIND=8), DIMENSION(:,:,:), POINTER :: pGc 
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhsR        
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: pGcX 
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: pGcY
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: pGcZ
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFX
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFXCell
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFY
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFYCell
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFZ
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFZCell
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFECell
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: PhiP
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: YTEMP
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: domgdx
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: domgdy
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: domgdz
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: drhodx
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: drhody
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: drhodz
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dpvxF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dpvyF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dpvzF
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: SDOX
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: SDOY
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: SDOZ
!---------------------------------------------------------------  
! Below added for pseudo turbulence
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JRSG11
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JRSG12
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JRSG13
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JRSG21
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JRSG22
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JRSG23
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JRSG31
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JRSG32
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JRSG33
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JTSG1
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JTSG2
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JTSG3
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: JRSGCell
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: JTSGCell
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: DivPhiRSG
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhog
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ugas
  INTEGER(KIND=4) :: j
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Qsg
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: DivPhiQsg
  !---------------------------------------------------------------  
  
  ! TLJ - added for Feedback term - 04/01/2025
  INTEGER, DIMENSION(:), ALLOCATABLE :: varInfoPicl
  INTEGER, DIMENSION(:), POINTER :: piclcvInfo
    REAL(KIND=8) :: dodx, dody, dodz,     &
                  omgx, omgy, omgz,     &
                  divu,                 &
                  dprdx, dprdy, dprdz,  &
                  dpdx, dpdy, dpdz,     &
                  phirho, ir, ir2 ,     &
                  dfxdx, dfxdy, dfxdz,  &
                  dfydx, dfydy, dfydz,  &
                  dfzdx, dfzdy, dfzdz   

#endif


   
!******************************************************************************

  RCSIdentString = '$RCSfile: PICL_TEMP_Runge.F90,v $ $Revision: 1.0 $'
 
  global => pRegion%global
  
  CALL RegisterFunction(global, 'PICL_TEMP_Runge',__FILE__ )



! Set pointers ----------------------------------------------------------------

    !pRegion => regions!pLevel%regions(iReg)
    pGrid   => pRegion%grid

!PPICLF Integration
#ifdef PICL

     piclIO = 100000000
     piclDtMin = REAL(global%dtMin,8)
     piclCurrentTime = REAL(global%currentTime,8)

     ! TLJ - 11/23/2024
     !     - This has now been removed
     doWrite = RFLU_DecideWrite(global)
     !Figure out piclIO call, might need to look into timestepping
     IF((doWrite .EQV. .TRUE.)) piclIO = 1


!PARTICLE stuff possbile needed
!    CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)


!allocate arrays to send to picl
    nCells = pRegion%grid%nCells
    ALLOCATE(rhoF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uxF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uyF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uzF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(csF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(tpF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error    

    ALLOCATE(ppF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error    

    ALLOCATE(vfP(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpxF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
    
    ALLOCATE(dpyF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpzF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRX(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRY(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRZ(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(rhsR(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcX(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcY(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcZ(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFX(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFXCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFY(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFYCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFZ(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFZCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFECell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFE(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(PhiP(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    IF(pRegion%mixtInput%axiFlag) THEN
      ALLOCATE(YTEMP(nCells),STAT=errorFlag)
      global%error = errorFlag
      IF(global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
      END IF ! global%error
    ENDIF

    ALLOCATE(domgdx(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
    
    ALLOCATE(domgdy(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(domgdz(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(drhodx(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(drhody(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(drhodz(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpvxF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpvyF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpvzF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDOX(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDOY(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDOZ(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

! Added for pseudo turbulence
!---------------------------------------------------------------  
    ALLOCATE(JRSG11(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JRSG12(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JRSG13(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JRSG21(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JRSG22(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JRSG23(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JRSG31(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JRSG32(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JRSG33(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JTSG1(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JTSG2(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JTSG3(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JRSGCell(9,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JTSGCell(3,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(DivPhiRSG(3,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(rhog(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(ugas(3,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(Qsg(3,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(DivPhiQsg(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

!---------------------------------------------------------------  

pGc => pRegion%mixt%gradCell
    ! 04/01/2025 - TLJ - we need feedback terms and their gradients to
    !       calculate the undisturbed torque component
    ! Internal definitions; some redundancy but just ignore
    ! We do not need energy, but might in the future
    DO i = 1,pRegion%grid%nCells
       JFXCell(i) = 0.0_RFREAL
       JFYCell(i) = 0.0_RFREAL
       JFZCell(i) = 0.0_RFREAL
       JFECell(i) = 0.0_RFREAL
       JFX(i) = 0.0_RFREAL
       JFY(i) = 0.0_RFREAL
       JFZ(i) = 0.0_RFREAL
       JFE(i) = 0.0_RFREAL
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JFX,JFX(i))  
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JFY,JFY(i))
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JFZ,JFZ(i))
       JFXCell(i) = JFX(i) 
       JFYCell(i) = JFY(i) 
       JFZCell(i) = JFZ(i) 
       pregion%mixt%piclFeedback(1,i) = JFXCell(i)
       pregion%mixt%piclFeedback(2,i) = JFYCell(i)
       pregion%mixt%piclFeedback(3,i) = JFZCell(i)
    END DO
    ! Now calculate the gradient of the feedback force
    ALLOCATE(varInfoPicl(3),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(piclcvInfo(3),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    varInfoPicl(1) = 1
    varInfoPicl(2) = 2
    varInfoPicl(3) = 3
    piclcvInfo = varInfoPicl
    CALL RFLU_ComputeGradCellsWrapper(pRegion,1,3,1,3,varInfoPicl, &
                                      pRegion%mixt%piclFeedback,&
                                      pRegion%mixt%piclgradFeedback)
    CALL RFLU_WENOGradCellsXYZWrapper(pRegion,1,3, &
                                      pRegion%mixt%piclgradFeedback)
    CALL RFLU_LimitGradCellsSimple(pRegion,1,3,1,3, &
                                   pRegion%mixt%piclFeedback,&
                                   piclcvInfo,&
                                   pRegion%mixt%piclgradFeedback)
    DEALLOCATE(varInfoPicl,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(piclcvInfo,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
  ! END - TLJ calculating gradient of feedback force

!Fill arrays for interp field
    DO i = 1,pRegion%grid%nCells
!Zero out phip
       PhiP(i) = 0.0_RFREAL

       ug(XCOORD) = pRegion%mixt%cv(CV_MIXT_XMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(YCOORD) = pRegion%mixt%cv(CV_MIXT_YMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(ZCOORD) = pRegion%mixt%cv(CV_MIXT_ZMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)
!------------------------------------------------------------------------------------- 
! Brad's Old formulation. Keep it here for reference
!         temp_drudtMixt = -pRegion%mixt%rhs(CV_MIXT_XMOM,i)/pRegion%grid%vol(i) &    
!                   +pRegion%mixt%cvOld(CV_MIXT_DENS,i)*DOT_PRODUCT(ug,pGc(:,2,i)) & 
!                   +ug(XCOORD)*DOT_PRODUCT(ug,pGc(:,1,i))                           
!                                                                                    
!         temp_drvdtMixt = -pRegion%mixt%rhs(CV_MIXT_YMOM,i)/pRegion%grid%vol(i) &    
!                   +pRegion%mixt%cvOld(CV_MIXT_DENS,i)*DOT_PRODUCT(ug,pGc(:,3,i))&  
!                   +ug(YCOORD)*DOT_PRODUCT(ug,pGc(:,1,i))                           
!                                                                                    
!         temp_drwdtMixt = -pRegion%mixt%rhs(CV_MIXT_ZMOM,i)/pRegion%grid%vol(i) &    
!                   +pRegion%mixt%cvOld(CV_MIXT_DENS,i)*DOT_PRODUCT(ug,pGc(:,4,i))&  
!                   +ug(ZCOORD)*DOT_PRODUCT(ug,pGc(:,1,i))                           
!------------------------------------------------------------------------------------- 
! IMPORTANT NOTE:
! 05/25/2025 - Thierry - rhs & diss arrays have to be divided by pRegion%grid%vol(i)
!                        when being used for interpolating values from Rocflu to ppiclF
!------------------------------------------------------------------------------------- 
      
       ! 03/11/2025 - Thierry - Du/Dt, Dv/Dt, Dw/Dt (not weighted by phi^g or rho^g)

       temp_dudtMixt  = (-pRegion%mixt%rhs(CV_MIXT_XMOM,i)/pRegion%grid%vol(i)& 
                         +ug(XCOORD)*pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i))&
                         /pRegion%mixt%cv(CV_MIXT_DENS,i)&
                         +DOT_PRODUCT(ug,pGc(:,2,i))

       temp_dvdtMixt  = (-pRegion%mixt%rhs(CV_MIXT_YMOM,i)/pRegion%grid%vol(i)& 
                         +ug(YCOORD)*pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i))&
                         /pRegion%mixt%cv(CV_MIXT_DENS,i)&
                         +DOT_PRODUCT(ug,pGc(:,3,i))
                         
       temp_dwdtMixt  = (-pRegion%mixt%rhs(CV_MIXT_ZMOM,i)/pRegion%grid%vol(i)& 
                         +ug(ZCOORD)*pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i))&
                         /pRegion%mixt%cv(CV_MIXT_DENS,i)&
                         +DOT_PRODUCT(ug,pGc(:,4,i))

       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JPHIP,vfP(i))
       vfP(i) = vfP(i)/pRegion%grid%vol(i)
       PhiP(i) = vfP(i)
       !VOL Frac cap
       IF(PhiP(i) .GT. 0.62) PhiP(i) = 0.62
       vfp(i) = PhiP(i)      

       ! TLJ - 02/07/2025 scaled conserved density by gas-phase volume fraction
       vFrac = 1.0_RFREAL - PhiP(i)!pRegion%mixt%piclVF(i)

       rhoF(i) = pRegion%mixt%cv(CV_MIXT_DENS,i) / vFrac
       uxF(i) =  pRegion%mixt%cv(CV_MIXT_XMOM,i) &
                /pRegion%mixt%cv(CV_MIXT_DENS,i)
       uyF(i) =  pRegion%mixt%cv(CV_MIXT_YMOM,i) &
                /pRegion%mixt%cv(CV_MIXT_DENS,i)
       uzF(i) =  pRegion%mixt%cv(CV_MIXT_ZMOM,i) &
                /pRegion%mixt%cv(CV_MIXT_DENS,i)

       csF(i) = pRegion%mixt%dv(DV_MIXT_SOUN,i)
       tpF(i) = pRegion%mixt%dv(DV_MIXT_TEMP,i) 
       ! Davin - added pressure to interpolation values 02/22/2025
       ppF(i) = pRegion%mixt%dv(DV_MIXT_PRES,i) 

       dpxF(i) = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_PRES,i) ! dp/dx
       dpyF(i) = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_PRES,i) ! dp/dy
       dpzF(i) = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_PRES,i) ! dp/dz

       dudx = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_XVEL,i)
       dudy = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_XVEL,i)
       dudz = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_XVEL,i)

       dvdx = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_YVEL,i)
       dvdy = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_YVEL,i)
       dvdz = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_YVEL,i)

       dwdx = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_ZVEL,i)
       dwdy = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_ZVEL,i)
       dwdz = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_ZVEL,i)

       domgdx(i) = dwdy - dvdz
       domgdy(i) = dudz - dwdx
       domgdz(i) = dvdx - dudy

       ! 04/01/2025 - TLJ - Calculate the substantial derivative of vorticity
       ! Internal definitions; some redundancy but just ignore
       dodx   = 0.0_RFREAL ! D(Omega_x)/DT
       dody   = 0.0_RFREAL ! D(Omega_y)/DT
       dodz   = 0.0_RFREAL ! D(Omega_z)/DT
       omgx   = dwdy - dvdz ! Omega_x
       omgy   = dudz - dwdx ! Omega_y
       omgz   = dvdx - dudy ! Omega_z
       divu   = dudx + dvdy + dwdz ! u_x+v_y+w_z; divergence of velocity
       dprdx  = pGc(XCOORD,1,i) ! d(rho phi)/dx
       dprdy  = pGc(YCOORD,1,i) ! d(rho phi)/dy
       dprdz  = pGc(ZCOORD,1,i) ! d(rho phi)/dz
       dpdx   = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_PRES,i) ! dp/dx
       dpdy   = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_PRES,i) ! dp/dy
       dpdz   = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_PRES,i) ! dp/dz
       dfxdx  = pRegion%mixt%piclgradFeedback(XCOORD,1,i) ! dFx/dx
       dfxdy  = pRegion%mixt%piclgradFeedback(YCOORD,1,i) ! dFx/dy
       dfxdz  = pRegion%mixt%piclgradFeedback(ZCOORD,1,i) ! dFx/dz
       dfydx  = pRegion%mixt%piclgradFeedback(XCOORD,2,i) ! dFy/dx
       dfydy  = pRegion%mixt%piclgradFeedback(YCOORD,2,i) ! dFy/dy
       dfydz  = pRegion%mixt%piclgradFeedback(ZCOORD,2,i) ! dFy/dz
       dfzdx  = pRegion%mixt%piclgradFeedback(XCOORD,3,i) ! dFz/dx
       dfzdy  = pRegion%mixt%piclgradFeedback(YCOORD,3,i) ! dFz/dy
       dfzdz  = pRegion%mixt%piclgradFeedback(ZCOORD,3,i) ! dFz/dz
       phirho = pRegion%mixt%cv(CV_MIXT_DENS,i) ! phi_g*rho_g
       ir     = 1.0_RFREAL / phirho
       ir2    = ir*ir
       ! 1. Vortex stretching
       dodx = omgx*dudx + omgy*dudy + omgz*dudz
       dody = omgx*dvdx + omgy*dvdy + omgz*dvdz
       dodz = omgx*dwdx + omgy*dwdy + omgz*dwdz
       ! 2. Vortex dilatation
       dodx = dodx - omgx*divu
       dody = dody - omgy*divu
       dodz = dodz - omgz*divu
       ! 3. Baroclinic
       dodx = dodx + (dprdy*dpdz - dprdz*dpdy)*ir2
       dody = dody + (dprdz*dpdx - dprdx*dpdz)*ir2
       dodz = dodz + (dprdx*dpdy - dprdy*dpdx)*ir2
       ! 4. Torque due to feedback force
       dodx = dodx + (dfzdy - dfydz)*ir
       dody = dody + (dfxdz - dfzdx)*ir
       dodz = dodz + (dfydx - dfxdy)*ir
       ! 5. Misalignment of phi*rho and feedback force
       dodx = dodx + (dprdy*JFZCell(i) - dprdz*JFYCell(i))*ir2
       dody = dody + (dprdz*JFXCell(i) - dprdx*JFZCell(i))*ir2
       dodz = dodz + (dprdx*JFYCell(i) - dprdy*JFXCell(i))*ir2
       ! 6. Add terms and store
       SDOX(i) = dodx
       SDOY(i) = dody
       SDOZ(i) = dodz
       ! End - TLJ - Calculate the substantial derivative of vorticity

       ! Substantial derivative of gas-phase velocity
       SDRX(i) = temp_dudtMixt ! Du/Dt
       SDRY(i) = temp_dvdtMixt ! Dv/Dt
       SDRZ(i) = temp_dwdtMixt ! Dw/Dt

       rhsR(i) = -pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i) ! \p(rho*phi)/\p(t)

       pGcX(i) = pGc(XCOORD,1,i) ! d(rho phi)/dx
       pGcY(i) = pGc(YCOORD,1,i) ! d(rho phi)/dy
       pGcz(i) = pGc(ZCOORD,1,i) ! d(rho phi)/dz

       ! Gradient of rho^g of mixture (not weighted by phi^g!)
       ! Using grad(rhog) directly
       drhodx(i) = pRegion%mixt%piclgradRhog(1,1,i) ! d(rho)/dx
       drhody(i) = pRegion%mixt%piclgradRhog(2,1,i) ! d(rho)/dy
       drhodz(i) = pRegion%mixt%piclgradRhog(3,1,i) ! d(rho)/dz

       ! Viscous term of pressure gradient (divergence of tau)
       dpvxF(i) = pRegion%mixt%diss(CV_MIXT_XMOM,i)/pRegion%grid%vol(i)
       dpvyF(i) = pRegion%mixt%diss(CV_MIXT_YMOM,i)/pRegion%grid%vol(i)
       dpvzF(i) = pRegion%mixt%diss(CV_MIXT_ZMOM,i)/pRegion%grid%vol(i)
     END DO

! Interp field calls
! TLJ - interpolates various fluid quantities onto the 
!       the ppiclf particle locations
! TLJ PPICLF_LRP_INT in PPICLF_USER.h must match the number
!     of calls to ppiclf_solve_InterpFieldUser
! Davin - added pressure 02/22/2025
      IF(PPICLF_LRP_INT .NE. 30) THEN
         WRITE(*,*) "Error: PPICLF_LRP_INT must be set to 30"
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:LRP_INT')
      END IF
 
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOF,rhoF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUX,uxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUY,uyF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUZ,uzF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDX,dpxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDY,dpyF)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDZ,dpzF)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JCS,csF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JT,tpF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JPHIP,vfP)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDRX,SDRX)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDRY,SDRY)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDRZ,SDRZ)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHSR,rhsR)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JPGCX,pGcX) 
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JPGCY,pGcY) 
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JPGCZ,pGcZ) 
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JXVOR,domgdx)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JYVOR,domgdy)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JZVOR,domgdz)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JP,ppF)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOGX,drhodx)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOGY,drhody)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOGZ,drhodz)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPVDX,dpvxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPVDY,dpvyF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPVDZ,dpvzF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDOX,SDOX)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDOY,SDOY)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JSDOZ,SDOZ)  

 ! Solve RK stage of time stepping particle solution
     CALL ppiclf_solve_IntegrateParticle(1,piclIO,piclDtMin,piclCurrentTime)

!FEED BACK TERMS
     JFXCell = 0.0_RFREAL
     JFYCell = 0.0_RFREAL
     JFZCell = 0.0_RFREAL
     JFECell = 0.0_RFREAL
     JRSGCell = 0.0_RFREAL
     DivPhiRSG = 0.0_RFREAL
     JTSGCell = 0.0_RFREAL
     Qsg = 0.0_RFREAL
     DivPhiQsg = 0.0_RFREAL

!Fill arrays for interp field
IF(global%piclFeedbackFlag == 1) THEN
    DO i = 1,pRegion%grid%nCells
       ug(XCOORD) = pRegion%mixt%cv(CV_MIXT_XMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(YCOORD) = pRegion%mixt%cv(CV_MIXT_YMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(ZCOORD) = pRegion%mixt%cv(CV_MIXT_ZMOM,i)&
                          /pRegion%mixt%cv(CV_MIXT_DENS,i)

       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JFX,JFX(i))  
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JFY,JFY(i))
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JFZ,JFZ(i))
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JE,JFE(i)) 
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JPHIP,vfP(i))
       PhiP(i) = vfP(i)/pRegion%grid%vol(i)

       !VOL Frac cap
       IF(PhiP(i) .GT. 0.62) PhiP(i) = 0.62
       vfp(i) = PhiP(i)      
   !---------------------------------------------------------------------------------------
       ! 07/21/2025 - Thierry - begins here - added for PseudoTurbulence
       if(global%piclPseudoTurbFlag .eq. 1) then
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG11,JRSG11(i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG12,JRSG12(i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG13,JRSG13(i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG21,JRSG21(i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG22,JRSG22(i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG23,JRSG23(i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG31,JRSG31(i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG32,JRSG32(i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JRSG33,JRSG33(i))

         call ppiclf_solve_GetProFld(i,PPICLF_P_JTSG1,JTSG1(i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JTSG2,JTSG2(i))
         call ppiclf_solve_GetProFld(i,PPICLF_P_JTSG3,JTSG3(i))
        
         JRSGCell(1,i) = JRSG11(i)
         JRSGCell(2,i) = JRSG12(i)
         JRSGCell(3,i) = JRSG13(i)
         JRSGCell(4,i) = JRSG21(i)
         JRSGCell(5,i) = JRSG22(i)
         JRSGCell(6,i) = JRSG23(i)
         JRSGCell(7,i) = JRSG31(i)
         JRSGCell(8,i) = JRSG32(i)
         JRSGCell(9,i) = JRSG33(i)
         JTSGCell(1,i) = JTSG1(i)
         JTSGCell(2,i) = JTSG2(i)
         JTSGCell(3,i) = JTSG3(i)

       endif ! piclPseudoTurbFlag
!---------------------------------------------------------------------------------------
       JFXCell(i) = JFX(i) 
       JFYCell(i) = JFY(i) 
       JFZCell(i) = JFZ(i) 
       !JE correction
       JFECell(i) = JFE(i) 
       energydotg = JFECell(i) ! includs KE feedback already

       IF(IsNan(JFXCell(i)) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-PX",i,JFXCell(i),ug(1),ug(2),ug(3)
         write(*,*) "JFY",i,JFYCell(i)
         write(*,*) "JFZ",i,JFZCell(i)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PX')
       END IF
       IF(IsNan(JFYCell(i)) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-PY",i,JFYCell(i),ug(1),ug(2),ug(3)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PY')
       END IF
       IF(IsNan(JFZCell(i)) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-PZ",i,JFZCell(i),ug(1),ug(2),ug(3)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PY')
       END IF
       IF(IsNan(energydotg) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-PE",energydotg,i,JFXCell(i),ug(1),JFYCell(i),ug(2),pregion%grid%vol(i)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PE')
       END IF
       IF (ANY(IsNan(JRSGCell(:,i))) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-RSG",i,JRSGCell(:,i)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken Reynolds SG')
       endif
       IF(ANY(IsNan(JTSGCell(:,i))) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-TSG",i,JTSGCell(:,i)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken TSG')
        endif

        pRegion%mixt%rhs(CV_MIXT_XMOM,i) &
                         = pRegion%mixt%rhs(CV_MIXT_XMOM,i) &
                         + JFXCell(i)
        
        pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                         = pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                         + JFYCell(i)

        pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                         = pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                         + JFZCell(i)

        pRegion%mixt%rhs(CV_MIXT_ENER,i) &
                         = pRegion%mixt%rhs(CV_MIXT_ENER,i) &
                         + energydotg
!---------------------------------------------------------------------------------------
! Compute the interpolated Reynolds Stress components before projecting them back.
!---------------------------------------------------------------------------------------

         if(global%piclPseudoTurbFlag .eq. 1) then
         ! Conservative to primitive variables
           rhog(i) = pRegion%mixt%cv(CV_MIXT_DENS,i) / (1.0_RFREAL - PhiP(i))
           ugas(XCOORD,i) = pRegion%mixt%cv(CV_MIXT_XMOM,i)&
                           /pRegion%mixt%cv(CV_MIXT_DENS,i)

           ugas(YCOORD,i) = pRegion%mixt%cv(CV_MIXT_YMOM,i)&
                            /pRegion%mixt%cv(CV_MIXT_DENS,i)

           ugas(ZCOORD,i) = pRegion%mixt%cv(CV_MIXT_ZMOM,i)&
                          /pRegion%mixt%cv(CV_MIXT_DENS,i)

           do j=1,9
             ! R_sg -> rho_g R_sg
             ! Multiply by density here as Osnes formulation does not include it
             JRSGCell(j,i) = JRSGCell(j,i) * rhog(i)
           end do

           ! K_sg = 1/(2*rhof) * tr(Rsg), dimension of (nCellsTot)
           ! K_sg is added to the Total Gas Energy term
           pRegion%mixt%piclKsg(i) = 1.0_RFREAL/(2.0_RFREAL*rhog(i)) &
                                  * (JRSGCell(1,i) + JRSGCell(5,i) + JRSGCell(9,i))

           ! Q_sg : Subgrid Energy Flux dimension (3, nCells)
           !Qsg = rhog*Tsg + ug.Rsg
           ! cp = pRegion%mixt%gv(GV_MIXT_CP,indCp*1:nCells)
           
           Qsg(XCOORD,i) = rhog(i)*JTSGCell(XCOORD,i) +   &
              (JRSGCell(1,i)*ugas(XCOORD,i) + JRSGCell(2,i)*ugas(YCOORD,i) + JRSGCell(3,i)*ugas(ZCOORD,i))

           Qsg(YCOORD,i) = rhog(i)*JTSGCell(YCOORD,i) +   & 
              (JRSGCell(4,i)*ugas(XCOORD,i) + JRSGCell(5,i)*ugas(YCOORD,i) + JRSGCell(6,i)*ugas(ZCOORD,i))

           Qsg(ZCOORD,:) = rhog(i)*JTSGCell(ZCOORD,i) +   &
              (JRSGCell(7,i)*ugas(XCOORD,i) + JRSGCell(8,i)*ugas(YCOORD,i) + JRSGCell(9,i)*ugas(ZCOORD,i))

           ! Storing for ParaView plotting
           do j=1,9
             ! \phi_g \rho_g R_sg
             pRegion%mixt%piclPhiRSG(j,i) = JRSGCell(j,i) * (1.0_RFREAL - PhiP(i))
           end do

           do j=1,3
             ! Q_sg -> \phi_g Q_sg
             pRegion%mixt%piclPhiQsg(j,i) = Qsg(j,i) * (1.0_RFREAL - PhiP(i))
           enddo

           IF (ANY(IsNan(pRegion%mixt%piclPhiQsg(:,i))) .EQV. .TRUE.) THEN
                   write(*,*) "BROKEN piclPhiQsg", i, pRegion%mixt%piclPhiQsg(:,i)
                   write(*,*) "rhog(i), ugas(1:3,i)", rhog(i), ugas(1:3,i)
                   write(*,*) "JTSGCell(1:3,i) ", JTSGCell(1:3,i)
                   write(*,*) "JRSGCell(1:9,i)", JRSGCell(:,i)
                   CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken PhiQSG')
           ENDIF
         
         endif ! piclPseudoTurbFlag

    END DO !nCells
       
       if(global%piclPseudoTurbFlag .eq. 1) then

         ALLOCATE(varInfoPicl(9),STAT=errorFlag)
         ALLOCATE(piclcvInfo(9),STAT=errorFlag)
         varInfoPicl = [(j, j=1, 9)]
         piclcvInfo = varInfoPicl
!
         CALL RFLU_ComputeGradCellsWrapper(pRegion,1,9,1,9,varInfoPicl, &   
                                           pRegion%mixt%piclPhiRsg,&                
                                           pRegion%mixt%piclGradPhiRsg)       
                                                                           
         CALL RFLU_WENOGradCellsXYZWrapper(pRegion,1,9, &                   
                                           pRegion%mixt%piclGradPhiRsg)       
                                                                           
         CALL RFLU_LimitGradCellsSimple(pRegion,1,9,1,9, &                  
                                        pRegion%mixt%piclPhiRsg,&                   
                                        piclcvInfo,&                        
                                        pRegion%mixt%piclGradPhiRsg) 

         DEALLOCATE(varInfoPicl,STAT=errorFlag)
         DEALLOCATE(piclcvInfo,STAT=errorFlag)

         ALLOCATE(varInfoPicl(3),STAT=errorFlag)
         ALLOCATE(piclcvInfo(3),STAT=errorFlag)
         varInfoPicl = [(j, j=1, 3)]
         piclcvInfo = varInfoPicl
!
         CALL RFLU_ComputeGradCellsWrapper(pRegion,1,3,1,3,varInfoPicl, &   
                                           pRegion%mixt%piclPhiQsg,&                
                                           pRegion%mixt%piclGradPhiQsg)       
                                                                           
         CALL RFLU_WENOGradCellsXYZWrapper(pRegion,1,3, &                   
                                           pRegion%mixt%piclGradPhiQsg)       
                                                                           
         CALL RFLU_LimitGradCellsSimple(pRegion,1,3,1,3, &                  
                                        pRegion%mixt%piclPhiQsg,&                   
                                        piclcvInfo,&                        
                                        pRegion%mixt%piclGradPhiQsg)          

         DEALLOCATE(varInfoPicl,STAT=errorFlag)
         DEALLOCATE(piclcvInfo,STAT=errorFlag)

!       ! Now compute Div(\phi_g R_sg)
        ! 07/23/2025 - Thierry Daoud 
        ! pRegion%mixt%piclgradPhiRSG dimension is (3,9,nCellsTot)
        ! nCellsTot = no. actual cells (nCells) + no. dummy (ghost?) cells
        ! When calculating the gradient, you need to use nCellsTot
        
! Div (\phi_g R_sg) - comma denotes partial derivative (,3 -> partial / partial x_3)
! x-direction: Div(\phi_g R_sg),x = (\phi R_11),1 + (\phi R_12),2 + (\phi R_13),3
! y-direction: Div(\phi_g R_sg),y = (\phi R_21),1 + (\phi R_22),2 + (\phi R_23),3
! z-direction: Div(\phi_g R_sg),z = (\phi R_31),1 + (\phi R_32),2 + (\phi R_33),3

    DO i = 1,pRegion%grid%nCells
         DivPhiRSG(XCOORD,i) = pregion%grid%vol(i)*(                   &
                               pRegion%mixt%piclGradPhiRsg(XCOORD,1,i) &
                             + pRegion%mixt%piclGradPhiRsg(YCOORD,2,i) &
                             + pRegion%mixt%piclGradPhiRsg(ZCOORD,3,i))

         DivPhiRSG(YCOORD,i) = pregion%grid%vol(i)*(            &
                               pRegion%mixt%piclGradPhiRsg(XCOORD,4,i) &
                             + pRegion%mixt%piclGradPhiRsg(YCOORD,5,i) &
                             + pRegion%mixt%piclGradPhiRsg(ZCOORD,6,i))

         DivPhiRSG(ZCOORD,i) = pregion%grid%vol(i)*(            &
                               pRegion%mixt%piclGradPhiRsg(XCOORD,7,i) &
                             + pRegion%mixt%piclGradPhiRsg(YCOORD,8,i) &
                             + pRegion%mixt%piclGradPhiRsg(ZCOORD,9,i))


! Div (\phi_g Q_sg) - comma denotes partial derivative (,3 -> partial / partial x_3)
! Scalar: Div(\phi_g Q_sg) = (\phi Q_1),1 + (\phi Q_2),2 + (\phi Q_3),3
         DivPhiQsg(i) =  pregion%grid%vol(i)*(                   &
                         pRegion%mixt%piclGradPhiQsg(XCOORD,1,i) &
                       + pRegion%mixt%piclGradPhiQsg(YCOORD,2,i) &
                       + pRegion%mixt%piclGradPhiQsg(ZCOORD,3,i))

 IF (IsNan(DivPhiQsg(i)) .EQV. .TRUE.) THEN
        write(*,*) "BROKEN- DivPhiQSG, i ", DivPhiQsg(i), i
        write(*,*) "pregion%grid%vol(i)", pregion%grid%vol(i)
        write(*,*) "pRegion%mixt%piclGradPhiQsg(XCOORD,1,i)", pRegion%mixt%piclGradPhiQsg(XCOORD,1,i)
        write(*,*) "pRegion%mixt%piclGradPhiQsg(YCOORD,2,i)", pRegion%mixt%piclGradPhiQsg(YCOORD,2,i)
        write(*,*) "pRegion%mixt%piclGradPhiQsg(ZCOORD,3,i)", pRegion%mixt%piclGradPhiQsg(ZCOORD,3,i)
        write(*,*) "rhog(i), ugas(1:3,i)", rhog(i), ugas(1:3,i)
        write(*,*) "JTSGCell(1:3,i) ", JTSGCell(1:3,i)
        write(*,*) "JRSGCell(1:9,i)", JRSGCell(:,i)
        CALL ErrorStop(global,ERR_INVALID_VALUE ,__LINE__,'PPICLF:Broken QSG')
endif
         
         ! Feedback Div(phi Rsg) to the Fluid Momentum Equations
         pRegion%mixt%rhs(CV_MIXT_XMOM,i) &
                          = pRegion%mixt%rhs(CV_MIXT_XMOM,i) &
                          + DivPhiRsg(XCOORD,i) 
         
         pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                          = pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                          + DivPhiRsg(YCOORD,i) 

         pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                          = pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                          + DivPhiRsg(ZCOORD,i)
                        
         ! Feedback Div(phi Qsg) to the Fluid Energy Equation
         pRegion%mixt%rhs(CV_MIXT_ENER,i) &
                          = pRegion%mixt%rhs(CV_MIXT_ENER,i) &
                          + DivPhiQsg(i)

        ENDDO

       endif ! piclPseudoTurbFlag
END IF ! global%piclFeedbackFlag

!Due to moving particle integration stuff stoping this for now
DO i = 1,pRegion%grid%nCells
!zero out PhiP
       PhiP(i) = 0.0D0
       CALL ppiclf_solve_GetProFld(i,PPICLF_P_JPHIP,vfP(i))
       vfP(i) = vfP(i)/pRegion%grid%vol(i)
       PhiP(i) = vfP(i)
!VOL Frac Cap
! Should we keep this???
       IF(PhiP(i) .GT. 0.62) PhiP(i) = 0.62
       pRegion%mixt%piclVF(i) = PhiP(i) 
END DO
!Deallocate arrays

    IF(pRegion%mixtInput%axiFlag) THEN
      DEALLOCATE(YTEMP,STAT=errorFlag)
      global%error = errorFlag
      IF(global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
      END IF ! global%error
    END IF

    DEALLOCATE(rhoF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uxF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uyF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uzF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(csF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(tpF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDRX,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDRY,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDRZ,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(rhsR,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(pGcX,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(pGcY,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF !global%error    

    DEALLOCATE(pGcZ,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF !global%error    

    DEALLOCATE(JFX,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFXCell,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFY,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFYCell,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFZ,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFZCell,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFE,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFECell,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error    

    DEALLOCATE(domgdx,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(domgdy,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(domgdz,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(drhodx,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(drhody,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(drhodz,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpvxF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpvyF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpvzF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDOX,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDOY,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDOZ,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(PhiP,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(ppF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(vfP,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpxF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpyF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
 
    DEALLOCATE(dpzF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

!---------------------------------------------------------------  
    DEALLOCATE(JRSG11,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JRSG12,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JRSG13,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JRSG21,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JRSG22,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JRSG23,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JRSG31,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JRSG32,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JRSG33,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
    DEALLOCATE(JTSG1,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JTSG2,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JTSG3,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JRSGCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JTSGCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(DivPhiRSG,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(rhog,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(ugas,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(Qsg,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(DivPhiQsg,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
!---------------------------------------------------------------  
#endif
!PPICLF Integration END

! finalize --------------------------------------------------------------------
  CALL DeregisterFunction(global )
END SUBROUTINE PICL_TEMP_Runge

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PICL_.F90,v $
!
!
!******************************************************************************

