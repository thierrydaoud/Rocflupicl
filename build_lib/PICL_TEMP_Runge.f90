










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
USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                             RFLU_ConvertCvPrim2Cons

 USE ModInterfaces, ONLY: RFLU_DecideWrite !BRAD added for picl
 



!DEC$ NOFREEFORM

! number of timesteps kept in history kernels
! maximum number of triangular patch boundaries

! y, y1, ydot, ydotc: 12

! rprop: 64

! map: 10






















!DEC$ FREEFORM


  IMPLICIT NONE


! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString


TYPE(t_global), POINTER :: global
TYPE(t_level), POINTER :: levels(:)
TYPE(t_region), POINTER :: pRegion
TYPE(t_grid), POINTER :: pGrid
!INTEGER :: errorFlag

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



   
!******************************************************************************

  RCSIdentString = '$RCSfile: PICL_TEMP_Runge.F90,v $ $Revision: 1.0 $'
 
  global => pRegion%global
  
  CALL RegisterFunction(global, 'PICL_TEMP_Runge',"../rocpicl/PICL_TEMP_Runge.F90" )



! Set pointers ----------------------------------------------------------------

    !pRegion => regions!pLevel%regions(iReg)
    pGrid   => pRegion%grid

!PPICLF Integration

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
      CALL ErrorStop(global,ERR_ALLOCATE,234,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uxF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,240,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uyF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,246,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uzF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,252,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(csF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,258,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(tpF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,264,'PPICLF:xGrid')
    END IF ! global%error    

    ALLOCATE(ppF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,270,'PPICLF:xGrid')
    END IF ! global%error    

    ALLOCATE(vfP(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,276,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpxF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,282,'PPICLF:xGrid')
    END IF ! global%error
    
    ALLOCATE(dpyF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,288,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpzF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,294,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRX(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,300,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRY(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,306,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRZ(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,312,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(rhsR(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,318,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcX(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,324,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcY(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,330,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcZ(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,336,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFX(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,342,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFXCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,348,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFY(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,354,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFYCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,360,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFZ(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,366,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFZCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,372,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFECell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,378,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFE(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,384,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(PhiP(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,390,'PPICLF:xGrid')
    END IF ! global%error

    IF(pRegion%mixtInput%axiFlag) THEN
      ALLOCATE(YTEMP(nCells),STAT=errorFlag)
      global%error = errorFlag
      IF(global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,397,'PPICLF:xGrid')
      END IF ! global%error
    ENDIF

    ALLOCATE(domgdx(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,404,'PPICLF:xGrid')
    END IF ! global%error
    
    ALLOCATE(domgdy(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,410,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(domgdz(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,416,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(drhodx(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,422,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(drhody(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,428,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(drhodz(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,434,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpvxF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,440,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpvyF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,446,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpvzF(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,452,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDOX(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,458,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDOY(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,464,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDOZ(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,470,'PPICLF:xGrid')
    END IF ! global%error


!Might need to update prim like plag does
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
 

       CALL ppiclf_solve_GetProFld(i,2,JFX(i))  
       CALL ppiclf_solve_GetProFld(i,3,JFY(i))
       CALL ppiclf_solve_GetProFld(i,4,JFZ(i))

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
      CALL ErrorStop(global,ERR_ALLOCATE,506,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(piclcvInfo(3),STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,512,'PPICLF:xGrid')
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
      CALL ErrorStop(global,ERR_ALLOCATE,531,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(piclcvInfo,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,537,'PPICLF:xGrid')
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

       CALL ppiclf_solve_GetProFld(i,1,vfP(i))
       PhiP(i) = vfP(i)/pRegion%grid%vol(i)

       ! TLJ - 02/07/2025 scaled conserved density by gas-phase volume fraction
       vFrac = 1.0_RFREAL - pRegion%mixt%piclVF(i)
       rhoF(i) = pRegion%mixt%cv(CV_MIXT_DENS,i) / vFrac
       uxF(i) = pRegion%mixt%cv(CV_MIXT_XMOM,i) &
                /pRegion%mixt%cv(CV_MIXT_DENS,i)
       uyF(i) = pRegion%mixt%cv(CV_MIXT_YMOM,i) &
                /pRegion%mixt%cv(CV_MIXT_DENS,i)
       uzF(i) = pRegion%mixt%cv(CV_MIXT_ZMOM,i) &
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

       
       !Dump back VolFrac
       !VOL Frac cap
       IF(PhiP(i) .GT. 0.62) PhiP(i) = 0.62
       vfp(i) = PhiP(i)      

       END DO 
! Interp field calls
! TLJ - interpolates various fluid quantities onto the 
!       the ppiclf particle locations
! TLJ 30 in PPICLF_USER.h must match the number
!     of calls to ppiclf_solve_InterpFieldUser
! Davin - added pressure 02/22/2025
      IF(30 .NE. 30) THEN
         WRITE(*,*) "Error: PPICLF_LRP_INT must be set to 30"
         CALL ErrorStop(global,ERR_INVALID_VALUE ,699,'PPICLF:LRP_INT')
      END IF
 
      CALL ppiclf_solve_InterpFieldUser(2,rhoF)
      CALL ppiclf_solve_InterpFieldUser(6,uxF)
      CALL ppiclf_solve_InterpFieldUser(7,uyF)
      CALL ppiclf_solve_InterpFieldUser(8,uzF)
      CALL ppiclf_solve_InterpFieldUser(10,dpxF)
      CALL ppiclf_solve_InterpFieldUser(11,dpyF)  
      CALL ppiclf_solve_InterpFieldUser(12,dpzF)  
      CALL ppiclf_solve_InterpFieldUser(9,csF)
      CALL ppiclf_solve_InterpFieldUser(24,tpF)
      CALL ppiclf_solve_InterpFieldUser(5,vfP)  
      CALL ppiclf_solve_InterpFieldUser(13,SDRX)
      CALL ppiclf_solve_InterpFieldUser(14,SDRY)  
      CALL ppiclf_solve_InterpFieldUser(15,SDRZ)  
      CALL ppiclf_solve_InterpFieldUser(16,rhsR)  
      CALL ppiclf_solve_InterpFieldUser(17,pGcX) 
      CALL ppiclf_solve_InterpFieldUser(18,pGcY) 
      CALL ppiclf_solve_InterpFieldUser(19,pGcZ) 
      CALL ppiclf_solve_InterpFieldUser(31,domgdx)
      CALL ppiclf_solve_InterpFieldUser(32,domgdy)  
      CALL ppiclf_solve_InterpFieldUser(33,domgdz)  
      CALL ppiclf_solve_InterpFieldUser(36,ppF)  
      CALL ppiclf_solve_InterpFieldUser(37,drhodx)
      CALL ppiclf_solve_InterpFieldUser(38,drhody)
      CALL ppiclf_solve_InterpFieldUser(39,drhodz)
      CALL ppiclf_solve_InterpFieldUser(40,dpvxF)
      CALL ppiclf_solve_InterpFieldUser(41,dpvyF)
      CALL ppiclf_solve_InterpFieldUser(42,dpvzF)
      CALL ppiclf_solve_InterpFieldUser(62,SDOX)  
      CALL ppiclf_solve_InterpFieldUser(63,SDOY)  
      CALL ppiclf_solve_InterpFieldUser(64,SDOZ)  
 
!FEED BACK TERM
!Fill arrays for interp field
IF(global%piclFeedbackFlag == 1) THEN
    DO i = 1,pRegion%grid%nCells
       ug(XCOORD) = pRegion%mixt%cv(CV_MIXT_XMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(YCOORD) = pRegion%mixt%cv(CV_MIXT_YMOM,i)&
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)

       ug(ZCOORD) = pRegion%mixt%cv(CV_MIXT_ZMOM,i)&
                          /pRegion%mixt%cv(CV_MIXT_DENS,i)

       JFXCell(i) = 0.0_RFREAL
       JFYCell(i) = 0.0_RFREAL
       JFZCell(i) = 0.0_RFREAL
       JFECell(i) = 0.0_RFREAL
       CALL ppiclf_solve_GetProFld(i,2,JFX(i))  
       CALL ppiclf_solve_GetProFld(i,3,JFY(i))
       CALL ppiclf_solve_GetProFld(i,4,JFZ(i))
       CALL ppiclf_solve_GetProFld(i,5,JFE(i)) 
  
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
         CALL ErrorStop(global,ERR_INVALID_VALUE ,767,'PPICLF:Broken PX')
       END IF
       IF(IsNan(JFYCell(i)) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-PY",i,JFYCell(i),ug(1),ug(2),ug(3)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,772,'PPICLF:Broken PY')
       END IF
       IF(IsNan(JFZCell(i)) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-PZ",i,JFZCell(i),ug(1),ug(2),ug(3)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,777,'PPICLF:Broken PY')
       END IF
       IF(IsNan(energydotg) .EQV. .TRUE.) THEN
         write(*,*) "BROKEN-PE",energydotg,i,JFXCell(i),ug(1),JFYCell(i),ug(2),pregion%grid%vol(i)
         write(*,*) "pregionvol", pregion%grid%vol(i)
         CALL ErrorStop(global,ERR_INVALID_VALUE ,782,'PPICLF:Broken PE')
       END IF

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
    END DO

END IF ! global%piclFeedbackFlag
!SOLVE
    CALL ppiclf_solve_IntegrateParticle(1,piclIO,piclDtMin,piclCurrentTime)
!Due to moving particle integration stuff stoping this for now
DO i = 1,pRegion%grid%nCells
!zero out PhiP
       PhiP(i) = 0.0D0
       CALL ppiclf_solve_GetProFld(i,1,vfP(i))
       vfP(i) = vfP(i)/pRegion%grid%vol(i)
       PhiP(i) = vfP(i)
!VOL Frac Cap
! Should we keep this???
       IF(PhiP(i) .GT. 0.62) THEN
         PhiP(i) = 0.62
       END IF
       pRegion%mixt%piclVF(i) = PhiP(i) 
END DO
!Deallocate arrays

    IF(pRegion%mixtInput%axiFlag) THEN
      DEALLOCATE(YTEMP,STAT=errorFlag)
      global%error = errorFlag
      IF(global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,825,'PPICLF:xGrid')
      END IF ! global%error
    END IF

    DEALLOCATE(rhoF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,832,'PPICLF:xGrid')
    END IF ! global%error
    DEALLOCATE(uxF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,837,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uyF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,843,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uzF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,849,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(csF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,855,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(tpF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,861,'PPICLF:xGrid')
    END IF ! global%error
    DEALLOCATE(SDRX,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,866,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDRY,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,872,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDRZ,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,878,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(rhsR,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,884,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(pGcX,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,890,'PPICLF:xGrid')
    END IF ! global%error
    DEALLOCATE(pGcY,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,895,'PPICLF:xGrid')
    END IF !global%error    

    DEALLOCATE(pGcZ,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,901,'PPICLF:xGrid')
    END IF !global%error    

    DEALLOCATE(JFX,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,907,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFXCell,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,913,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFY,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,919,'PPICLF:xGrid')
    END IF ! global%error


    DEALLOCATE(JFYCell,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,926,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFZ,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,932,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFZCell,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,938,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFE,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,944,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFECell,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,950,'PPICLF:xGrid')
    END IF ! global%error    

    DEALLOCATE(domgdx,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,956,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(domgdy,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,962,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(domgdz,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,968,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(drhodx,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,974,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(drhody,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,980,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(drhodz,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,986,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpvxF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,992,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpvyF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,998,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpvzF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,1004,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDOX,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,1010,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDOY,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,1016,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDOZ,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,1022,'PPICLF:xGrid')
    END IF ! global%error


    DEALLOCATE(PhiP,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,1029,'PPICLF:xGrid')
    END IF ! global%error


    DEALLOCATE(ppF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,1036,'PPICLF:xGrid')
    END IF ! global%error
    DEALLOCATE(vfP,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,1041,'PPICLF:xGrid')
    END IF ! global%error
    DEALLOCATE(dpxF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,1046,'PPICLF:xGrid')
    END IF ! global%error
    DEALLOCATE(dpyF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,1051,'PPICLF:xGrid')
    END IF ! global%error
 
    DEALLOCATE(dpzF,STAT=errorFlag)
    global%error = errorFlag
    IF(global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,1057,'PPICLF:xGrid')
    END IF ! global%error
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

