










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
! Purpose: Shut down Rocflu-MP.
!
! Description: None.
!
! Input: 
!   levels      Data associated with levels
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_EndFlowSolver.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_EndFlowSolver(levels)

  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModDataTypes
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt_input   
  USE ModBndPatch, ONLY: t_patch
  USE ModError
  USE ModMPI  

  USE RFLU_ModABC 
  USE RFLU_ModAxisymmetry 
  USE RFLU_ModBFaceGradAccessList
  USE RFLU_ModBoundLists 
  USE RFLU_ModBoundXvUtils
  USE RFLU_ModNSCBC, ONLY: RFLU_NSCBC_DecideHaveNSCBC
  USE RFLU_ModCellMapping 
  USE RFLU_ModCommLists 
  USE RFLU_ModDimensions
  USE RFLU_ModEdgeList   
  USE RFLU_ModFaceList  
  USE RFLU_ModForcesMoments
  USE RFLU_ModGeometry      
  USE RFLU_ModGFM
  USE RFLU_ModGlobalIds
  USE RFLU_ModHypre
  USE RFLU_ModMPI
  USE RFLU_ModMovingFrame, ONLY: RFLU_MVF_DestroyPatchVelAccel
  USE RFLU_ModOLES
  USE RFLU_ModPatchCoeffs
  USE RFLU_ModProbes
  USE RFLU_ModReadWriteAuxVars
  USE RFLU_ModReadWriteFlow
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModReadWriteGridSpeeds  
  USE RFLU_ModRenumberings
  USE RFLU_ModStencilsBFaces
  USE RFLU_ModStencilsCells
  USE RFLU_ModStencilsFaces
  USE RFLU_ModWeights
  
  
  
  USE ModInterfaces, ONLY: RFLU_DeallocateMemoryWrapper, &
                           RFLU_DecideNeedBGradFace, &
                           RFLU_DecideNeedStencils, &
                           RFLU_DecideNeedWeights, &
                           RFLU_DestroyGrid, & 
                           RFLU_PrintFlowInfoWrapper, & 
                           RFLU_PrintGridInfo, & 
                           RFLU_PrintWarnInfo, &  
                           RFLU_WriteRestartInfo
    

  USE PICL_ModInterfaces, ONLY: PICL_TEMP_WriteVTU


  IMPLICIT NONE

! ******************************************************************************
! Arguments
! ******************************************************************************

  TYPE(t_level), POINTER :: levels(:)

! ******************************************************************************
! Locals
! ******************************************************************************

  CHARACTER(CHRLEN) :: RCSIdentString
  LOGICAL :: moveGrid  
  INTEGER :: errorFlag,iPatch,iProbe,iReg
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput 
  TYPE(t_patch), POINTER :: pPatch   
  TYPE(t_region), POINTER :: pRegion,pRegionSerial  

  INTEGER :: icg
  REAL(KIND=8) :: vFrac, ir

  REAL(KIND=8) :: timerStart
  REAL(KIND=8) :: timerEnd
  REAL(KIND=8) :: elapsedtime
  REAL(KIND=8) :: elapsedtime_hour

! ******************************************************************************
! Start, initialize some variables
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_EndFlowSolver.F90,v $ $Revision: 1.1.1.1 $'

  moveGrid = .FALSE. 

! ******************************************************************************
! Set global pointer and initialize global type, register function
! ****************************************************************************** 

  global => levels(1)%regions(1)%global

  CALL RegisterFunction(global,'RFLU_EndFlowSolver',"../rocflu/RFLU_EndFlowSolver.F90")

! ******************************************************************************
! Determine whether to print final soln or final CFL timestep
! ******************************************************************************
  !
  ! TLJ - 11/23/2024 - Modified final write solution files
  !       if printEndTime = 1, print solution at t=RUNTIME
  !                       = 2, print solution at final CFL timestep
  !
  IF (global%printEndTime==1) THEN

! ******************************************************************************
! Determine whether have moving grid
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => levels(1)%regions(iReg)
  
    IF ( pRegion%mixtInput%moveGrid .EQV. .TRUE. ) THEN 
      moveGrid = .TRUE.
      EXIT
    END IF ! regions
  END DO ! iReg

! ******************************************************************************
! Write grid file (if necessary) and flow file. Write restart info file after 
! flow (and grid) file so that incomplete flow (and grid) files due to
! exceeding time limit do not show up as iteration or time stamp in restart 
! info file.
! ******************************************************************************


  DO iReg = 1,global%nRegionsLocal
    pRegion => levels(1)%regions(iReg) ! single-level grids for now
    
    CALL RFLU_WriteDimensionsWrapper(pRegion,WRITE_DIMENS_MODE_MAYBE)    

    IF ( moveGrid .EQV. .TRUE. ) THEN
      CALL RFLU_WriteGridWrapper(pRegion)
      CALL RFLU_WriteGridSpeedsWrapper(pRegion)
    END IF ! moveGrid    
    
        IF ( global%piclUsed .EQV. .TRUE. ) THEN
        ! TLJ added to plot primitive variables - 02/19/2025
        DO icg = 1,pGrid%nCellsTot
           vFrac = 1.0_RFREAL - pRegion%mixt%piclVF(icg)
           ir = 1.0_RFREAL/pRegion%mixt%cv(CV_MIXT_DENS,icg)
           pRegion%mixt%cv(CV_MIXT_DENS,icg) = pRegion%mixt%cv(CV_MIXT_DENS,icg)/vFrac
           pRegion%mixt%cv(CV_MIXT_XMOM,icg) = ir*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
           pRegion%mixt%cv(CV_MIXT_YMOM,icg) = ir*pRegion%mixt%cv(CV_MIXT_YMOM,icg)
           pRegion%mixt%cv(CV_MIXT_ZMOM,icg) = ir*pRegion%mixt%cv(CV_MIXT_ZMOM,icg)
           pRegion%mixt%cv(CV_MIXT_ENER,icg) = ir*pRegion%mixt%cv(CV_MIXT_ENER,icg)
        ENDDO
        ENDIF

        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( global%piclUsed .EQV. .TRUE. ) THEN
        DO icg = 1,pGrid%nCellsTot
           vFrac = 1.0_RFREAL - pRegion%mixt%piclVF(icg)
           ir = vFrac*pRegion%mixt%cv(CV_MIXT_DENS,icg)
           pRegion%mixt%cv(CV_MIXT_DENS,icg) = ir
           pRegion%mixt%cv(CV_MIXT_XMOM,icg) = ir*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
           pRegion%mixt%cv(CV_MIXT_YMOM,icg) = ir*pRegion%mixt%cv(CV_MIXT_YMOM,icg)
           pRegion%mixt%cv(CV_MIXT_ZMOM,icg) = ir*pRegion%mixt%cv(CV_MIXT_ZMOM,icg)
           pRegion%mixt%cv(CV_MIXT_ENER,icg) = ir*pRegion%mixt%cv(CV_MIXT_ENER,icg)
        ENDDO
        ENDIF

    IF ( global%piclUsed .EQV. .TRUE. ) THEN
       CALL PICL_TEMP_WriteVTU(pRegion)
    END IF

    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
      CALL RFLU_WriteAuxVarsWrapper(pRegion)
    END IF ! solverType

    CALL RFLU_BXV_WriteVarsWrapper(pRegion)
    
    IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN
      CALL RFLU_WritePatchCoeffsWrapper(pRegion)
    END IF ! global%patchCoeffFlag
    
    
  END DO ! iReg

  CALL RFLU_WriteRestartInfo(global)
  

! ******************************************************************************
! Print (if necessary, grid, and) flow information 
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    IF ( moveGrid .EQV. .TRUE. ) THEN
      DO iReg = 1,global%nRegionsLocal
        pRegion => levels(1)%regions(iReg) ! single-level grids for now
        CALL RFLU_PrintGridInfo(pRegion)
      END DO ! iReg     
    END IF ! moveGrid       
       
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg) ! single-level grids for now
      CALL RFLU_PrintFlowInfoWrapper(pRegion)

    END DO ! iReg    
  END IF ! global%verbLevel


  ENDIF ! end global%printEndTime

! ******************************************************************************
! Deallocate memory. NOTE must be done before calling deallocation wrapper. 
! NOTE weights must be deallocated before stencils.
! ******************************************************************************

! ==============================================================================
! Weights for cell and face gradients
! ==============================================================================

  IF ( RFLU_DecideNeedWeights(pRegion) .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal 
      pRegion => levels(1)%regions(iReg)
      pMixtInput => pRegion%mixtInput   

      IF ( pMixtInput%spaceOrder > 1 ) THEN 
        CALL RFLU_DestroyWtsC2CWrapper(pRegion)
      END IF ! pMixtInput%spaceOrder      

      IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN      
        CALL RFLU_DestroyWtsF2CWrapper(pRegion)
      END IF ! pMixtInput%flowModel 
      
      DO iPatch = 1,pRegion%grid%nPatches
        pPatch => pRegion%patches(iPatch)

        IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
          CALL RFLU_DestroyWtsBF2CWrapper(pRegion,pPatch)      
        END IF ! RFLU_DecideNeedBGradFace
      END DO ! iPatch
    END DO ! iReg   
  END IF ! RFLU_DecideNeedWeights

! ==============================================================================
! Weights for optimal LES approach
! ==============================================================================

  IF ( RFLU_DecideNeedWeights(pRegion) .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg) 
      pMixtInput => pRegion%mixtInput 
         
      IF ( pMixtInput%spaceDiscr == DISCR_OPT_LES ) THEN 
        CALL RFLU_DestroyStencilsWeightsOLES(pRegion)
      END IF ! pMixtInput
    END DO ! iReg
  END IF ! RFLU_DecideNeedWeights

! ******************************************************************************
! Descale geometry used in axisymmtry flow
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => levels(1)%regions(iReg)

    IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
      CALL RFLU_AXI_DescaleGeometry(pRegion)
    END IF ! pRegion%mixtInput%axiFlag
  END DO ! iReg
 
! ******************************************************************************
! Deallocate memory for stencils
! ******************************************************************************

  IF ( RFLU_DecideNeedStencils(pRegion) .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal 
      pRegion => levels(1)%regions(iReg)
      pMixtInput => pRegion%mixtInput 
    
      IF ( pMixtInput%spaceOrder > 1 ) THEN 
        CALL RFLU_DestroyC2CStencilWrapper(pRegion)
        CALL RFLU_DestroyListCC2CStencil(pRegion)  
      END IF ! pMixtInput%spaceOrder
    
      IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN      
        CALL RFLU_DestroyF2CStencilWrapper(pRegion) 
        CALL RFLU_DestroyListCF2CStencil(pRegion)
      END IF ! pMixtInput%flowModel

      DO iPatch = 1,pRegion%grid%nPatches
        pPatch => pRegion%patches(iPatch)
    
        IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
          CALL RFLU_DestroyBF2CStencilWrapper(pRegion,pPatch)
        END IF ! RFLU_DecideNeedBGradFace
      END DO ! iPatch
    END DO ! iReg
  END IF ! RFLU_DecideNeedStencils

! ******************************************************************************
! Destroy boundary-face gradient access list for viscous flows
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal 
    pRegion => levels(1)%regions(iReg) 
    
    IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN          
      CALL RFLU_DestroyBFaceGradAccessList(pRegion)
    END IF ! pRegion%mixtInput
  END DO ! iReg  

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => levels(1)%regions(iReg) ! single-level grids for now
    pGrid   => pRegion%grid
        
    CALL RFLU_DestroyPatchCoeffs(pRegion)

    IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
      CALL RFLU_BXV_DestroyVarsCv(pRegion)
      CALL RFLU_BXV_DestroyVarsDv(pRegion)
      CALL RFLU_BXV_DestroyVarsTStep(pRegion)
    END IF ! RFLU_NSCBC_DecideHaveNSCBC(pRegion)   

    IF ( global%forceFlag .EQV. .TRUE. ) THEN 
      CALL RFLU_DestroyForcesMoments(pRegion)   
      CALL RFLU_DestroyGlobalThrustFlags(pRegion)   
    END IF ! global%forceFlag

    IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
      CALL RFLU_MVF_DestroyPatchVelAccel(pRegion)
    END IF ! global%mvFrameFlag

    IF ( global%gfmFlag .EQV. .TRUE. ) THEN
      IF ( global%nRegionsLocal > 1 ) THEN
        CALL RFLU_RNMB_DestroyPC2SCMap(pRegion)
        CALL RFLU_RNMB_DestroySC2PCMap(pRegion)
      END IF ! global%nRegionsLocal > 1

      CALL RFLU_GFM_DestroyLevelSet(pRegion)
    END IF ! global%gfmFlag

            
    IF ( pGrid%nBorders > 0 ) THEN 
      CALL RFLU_MPI_DestroyBuffersWrapper(pRegion)


      CALL RFLU_COMM_DestroyCommLists(pRegion)
      CALL RFLU_COMM_DestroyBorders(pRegion)
    END IF ! pGrid%nBorders 
        
    CALL RFLU_DeallocateMemoryWrapper(pRegion)

    CALL RFLU_DestroyCellMapping(pRegion)
    CALL RFLU_DestroyFaceList(pRegion)
    
    
    IF ( pRegion%mixtInput%movegrid .EQV. .TRUE. ) THEN
      CALL RFLU_DestroyEdgeList(pRegion)
    END IF ! pRegion 
    
    CALL RFLU_DestroyGeometry(pRegion)
    CALL RFLU_DestroyGrid(pRegion)           
  END DO ! iReg

! ******************************************************************************
! Deallocate memory for absorbing boundary condition
! ******************************************************************************

  IF ( global%abcFlag .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg) ! single-level grids for now
      pMixtInput => pRegion%mixtInput

      IF ( global%abcKind == 0 ) THEN
        CALL RFLU_ABC_DestroySigma(pRegion)
      END IF ! global%abcKind
    END DO ! iReg
  END IF ! global%abcFlag

! ******************************************************************************
! Deatroy non-dissipative solver specific arrays
! ******************************************************************************
  
  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN

! ==============================================================================
! Destroy cell number offsets in each region
! ==============================================================================

    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg) ! single-level grids for now

      CALL RFLU_GID_DestroynCellsOffset(pRegion)
    END DO ! iReg

! ==============================================================================
! Destroy global cell numbers of virtual cells.
! ==============================================================================
 
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg) ! single-level grids for now
      pGrid   => pRegion%grid
  
      IF ( pGrid%nBorders > 0 ) THEN
        CALL RFLU_GID_DestroyGlobalIds(pRegion)
      END IF ! pRegion%grid%nBorders
    END DO ! iReg

! ******************************************************************************
! Destroy Hypre objects neeeded for SOLV_IMPLICIT_HM
! ******************************************************************************

    CALL RFLU_HYPRE_DestroyObjects(levels(1)%regions)
  END IF ! global%solverType

! ******************************************************************************
! Close convergence, mass, and probe files
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC ) THEN
    CLOSE(IF_CONVER,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= 0 ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,579)
    END IF ! global%error
  END IF ! global%myProcid

  IF ( (global%myProcid == MASTERPROC) .AND. (moveGrid .EQV. .TRUE.) ) THEN
    CLOSE(IF_MASS,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= 0 ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,587)
    END IF ! global%error
  END IF ! global%myProcid  

  IF ( global%nProbes > 0 ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg) ! single-level grids for now 
         
      CALL RFLU_CloseProbeFiles(pRegion)
    END DO ! iReg
  END IF ! global%nProbes

!begin BBR
    IF ( global%myProcid == MASTERPROC ) THEN
    CLOSE(IF_PM,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= 0 ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,604)
    END IF ! global%error
  END IF ! global%myProcid

!  IF ( global%myProcid == MASTERPROC ) THEN
!    CLOSE(IF_INTEG,IOSTAT=errorFlag)
!    global%error = errorFlag
!    IF ( global%error /= 0 ) THEN
!      CALL ErrorStop(global,ERR_FILE_CLOSE,612)
!    END IF ! global%error
!  END IF ! global%myProcid
!end BBR

! ******************************************************************************
! Print Simulation (rflump only) Run Time
! ******************************************************************************

  IF(global%myProcid == MASTERPROC) THEN
     timerStart  = global%timingSubRout(1)
     timerEnd    = MPI_Wtime()
     elapsedtime = (timerEnd - timerStart) / 60.0_RFREAL
     elapsedtime_hour = elapsedtime / 60.0_RFREAL
     print*,"*** Total rflump timing = ",elapsedtime," minutes"
     print*,"*** Total rflump timing = ",elapsedtime_hour," hours"
  ENDIF

! ******************************************************************************
! Print info about warnings
! ******************************************************************************

  CALL RFLU_PrintWarnInfo(global)

! ******************************************************************************
! Deallocate PETSc memory & finalize PETSc
! ******************************************************************************


! ******************************************************************************
! Finalize Rocprof
! ****************************************************************************** 


! ******************************************************************************
! Finalize MPI
! ****************************************************************************** 

  CALL MPI_Finalize(errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_OUTPUT,664)
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel /= VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A)') SOLVER_NAME    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Finalization done.'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Program finished.'
    WRITE(STDOUT,'(A)') SOLVER_NAME         
  END IF ! global%myProcid 


  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_EndFlowSolver

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_EndFlowSolver.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.10  2009/07/08 20:53:54  mparmar
! Removed RFLU_ModHouMaheshBoundCond
!
! Revision 1.9  2009/07/08 19:12:13  mparmar
! Added deallocation for absorbing layer
!
! Revision 1.8  2008/12/06 08:43:48  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2008/03/27 12:13:25  haselbac
! Added axisymmetry capability
!
! Revision 1.5  2008/01/19 20:19:45  haselbac
! Added calls to PLAG_DecideHaveSurfStats
!
! Revision 1.4  2007/12/03 16:34:47  mparmar
! Removed RFLU_DestroyPatchVelocity
!
! Revision 1.3  2007/11/28 23:05:28  mparmar
! Deallocating SOLV_IMPLICIT_HM related arrays
!
! Revision 1.2  2007/06/18 18:09:08  mparmar
! Added closing of moving reference frame
!
! Revision 1.1  2007/04/09 18:49:57  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.60  2007/03/31 23:53:39  haselbac
! Added calls to determine, write, and print nPclsGlobal
!
! Revision 1.59  2006/10/20 21:32:08  mparmar
! Added call to RFLU_DestroyGlobalThrustFlags
!
! Revision 1.58  2006/08/19 15:46:41  mparmar
! Changed logic for NSCBC and added calls to deallocate patch arrays
!
! Revision 1.57  2006/08/18 14:04:02  haselbac
! Added call to destroy AVFace2Patch list
!
! Revision 1.56  2006/04/07 16:04:03  haselbac
! Adapted to changes in bf2c wts computation
!
! Revision 1.55  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.54  2006/04/07 14:53:11  haselbac
! Adapted to changes in bface stencil routines
!
! Revision 1.53  2006/03/09 14:09:49  haselbac
! Now call wrapper routines for stencils
!
! Revision 1.52  2006/01/06 22:15:07  haselbac
! Adapted to name changes
!
! Revision 1.51  2005/11/10 16:51:29  fnajjar
! Added plagUsed IF statement around PLAG routines
!
! Revision 1.50  2005/10/27 19:20:06  haselbac
! Adapted to changes in stencil routine names
!
! Revision 1.49  2005/10/25 19:39:23  haselbac
! Added IF on forceFlag
!
! Revision 1.48  2005/10/05 14:18:41  haselbac
! Adapted to changes in stencil modules, added call to destroy bface wts
!
! Revision 1.47  2005/09/14 15:59:33  haselbac
! Minor clean-up
!
! Revision 1.46  2005/09/13 20:44:06  mtcampbe
! Added Rocprof finalization
!
! Revision 1.45  2005/08/09 00:59:42  haselbac
! Enclosed writing of patch coeffs within IF (patchCoeffFlag)
!
! Revision 1.44  2005/08/03 18:28:11  hdewey2
! Enclosed PETSc deallocation calls within IF
!
! Revision 1.43  2005/08/02 18:24:32  hdewey2
! Added PETSc support
!
! Revision 1.42  2005/05/18 22:12:55  fnajjar
! ACH: Added destruction of iPclSend buffers, now use nFacesAV
!
! Revision 1.41  2005/04/29 23:03:22  haselbac
! Added destruction of avf2b list
!
! Revision 1.40  2005/04/29 12:48:23  haselbac
! Updated closing of probe files to changes in probe handling
!
! Revision 1.39  2005/04/15 16:31:18  haselbac
! Removed calls to XyzEdge2RegionDegrList routines
!
! Revision 1.38  2005/04/15 15:07:15  haselbac
! Converted to MPI
!
! Revision 1.37  2005/01/18 15:18:19  haselbac
! Commented out COMM calls for now
!
! Revision 1.36  2005/01/14 21:33:56  haselbac
! Added calls to destroy comm lists and borders
!
! Revision 1.35  2005/01/03 16:14:57  haselbac
! Added call to destroy bface stencils
!
! Revision 1.34  2004/12/21 15:05:11  fnajjar
! Included calls for PLAG surface statistics
!
! Revision 1.33  2004/10/19 19:29:15  haselbac
! Adapted to GENX changes
!
! Revision 1.32  2004/07/06 15:14:50  haselbac
! Adapted to changes in libflu and modflu, cosmetics
!
! Revision 1.31  2004/06/16 20:01:08  haselbac
! Added writing of patch coeffs, destruction of memory
!
! Revision 1.30  2004/03/17 04:28:26  haselbac
! Adapted call to RFLU_WriteDimensionsWrapper
!
! Revision 1.29  2004/03/11 16:33:22  fnajjar
! ACH: Moved call to RFLU_WriteDimWrapper outside if bcos of Rocpart
!
! Revision 1.28  2004/03/08 22:01:43  fnajjar
! ACH: Changed call to RFLU_WriteDimensionsWrapper so PLAG data also written
!
! Revision 1.27  2004/01/29 22:59:19  haselbac
! Removed hardcoded error computation for supersonic vortex
!
! Revision 1.26  2003/12/04 03:30:01  haselbac
! Added destruction calls for gradients, cleaned up
!
! Revision 1.25  2003/11/25 21:04:40  haselbac
! Added call to RFLU_PrintFlowInfoWrapper, cosmetic changes
!
! Revision 1.24  2003/11/03 03:51:17  haselbac
! Added call to destroy boundary-face gradient access list
!
! Revision 1.23  2003/08/07 15:33:30  haselbac
! Added call to RFLU_PrintWarnInfo
!
! Revision 1.22  2003/07/22 02:07:00  haselbac
! Added writing out of warnings
!
! Revision 1.21  2003/06/20 22:35:26  haselbac
! Added call to RFLU_WriteRestartInfo
!
! Revision 1.20  2003/01/28 14:36:53  haselbac
! Added deallocation, merged printing of grid and flow info, cosmetics
!
! Revision 1.19  2002/12/20 23:20:06  haselbac
! Fixed output bug: no output for verbosity=0
!
! Revision 1.18  2002/11/08 21:31:32  haselbac
! Added closing of total-mass file
!
! Revision 1.17  2002/11/02 02:03:20  wasistho
! Added TURB statistics
!
! Revision 1.16  2002/10/27 19:12:09  haselbac
! Added writing of grid for moving grid calcs
!
! Revision 1.15  2002/10/19 16:12:16  haselbac
! Cosmetic changes to output
!
! Revision 1.14  2002/10/16 21:16:06  haselbac
! Added writing of header when running
!
! Revision 1.13  2002/10/12 14:57:37  haselbac
! Enclosed RFLU_WriteFlowWrapper between ifndef GENX
!
! Revision 1.12  2002/10/08 15:49:29  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.11  2002/10/05 19:20:55  haselbac
! GENX integration, close probe file, use flow wrapper routine
!
! Revision 1.10  2002/09/09 15:49:26  haselbac
! global and mixtInput now under region
!
! Revision 1.9  2002/07/25 14:25:07  haselbac
! No longer called finalize for CHARM=1
!
! Revision 1.8  2002/06/27 15:27:11  haselbac
! Change name if CHARM defined, comment out destruction for now (crashes with CHARM)
!
! Revision 1.7  2002/06/17 13:34:12  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.6  2002/06/14 22:26:08  wasistho
! update statistics
!
! Revision 1.5  2002/06/14 21:54:35  wasistho
! Added time avg statistics
!
! Revision 1.4  2002/06/14 20:20:43  haselbac
! Deleted ModLocal, changed local%nRegions to global%nRegionsLocal, added destroy flag
!
! Revision 1.3  2002/05/04 17:08:13  haselbac
! Close convergence file and print flow information
!
! Revision 1.2  2002/04/11 19:03:22  haselbac
! Added calls and cosmetic changes
!
! Revision 1.1  2002/03/14 19:11:00  haselbac
! Initial revision
!
! ******************************************************************************

