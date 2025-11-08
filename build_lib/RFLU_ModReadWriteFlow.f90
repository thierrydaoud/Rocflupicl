










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
! Purpose: Suite of routines to read and write mixture solution files.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModReadWriteFlow.F90,v 1.6 2016/01/31 04:54:40 rahul Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModReadWriteFlow

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMPI
  USE ModGrid, ONLY: t_grid

  USE ModBuildFileNames, ONLY: BuildFileNameSteady, &
                               BuildFileNameUnsteady, &
                               BuildFileNameSteadyVTK, &
                               BuildFileNameUnsteadyVTK, &
                               BuildFileNameSteadyPVTK, &
                               BuildFileNameUnsteadyPVTK

  ! TLJ - added for species plotting
  USE RFLU_ModConvertCv,   ONLY: RFLU_ConvertCvCons2Prim, &
                                 RFLU_ConvertCvPrim2Cons, &
                                 RFLU_ScalarConvertCvCons2Prim, &
                                 RFLU_ScalarConvertCvPrim2Cons
  USE RFLU_ModJWL

  USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
  USE ModSpecies, ONLY: t_spec_input, &
                        t_spec_type
 !Fred - Adding Paraview output for species fractions - 11/1/17


  USE ModInterfaces, ONLY: RFLU_GetCvLoc

  USE RFLU_ModPlottingVars

  USE RFLU_ModIRPrecision
  USE RFLU_ModLibVTKIO
!  USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout=>OUTPUT_UNIT, stderr=>ERROR_UNIT

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_ReadFlowWrapper, &
            RFLU_WriteFlowWrapper, &
            !begin BBR
            RFLU_WriteVTK , &
            RFLU_WriteSlice, &
            RFLU_WriteRProf
            !end BBR

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN) :: RCSIdentString = &
    '$RCSfile: RFLU_ModReadWriteFlow.F90,v $ $Revision: 1.6 $'

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Read flow dataset for mixture in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iFile       File number
!   nVars       Number of variables to be read
!
! Output: None.
!
! Notes:
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadFlowV1ASCII(pRegion,iFile,nVars)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iFile,nVars
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString,sectionString
    INTEGER :: cvMixtDens,cvMixtEner,cvMixtPres,cvMixtXMom,cvMixtXVel, &
               cvMixtYMom,cvMixtYVel,cvMixtZMom,cvMixtZVel,errorFlag,i, &
               iVars,j,loopCounter
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, allocate pointers
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadFlowV1ASCII',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII flow dataset...'
    END IF ! global%verbLevel

    pGrid => pRegion%grid

! ******************************************************************************
!   File is already open, just read the dataset
! ******************************************************************************

    iVars       = 0
    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile,'(A)') sectionString

      SELECT CASE ( TRIM(sectionString) )

! ------------------------------------------------------------------------------
!       Mixture density - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture density' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture density...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtDens = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_DENS)

          READ(iFile,'(5(E23.16))') (pCv(cvMixtDens,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture x-momentum - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture x-momentum' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture x-momentum...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtXMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_XMOM)

          READ(iFile,'(5(E23.16))') (pCv(cvMixtXMom,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture x-velocity - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture x-velocity' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture x-velocity...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtXVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_XVEL)

          READ(iFile,'(5(E23.16))') (pCv(cvMixtXVel,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture y-momentum - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture y-momentum' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture y-momentum...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtYMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_YMOM)

          READ(iFile,'(5(E23.16))') (pCv(cvMixtYMom,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture y-velocity - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture y-velocity' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture y-velocity...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtYVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_YVEL)

          READ(iFile,'(5(E23.16))') (pCv(cvMixtYVel,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture z-momentum - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture z-momentum' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture z-momentum...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtZMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_ZMOM)

          READ(iFile,'(5(E23.16))') (pCv(cvMixtZMom,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture z-velocity - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture z-velocity' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture z-velocity...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtZVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_ZVEL)

          READ(iFile,'(5(E23.16))') (pCv(cvMixtZVel,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture total internal energy - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture total internal energy' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture total internal '// &
                                     'energy...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtEner = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_ENER)

          READ(iFile,'(5(E23.16))') (pCv(cvMixtEner,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture pressure - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture pressure' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture pressure...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtPres = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_PRES)

          READ(iFile,'(5(E23.16))') (pCv(cvMixtPres,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------

        CASE ( '# End' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel

          EXIT

! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------

        CASE DEFAULT
          IF ( global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
          END IF ! verbosityLevel

          CALL ErrorStop(global,ERR_INVALID_MARKER,396,sectionString)

      END SELECT ! TRIM

! ------------------------------------------------------------------------------
!     Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN
        CALL ErrorStop(global,ERR_INFINITE_LOOP,405)
      END IF ! loopCounter
    END DO ! <empty>

! ==============================================================================
!   Check and information about number of variables read
! ==============================================================================

    IF ( iVars /= nVars ) THEN
      CALL ErrorStop(global,ERR_INVALID_NVARS,414)
    END IF ! iVar

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII flow dataset done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadFlowV1ASCII






! ******************************************************************************
!
! Purpose: Read flow dataset for mixture in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iFile       File number
!   nVars       Number of variables to be read
!
! Output: None.
!
! Notes:
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadFlowV1Binary(pRegion,iFile,nVars)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iFile,nVars
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString,sectionString
    INTEGER :: cvMixtDens,cvMixtEner,cvMixtPres,cvMixtXMom,cvMixtXVel, &
               cvMixtYMom,cvMixtYVel,cvMixtZMom,cvMixtZVel,errorFlag,i, &
               iVars,j,loopCounter
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, allocate pointers
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadFlowV1Binary',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary flow dataset...'
    END IF ! global%verbLevel

    pGrid => pRegion%grid

! ******************************************************************************
!   File is already open, just read the dataset
! ******************************************************************************

    iVars       = 0
    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile) sectionString

      SELECT CASE ( TRIM(sectionString) )

! ------------------------------------------------------------------------------
!       Mixture density - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture density' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture density...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtDens = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_DENS)

          READ(iFile) (pCv(cvMixtDens,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture x-momentum - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture x-momentum' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture x-momentum...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtXMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_XMOM)

          READ(iFile) (pCv(cvMixtXMom,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture x-velocity - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture x-velocity' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture x-velocity...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtXVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_XVEL)

          READ(iFile) (pCv(cvMixtXVel,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture y-momentum - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture y-momentum' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture y-momentum...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtYMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_YMOM)

          READ(iFile) (pCv(cvMixtYMom,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture y-velocity - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture y-velocity' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture y-velocity...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtYVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_YVEL)

          READ(iFile) (pCv(cvMixtYVel,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture z-momentum - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture z-momentum' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture z-momentum...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtZMom = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_ZMOM)

          READ(iFile) (pCv(cvMixtZMom,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture z-velocity - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture z-velocity' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture z-velocity...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtZVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_ZVEL)

          READ(iFile) (pCv(cvMixtZVel,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture total internal energy - compressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture total internal energy' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture total internal '// &
                                     'energy...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtEner = RFLU_GetCvLoc(global,FLUID_MODEL_COMP,CV_MIXT_ENER)

          READ(iFile) (pCv(cvMixtEner,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       Mixture pressure - incompressible solver
! ------------------------------------------------------------------------------

        CASE ( '# Mixture pressure' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture pressure...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          cvMixtPres = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_PRES)

          READ(iFile) (pCv(cvMixtPres,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------

        CASE ( '# End' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel

          EXIT

! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------

        CASE DEFAULT
          IF ( global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
          END IF ! verbosityLevel

          CALL ErrorStop(global,ERR_INVALID_MARKER,692,sectionString)

      END SELECT ! TRIM

! ------------------------------------------------------------------------------
!     Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN
        CALL ErrorStop(global,ERR_INFINITE_LOOP,701)
      END IF ! loopCounter

    END DO ! <empty>

! ==============================================================================
!   Check and information about number of variables read
! ==============================================================================

    IF ( iVars /= nVars ) THEN
      CALL ErrorStop(global,ERR_INVALID_NVARS,711)
    END IF ! iVar

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary flow dataset done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadFlowV1Binary







! ******************************************************************************
!
! Purpose: Read flow dataset (v2) for mixture in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iFile       File number
!   nVars       Number of variables to be read
!
! Output: None.
!
! Notes:
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadFlowV2ASCII(pRegion,iFile,nVars)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iFile,nVars
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString,sectionString
    INTEGER :: errorFlag,i,iLoc,iVars,j,loopCounter
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, allocate pointers
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadFlowV2ASCII',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII flow dataset ' // &
                                           '(v2)...'
    END IF ! global%verbLevel

    pGrid => pRegion%grid

! ******************************************************************************
!   File is already open, just read the dataset
! ******************************************************************************

    iVars       = 0
    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile,'(A)') sectionString

! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------

      IF ( sectionString(1:2) /= '# ' ) THEN
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel

        CALL ErrorStop(global,ERR_INVALID_MARKER,812,sectionString)
      END IF ! sectionString(1:2)

      sectionString = sectionString(3:)

      SELECT CASE ( TRIM(sectionString) )

! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------

        CASE ( 'End' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel

          EXIT

! ------------------------------------------------------------------------------
!       Extract dataset location in cv from section string
! ------------------------------------------------------------------------------

        CASE DEFAULT
          READ(sectionString,'(I2)') iLoc

          IF ( iLoc < 1 .AND. iLoc > pRegion%mixtInput%nCv ) THEN
            IF ( global%verbLevel > VERBOSE_LOW ) THEN
              sectionString = '# ' // sectionString
              WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
            END IF ! verbosityLevel

            CALL ErrorStop(global,ERR_INVALID_CVILOC,844)
          END IF ! iLoc

          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A,I2,A)') SOLVER_NAME,'Cv location ',iLoc,'...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          READ(iFile,'(5(E23.16))') (pCv(iLoc,j),j=1,pGrid%nCellsTot)
      END SELECT ! TRIM

! ------------------------------------------------------------------------------
!     Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN
        CALL ErrorStop(global,ERR_INFINITE_LOOP,864)
      END IF ! loopCounter
    END DO ! <empty>

! ==============================================================================
!   Check and information about number of variables read
! ==============================================================================

    IF ( iVars /= nVars ) THEN
      CALL ErrorStop(global,ERR_INVALID_NVARS,873)
    END IF ! iVar

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII flow dataset ' // &
                                           '(v2) done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadFlowV2ASCII






! ******************************************************************************
!
! Purpose: Read flow dataset (v2) for mixture in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iFile       File number
!   nVars       Number of variables to be read
!
! Output: None.
!
! Notes:
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadFlowV2Binary(pRegion,iFile,nVars)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iFile,nVars
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString,sectionString
    INTEGER :: errorFlag,i,iLoc,iVars,j,loopCounter
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, allocate pointers
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadFlowV2Binary',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary flow dataset ' // &
                                           '(v2)...'
    END IF ! global%verbLevel

    pGrid => pRegion%grid

! ******************************************************************************
!   File is already open, just read the dataset
! ******************************************************************************

    iVars       = 0
    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile) sectionString

! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------

      IF ( sectionString(1:2) /= '# ' ) THEN
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel

        CALL ErrorStop(global,ERR_INVALID_MARKER,974,sectionString)
      END IF ! sectionString(1:2)

      sectionString = sectionString(3:)

      SELECT CASE ( TRIM(sectionString) )

! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------

        CASE ( 'End' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel

          EXIT

! ------------------------------------------------------------------------------
!       Extract dataset location in cv from section string
! ------------------------------------------------------------------------------

        CASE DEFAULT
          READ(sectionString,'(I2)') iLoc

          IF ( iLoc < 1 .AND. iLoc > pRegion%mixtInput%nCv ) THEN
            IF ( global%verbLevel > VERBOSE_LOW ) THEN
              sectionString = '# ' // sectionString
              WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
            END IF ! verbosityLevel

            CALL ErrorStop(global,ERR_INVALID_CVILOC,1006)
          END IF ! iLoc

          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A,I2,A)') SOLVER_NAME,'Cv location ',iLoc,'...'
          END IF ! global%verbLevel

          pCv => pRegion%mixt%cv

          iVars = iVars + 1

          READ(iFile) (pCv(iLoc,j),j=1,pGrid%nCellsTot)
      END SELECT ! TRIM

! ------------------------------------------------------------------------------
!     Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN
        CALL ErrorStop(global,ERR_INFINITE_LOOP,1026)
      END IF ! loopCounter

    END DO ! <empty>

! ==============================================================================
!   Check and information about number of variables read
! ==============================================================================

    IF ( iVars /= nVars ) THEN
      CALL ErrorStop(global,ERR_INVALID_NVARS,1036)
    END IF ! iVar

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary flow dataset ' // &
                                           '(v2) done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadFlowV2Binary








! ******************************************************************************
!
! Purpose: Read flow file for mixture in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Read initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadFlowASCII(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    LOGICAL :: fileExists
    CHARACTER(CHRLEN) :: errorString,iFileName,iFileNameOld,sectionString, &
                         timeString1,timeString2
    INTEGER :: errorFlag,i,iDataSet,iFile,j,nCellsTot,nCellsExpected,nVars, &
               nVarsExpected,precActual,precExpected,rangeActual,rangeExpected
    REAL(RFREAL) :: currentTime
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadFlowASCII',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII flow file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.mixt.cva', & 
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName) 

      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.floa', & 
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileNameOld)
                                
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.mixt.cva', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)      

      CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.floa', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileNameOld)                             

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT
       
    INQUIRE(FILE=iFileName,EXIST=fileExists)
    
    IF ( fileExists .EQV. .TRUE. ) THEN 
      OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", &
           IOSTAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_FILE_OPEN,1162,iFileName)
      END IF ! global%error
    ELSE 
      OPEN(iFile,FILE=iFileNameOld,FORM="FORMATTED",STATUS="OLD", &
           IOSTAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_FILE_OPEN,1169,iFileNameOld)
      END IF ! global%error    
    END IF ! fileExists

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) == '# ROCFLU flow file' ) THEN
      iDataSet = CV_DATASET_OLD
    ELSE IF ( TRIM(sectionString) == '# ROCFLU flow file (v2)' ) THEN
      iDataSet = CV_DATASET_NEW
    ELSE
      CALL ErrorStop(global,ERR_INVALID_MARKER,1188,sectionString)
    END IF ! TRIM

! -----------------------------------------------------------------------------
!   Precision and range
! -----------------------------------------------------------------------------

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Precision and range' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,1197,sectionString)
    END IF ! TRIM

    precExpected  = PRECISION(1.0_RFREAL)
    rangeExpected = RANGE(1.0_RFREAL)

    READ(iFile,'(2(I8))') precActual,rangeActual
    IF ( precActual < precExpected .OR. rangeActual < rangeExpected ) THEN
      CALL ErrorStop(global,ERR_PREC_RANGE,1205)
    END IF ! precActual

! -----------------------------------------------------------------------------
!   Initial residual and physical time
! -----------------------------------------------------------------------------

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Initial residual' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,1214,iFileName)
    END IF ! TRIM

    READ(iFile,'(E23.16)') global%resInit

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Physical time' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,1221,iFileName)
    END IF ! TRIM

    READ(iFile,'(E23.16)') currentTime

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      IF ( global%currentTime < 0.0_RFREAL ) THEN
        global%currentTime = currentTime
      ELSE
        WRITE(timeString1,'(1PE11.5)') global%currentTime
        WRITE(timeString2,'(1PE11.5)') currentTime
        
        IF ( TRIM(timeString1) /= TRIM(timeString2) ) THEN
          WRITE(STDOUT,'(A,4(1X,A))') SOLVER_NAME, & 
                                      '*** WARNING *** Time mismatch:', & 
                                      TRIM(timeString1),'vs.',TRIM(timeString2)
                                      
          global%warnCounter = global%warnCounter + 1
        END IF ! global%currentTime
      END IF ! global%currentTime
    END IF ! global%flowType  

! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    nVarsExpected  = pRegion%mixtInput%nCv
    nCellsExpected = pGrid%nCellsTot

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,1254,sectionString)
    END IF ! TRIM
  
    READ(iFile,'(2(I16))') nCellsTot,nVars
    IF ( nCellsTot /= nCellsExpected ) THEN
      WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nCellsTot, &
                                                'but expected:',nCellsExpected
      CALL ErrorStop(global,ERR_INVALID_NCELLS,1261,errorString)
    END IF ! nCellsExpected

    IF ( nVars /= nVarsExpected ) THEN
      WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nVars, &
                                                'but expected:',nVarsExpected
      CALL ErrorStop(global,ERR_INVALID_NVARS,1267)
    END IF ! nVarsExpected

! ==============================================================================
!   Rest of file
! ==============================================================================

    IF ( iDataSet == CV_DATASET_OLD ) THEN
      CALL RFLU_ReadFlowV1ASCII(pRegion,iFile,nVars)
    ELSE
      CALL RFLU_ReadFlowV2ASCII(pRegion,iFile,nVars)
    END IF ! iDataSet

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,1287,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII flow file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadFlowASCII







! ******************************************************************************
!
! Purpose: Read flow file for mixture in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Read initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadFlowBinary(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    LOGICAL :: fileExists
    CHARACTER(CHRLEN) :: errorString,iFileName,iFileNameOld,sectionString, &
                         timeString1,timeString2
    INTEGER :: errorFlag,i,iDataSet,iFile,j,nCellsTot,nCellsExpected,nVars, &
               nVarsExpected,precActual,precExpected,rangeActual,rangeExpected
    REAL(RFREAL) :: currentTime
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadFlowBinary',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary flow file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.mixt.cv', & 
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.flo', & 
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileNameOld)                                 

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.mixt.cv', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)
      CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.flo', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileNameOld)                               

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT
       
    INQUIRE(FILE=iFileName,EXIST=fileExists)
    
    IF ( fileExists .EQV. .TRUE. ) THEN 
!      OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
!           IOSTAT=errorFlag)
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
      OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
           IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_L )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
         ACCESS="SEQUENTIAL",CONVERT="LITTLE_ENDIAN",IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_B )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
         ACCESS="SEQUENTIAL",CONVERT="BIG_ENDIAN",IOSTAT=errorFlag)
    END IF
! BBR - end
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_FILE_OPEN,1421,iFileName)
      END IF ! global%error
    ELSE 
!      OPEN(iFile,FILE=iFileNameOld,FORM="UNFORMATTED",STATUS="OLD", &
!           IOSTAT=errorFlag)
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
      OPEN(iFile,FILE=iFileNameOld,FORM="UNFORMATTED",STATUS="OLD", &
           IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_L )THEN
      OPEN(iFile,FILE=iFileNameOld,FORM="UNFORMATTED",STATUS="OLD", &
           ACCESS="SEQUENTIAL",CONVERT="LITTLE_ENDIAN",IOSTAT=errorFlag)    
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_B )THEN
      OPEN(iFile,FILE=iFileNameOld,FORM="UNFORMATTED",STATUS="OLD", &
           ACCESS="SEQUENTIAL",CONVERT="BIG_ENDIAN",IOSTAT=errorFlag)
    END IF
! BBR - end
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_FILE_OPEN,1440,iFileNameOld)
      END IF ! global%error    
    END IF ! fileExists

! ==============================================================================
! Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile) sectionString
    IF ( TRIM(sectionString) == '# ROCFLU flow file' ) THEN
      iDataSet = CV_DATASET_OLD
    ELSE IF ( TRIM(sectionString) == '# ROCFLU flow file (v2)' ) THEN
      iDataSet = CV_DATASET_NEW
    ELSE
      CALL ErrorStop(global,ERR_INVALID_MARKER,1459,sectionString)
    END IF ! TRIM

! -----------------------------------------------------------------------------
!   Precision and range
! -----------------------------------------------------------------------------

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Precision and range' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,1468,sectionString)
    END IF ! TRIM

    precExpected  = PRECISION(1.0_RFREAL)
    rangeExpected = RANGE(1.0_RFREAL)

    READ(iFile) precActual,rangeActual
    IF ( precActual < precExpected .OR. rangeActual < rangeExpected ) THEN
      CALL ErrorStop(global,ERR_PREC_RANGE,1476)
    END IF ! precActual

! -----------------------------------------------------------------------------
!   Initial residual and physical time
! -----------------------------------------------------------------------------

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Initial residual' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,1485,iFileName)
    END IF ! TRIM

    READ(iFile) global%resInit

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Physical time' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,1492,iFileName)
    END IF ! TRIM

    READ(iFile) currentTime
      
    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      IF ( global%currentTime < 0.0_RFREAL ) THEN
        global%currentTime = currentTime
      ELSE
        WRITE(timeString1,'(1PE11.5)') global%currentTime
        WRITE(timeString2,'(1PE11.5)') currentTime

        IF ( TRIM(timeString1) /= TRIM(timeString2) ) THEN
          WRITE(STDOUT,'(A,4(1X,A))') SOLVER_NAME, & 
                                      '*** WARNING *** Time mismatch:', & 
                                      TRIM(timeString1),'vs.',TRIM(timeString2)
                                      
          global%warnCounter = global%warnCounter + 1
        END IF ! global%currentTime
      END IF ! global%currentTime
    END IF ! global%flowType  

! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    nVarsExpected  = pRegion%mixtInput%nCv
    nCellsExpected = pGrid%nCellsTot

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,1525,sectionString)
    END IF ! TRIM

    READ(iFile) nCellsTot,nVars

    IF ( nCellsTot /= nCellsExpected ) THEN
      WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nCellsTot, &
                                                'but expected:',nCellsExpected
      CALL ErrorStop(global,ERR_INVALID_NCELLS,1533,errorString)
    END IF ! nCellsExpected

    IF ( nVars /= nVarsExpected ) THEN
      WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nVars, &
                                                'but expected:',nVarsExpected
      CALL ErrorStop(global,ERR_INVALID_NVARS,1539)
    END IF ! nVarsExpected

! ==============================================================================
!   Rest of file
! ==============================================================================

    IF ( iDataSet == CV_DATASET_OLD ) THEN
      CALL RFLU_ReadFlowV1Binary(pRegion,iFile,nVars)
    ELSE
      CALL RFLU_ReadFlowV2Binary(pRegion,iFile,nVars)
    END IF ! iDataSet

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,1559,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary flow file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadFlowBinary


! ******************************************************************************
!
! Purpose: Read flow file for mixture in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Read initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_PICL_ReadFlowBinary(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    LOGICAL :: fileExists
    CHARACTER(CHRLEN) :: errorString,iFileName,iFileNameOld,sectionString, &
                         timeString1,timeString2
    INTEGER :: errorFlag,i,iDataSet,iFile,j,nCellsTot,nCellsExpected,nVars, &
               nVarsExpected,precActual,precExpected,rangeActual,rangeExpected
    REAL(RFREAL) :: currentTime
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PICL_ReadFlowBinary',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading picl binary flow file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.picl.vf', &
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.picl.vf', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT

    INQUIRE(FILE=iFileName,EXIST=fileExists)

    IF ( fileExists .EQV. .TRUE. ) THEN
!      OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
!           IOSTAT=errorFlag)
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
      OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
           IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_L )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
         ACCESS="SEQUENTIAL",CONVERT="LITTLE_ENDIAN",IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_B )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
         ACCESS="SEQUENTIAL",CONVERT="BIG_ENDIAN",IOSTAT=errorFlag)
    END IF
! BBR - end
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_FILE_OPEN,1682,iFileName)
      END IF ! global%error
    ELSE
!      OPEN(iFile,FILE=iFileNameOld,FORM="UNFORMATTED",STATUS="OLD", &
!           IOSTAT=errorFlag)
! BBR - begin
    !END IF
! BBR - end
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_FILE_OPEN,1692,iFileNameOld)
      END IF ! global%error
    END IF ! fileExists

! ==============================================================================
! Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile) sectionString
    IF ( TRIM(sectionString) == '# ROCFLU flow file' ) THEN
      iDataSet = CV_DATASET_OLD
    ELSE IF ( TRIM(sectionString) == '# ROCFLU flow file (v2)' ) THEN
      iDataSet = CV_DATASET_NEW
    ELSE
      CALL ErrorStop(global,ERR_INVALID_MARKER,1711,sectionString)
    END IF ! TRIM

! -----------------------------------------------------------------------------
!   Precision and range
! -----------------------------------------------------------------------------

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Precision and range' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,1720,sectionString)
    END IF ! TRIM

    precExpected  = PRECISION(1.0_RFREAL)
    rangeExpected = RANGE(1.0_RFREAL)

    READ(iFile) precActual,rangeActual
    IF ( precActual < precExpected .OR. rangeActual < rangeExpected ) THEN
      CALL ErrorStop(global,ERR_PREC_RANGE,1728)
    END IF ! precActual

! -----------------------------------------------------------------------------
!   Initial residual and physical time
! -----------------------------------------------------------------------------

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Initial residual' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,1737,iFileName)
    END IF ! TRIM

    READ(iFile) global%resInit

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Physical time' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,1744,iFileName)
    END IF ! TRIM

    READ(iFile) currentTime

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      IF ( global%currentTime < 0.0_RFREAL ) THEN
        global%currentTime = currentTime
      ELSE
        WRITE(timeString1,'(1PE11.5)') global%currentTime
        WRITE(timeString2,'(1PE11.5)') currentTime

        IF ( TRIM(timeString1) /= TRIM(timeString2) ) THEN
          WRITE(STDOUT,'(A,4(1X,A))') SOLVER_NAME, &
                                      '*** WARNING *** Time mismatch:', &
                                      TRIM(timeString1),'vs.',TRIM(timeString2)

          global%warnCounter = global%warnCounter + 1
        END IF ! global%currentTime
      END IF ! global%currentTime
    END IF ! global%flowType

! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    nVarsExpected  = pRegion%mixtInput%nCv
    nCellsExpected = pGrid%nCellsTot

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,1777,sectionString)
    END IF ! TRIM

    READ(iFile) nCellsTot,nVars

    IF ( nCellsTot /= nCellsExpected ) THEN
      WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nCellsTot, &
                                                'but expected:',nCellsExpected
      CALL ErrorStop(global,ERR_INVALID_NCELLS,1785,errorString)
    END IF ! nCellsExpected

    IF ( nVars /= nVarsExpected ) THEN
      WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nVars, &
                                                'but expected:',nVarsExpected
      CALL ErrorStop(global,ERR_INVALID_NVARS,1791)
    END IF ! nVarsExpected

! ==============================================================================
!   Rest of file
! ==============================================================================

    IF ( iDataSet == CV_DATASET_OLD ) THEN
      !CALL RFLU_ReadFlowV1Binary(pRegion,iFile,nVars)
    ELSE
      CALL RFLU_PICL_ReadFlowV2Binary(pRegion,iFile,nVars)
    END IF ! iDataSet

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,1811,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary flow file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PICL_ReadFlowBinary

! ******************************************************************************
!
! Purpose: Read flow dataset (v2) for mixture in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iFile       File number
!   nVars       Number of variables to be read
!
! Output: None.
!
! Notes:
!
! ******************************************************************************

  SUBROUTINE RFLU_PICL_ReadFlowV2Binary(pRegion,iFile,nVars)
   
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iFile,nVars
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString,sectionString
    INTEGER :: errorFlag,i,iLoc,iVars,j,loopCounter
    !REAL(RFREAL), DIMENSION(:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, allocate pointers
! ******************************************************************************
    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PICL_ReadFlowV2Binary',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading PICL binary flow dataset ' // &
                                           '(v2)...'
    END IF ! global%verbLevel
    nVars = 1
    pGrid => pRegion%grid

! ******************************************************************************
!   File is already open, just read the dataset
! ******************************************************************************

    iVars       = 0
    loopCounter = 0

!        ALLOCATE(pRegion%mixt%piclVF(pGrid%nCellsTot),STAT=errorFlag)
!        global%error = errorFlag
!        IF (global%error /= ERR_NONE) THEN
!        CALL ErrorStop(global,ERR_ALLOCATE,1894,'pRegion%mixt%piclfVF')
!        END IF ! global%error


    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile) sectionString

! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------

      IF ( sectionString(1:2) /= '# ' ) THEN
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel

        CALL ErrorStop(global,ERR_INVALID_MARKER,1912,sectionString)
      END IF ! sectionString(1:2)

      sectionString = sectionString(3:)

      SELECT CASE ( TRIM(sectionString) )

! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------

        CASE ( 'End' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel

          EXIT

! ------------------------------------------------------------------------------
!       Extract dataset location in cv from section string
! ------------------------------------------------------------------------------

        CASE DEFAULT
          READ(sectionString,'(I2)') iLoc

          IF ( iLoc < 1 .AND. iLoc > pRegion%mixtInput%nCv ) THEN
            IF ( global%verbLevel > VERBOSE_LOW ) THEN
              sectionString = '# ' // sectionString
              WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
            END IF ! verbosityLevel

            CALL ErrorStop(global,ERR_INVALID_CVILOC,1944)
          END IF ! iLoc

          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A,I2,A)') SOLVER_NAME,'Picl vf location ',iLoc,'...'
          END IF ! global%verbLevel

          !pCv => pRegion%mixt%piclVF!pRegion%mixt%cv

           iVars = iVars + 1
        !   write(*,*) "Pcv",pCv(1)
           write(*,*) "DEBUG BEFORE PICL READ"
          ! do j=1,pGrid%nCellsTot
           !     READ(iFile) pCv(j)
           !     write(*,*) "j vf", j,pCv(j)
           !end do
  IF ( global%piclUsed) THEN
          READ(iFile) (pRegion%mixt%piclVF(j),j=1,pGrid%nCellsTot)
  END IF
!(pCv(j),j=1,pGrid%nCellsTot)
           do j=1,pGrid%nCellsTot
                !READ(iFile) pCv(j)
                !write(*,*) "j vf", j,pRegion%mixt%piclVF(j)
           end do
           write(*,*) "DEBUG AFTER PICL READ"
      END SELECT ! TRIM

! ------------------------------------------------------------------------------
!     Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN
        CALL ErrorStop(global,ERR_INFINITE_LOOP,1979)
      END IF ! loopCounter

    END DO ! <empty>

! ==============================================================================
!   Check and information about number of variables read
! ==============================================================================

    IF ( iVars /= nVars ) THEN
      CALL ErrorStop(global,ERR_INVALID_NVARS,1989)
    END IF ! iVar

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary flow dataset ' // &
                                           '(v2) done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PICL_ReadFlowV2Binary


! ******************************************************************************
!
! Purpose: Wrapper for reading of flow files in ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadFlowWrapper(pRegion)



    USE SPEC_RFLU_ModReadWriteVars, ONLY: SPEC_RFLU_ReadCvASCII, &
                                          SPEC_RFLU_ReadCvBinary, &     
                                          SPEC_RFLU_ReadEEvASCII, & 
                                          SPEC_RFLU_ReadEEvBinary

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadFlowWrapper',"../modflu/RFLU_ModReadWriteFlow.F90")

! ******************************************************************************
!   Read mixture solution files
! ******************************************************************************

      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL RFLU_ReadFlowASCII(pRegion)
!      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
! BBR - begin
       ELSE IF ( global%solutFormat == FORMAT_BINARY .OR. & 
                global%solutFormat == FORMAT_BINARY_L .OR. &
                global%solutFormat == FORMAT_BINARY_B ) THEN
! BBR - end
        CALL RFLU_ReadFlowBinary(pRegion)
  IF ( global%piclUsed) THEN
!        if ( global%moduleType == MODULE_TYPE_POSTPROC ) THEN
                CALL RFLU_PICL_ReadFlowBinary(pRegion)
!        end if
  END IF
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,2093)
      END IF ! global%solutFormat

! ******************************************************************************
!   Read physical module solution files
! ******************************************************************************


! ==============================================================================
!   Species
! ==============================================================================

    IF ( global%specUsed .EQV. .TRUE. ) THEN
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL SPEC_RFLU_ReadCvASCII(pRegion)
                
        IF ( global%moduleType == MODULE_TYPE_POSTPROC ) THEN
          IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN  
            CALL SPEC_RFLU_ReadEEvASCII(pRegion)
          END IF ! pRegion%specInput%nSpeciesEE
        END IF ! global%moduleType 
!      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
! BBR - begin
       ELSE IF ( global%solutFormat == FORMAT_BINARY .OR. & 
                global%solutFormat == FORMAT_BINARY_L .OR. &
                global%solutFormat == FORMAT_BINARY_B ) THEN
! BBR - end
        CALL SPEC_RFLU_ReadCvBinary(pRegion)
        
        IF ( global%moduleType == MODULE_TYPE_POSTPROC ) THEN
          IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN  
            CALL SPEC_RFLU_ReadEEvBinary(pRegion)
          END IF ! pRegion%specInput%nSpeciesEE
        END IF ! global%moduleType        
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,2156)
      END IF ! global%solutFormat
    END IF ! global%specUsed

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadFlowWrapper

! ******************************************************************************
!
! Purpose: Write flow file for mixture in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Write initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteFlowASCII(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString,iFileName1
    INTEGER :: errorFlag,i,iFile,iLoc,j,iFile1
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv, pDv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteFlowASCII',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII flow file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.mixt.cva', &
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.mixt.dva', &
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName1)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.mixt.cva', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT      
    iFile1 = IF_SOLUT+1
    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", &
         IOSTAT=errorFlag)
    OPEN(iFile1,FILE=iFileName1,FORM="FORMATTED",STATUS="UNKNOWN", &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,2261,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU flow file (v2)'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile1,'(A)') sectionString

    sectionString = '# Precision and range'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(2(I8))') PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)
    WRITE(iFile1,'(A)') sectionString
    WRITE(iFile1,'(2(I8))') PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

    sectionString = '# Initial residual'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(E23.16)') global%resInit
    WRITE(iFile1,'(A)') sectionString
    WRITE(iFile1,'(E23.16)') global%resInit

    sectionString = '# Physical time'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(E23.16)') global%currentTime
    WRITE(iFile1,'(A)') sectionString
    WRITE(iFile1,'(E23.16)') global%currentTime
    
! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    sectionString = '# Dimensions'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(2(I16))') pGrid%nCellsTot,pRegion%mixtInput%nCv
    WRITE(iFile1,'(A)') sectionString
    WRITE(iFile1,'(2(I16))') pGrid%nCellsTot,4 !pRegion%mixtInput%nCv

! ==============================================================================
!   Data
! ==============================================================================

    pCv => pRegion%mixt%cv
    pDv => pRegion%mixt%dv

    DO iLoc=1,pRegion%mixtInput%nCv
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_LOW ) THEN
        WRITE(STDOUT,'(A,3X,A,I2,A)') SOLVER_NAME,'Cv location ',iLoc,'...'
      END IF ! global%verbLevel

      WRITE(sectionString,'(A,I2)') '# ',iLoc
      WRITE(iFile,'(A)') sectionString
      WRITE(iFile,'(5(E23.16))') (pCv(iLoc,j),j=1,pGrid%nCellsTot)
    END DO ! iLoc

    DO iLoc = 1,4
      WRITE(sectionString,'(A,I2)') '# ',iLoc
      WRITE(iFile1,'(A)') sectionString
      WRITE(iFile1,'(5(E23.16))') (pDv(iLoc,j),j=1,pGrid%nCellsTot)
    END DO ! iLoc


! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile1,'(A)') sectionString

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    CLOSE(iFile1,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,2353,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII flow file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteFlowASCII






! ******************************************************************************
!
! Purpose: Write flow file for mixture in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Write initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteFlowBinary(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,i,iFile,iLoc,j
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteFlowBinary',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary flow file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.mixt.cv', &
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)    

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.mixt.cv', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT
!    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
!         IOSTAT=errorFlag)
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_L )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         ACCESS="SEQUENTIAL",CONVERT="LITTLE_ENDIAN",IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_B )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         ACCESS="SEQUENTIAL",CONVERT="BIG_ENDIAN",IOSTAT=errorFlag)
    END IF
! BBR - end
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,2472,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU flow file (v2)'
    WRITE(iFile) sectionString

    sectionString = '# Precision and range'
    WRITE(iFile) sectionString
    WRITE(iFile) PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

    sectionString = '# Initial residual'
    WRITE(iFile) sectionString
    WRITE(iFile) global%resInit

    sectionString = '# Physical time'
    WRITE(iFile) sectionString
    WRITE(iFile) global%currentTime
    
! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    sectionString = '# Dimensions'
    WRITE(iFile) sectionString
    WRITE(iFile) pGrid%nCellsTot,pRegion%mixtInput%nCv

! ==============================================================================
!   Data
! ==============================================================================

    pCv => pRegion%mixt%cv

    DO iLoc=1,pRegion%mixtInput%nCv
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_LOW ) THEN
        WRITE(STDOUT,'(A,3X,A,I2,A)') SOLVER_NAME,'Cv location ',iLoc,'...'
      END IF ! global%verbLevel

      WRITE(sectionString,'(A,I2)') '# ',iLoc
      WRITE(iFile) sectionString
      WRITE(iFile) (pCv(iLoc,j),j=1,pGrid%nCellsTot)
    END DO ! iLoc

! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile) sectionString

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,2545,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary flow file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteFlowBinary


! ******************************************************************************
!
! Purpose: Write particle cell properties file in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Write initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_PICL_WriteFlowBinary(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,i,iFile,iLoc,j
    !REAL(RFREAL), DIMENSION(:), POINTER :: pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PICL_WriteFlowBinary',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing PICL binary flow file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.picl.vf', &
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.picl.vf', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT
!    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
!         IOSTAT=errorFlag)
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_L )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         ACCESS="SEQUENTIAL",CONVERT="LITTLE_ENDIAN",IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_B )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         ACCESS="SEQUENTIAL",CONVERT="BIG_ENDIAN",IOSTAT=errorFlag)
    END IF
! BBR - end
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,2660,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU flow file (v2)'
    WRITE(iFile) sectionString

    sectionString = '# Precision and range'
    WRITE(iFile) sectionString
    WRITE(iFile) PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

    sectionString = '# Initial residual'
    WRITE(iFile) sectionString
    WRITE(iFile) global%resInit

    sectionString = '# Physical time'
    WRITE(iFile) sectionString
    WRITE(iFile) global%currentTime

! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    sectionString = '# Dimensions'
    WRITE(iFile) sectionString
    WRITE(iFile) pGrid%nCellsTot,pRegion%mixtInput%nCv

! ==============================================================================
!   Data
! ==============================================================================

    !pCv => pRegion%mixt%piclVF

!    DO iLoc=1,pRegion%mixtInput%nCv
     iLoc=1
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_LOW ) THEN
        WRITE(STDOUT,'(A,3X,A,I2,A)') SOLVER_NAME,'VF location ',iLoc,'...'
      END IF ! global%verbLevel

      WRITE(sectionString,'(A,I2)') '# ',iLoc
      WRITE(iFile) sectionString
  IF ( global%piclUsed) THEN
      WRITE(iFile) (pRegion%mixt%piclVF(j),j=1,pGrid%nCellsTot)
      ! 03/20/2025 - Thierry - begins here
      WRITE(iFile) (pRegion%mixt%piclgradRhog(:,1,j),j=1,pGrid%nCellsTot)
      ! 03/20/2025 - Thierry - ends here

  END IF
!(pCv(j),j=1,pGrid%nCellsTot)
      !if (pRegion%iRegionGlobal .eq. 1) then
      !do j=1,pGrid%nCellsTot
       ! write(*,*)"pCv",j, pCv(j),pRegion%mixt%piclVF(j)
      !end do
      !end if
!    END DO ! iLoc

! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile) sectionString

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,2748,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary piclVF file done.'
      ! 03/20/2025 - Thierry - begins here
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary piclgradRhog file&
                                            done.'
      ! 03/20/2025 - Thierry - ends here
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PICL_WriteFlowBinary




! ******************************************************************************
!
! Purpose: Write flow file for mixture in binary PARAVIEW format.
!
! Description: Theses files are meant for parallel read from PARAVIEW
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Write initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa
!
!   2. The species output is currently written assuming a 2 species fluid.
!      Routine will have to be modified for 3+ species.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteVTK(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************


! number of timesteps kept in history kernels
! maximum number of triangular patch boundaries

! y, y1, ydot, ydotc: 12

! rprop: 53

! rprop4: PPICLF_LRP4 - Reynolds Subgrid Stress Components

! rprop5: PPICLF_LRP5 - Storing Force Models
!#define PPICLF_R_FQSX 1
!#define PPICLF_R_FQSY 2
!#define PPICLF_R_FQSZ 3
!#define PPICLF_R_FAMX 4
!#define PPICLF_R_FAMY 5
!#define PPICLF_R_FAMZ 6
!#define PPICLF_R_FAMBX 7
!#define PPICLF_R_FAMBY 8
!#define PPICLF_R_FAMBZ 9
!#define PPICLF_R_FCX 10
!#define PPICLF_R_FCY 11
!#define PPICLF_R_FCZ 12
!#define 49 13
!#define 50 14
!#define 51 15
!#define 52 16
!#define PPICLF_R_FPGX 17 
!#define PPICLF_R_FPGY 18 
!#define PPICLF_R_FPGZ 19 

! map: 22
!--- x,y,z Forces Fedback to Rocflu
!---
!--- Add comment about these terms 
!--- Reynolds Subgrid Stress Tensor
!--- Pseudo Turbulent Kinetic Energy























! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,i,iFile,iLoc,j,regId,icl,icg,ivl,ivg,p
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pPv,pCv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    !TYPE(t_region), POINTER :: pRegion

    CHARACTER(CHRLEN) :: cnum,dnum
    
    INTEGER(I4P)::D_IO,E_IO,Nn,Ne,offset_type,celltype,indx
    
    INTEGER(I1P), ALLOCATABLE, DIMENSION(:) :: cell_type
    INTEGER(I4P), ALLOCATABLE, DIMENSION(:) :: offset,connect
    REAL(R4P), ALLOCATABLE, DIMENSION(:,:) :: points
    REAL(R8P), ALLOCATABLE, DIMENSION(:) :: temp

    INTEGER(I4p), DIMENSION(:,:), POINTER :: cell2v


  LOGICAL :: scalarConvFlag
  INTEGER :: iCvSpecAir,iCvSpecProducts,iSpec
  INTEGER :: iCvSpecArgon
  REAL(RFREAL) :: YProducts
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
  TYPE(t_spec_input), POINTER :: pSpecInput
  TYPE(t_spec_type), POINTER :: pSpecType

!Fred - For Species plotting capability
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: phiP
   REAL(KIND=8),DIMENSION(:,:,:,:), ALLOCATABLE :: vfP,vfp2
   REAL(RFREAL) :: tester, total_vol,tester2
   integer :: lx,ly,lz

   INTEGER :: gasModel
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start, open vtu file
! ******************************************************************************

    global => pRegion%global 
    gasModel = pRegion%mixtInput%gasModel



    IF ( global%specUsed .EQV. .TRUE. ) THEN
      pSpecInput => pRegion%specInput
      pCvSpec => pRegion%spec%cv

      ! TLJ - 12/14/2024
      IF (gasModel == GAS_MODEL_MIXT_JWL) THEN
         iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                   'PRODUCTS')
      ELSE
         iCvSpecArgon = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                   'ARGON')
      ENDIF
      iCvSpecAir = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                   'AIR')
    END IF
!Fred - Adding species plotting capabilities

    CALL RegisterFunction(global,'RFLU_WriteVTK',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing PARAVIEW file...'
      IF (global%nSpecies > 0) THEN
         WRITE(STDOUT,'(A,2x,2I4)') '  TLJ: iCvSpecAir, iCvSpecProducts = ', &
              iCvSpecAir, iCvSpecProducts
      ELSE
         WRITE(STDOUT,'(A,2x,2I4)') '  TLJ: Single-species Gas'
      END IF
    END IF ! global%verbLevel
    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteadyVTK(global,FILEDEST_OUTDIR,'.vtu', &
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteadyVTK(global,FILEDEST_OUTDIR,'.vtu', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

! ==============================================================================
!   Write pvtu file content for tet, hex and prism elements
! ==============================================================================
   Nn   = pRegion%grid%nVert
   Ne   = pRegion%grid%nCells

   ALLOCATE(points(XCOORD:ZCOORD,1:Nn))

   IF ( pRegion%grid%nTets > 0 ) THEN
     indx = 4
     offset_type = 4_I4P
     celltype = 10_I1P
     cell2v => pRegion%grid%tet2v
   ELSE IF ( pRegion%grid%nHexs > 0 ) THEN
     indx = 8
     offset_type = 8_I4P
     celltype = 12_I1P
     cell2v => pRegion%grid%hex2v
   ELSE IF ( pRegion%grid%nPris > 0 ) THEN
     indx = 6
     offset_type = 6_I4P
     celltype = 13_I1P
     cell2v => pRegion%grid%pri2v

   END IF

   DO ivl = 1,Nn
     points(XCOORD,ivl) = pRegion%grid%xyz(XCOORD,ivl)
     points(YCOORD,ivl) = pRegion%grid%xyz(YCOORD,ivl)
     points(ZCOORD,ivl) = pRegion%grid%xyz(ZCOORD,ivl)
   END DO
 
   ALLOCATE(offset(1:Ne))
   DO j = 1, Ne
     offset(j) = j*offset_type
   END DO
 
   ALLOCATE(cell_type(1:Ne))
   DO j = 1, Ne
     cell_type(j) = celltype
   END DO

   ALLOCATE(connect(1:indx*Ne))
   p=0
   DO j=1,Ne
      DO i=1,indx
       p=p+1
       connect(p)=cell2v(i,j)-1_I4P
      END DO
   END DO

!! Rahul- The code block commented below works only for hexahedral elements.
!! It is now updated to accommodate tetrahedral elements too.
!  DO ivl = 1,Nn
!    points(XCOORD,ivl) = pRegion%grid%xyz(XCOORD,ivl) 
!    points(YCOORD,ivl) = pRegion%grid%xyz(YCOORD,ivl)
!    points(ZCOORD,ivl) = pRegion%grid%xyz(ZCOORD,ivl)
!  END DO

!  ALLOCATE(offset(1:Ne))
!  DO j = 1, Ne
!    offset(j) = j*8_I4P
!  END DO

!  ALLOCATE(cell_type(1:Ne))
!  DO j = 1, Ne
!    cell_type(j) = 12_I1P
!  END DO

!  ALLOCATE(connect(1:8*Ne))
!  p=0
!  DO j=1,Ne
!     DO i=1,8
!     p=p+1
!     connect(p)=pRegion%grid%hex2v(i,j)-1_I4P
!     END DO
!  END DO 
!! End rahul

  ALLOCATE(temp(1:Ne))
  E_IO = VTK_INI_XML(output_format = 'binary', filename = iFilename, &
         mesh_topology = 'UnstructuredGrid')
  E_IO = VTK_GEO_XML(NN = Nn, NC = Ne, &
         X = points(XCOORD,1:Nn), Y = points(YCOORD,1:Nn), Z = points(ZCOORD,1:Nn))
  E_IO = VTK_CON_XML(NC = Ne, connect = connect, offset = offset &
         , cell_type = cell_type )
  E_IO = VTK_DAT_XML(var_location = 'cell', var_block_action = 'open')

  ! Gas-phase variables
  ! TLJ - I changed these to be primitive variables
  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Gas Density', &
         var = pRegion%mixt%cv(CV_MIXT_DENS,1:Ne))
  !E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Momentum', &
  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Velocity', &
         varX=pRegion%mixt%cv(CV_MIXT_XMOM,1:Ne), &
         varY=pRegion%mixt%cv(CV_MIXT_YMOM,1:Ne), &
         varZ=pRegion%mixt%cv(CV_MIXT_ZMOM,1:Ne))

  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Energy', &
         var=pRegion%mixt%cv(CV_MIXT_ENER,1:Ne))
  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Pressure', &
         var=pRegion%mixt%dv(DV_MIXT_PRES,1:Ne))
  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Temperature', &
         var=pRegion%mixt%dv(DV_MIXT_TEMP,1:Ne))
  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Speed of Sound', &
         var=pRegion%mixt%dv(DV_MIXT_SOUN,1:Ne))

  ! TLJ - 12/14/2024
  IF (gasModel == GAS_MODEL_MIXT_JWL) THEN
    E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Products Energy', &
           var=pRegion%mixt%dv(DV_MIXT_EJWL,1:Ne))
  ENDIF

  IF ( global%piclUsed) THEN
    !E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'phi_p', &
    !     var=pRegion%mixt%piclVF(1:Ne))
    E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Particle Volume Fraction', &
         var=pRegion%mixt%piclVF(1:Ne))
      ! 03/20/2025 - Thierry - begins here
    E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Gradient Gas Density', &
         varX=pRegion%mixt%piclgradRhog(1,1,1:Ne), &
         varY=pRegion%mixt%piclgradRhog(2,1,1:Ne), &
         varZ=pRegion%mixt%piclgradRhog(3,1,1:Ne))
      ! 03/20/2025 - Thierry - ends here
  END IF


!Fred - Adding species plotting options for PARAVIEW files - 11/1/2017
!TLJ  - Added conversion from conservative \rho*Y to Y - 12/14/2024
  IF ( global%specUsed .EQV. .TRUE. ) THEN
     IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
        scalarConvFlag = .FALSE.
     ELSE
        scalarConvFlag = .TRUE.
        CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv, &
                     pRegion%spec%cvState)
     END IF

     ! TLJ - 12/14/2024
     IF (gasModel == GAS_MODEL_MIXT_JWL) THEN
        E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Explosive Mass Fraction', &
               var = pCvSpec(iCvSpecProducts,1:Ne))
     ELSE
        E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Argon Mass Fraction', &
               var = pCvSpec(iCvSpecArgon,1:Ne))
     ENDIF
     E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Air Mass Fraction', &
            var = pCvSpec(iCvSpecAir,1:Ne))

     IF (scalarConvFlag .EQV. .TRUE.) THEN
        CALL RFLU_ScalarConvertcvPrim2Cons(pRegion,pRegion%spec%cv, &
                      pRegion%spec%cvState)
     END IF
  END IF
  
  E_IO = VTK_DAT_XML(var_location = 'cell', var_block_action = 'close')
  E_IO = VTK_GEO_XML()
  E_IO = VTK_END_XML()

  DEALLOCATE(points)
  DEALLOCATE(offset)
  DEALLOCATE(cell_type)
  DEALLOCATE(connect)
! ==============================================================================
!   Close vtu file
! ==============================================================================


! ******************************************************************************
!   Start, open pvtu file
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC) THEN

!      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Opening pvtu files...'

      IF ( global%flowType == FLOW_UNSTEADY ) THEN
        CALL BuildFileNameUnsteadyPVTK(global,FILEDEST_OUTDIR,'.pvtu', &
                           global%currentTime,iFileName)
      ELSE
        CALL BuildFileNameSteadyPVTK(global,FILEDEST_OUTDIR,'.pvtu', &
                           global%currentIter,iFileName)    
      END IF

      E_IO = PVTK_INI_XML(filename = iFilename,  &
             mesh_topology = 'PUnstructuredGrid', tp='Float32')      

      DO regId = 1, global%nRegions 
        IF ( global%flowType == FLOW_UNSTEADY ) THEN
!          CALL BuildFileNameUnsteadyVTK(global,FILEDEST_OUTDIR,'.vtu', &
!                                   regId,global%currentTime, &
!                                   iFileName)
        WRITE(iFileName,'(A,1PE11.5,A,I5.5,A)') &
                          TRIM('./')//TRIM(global%caseName)// &
                          '_',global%currentTime,TRIM('_') &
                          ,regId,TRIM('.vtu')       
        ELSE
!          CALL BuildFileNameSteadyVTK(global,FILEDEST_OUTDIR,'.vtu', &
!                                 regId,global%currentIter, &
!                                 iFileName)
        WRITE(iFileName,'(A,1PE11.5,A,I6.6,A)') &
                          TRIM('./')//TRIM(global%caseName)// &
                          '_',global%currentIter,TRIM('_') &
                          ,regId,TRIM('.vtu')
        ENDIF ! global%flowType
        E_IO = PVTK_GEO_XML(source=iFilename)
      END DO

      E_IO = PVTK_DAT_XML(var_location = 'cell', var_block_action = 'OPEN')
      E_IO = PVTK_VAR_XML(varname = 'Gas Density', tp='Float64')
      !E_IO = PVTK_VAR_XML(Nc = 3, varname = 'Momentum', tp='Float64' )
      E_IO = PVTK_VAR_XML(Nc = 3, varname = 'Velocity', tp='Float64' )
      E_IO = PVTK_VAR_XML(varname = 'Energy', tp='Float64')
      E_IO = PVTK_VAR_XML(varname = 'Pressure', tp='Float64')
      E_IO = PVTK_VAR_XML(varname = 'Temperature', tp='Float64')
      E_IO = PVTK_VAR_XML(varname = 'Speed of Sound', tp='Float64')
      IF (gasModel == GAS_MODEL_MIXT_JWL) THEN
         E_IO = PVTK_VAR_XML(varname = 'Products Energy', tp='Float64')
      ENDIF

  IF ( global%piclUsed) THEN
      E_IO = PVTK_VAR_XML(varname = 'Particle Volume Fraction', tp='Float64')
      ! 03/20/2025 - Thierry - begins here
      E_IO = PVTK_VAR_XML(Nc = 3, varname = 'Gradient Gas Density',& 
                                                               tp='Float64' )
      ! 03/20/2025 - Thierry - ends here
  END IF


IF ( global%specUsed .EQV. .TRUE. ) THEN
     ! TLJ - 12/14/2024
     IF (gasModel == GAS_MODEL_MIXT_JWL) THEN
        E_IO = PVTK_VAR_XML(varname = 'Explosive Mass Fraction', tp='Float64')    
     ELSE
        E_IO = PVTK_VAR_XML(varname = 'Argon Mass Fraction', tp='Float64')    
     ENDIF
     E_IO = PVTK_VAR_XML(varname = 'Air Mass Fraction', tp='Float64')    
END IF

      E_IO = PVTK_DAT_XML(var_location = 'cell', var_block_action = 'Close') 
      E_IO = PVTK_END_XML()

    END IF ! global%myProcid

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing PARAVIEW file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteVTK


! ******************************************************************************
!
! Purpose: Write midsection slice of solution in binary PARAVIEW format.
!
! Description: Theses files are meant for parallel read from PARAVIEW
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Write initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteSlice(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,i,iFile,iLoc,j,regId,sLoc,ivl,ivg,p
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pPv
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

    CHARACTER(CHRLEN) :: cnum,dnum

    INTEGER(I4P)::D_IO,E_IO,Nn,Ne,Nz,islcn,islce

    INTEGER(I1P), ALLOCATABLE, DIMENSION(:) :: cell_type
    INTEGER(I4P), ALLOCATABLE, DIMENSION(:) :: offset,connect
    REAL(R4P), ALLOCATABLE, DIMENSION(:,:) :: points
    REAL(R8P), ALLOCATABLE, DIMENSION(:) :: temp



! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start, open slice vtu file
! ******************************************************************************

    global => pRegion%global


    CALL RegisterFunction(global,'RFLU_WriteSlice',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A)') '************************************************'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing Slice PARAVIEW file...'
      WRITE(STDOUT,'(A)') '************************************************'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteadyVTK(global,FILEDEST_OUTDIR,'.slice.vtu', &
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)
    ELSE
      CALL BuildFileNameSteadyVTK(global,FILEDEST_OUTDIR,'.slice.vtu', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)
    ENDIF ! global%flowType

! ==============================================================================
!   Write slice vtu file content
! ==============================================================================


    Nz = pRegion%grid%nCells/pRegion%patches(2)%nBQuads

    Nn = pRegion%grid%nVert / Nz
    Ne = pRegion%grid%nCells / Nz 

    islcn = INT(Nz/2)*Nn 
    islce = INT(Nz/2)*Ne

  ALLOCATE(points(XCOORD:ZCOORD,1:Nn))
   
  p = 0
!  DO ivl = 1,Nn
  DO ivl = islcn+1,islcn+Nn
    p = p + 1
   ! points(XCOORD,ivl) = pRegion%grid%xyz(XCOORD,ivl)
   ! points(YCOORD,ivl) = pRegion%grid%xyz(YCOORD,ivl)
   ! points(ZCOORD,ivl) = pRegion%grid%xyz(ZCOORD,ivl)
    points(XCOORD,p) = pRegion%grid%xyz(XCOORD,ivl)
    points(YCOORD,p) = pRegion%grid%xyz(YCOORD,ivl)
    points(ZCOORD,p) = pRegion%grid%xyz(ZCOORD,ivl)
  END DO

  ALLOCATE(offset(1:Ne))
  DO j = 1, Ne
    offset(j) = j*8_I4P
  END DO

  ALLOCATE(cell_type(1:Ne))
  DO j = 1, Ne
    cell_type(j) = 12_I1P
  END DO

  ALLOCATE(connect(1:8*Ne))
  p=0
  DO j=1,Ne
     DO i=1,8
     p=p+1
     connect(p)=pRegion%grid%hex2v(i,j)-1
     END DO
  END DO

  E_IO = VTK_INI_XML(output_format = 'binary', filename = iFilename, &
         mesh_topology = 'UnstructuredGrid')
!  E_IO = VTK_GEO_XML(NN = Nn, NC = Ne, &
!         X = points(XCOORD,1:Nn), Y = points(YCOORD,1:Nn), Z = points(ZCOORD,1:Nn))
!  E_IO = VTK_CON_XML(NC = Ne, connect = connect, offset = offset &
!         , cell_type = cell_type )
!  E_IO = VTK_DAT_XML(var_location = 'cell', var_block_action = 'open')
!  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Gas Density', &
!         var = pRegion%mixt%cv(CV_MIXT_DENS,1:Ne))
!  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Momentum', &
!         varX=pRegion%mixt%cv(CV_MIXT_XMOM,1:Ne), &
!         varY=pRegion%mixt%cv(CV_MIXT_YMOM,1:Ne), &
!         varZ=pRegion%mixt%cv(CV_MIXT_ZMOM,1:Ne))
!  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Energy', &
!         var = pRegion%mixt%cv(CV_MIXT_ENER,1:Ne))
!  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Pressure', &
!         var = pRegion%mixt%cv(CV_MIXT_PRES,1:Ne))
!  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Temperature', &
!         var = pRegion%mixt%dv(DV_MIXT_TEMP,1:Ne))
!  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Speed of Sound', &
!         var = pRegion%mixt%dv(DV_MIXT_SOUN,1:Ne))

  E_IO = VTK_GEO_XML(NN = Nn, NC = Ne, &
         X = points(XCOORD,islcn+1:islcn+Nn), &
         Y = points(YCOORD,islcn+1:islcn+Nn), &
         Z = points(ZCOORD,islcn+1:islcn+Nn))
  E_IO = VTK_CON_XML(NC = Ne, connect = connect, offset = offset &
         , cell_type = cell_type )
  E_IO = VTK_DAT_XML(var_location = 'cell', var_block_action = 'open')
  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Gas Density', &
         var = pRegion%mixt%cv(CV_MIXT_DENS,islce+1:islce+Ne))
  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Momentum', &
         varX=pRegion%mixt%cv(CV_MIXT_XMOM,islce+1:islce+Ne), &
         varY=pRegion%mixt%cv(CV_MIXT_YMOM,islce+1:islce+Ne), &
         varZ=pRegion%mixt%cv(CV_MIXT_ZMOM,islce+1:islce+Ne))
  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Energy', &
         var = pRegion%mixt%cv(CV_MIXT_ENER,islce+1:islce+Ne))
  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Pressure', &
         var = pRegion%mixt%dv(DV_MIXT_PRES,islce+1:islce+Ne))
  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Temperature', &
         var = pRegion%mixt%dv(DV_MIXT_TEMP,islce+1:islce+Ne))
  E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Speed of Sound', &
         var = pRegion%mixt%dv(DV_MIXT_SOUN,islce+1:islce+Ne))
  !E_IO = VTK_VAR_XML(NC_NN = Ne, varname = 'Products Energy', &
  !       var = pRegion%mixt%dv(DV_MIXT_EJWL,islce+1:islce+Ne))


  E_IO = VTK_DAT_XML(var_location = 'cell', var_block_action = 'close')
  E_IO = VTK_GEO_XML()
  E_IO = VTK_END_XML()

  DEALLOCATE(points)
  DEALLOCATE(offset)
  DEALLOCATE(cell_type)
  DEALLOCATE(connect)

! ==============================================================================
!   Close slice vtu file
! ==============================================================================


! ******************************************************************************
!   Start, open slice pvtu file
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC) THEN

      IF ( global%flowType == FLOW_UNSTEADY ) THEN
        CALL BuildFileNameUnsteadyPVTK(global,FILEDEST_OUTDIR,'.slice.pvtu', &
                           global%currentTime,iFileName)
      ELSE
        CALL BuildFileNameSteadyPVTK(global,FILEDEST_OUTDIR,'.slice.pvtu', &
                           global%currentIter,iFileName)
      END IF

      E_IO = PVTK_INI_XML(filename = iFilename,  &
             mesh_topology = 'PUnstructuredGrid', tp='Float32')

      DO regId = 1, global%nRegions
        IF ( global%flowType == FLOW_UNSTEADY ) THEN
          CALL BuildFileNameUnsteadyVTK(global,FILEDEST_OUTDIR,'.slice.vtu', &
                                   regId,global%currentTime, &
                                   iFileName)
        ELSE
          CALL BuildFileNameSteadyVTK(global,FILEDEST_OUTDIR,'.slice.vtu', &
                                 regId,global%currentIter, &
                                 iFileName)
        ENDIF ! global%flowType
        E_IO = PVTK_GEO_XML(source=iFilename)
      END DO

      E_IO = PVTK_DAT_XML(var_location = 'cell', var_block_action = 'open')
      E_IO = PVTK_VAR_XML(varname = 'Gas Density', tp='Float64')
      E_IO = PVTK_VAR_XML(Nc = 3, varname = 'Momentum', tp='Float64' )
      E_IO = PVTK_VAR_XML(varname = 'Energy', tp='Float64')
      E_IO = PVTK_VAR_XML(varname = 'Pressure', tp='Float64')
      E_IO = PVTK_VAR_XML(varname = 'Temperature', tp='Float64')
      E_IO = PVTK_VAR_XML(varname = 'Speed of Sound', tp='Float64')
      !E_IO = PVTK_VAR_XML(varname = 'Products Energy', tp='Float64')
      E_IO = PVTK_DAT_XML(var_location = 'cell', var_block_action = 'Close')
      E_IO = PVTK_END_XML()

    END IF ! global%myProcid

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A)') '**************************************************'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing Slice PARAVIEW file done.'
      WRITE(STDOUT,'(A)') '**************************************************'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteSlice

!******************************************************************************
!
! Purpose: Writes the Paraview-Readable Particles files.
!
! Description: Theses files are meant for parallel read from PARAVIEW
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Write initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_WritePcls(pRegion)

    USE ModDataStruct, ONLY: t_level

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,i,iFile,iLoc,j,regId,iPcl,p
    TYPE(t_global), POINTER :: global

    INTEGER(I4P)::E_IO,Np

    INTEGER(I1P), ALLOCATABLE, DIMENSION(:) :: cell_type
    INTEGER(I4P), ALLOCATABLE, DIMENSION(:) :: offset,connect
    REAL(R4P), ALLOCATABLE, DIMENSION(:,:) :: points
    
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion,pLocal
    TYPE(t_level), POINTER :: levels(:)

! ******************************************************************************
!   Start, open vtu file
! ******************************************************************************

    global => pRegion%global


    CALL RegisterFunction(global,'RFLU_WritePcls',"../modflu/RFLU_ModReadWriteFlow.F90")


  CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WritePcls


!******************************************************************************
!
! Purpose: Write theta- and z-averaged variable of interest (\rho, p, vfp).
!
! Description: Theses files are meant for parallel read from PARAVIEW
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Write initial residual and physical time for both steady and unsteady
!      flows so that could use steady solution as input for unsteady run and
!      vice versa.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteRProf(regions)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,i,iFile,iLoc,j
 
    INTEGER :: Nsmpl,Nqty,iReg
    INTEGER, PARAMETER :: i_p=1,i_r=2,i_vf=3 
    REAL(RFREAL) :: maxrad_l,maxrad
    REAL(RFREAL), ALLOCATABLE, DIMENSION(:,:) :: avqty_l,avqty
    REAL(RFREAL), ALLOCATABLE, DIMENSION(:) :: rcnt_l,rcnt

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: regions(:)
    TYPE(t_region), POINTER :: pRegion
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open r-profile file
! ******************************************************************************

    global => regions(1)%global
    !pGrid => pRegion%grid

    CALL RegisterFunction(global,'RFLU_WriteRProf',"../modflu/RFLU_ModReadWriteFlow.F90")

    IF ( global%myProcid == MASTERPROC) THEN

      WRITE(STDOUT,'(A)') '************************************************'
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing r-profile file...'
      WRITE(STDOUT,'(A)') '************************************************'

      IF ( global%flowType == FLOW_UNSTEADY ) THEN
        CALL BuildFileNameUnsteadyPVTK(global,FILEDEST_OUTDIR,'.rprof.dat', &
                           global%currentTime,iFileName)
      ELSE
        CALL BuildFileNameSteadyPVTK(global,FILEDEST_OUTDIR,'.rprof.dat', &
                           global%currentIter,iFileName)
      END IF

    iFile = IF_SOLUT
    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,3935,iFileName)
    END IF ! global%error

    END IF

! ==============================================================================
!   Compute and Write Data
! ==============================================================================
  
    Nqty = 3

!    maxrad_l = MAXVAL(SQRT(pGrid%cofg(XCOORD,:)**2.0_RFREAL + &
!             pGrid%cofg(YCOORD,:)**2.0_RFREAL))

!    CALL MPI_Allreduce(maxrad_l,maxrad,1,MPI_INTEGER,MPI_MAX, &
!                       global%mpiComm,global%mpierr )

!    Nsmpl = NINT(maxrad)*1000
!    Nsmpl = Nsmpl + 1
    Nsmpl = 151   ! BBR -WARNING- value true only for case 1
                  ! need a clean way to get max radius

    ALLOCATE(avqty(1:Nsmpl,1:Nqty))
    ALLOCATE(avqty_l(1:Nsmpl,1:Nqty))
    ALLOCATE(rcnt_l(1:Nsmpl))
    ALLOCATE(rcnt(1:Nsmpl))
    avqty_l = 0.0_RFREAL;avqty = 0.0_RFREAL;
    rcnt_l = 0.0_RFREAL; rcnt = 0.0_RFREAL;

  DO iReg = 1,global%nRegionsLocal
    pRegion => regions(iReg)


!    Nqty = 3
 
!!    maxrad_l = MAXVAL(SQRT(pGrid%cofg(XCOORD,:)**2.0_RFREAL + &
!!             pGrid%cofg(YCOORD,:)**2.0_RFREAL))

!!    CALL MPI_Allreduce(maxrad_l,maxrad,1,MPI_INTEGER,MPI_MAX, &
!!                       global%mpiComm,global%mpierr )

!!    Nsmpl = NINT(maxrad)*1000
!!    Nsmpl = Nsmpl + 1
!    Nsmpl = 301   ! BBR -WARNING- value true only for case 1
                  ! need a clean way to get max radius

!    ALLOCATE(avqty(1:Nsmpl,1:Nqty))
!    ALLOCATE(avqty_l(1:Nsmpl,1:Nqty))
!    ALLOCATE(rcnt_l(1:Nsmpl))
!    ALLOCATE(rcnt(1:Nsmpl))
!    avqty_l = 0.0_RFREAL;avqty = 0.0_RFREAL;  
!    rcnt_l = 0.0_RFREAL; rcnt = 0.0_RFREAL;

    DO i = 1,pRegion%grid%nCells
       j = NINT(0.5*1000.*SQRT(pRegion%grid%xyz(XCOORD,i)**2 &
                        + pRegion%grid%xyz(YCOORD,i)**2)) + 1
       avqty_l(j,i_r) = avqty_l(j,i_r) + pRegion%mixt%cv(CV_MIXT_DENS,i) 
       avqty_l(j,i_p) = avqty_l(j,i_p) + pRegion%mixt%dv(DV_MIXT_PRES,i)
       rcnt_l(j) = rcnt_l(j) + 1.0_RFREAL
    END DO


  END DO !iReg

    CALL MPI_ALLREDUCE(avqty_l(:,i_r),avqty(:,i_r),Nsmpl,MPI_RFREAL, &
                       MPI_SUM,global%mpiComm,global%mpierr)
    CALL MPI_ALLREDUCE(avqty_l(:,i_p),avqty(:,i_p),Nsmpl,MPI_RFREAL, &
                       MPI_SUM,global%mpiComm,global%mpierr)
    CALL MPI_ALLREDUCE(avqty_l(:,i_vf),avqty(:,i_vf),Nsmpl,MPI_RFREAL, &
                       MPI_SUM,global%mpiComm,global%mpierr)
    CALL MPI_Allreduce(rcnt_l,rcnt,Nsmpl,MPI_RFREAL,MPI_SUM, &
                       global%mpiComm,global%mpierr )

    DEALLOCATE(avqty_l,rcnt_l)

    IF ( global%myProcid == MASTERPROC) THEN    

      DO i = 1,Nsmpl
        WRITE(iFile,'(4(1X,E13.6))') DBLE(i-1)/500. &
               ,avqty(i,i_p)/MAX(rcnt(i),1.0_RFREAL) &
               ,avqty(i,i_r)/MAX(rcnt(i),1.0_RFREAL) &
               ,avqty(i,i_vf)/MAX(rcnt(i),1.0_RFREAL) 
      END DO

      CLOSE(iFile,IOSTAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_FILE_CLOSE,4056,iFileName)
      END IF ! global%error
    END IF ! masterproc

    DEALLOCATE(avqty)
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteRProf





! ==============================================================================
! Purpose: Wrapper for writing of flow files in ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteFlowWrapper(pRegion)



    USE SPEC_RFLU_ModReadWriteVars, ONLY: SPEC_RFLU_WriteCvASCII, &
                                          SPEC_RFLU_WriteCvBinary, & 
                                          SPEC_RFLU_WriteEEvASCII, &
                                          SPEC_RFLU_WriteEEvBinary

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    TYPE(t_global), POINTER :: global
    CHARACTER(CHRLEN) :: NameSolDir
    LOGICAL :: flag
    INTEGER :: errorFlag

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteFlowWrapper',"../modflu/RFLU_ModReadWriteFlow.F90")

! ******************************************************************************
!   Write mixture solution files
! ******************************************************************************


! TLJ - this block sometimes gives me trouble
      IF ( global%flowType == FLOW_UNSTEADY ) THEN
      WRITE(NameSolDir,'(1PE11.5)') global%currentTime
      ELSEIF ( global%flowType == FLOW_STEADY ) THEN
      WRITE(NameSolDir,'(i10)') global%currentIter
      END IF
      IF(global%myProcid==MASTERPROC) THEN
        CALL system( 'mkdir -p SOL_'//NameSolDir )
        CALL system( 'mkdir -p PARAVIEW_'//NameSolDir )
      END IF

      ! Need to check if using MPI or not since this routine is called by the
      ! initializer, which doesn't use MPI
      CALL MPI_INITIALIZED(flag, errorFlag)
      IF (flag) CALL MPI_Barrier(global%mpiComm,errorFlag)
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL RFLU_WriteFlowASCII(pRegion)
!      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
! BBR - begin
       ELSE IF ( global%solutFormat == FORMAT_BINARY .OR. & 
                global%solutFormat == FORMAT_BINARY_L .OR. &
                global%solutFormat == FORMAT_BINARY_B ) THEN
! BBR - end
        CALL RFLU_WriteFlowBinary(pRegion)
        IF ( global%piclUsed .EQV. .TRUE. ) THEN
         CALL RFLU_PICL_WriteFlowBinary(pRegion)
        END IF
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,4178)
      END IF ! global%solutFormat

! ******************************************************************************
!   Write physical module solution files
! ******************************************************************************


! ==============================================================================
!   Species
! ==============================================================================

    IF ( global%specUsed .EQV. .TRUE. ) THEN
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL SPEC_RFLU_WriteCvASCII(pRegion)
                       
        IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 
          CALL SPEC_RFLU_WriteEEvASCII(pRegion)
        END IF ! pRegion%specInput%nSpeciesEE
!      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
! BBR - begin
       ELSE IF ( global%solutFormat == FORMAT_BINARY .OR. & 
                global%solutFormat == FORMAT_BINARY_L .OR. &
                global%solutFormat == FORMAT_BINARY_B ) THEN
! BBR - end
        CALL SPEC_RFLU_WriteCvBinary(pRegion)

        IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 
          CALL SPEC_RFLU_WriteEEvBinary(pRegion)
        END IF ! pRegion%specInput%nSpeciesEE
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,4238)
      END IF ! global%solutFormat
    END IF ! global%specUsed

! begin BBR - Write VTK
! TLJ added - sqdet (Sept 30, 2024)
! TLJ added - barrelExp (Sept 30, 2024)
! TLJ added - wedge (Nov 26, 2024)
! TLJ added - rectshktb (Jan 11, 2025)

IF ( global%casename .EQ. "cyldet" .OR. &
     global%casename .EQ. "shktb"  .OR. &
     global%casename .EQ. "cylds"  .OR. &
     global%casename .EQ. "rectshktb"  .OR. &
     global%casename .EQ. "barrelExp"  .OR. &
     global%casename .EQ. "wedge"  .OR. &
     global%casename .EQ. "sqdet"  .OR. &
     global%casename .EQ. "hass"   .OR. &
     global%casename .EQ. "xdsp"    ) THEN

  IF (global%moduleType .EQ. MODULE_TYPE_SOLVER)THEN
    CALL RFLU_WriteVTK(pRegion)
  END IF ! global%moduleType

END IF ! "cyldet"

! end BBR - Write VTK

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteFlowWrapper







! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModReadWriteFlow


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModReadWriteFlow.F90,v $
! Revision 1.6  2016/01/31 04:54:40  rahul
! Fixed a bug from previous check-in.
!
! Revision 1.5  2015/12/30 07:52:47  rahul
! Added capability to write .pvtu files for tetrahedral and prisms. Cylds and
! hass test-cases have been added to write the .pvtu files.
!
! Revision 1.4  2015/08/12 19:41:42  brollin
! Updating module declaration in rfluinit.F90
!
! Revision 1.3  2015/07/27 04:45:42  brollin
! 1) Corrected bug in RFLUCONV where global%gridFormat was used instead of global%gridSrcFormat
! 2) Implemented new subroutine for shock tube problems (Shktb)
!
! Revision 1.2  2015/07/23 23:11:18  brollin
! 1) The pressure coefficient of the  collision model has been changed back to its original form
! 2) New options in the format of the grid and solutions have been added. Now the user can choose the endianness, and convert from one to the over in rfluconv.
! 3) The solutions are now stored in folders named by timestamp or iteration number.
! 4) The address enty in the hashtable has been changed to an integer(8) for cases when the grid becomes very large.
! 5) RFLU_WritePM can now compute PM2 on the fly for the Macroscale problem
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:36  brollin
! New Stable version
!
! Revision 1.5  2008/12/06 08:43:44  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:16:57  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/11/30 05:53:48  mparmar
! Fixed a bug in writing sectionString in binary file
!
! Revision 1.2  2007/11/28 23:05:26  mparmar
! Added reading/writing of data set v2 while maintaining backward compatibility
!
! Revision 1.1  2007/04/09 18:49:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.11  2006/12/15 13:25:27  haselbac
! Fixed bug in format statement, found by ifort
!
! Revision 1.10  2006/03/26 20:22:07  haselbac
! Removed error traps for GL model
!
! Revision 1.9  2006/01/12 09:40:50  wasistho
! timeStamp to currentTime in turb readFlow
!
! Revision 1.8  2006/01/10 05:04:23  wasistho
! Get turbulence data from Genx
!
! Revision 1.7  2005/12/29 19:54:51  wasistho
! modified Rocturb part in ReadFlowWrapper
!
! Revision 1.6  2005/11/27 01:51:39  haselbac
! Added calls to EEv routines, changed extensions in backw-compatible way
!
! Revision 1.5  2004/11/06 03:18:26  haselbac
! Substantial additions to allow reading/writing of data for other fluid models
!
! Revision 1.4  2004/10/19 19:28:24  haselbac
! Adapted to changes in GENX logic
!
! Revision 1.3  2004/08/23 23:08:43  fnajjar
! Activated binary IO routine calls
!
! Revision 1.2  2004/07/28 15:29:20  jferry
! created global variable for spec use
!
! Revision 1.1  2004/07/06 15:14:31  haselbac
! Initial revision
!
! ******************************************************************************

