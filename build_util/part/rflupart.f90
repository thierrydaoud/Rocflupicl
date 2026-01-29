










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
! Purpose: Driver routine for rflupart. 
!
! Description: None.
!
! Input: 
!   caseString	String with casename
!   verbLevel	Verbosity level
!
! Output: None.
!
! Notes:
!   1. There is no call to RFLU_AllocateMemoryWrapper and the corresponding 
!      call to RFLU_DeallocateMemoryWrapper because the grid quantities are 
!      determined in the grid conversion routines. 
!
! ******************************************************************************
!
! $Id: rflupart.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE rflupart(caseString,verbLevel)

  USE ModError
  USE ModDataTypes
  USE ModParameters  
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModMPI
  
  USE RFLU_ModAllocateMemory
  USE RFLU_ModBoundLists
  USE RFLU_ModCellMapping
  USE RFLU_ModColoring
  USE RFLU_ModCommLists
  USE RFLU_ModDeallocateMemory
  USE RFLU_ModDimensions
  USE RFLU_ModFaceList
  USE RFLU_ModGeometry
  USE RFLU_ModGridSpeedUtils
  USE RFLU_ModGridUtils
  USE RFLU_ModPartitionRegion
  USE RFLU_ModPatchUtils
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModReadWriteGridSpeeds
  USE RFLU_ModRegionMapping
  USE RFLU_ModRenumberings
  USE RFLU_ModStencilsCells
  USE RFLU_ModStencilsUtils
  USE RFLU_ModSymmetryPeriodic
  USE RFLU_ModVertexLists


  USE ModInterfaces, ONLY: RFLU_BuildDataStruct, &   
                           RFLU_CreateGrid, &
                           RFLU_DestroyGrid, &
                           RFLU_GetUserInput, &
                           RFLU_InitGlobal, & 
                           RFLU_PrintGridInfo, & 
                           RFLU_PrintHeader, & 
                           RFLU_PrintWarnInfo, &
                           RFLU_RandomInit, &
                           RFLU_ReadConvGridWrapper, &
                           RFLU_ReadRestartInfo, & 
                           RFLU_SetModuleType, &  
                           RFLU_SetRestartTimeFlag, &
                           RFLU_USER_EnforcePatchCoords, & 
                           RFLU_WriteRestartInfo, & 
                           RFLU_WriteVersionString, &
                           ScaleRotateVector
                           
  IMPLICIT NONE


! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*) :: caseString
  INTEGER, INTENT(IN) :: verbLevel

!BAD 2/10/15 for parallel first attempt variable init
  INTEGER :: ERROR_M, world_size, world_rank,chunk_size,Reg_start,Reg_end,left_overs,iVar
!TYPE(t_region), DIMENSION(:), ALLOCATABLE :: regionDataBuff
!TYPE(t_region), DIMENSION(:),ALLOCATABLE :: regionDataBuff_get

!integer, DIMENSION(:), ALLOCATABLE :: BorderCntrDataBuff
!integer, DIMENSION(:),ALLOCATABLE :: BorderCntrDataBuff_get
!END BAD End variable init

!BAD TIMING VARAIBLE DECLARE

   !INTEGER :: clock_rate, clock_max, timeroutineat, dummy

   !REAL(RFREAL) :: elapsedtime_serial, elapsedtime_total


!END BAD TIMING VARAIBLE



! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: casename
  INTEGER :: errorFlag,iLev,iReg
  TYPE(t_region), POINTER :: pRegion,pRegionSerial
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER :: levels(:)

! TEMPORARY
!  INTEGER :: order
! END TEMPORARY

  ! BBR
  CHARACTER(CHRLEN) :: ext

! Sam - logical to switch between parallel and serial partitioning
  LOGICAL ppart
! Sam - end

! ******************************************************************************
! Start, initialize global data
! ******************************************************************************  

  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR - Pointer allocation failed.'
    STOP
  END IF ! errorFlag 

  casename = caseString(1:LEN(caseString))

  CALL RFLU_InitGlobal(casename,verbLevel,MPI_COMM_WORLD,global)

  CALL RegisterFunction(global,'rflupart', & 
                        "../../utilities/part/rflupart.F90")

  CALL RFLU_SetModuleType(global,MODULE_TYPE_PART)

!BAD TIMING

!Call system_clock(dummy, clock_rate, clock_max)

!END BAD TIMING

  CALL MPI_Init(ERROR_M)
  CALL MPI_Comm_size(MPI_COMM_WORLD, world_size, ERROR_M)
  CALL MPI_COMM_rank(MPI_COMM_WORLD, world_rank, ERROR_M)

  ppart = (world_size .gt. 1) ! parallel partitioning if world_size > 1

!To prevent printing bottleneck
  IF (world_rank > 0) then

    global%myProcid = world_rank

  END IF

! ******************************************************************************
! Print header and write version string
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC ) THEN
    CALL RFLU_WriteVersionString(global)     
    IF ( global%verbLevel /= VERBOSE_NONE ) THEN
      CALL RFLU_PrintHeader(global)
    END IF ! global%verbLevel
  END IF ! global%myProcid

! ******************************************************************************
! Read mapping file, impose serial mapping, and build basic data structure
! ******************************************************************************

  CALL RFLU_ReadRegionMappingFile(global,MAPFILE_READMODE_PEEK,global%myProcId)
  CALL RFLU_SetRegionMappingSerial(global)  
  CALL RFLU_CreateRegionMapping(global,MAPTYPE_REG)
  CALL RFLU_ImposeRegionMappingSerial(global)

  CALL RFLU_BuildDataStruct(global,levels) 
  CALL RFLU_ApplyRegionMapping(global,levels)
  CALL RFLU_DestroyRegionMapping(global,MAPTYPE_REG)  


! ******************************************************************************
! Initialize random number generator. NOTE needed in order to write sensible 
! data when writing Rocpart solution files. 
! ******************************************************************************

  CALL RFLU_RandomInit(levels(1)%regions)

! ******************************************************************************
! Read input file and restart info. NOTE need restart info for GENX runs to 
! determine whether have a restart. 
! ******************************************************************************

  CALL RFLU_GetUserInput(levels(1)%regions,.TRUE.) 
  IF ( global%myProcid == MASTERPROC ) THEN !Added to prevent mpi hanging BAD
    CALL RFLU_ReadRestartInfo(global)
    CALL RFLU_SetRestartTimeFlag(global)
  END IF ! global%myProcid
 ! BBR - begin - creating folder SOL_0.00000E+00 and PARAVIEW_0.00000E+00 
 IF ( global%myProcid == MASTERPROC ) THEN
   CALL system( 'mkdir GRIDDATA' )
   IF ( global%flowtype == FLOW_UNSTEADY ) THEN  
     WRITE(ext,'(1PE11.5)')global%currentTime
     CALL system( 'mkdir SOL_'//TRIM(ext) )
     CALL system( 'mkdir PARAVIEW_'//TRIM(ext) )
   ELSE
     WRITE(ext,'(I6.6)')global%currentIter
     CALL system( 'mkdir SOL_'//TRIM(ext) )
     CALL system( 'mkdir PARAVIEW_'//TRIM(ext) )
   END IF 
 END IF ! global%myProcid
 ! BBR - end 

! ******************************************************************************
! Read and convert grid
! ******************************************************************************

  pRegionSerial => levels(1)%regions(0) ! NOTE must set otherwise get core dump 
    
  CALL RFLU_ReadConvGridWrapper(pRegionSerial)

! ******************************************************************************
! Read boundary condition file
! ******************************************************************************
    
  IF ( pRegionSerial%grid%nPatches > 0 ) THEN        
    CALL RFLU_ReadBCInputFileWrapper(pRegionSerial)    
  END IF ! pRegionSerial%grid%nPatches 
  
! ******************************************************************************
! Set maximum dimensions
! ******************************************************************************

  CALL RFLU_SetMaxDimensions(pRegionSerial)

! ******************************************************************************
! Impose patch coordinates, scale and rotate if desired
! ******************************************************************************

  IF ( global%enforceFlag .EQV. .TRUE. ) THEN 
    CALL RFLU_USER_EnforcePatchCoords(pRegionSerial)
  END IF ! global%enforceFlag
  
  IF ( global%transformFlag .EQV. .TRUE. ) THEN 
    CALL ScaleRotateVector(global,pRegionSerial%grid%xyz)
  END IF ! global%transformFlag

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    CALL RFLU_PrintGridInfo(pRegionSerial)
  END IF ! global%verbLevel

! ******************************************************************************
! Build data structures
! ******************************************************************************

  CALL RFLU_CreateCellMapping(pRegionSerial)
  CALL RFLU_BuildLoc2GlobCellMapping(pRegionSerial)
  CALL RFLU_BuildGlob2LocCellMapping(pRegionSerial)

  CALL RFLU_CreateBVertexLists(pRegionSerial)
  CALL RFLU_BuildBVertexLists(pRegionSerial)

  CALL RFLU_CreateFaceList(pRegionSerial)
  CALL RFLU_BuildFaceList(pRegionSerial)
  CALL RFLU_RenumberBFaceLists(pRegionSerial)

! ******************************************************************************
! Distort grid (NOTE needs to be done after having boundary vertex lists)
! ******************************************************************************

  IF ( global%distortFlag .EQV. .TRUE. ) THEN
    CALL RFLU_DistortGrid(pRegionSerial)
  END IF ! global%distortFlag

! ******************************************************************************
! Add virtual cells on symmetry and periodic boundaries
! ******************************************************************************

  IF ( RFLU_SYPE_HaveSyPePatches(pRegionSerial) .EQV. .TRUE. ) THEN
    CALL RFLU_CreatePatchNeighborMaps(pRegionSerial)
    CALL RFLU_BuildPatchNeighborMaps(pRegionSerial)
  
    CALL RFLU_COMM_CountBordersSerial(pRegionSerial)
    CALL RFLU_COMM_CreateBorders(pRegionSerial)
    
    CALL RFLU_CreateGeometry(pRegionSerial)
    CALL RFLU_BuildGeometry(pRegionSerial,.FALSE.) 
   
    CALL RFLU_ComputePatchNormalsLocal(pRegionSerial)
    CALL RFLU_CheckPatchBcConsistency(pRegionSerial)
    CALL RFLU_SYPE_BuildTransforms(pRegionSerial)
    CALL RFLU_SYPE_WriteTransforms(pRegionSerial)    

    CALL RFLU_SYPE_CreateVertexMaps(pRegionSerial)
    CALL RFLU_SYPE_BuildVertexMaps(pRegionSerial)

    CALL RFLU_CreateVert2CellList(pRegionSerial)
    CALL RFLU_BuildVert2CellList(pRegionSerial)

    CALL RFLU_SYPE_AddVirtualCells(pRegionSerial)
    CALL RFLU_SYPE_BuildP2VCListSerial(pRegionSerial)    
    
    CALL RFLU_DestroyVert2CellList(pRegionSerial)
    CALL RFLU_DestroyGeometry(pRegionSerial)
    
    CALL RFLU_DestroyPatchNeighborMaps(pRegionSerial)    
    
    CALL RFLU_PrintGridInfo(pRegionSerial)
    
    CALL RFLU_COMM_WriteCommLists(pRegionSerial)
    
    CALL RFLU_DestroyBVertexLists(pRegionSerial)
    CALL RFLU_DestroyFaceList(pRegionSerial)
    
    CALL RFLU_CreateBVertexLists(pRegionSerial)
    CALL RFLU_BuildBVertexLists(pRegionSerial)

    CALL RFLU_CreateFaceList(pRegionSerial)
    CALL RFLU_BuildFaceList(pRegionSerial)
    CALL RFLU_RenumberBFaceLists(pRegionSerial)    
  END IF ! RFLU_SYPE_HaveSyPePatches  

! ******************************************************************************
! Allocate memory 
! ******************************************************************************
 
  IF ( RFLU_DecideNeedGridSpeeds(pRegionSerial) .EQV. .TRUE. ) THEN 
    CALL RFLU_AllocateMemoryGSpeeds(pRegionSerial)
  END IF ! RFLU_DecideNeedGridSpeeds

! ******************************************************************************
! Write dimensions. Also write cell mapping for serial grid, which is needed by 
! post-processor if the regions are merged. Write grid if running in serial or
! if running with particles. In the latter case, the serial grid is needed to 
! initialize the particle field.
! ******************************************************************************
  
  CALL RFLU_WriteDimensions(pRegionSerial)  
  CALL RFLU_WriteLoc2GlobCellMapping(pRegionSerial)

  IF ( (global%nRegionsLocal == 1) .OR. (global%plagUsed .EQV. .TRUE.) ) THEN 
    CALL RFLU_WriteGridWrapper(pRegionSerial)
  END IF ! global%nRegionsLocal

! ******************************************************************************
! Build coloring for implicit solver, write for serial runs
! ******************************************************************************

! TO DO
! Need to integrate this better with partitioning so that there is no need to
! build geometry twice. This probably means that need to have two calls to 
! build coloring, one for purely serial runs and another for parallel runs.
! END TO DO

  IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN 
    IF ( pRegionSerial%mixtInput%spaceOrder > 1 ) THEN
      CALL RFLU_CreateVert2CellList(pRegionSerial)
      CALL RFLU_BuildVert2CellList(pRegionSerial)

      CALL RFLU_CreateGeometry(pRegionSerial)
      CALL RFLU_BuildGeometry(pRegionSerial)

      CALL RFLU_SetInfoC2CStencilWrapper(pRegionSerial, & 
                                         pRegionSerial%mixtInput%spaceOrder-1)
      CALL RFLU_CreateC2CStencilWrapper(pRegionSerial)
      CALL RFLU_BuildC2CStencilWrapper(pRegionSerial,constrInput=CONSTR_NONE)

      CALL RFLU_DestroyGeometry(pRegionSerial)
      CALL RFLU_DestroyVert2CellList(pRegionSerial)
    END IF ! pRegionSerial%mixtInput%spaceOrder

    CALL RFLU_CreateCell2FaceList(pRegionSerial)  
    CALL RFLU_BuildCell2FaceList(pRegionSerial)
  
    CALL RFLU_COL_CreateColoring(pRegionSerial)
    CALL RFLU_COL_BuildColoring(pRegionSerial)

    IF ( pRegionSerial%mixtInput%spaceOrder > 1 ) THEN
      CALL RFLU_DestroyC2CStencilWrapper(pRegionSerial)
    ELSE
      CALL RFLU_DestroyCell2FaceList(pRegionSerial)      
    END IF ! pRegionSerial%mixtInput%spaceOrder

    IF ( global%nRegionsLocal == 1 ) THEN 
      CALL RFLU_COL_WriteColoring(pRegionSerial)
    END IF ! global%nRegionsLocal
  END IF ! global%solverType

! ******************************************************************************
! Partition
! ******************************************************************************

  IF ( global%nRegionsLocal > 1 ) THEN 

! ==============================================================================
!   Partition region and build lists
! ==============================================================================

    CALL RFLU_RNMB_CreateSC2RMap(pRegionSerial)    
    CALL RFLU_PART_PartitionRegion(pRegionSerial)  
    CALL RFLU_RNMB_WriteSC2RMap(pRegionSerial)
    
    CALL RFLU_PART_CreateReg2CellMap(pRegionSerial)    
    CALL RFLU_PART_BuildReg2CellMap(pRegionSerial)
    
    CALL RFLU_PART_BuildBorderFaceList(pRegionSerial)

    CALL RFLU_CreateVert2CellList(pRegionSerial)
    CALL RFLU_BuildVert2CellList(pRegionSerial)

    CALL RFLU_CreateBCellMList(pRegionSerial)
    CALL RFLU_BuildBCellMList(pRegionSerial) 

    IF ( pRegionSerial%mixtInput%spaceOrder > 1 ) THEN 
      CALL RFLU_CreateGeometry(pRegionSerial)
      CALL RFLU_BuildGeometry(pRegionSerial) 
    END IF ! pRegionSerial
         
! ==============================================================================
!   Write out grid, solution, and renumbering lists for each region
! ==============================================================================

!BAD 2/10/15 first attempt at parallel part
!"Hold on to your butts" - Sam Jackson

!write(*,*)'**********************************************************************'
!write(*,*) 'BAD P TEST Part 0'


!6-3-15 the MPI init size and rank moved to top to stop print bottleneck
!CALL MPI_Init(ERROR_M)
!CALL MPI_Comm_size(MPI_COMM_WORLD, world_size, ERROR_M)
!CALL MPI_COMM_rank(MPI_COMM_WORLD, world_rank, ERROR_M)

!BAD TIMING

    IF (ppart) then
    
      CALL MPI_BARRIER(MPI_COMM_WORLD, ERROR_M)
      
      !IF (world_rank == 0) then
      !
      !        Call SYSTEM_CLOCK(timeroutineat, clock_rate, clock_max)
      !
      !        elapsedtime_serial = REAL(timeroutineat - dummy) / Real(clock_rate)
      !
      !        open (998,file = 'BAD_RUN_TIMES_serial.tim', status = 'unknown')
      !
      !        write (998,*) 'Serial time = ',elapsedtime_serial
      !
      !        close(998)
      !
      !
      !
      !END IF
      !END BAD TIMING
      
      !write(*,*)'**********************************************************************'
      !write(*,*) 'BAD P TEST Part 1'
      
      chunk_size  = global%nRegionsLocal / world_size
      Reg_start = (chunk_size * world_rank) + 1
      Reg_end = (world_rank+1) * chunk_size
      
      !write (*,*) world_rank,world_size,chunk_size,Reg_start,Reg_end
      
      left_overs = mod(global%nRegionsLocal, world_size)
      
      IF (world_rank == world_size -1) then
              Reg_end = Reg_end + left_overs
      END IF 
    ELSE
      Reg_start = 1
      Reg_end = global%nRegionsLocal
    END IF

!write(*,*) world_rank,Reg_end

!    DO iReg = 1,global%nRegionsLocal OLD CODE
    DO iReg = Reg_start,Reg_end       
      pRegion => levels(1)%regions(iReg)

      CALL RFLU_PART_CreateCellLists(pRegion,pRegionSerial)
      CALL RFLU_RNMB_CreatePC2SCMap(pRegion)

      CALL RFLU_CreateCellMapping(pRegion)
      CALL RFLU_BuildLoc2GlobCellMapping(pRegion)
      CALL RFLU_BuildGlob2LocCellMapping(pRegion)      
      
      CALL RFLU_PART_BuildCellLists(pRegion,pRegionSerial)
      CALL RFLU_PART_AddVirtualCells(pRegion,pRegionSerial)      
      CALL RFLU_RNMB_BuildSC2PCMap(pRegion)

      CALL RFLU_PART_BuildVertexLists(pRegion,pRegionSerial)
      CALL RFLU_RNMB_BuildSV2PVMap(pRegion)
      
      CALL RFLU_RNMB_BuildSBC2PCMap(pRegion,pRegionSerial)      
      CALL RFLU_PART_CreatePatchLists(pRegion,pRegionSerial)
      CALL RFLU_RNMB_CreatePBF2SBFMap(pRegion)
      CALL RFLU_PART_BuildPatchLists(pRegion,pRegionSerial)
      CALL RFLU_RNMB_DestroySBC2PCMap(pRegion)
      CALL RFLU_PART_RenumberVertexLists(pRegion)
      CALL RFLU_PART_BuildVertexData(pRegion,pRegionSerial) 
      
      CALL RFLU_SetMaxDimensions(pRegion)
      CALL RFLU_WriteDimensions(pRegion)
      CALL RFLU_WriteLoc2GlobCellMapping(pRegion)


      CALL RFLU_WriteGridWrapper(pRegion)

      
      CALL RFLU_RNMB_WritePxx2SxxMaps(pRegion)
      
      IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN 
        CALL RFLU_COL_CreateColoring(pRegion)
        CALL RFLU_COL_BuildColoring(pRegion,pRegionSerial)
        CALL RFLU_COL_WriteColoring(pRegion)
        CALL RFLU_COL_DestroyColoring(pRegion)
      END IF ! global%solverType      
      
      CALL RFLU_DestroyCellMapping(pRegion)
      
      CALL RFLU_PART_DestroyVertexData(pRegion)
                        
      CALL RFLU_RNMB_DestroyPC2SCMap(pRegion)
      CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegion)
      CALL RFLU_RNMB_DestroyPV2SVMap(pRegion)      
      
      CALL RFLU_RNMB_DestroySC2PCMap(pRegion)
      CALL RFLU_RNMB_DestroySV2PVMap(pRegion)

      CALL RFLU_PART_DestroyPatchLists(pRegion)      
      CALL RFLU_PART_DestroyCellLists(pRegion)       
    END DO ! iReg

!CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR_M)
!write(*,*)'**********************************************************************'
!write(*,*) 'BAD P TEST Part 2'

    CALL RFLU_DestroyBCellMList(pRegionSerial)        
    CALL RFLU_DestroyVert2CellList(pRegionSerial)   
    CALL RFLU_PART_BuildBorderFaceList(pRegionSerial)

    IF ( pRegionSerial%mixtInput%spaceOrder > 1 ) THEN 
      CALL RFLU_DestroyGeometry(pRegionSerial)
    END IF ! pRegionSerial

! ==============================================================================
!   Build communication lists
! ==============================================================================

! ------------------------------------------------------------------------------
!   Read dimensions and renumbering lists. NOTE need to read bc file because of
!   sype cases, in which mapping from patches to virtual cells is constructed
!   below.
! ------------------------------------------------------------------------------

!    DO iReg = 1,global%nRegionsLocal
    DO iReg = Reg_start,Reg_end    
      pRegion => levels(1)%regions(iReg)

      CALL RFLU_ReadDimensions(pRegion) 
! TEMPORARY - THIS IS A BUG - MUST NOT BE CALLED BEFORE CREATEGRID BCOS 
!             PATCHES NOT YET CREATED - MAX DIMS READ FROM FILE ANYWAY
!      CALL RFLU_SetMaxDimensions(pRegion)
! END TEMPORARY              
              
      CALL RFLU_CreateGrid(pRegion)
      
      IF ( pRegion%grid%nPatches > 0 ) THEN        
        CALL RFLU_ReadBCInputFileWrapper(pRegion)    
      END IF ! pRegion%grid%nPatches       


      CALL RFLU_ReadGridWrapper(pRegion)
     
      CALL RFLU_CreateCellMapping(pRegion)
      CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
      CALL RFLU_BuildGlob2LocCellMapping(pRegion)
     
      CALL RFLU_RNMB_CreatePC2SCMap(pRegion)    
      CALL RFLU_RNMB_CreatePV2SVMap(pRegion)
      CALL RFLU_RNMB_CreatePBF2SBFMap(pRegion)
      
      CALL RFLU_RNMB_ReadPxx2SxxMaps(pRegion)
            
      CALL RFLU_RNMB_BuildSC2PCMap(pRegion)
      CALL RFLU_RNMB_BuildSV2PVMap(pRegion)
      
! TO DO DELETE WINDOW
!
! END TO DO DELETE WINDOW
    END DO ! iReg
!CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR_M)
!write(*,*)'**********************************************************************'
!write(*,*) 'BAD P TEST Part 3'

! ------------------------------------------------------------------------------
!   Count number of borders, create borders, and build communication lists
! ------------------------------------------------------------------------------

!    DO iReg = 1,global%nRegionsLocal
    Do iReg = Reg_start,Reg_end
      pRegion => levels(1)%regions(iReg)
    
      CALL RFLU_COMM_CreateBorderCntr(pRegion)
      CALL RFLU_COMM_CountBorders(pRegion,pRegionSerial)      
    END DO ! iReg
!CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR_M)
!write(*,*)'**********************************************************************'
!write(*,*) 'BAD P TEST Part 4'

!Need to compile all the seperate regions(iReg) into one "levels(1)%regions" file 

! Variable statment needed, may not need INTENT, CALL it regionDataBuff?
!REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: cellDataBuff
!ALLOCATE (BorderCntrDataBuff(Reg_end-Reg_start+1) )
!ALLOCATE (BorderCntrDataBuff_get( global%nRegionsLocal))
!iReg = Reg_start
!do iVar = 1,(Reg_end-Reg_start+1)
!        BorderCntrDataBuffDataBuff(iVar)=&
!levels(1)%regions(iReg)%grid%bordercntr
!        write(*,*) 'iRegionGlobal =',regionDataBuff(iVar)%iRegionGlobal        
!        iReg = iReg + 1
!write(*,*) ' Iteration done'
!end do
!write(*,*) ' Compiled worked'
!CALL MPI_AllGather(levels(1)%regions(1)%iRegionGlobal,(1), MPI_INTEGER,&
!              levels(1)%regions(5)%iRegionGlobal, (1), MPI_INTEGER,MPI_COMM_WORLD,&
!ERROR_M)
!write(*,*) 'All Gather worked'
!write(*,*) 'PASS TEST =',levels(1)%regions(1)%iRegionGlobal
!Unpack the data
!do iReg = 1,global%nRegionsLocal
!        levels(1)%regions(iReg) = regionDataBuff_get(iReg)
!        write(*,*) 'iRegionGlobal,get =',regionDataBuff_get(iReg)%iRegionGlobal
!        write(*,*) 'iRegionGlobal,final =',levels(1)%regions(iReg)%iRegionGlobal
!end do
!write(*,*) 'Decompiled worked'

    CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR_M)!might need for consistency

    CALL RFLU_COMM_CheckCountBorders(levels(1)%regions)

    CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR_M)!might need for consistency

!write(*,*)'**********************************************************************'
!write(*,*) 'BAD P TEST Part 4.5'
!    DO iReg = 1,global%nRegionsLocal
    Do iReg = Reg_start,Reg_end
      pRegion => levels(1)%regions(iReg)
    
      CALL RFLU_COMM_CreateBorders(pRegion,CREATE_BORDERS_MODE_DIM_UNKNOWN)    
    END DO ! iReg

!CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR_M)
!write(*,*)'**********************************************************************'
!write(*,*) 'BAD P TEST Part 5'
    CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR_M)!might need for consistency
    CALL RFLU_COMM_BuildCommLists(levels(1)%regions)      
    CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR_M)!might need for consistency

!    DO iReg = 1,global%nRegionsLocal
    DO iReg = Reg_start,Reg_end
      pRegion => levels(1)%regions(iReg)
    
      CALL RFLU_COMM_DestroyBorderCntr(pRegion)    
    END DO ! iReg
!CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR_M)
!write(*,*)'**********************************************************************'
!write(*,*) 'BAD P TEST Part 6'
! ------------------------------------------------------------------------------
!   Write communication lists, destroy renumbering and communication lists. 
!   NOTE need to write dimensions again because now also have information about
!   number of borders.
! ------------------------------------------------------------------------------

!    DO iReg = 1,global%nRegionsLocal
    DO iReg = Reg_start,Reg_end    
      pRegion => levels(1)%regions(iReg)
     
      CALL RFLU_WriteDimensionsBorders(pRegion)
      
      CALL RFLU_COMM_WriteCommLists(pRegion)

      IF ( RFLU_SYPE_HaveSyPePatches(pRegion) .EQV. .TRUE. ) THEN 
        CALL RFLU_SYPE_BuildP2VCList(pRegion,pRegionSerial)
        CALL RFLU_WriteDimensions(pRegion) 
        CALL RFLU_WriteGridWrapper(pRegion)  
        CALL RFLU_SYPE_DestroyP2VCList(pRegion)                 
      END IF ! RFLU_SYPE_HaveSyPePatches

     
      CALL RFLU_DestroyCellMapping(pRegion)
     
      CALL RFLU_RNMB_DestroyPC2SCMap(pRegion)    
      CALL RFLU_RNMB_DestroyPV2SVMap(pRegion)
      CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegion)   

      CALL RFLU_RNMB_DestroySC2PCMap(pRegion)
      CALL RFLU_RNMB_DestroySV2PVMap(pRegion)      
   
      CALL RFLU_DestroyGrid(pRegion)
    END DO ! iReg
!CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR_M)
!write(*,*)'**********************************************************************'
!write(*,*) 'BAD P TEST Part 7'

! ==============================================================================
!   Destroy partition mapping
! ==============================================================================

    CALL RFLU_RNMB_DestroySC2RMap(pRegionSerial) 
  END IF ! global%nRegionsLocal                 

! ******************************************************************************
! Write restart info file
! ******************************************************************************
  IF (world_rank == 0) then
    CALL RFLU_WriteRestartInfo(global) 
  END IF
! ******************************************************************************
! Deallocate memory 
! ******************************************************************************
    
  IF ( RFLU_SYPE_HaveSyPePatches(pRegionSerial) .EQV. .TRUE. ) THEN 
    CALL RFLU_SYPE_DestroyP2VCList(pRegionSerial)                 
  END IF ! RFLU_SYPE_HaveSyPePatches    
    
  IF ( RFLU_DecideNeedGridSpeeds(pRegionSerial) .EQV. .TRUE. ) THEN   
    CALL RFLU_DeallocateMemoryGSpeeds(pRegionSerial)
  END IF ! RFLU_DecideNeedGridSpeeds
  
  IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN 
    CALL RFLU_COL_DestroyColoring(pRegionSerial)
  END IF ! global%solverType  
  
  CALL RFLU_DestroyFaceList(pRegionSerial)  
  CALL RFLU_DestroyBVertexLists(pRegionSerial)
  CALL RFLU_DestroyCellMapping(pRegionSerial)

! ******************************************************************************
! Print info about warnings
! ******************************************************************************
  IF (world_rank == 0) then
    CALL RFLU_PrintWarnInfo(global)
  END IF

!BAD TIMING

  CALL MPI_BARRIER(MPI_COMM_WORLD, ERROR_M)

!IF (world_rank == 0) then
!
!        Call SYSTEM_CLOCK(timeroutineat, clock_rate, clock_max)
!
!        elapsedtime_total = REAL(timeroutineat - dummy) / Real(clock_rate)
!
!        open (999,file = 'BAD_RUN_TIMES.tim', status = 'unknown')
!
!        write (999,*) 'Serial time = ',elapsedtime_serial, ' Total_time = ', &
!                        elapsedtime_total
! 
!        close(999)
!END IF


!END BAD TIMING



! ******************************************************************************
! End
! ******************************************************************************
  CALL MPI_Finalize(ERROR_M)
!!!!!!!!!!END BAD  
  CALL DeregisterFunction(global)

END SUBROUTINE rflupart

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: rflupart.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:56  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:11  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:56:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.20  2006/08/18 14:05:26  haselbac
! Added call to write transforms
!
! Revision 1.19  2006/08/04 03:06:52  haselbac
! Added grid distortion capability
!
! Revision 1.18  2006/04/17 19:57:50  haselbac
! Bug fix: Must read bc file before building patch-to-virtual cell mapping
!
! Revision 1.17  2006/04/15 17:02:12  haselbac
! Bug fix: bc file must be read before face list construction
!
! Revision 1.16  2006/03/25 22:06:23  haselbac
! Added capability of dealing with sype patches
!
! Revision 1.15  2006/02/06 23:55:55  haselbac
! Added comm argument to RFLU_InitGlobal
!
! Revision 1.14  2006/01/06 22:17:38  haselbac
! Adapted to name changes
!
! Revision 1.13  2005/10/27 19:21:39  haselbac
! Adapted to changes in stencil routine names
!
! Revision 1.12  2005/10/05 14:25:51  haselbac
! Adapted to changes in stencil modules, added use of vertex list module
!
! Revision 1.11  2005/09/22 17:13:59  hdewey2
! Moved stencil building calls for the implicit solver so the Jacobian can always be colored based on the 1st order stencils.
!
! Revision 1.10  2005/09/19 18:41:24  haselbac
! Bug fix for asymmetric borders: Now check and correct
!
! Revision 1.9  2005/08/24 01:39:28  haselbac
! Added capability to build coloring for second-order runs
!
! Revision 1.8  2005/08/19 16:37:53  haselbac
! Bug fix: Added calls to build cell-to-face list before coloring
!
! Revision 1.7  2005/08/19 02:36:54  haselbac
! Added calls to coloring routines
!
! Revision 1.6  2005/07/01 16:19:47  haselbac
! Set max dims bfor writing to reduce mem in rfluinit and rflump
!
! Revision 1.5  2005/05/09 17:44:25  haselbac
! Bug fix: Call RFLU_ReadDimensions, not wrapper
!
! Revision 1.4  2005/05/05 18:38:39  haselbac
! Removed MPI calls after bug in Rocin/out fixed
!
! Revision 1.3  2005/05/04 03:36:46  haselbac
! Added init and finalize of MPI when running within GENX
!
! Revision 1.2  2005/05/03 03:11:31  haselbac
! Converted to C++ reading of command-line
!
! Revision 1.1  2005/04/15 15:09:15  haselbac
! Initial revision
!
! ******************************************************************************

