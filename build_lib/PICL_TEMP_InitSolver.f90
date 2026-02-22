










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

SUBROUTINE PICL_TEMP_InitSolver( pRegion)

!  USE 

!USE ModInterfaces, ONLY: 
 USE ModDataTypes
  USE ModDataStruct, ONLY : t_level,t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMaterials, ONLY  : t_material
  USE ModError
  USE ModParameters
  USE ModMPI      
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt

!1
  use ppiclf_solve, only: ppiclf_solve_InitParticle, ppiclf_solve_Initialize
  use ppiclf_m_comm, only: ppiclf_comm_InitOverlapGrid
  USE ModRandom, ONLY: Rand1Uniform,Rand1Normal
  USE RFLU_ModInCellTest

  IMPLICIT NONE
!DEC$ NOFREEFORM

! number of timesteps kept in history kernels

!Change here when viscous unsteady on
!#define PPICLF_VU 0
!#define PPICLF_LRP3 6*PPICLF_VU

! maximum number of triangular patch boundaries

! y, y1, ydot, ydotc: 12

! rprop: 48

! rprop5: 0 - Storing Force Models

! map: 22
!--- Particle Volume Fraction Feedback
!--- x,y,z Forces Feedback
!---Energy Feedback
!--- More VF quanities. ***NEED TO CONFIRM THEY ARE USED ***
!--- Reynolds Subgrid Stress Tensor
!--- Pseudo Turbulent Kinetic Energy


























!DEC$ FREEFORM



! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

!TYPE(t_region), DIMENSION(:), POINTER :: regions
TYPE(t_global), POINTER :: global
!TYPE(t_level), POINTER :: levels(:)
TYPE(t_region), POINTER :: pRegion
TYPE(t_grid), POINTER :: pGrid
TYPE(t_material), POINTER :: material
INTEGER :: errorFlag,icg

   CHARACTER(CHRLEN) :: endString,iFileName,matName,comment
   CHARACTER(12) :: vtuFile,vtuFile1
   LOGICAL :: notfoundFlag, pf_restart,pf_rpInit,pf_settle,&
              wall_exists, fexists, foundMat
   INTEGER :: i,npart,nCells,vi,vii,ii,jj,kk,loopCounter,ipart,icl
   INTEGER :: PPC,numPclCells,npart_local,i_global,i_global_min,i_global_max,&
              iFile,iMat, k, j, l, m
   REAL(RFREAL) :: dp_min,dp_max,rhop,tester,ratio,total_vol,xMinCurt,&
                   xMaxCurt,yMinCurt,yMaxCurt,xMinCell,xMaxCell,yMinCell,&
                   yMaxCell,zMinCell,zMaxCell,x,vFrac,volpclsum,xLoc,yLoc,zLoc,yL, &
                   zpf_factor,xpf_factor,dp,neighborWidth,xp_min,xp_max, &
                   yp_min, yp_max, zp_min, zp_max, MinFluidCells
   REAL(RFREAL) :: y(12, 80000), &
                   rprop(48, 80000)
   REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: rocGrid 
   REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: xGrid, yGrid, zGrid
   REAL(RFREAL),ALLOCATABLE,DIMENSION(:) :: xData,yData,zData,rData,dumData    
   REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: volp,SPL 
   REAL(RFREAL), DIMENSION(3) :: tpw1,tpw2,tpw3         
   REAL(RFREAL) :: xin, wout, pi
   REAL(RFREAL) :: rmass

   INTEGER :: seed(33), isize, CellVertices

   INTEGER :: stationary, qs_flag, am_flag, pg_flag, &
        collisional_flag, heattransfer_flag, feedback_flag, &
        qs_fluct_flag, ppiclf_debug, rmu_flag, &
        rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag, &
        qs_fluct_filter_adapt_flag, &
        ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH, &
        sbNearest_flag, burnrate_flag, flow_model, pseudoTurb_flag
   REAL(RFREAL) :: rmu_ref, tref, suth, ksp, erest
   COMMON /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag, &
        collisional_flag, heattransfer_flag, feedback_flag, &
        qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth, &
        rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag, &
        qs_fluct_filter_adapt_flag, ksp, erest, &
        ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH, &
        sbNearest_flag, burnrate_flag, flow_model, pseudoTurb_flag
   REAL(RFREAL) :: ppiclf_rcp_part
   CHARACTER(12) :: ppiclf_matname
   COMMON /RFLU_ppiclf_misc01/ ppiclf_rcp_part
   COMMON /RFLU_ppiclf_misc02/ ppiclf_matname

   ! 08/19/24 - Thierry - added for Periodicity - begins here
   INTEGER :: x_per_flag, y_per_flag, z_per_flag, ang_per_flag 
   REAL(RFREAL) :: x_per_min, x_per_max, & 
                   y_per_min, y_per_max, & 
                   z_per_min, z_per_max, & 
                   ang_per_angle, ang_per_xangle, &
                   ang_per_rin, ang_per_rout
   ! 08/19/24 - Thierry - added for Periodicity - ends here

   REAL(RFREAL), DIMENSION(3) :: MaxPoint, MinPoint, CellLen, Max_CellLen, &
                                 filter, filter_local
   
   ! 04/04/2025 - TLJ - added min/max grid for periodicity
   INTEGER :: errorFrag
   REAL(RFREAL) :: gridmin,gridmax

   
!******************************************************************************

  RCSIdentString = '$RCSfile: PICL_TEMP_InitSolver.F90,v $ $Revision: 1.1.1.1 $'
  
  global => pRegion%global !pRegion%global

  CALL RegisterFunction( global,'PICL_TEMP_InitSolver',"../rocpicl/PICL_TEMP_InitSolver.F90" )


! Set pointers ----------------------------------------------------------------

    !levels(0)=>pLevel    
  !  pRegion => regions(1)         !pLevel%regions(iReg)
    pGrid   => pRegion%grid!pRegion%grid

!MOVING 1 HERE TO AVOID SEG FAULT OF ROCFLU STORED VF

IF (global%rkscheme /= RK_SCHEME_3_WRAY) THEN
  CALL ErrorStop(global,ERR_PICL_WRONG_RK,196,'Wrong RK Scheme for ppiclf. Needs RK3')
END IF

! Setting flags from input file
! ************************************************************
stationary = global%piclStationaryFlag
qs_flag = global%piclQsFlag
am_flag = global%piclAmFlag
pg_flag = global%piclPgFlag
collisional_flag = global%piclCollisionFlag
ViscousUnsteady_flag = global%piclViscousUnsteady
heattransfer_flag = global%piclHeatTransferFlag
feedback_flag = global%piclFeedbackFlag
qs_fluct_flag = global%piclQsFluctFlag
ppiclf_debug = global%piclDebug
sbNearest_flag = global%piclSBNearFlag
burnrate_flag = global%piclBurnRateFlag

rmu_flag = pRegion%mixtInput%viscModel
rmu_ref = pRegion%mixtInput%refVisc
tref = pRegion%mixtInput%refViscTemp
suth = pRegion%mixtInput%suthCoef
rmu_suth_param = VISC_SUTHR
rmu_fixed_param = VISC_FIXED
flow_model = int(pRegion%mixtInput%flowModel)

ksp = global%piclKsp
erest = global%piclERest

qs_fluct_filter_flag = global%piclQsFluctFilterFlag
qs_fluct_filter_adapt_flag = global%piclQsFluctFilterAdaptFlag

x_per_flag = global%piclPeriodicXFlag 
! Find min/max grid coordinates across MPI ranks
gridmin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
gridmax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
CALL MPI_AllReduce(gridmin,x_per_min,1,MPI_RFREAL,MPI_MIN, &
        global%mpiComm,errorFlag)
CALL MPI_AllReduce(gridmax,x_per_max,1,MPI_RFREAL,MPI_MAX, &
        global%mpiComm,errorFlag)

y_per_flag = global%piclPeriodicYFlag    
! Find min/max grid coordinates across MPI ranks
gridmin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
gridmax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
CALL MPI_AllReduce(gridmin,y_per_min,1,MPI_RFREAL,MPI_MIN, &
        global%mpiComm,errorFlag)
CALL MPI_AllReduce(gridmax,y_per_max,1,MPI_RFREAL,MPI_MAX, &
        global%mpiComm,errorFlag)

z_per_flag = global%piclPeriodicZFlag 
! Find min/max grid coordinates across MPI ranks
gridmin = MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))
gridmax = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))
CALL MPI_AllReduce(gridmin,z_per_min,1,MPI_RFREAL,MPI_MIN, &
        global%mpiComm,errorFlag)
CALL MPI_AllReduce(gridmax,z_per_max,1,MPI_RFREAL,MPI_MAX, &
        global%mpiComm,errorFlag)

ang_per_flag   = global%piclAngularPeriodicFlag
ang_per_angle  = global%piclAngularPeriodicAngle
ang_per_xangle = global%piclAngularPeriodicXAngle
ang_per_rin    = global%piclAngularPeriodicRin
ang_per_rout   = global%piclAngularPeriodicRout

pi = acos(-1.0D0)
ang_per_angle  = global%piclAngularPeriodicAngle * pi / 180.0D0
ang_per_xangle = global%piclAngularPeriodicXAngle * pi / 180.0D0


! Sanity check for viscosity
IF(rmu_ref .LT. 0.0d0) THEN
    CALL ErrorStop(global,ERR_PICL_INVALID_VISC,268,&
        'Negative viscosity for ppiclF')
END IF

 ! Initialization for viscous unsteady term
 ppiclf_nTimeBH = 1
 ppiclf_nUnsteadyData = 0
! ************************************************************

! Josh Gillis - Fixed restart probelm
! TLJ - we should probably use the .rin file instead
global%restartFromScratch = .true.
vtuFile1 = 'par00002.vtu'
INQUIRE(FILE=trim(vtuFile1), EXIST=fexists)
IF( global%myProcid == MASTERPROC) then
   write(*,*) fexists, "fexists"
ENDIF
IF(fexists) THEN
   global%restartFromScratch = .false.
   IF ( global%myProcid == MASTERPROC) THEN
      PRINT*, " "
      PRINT*, " ======================================="
      PRINT*, " "
      WRITE(*,*) 'Starting PPICLF Restart'
   END IF
END IF

IF(global%restartFromScratch) THEN
   IF ( global%myProcid == MASTERPROC) then
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading points.dat file...'
   END IF
  
   iFileName = 'points.dat'
   iFile = 0

   ! open data file
   OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
   global%error = errorFlag   
   IF( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,307,iFileName)
   END IF

   ! check for comments at beginning of file
   ! No comments allowed after beginning of file
   comment = '!'
   DO WHILE (comment(1:1) == '!')
      READ(iFile, '(A)') comment
      comment = ADJUSTL(TRIM(comment))
   END DO

   BACKSPACE(iFile, IOSTAT=ErrorFlag)
  
   READ(iFile,*) npart ! global number of particles
   IF (npart .gt. 80000*global%nProcs) THEN
      CALL ErrorStop(global,ERR_ILLEGAL_VALUE,322,'PPICLF:too &
        many particles to initialize')
   END IF
  
   npart_local = CEILING(REAL(npart)/(REAL(global%nProcs)))
   i = 1
   i_global = 1
   i_global_min = npart_local*global%myProcid
   i_global_max = npart_local*(global%myProcid+1)
   IF(i_global_max > npart) i_global_max = npart

   rprop(1:48,1:80000) = 0.0D0
   dp_max = 0.0D0
   xp_min =  17400000.0
   yp_min =  17400000.0
   zp_min =  17400000.0
   xp_max = -17400000.0
   yp_max = -17400000.0
   zp_max = -17400000.0
   DO i_global=1,npart
      READ(iFile,*) matName, (y(ii,i),ii=1,3),dp !points.dat not formated

      ! These are global max/min's since in 1:npart loop
      dp_max = max(dp_max,dp)
      xp_min = min(xp_min,y(1,i)-dp/2.0)
      xp_max = max(xp_max,y(1,i)+dp/2.0)
      yp_min = min(yp_min,y(2,i)-dp/2.0)
      yp_max = max(yp_max,y(2,i)+dp/2.0)
      zp_min = min(zp_min,y(3,i)-dp/2.0)
      zp_max = max(zp_max,y(3,i)+dp/2.0)
  
      ! if in range for this processor set all the other properties and increment i
      IF((i_global .GT. i_global_min) .AND. (i_global .LE. i_global_max)) THEN
         y(4,i) = 0.0D0
         y(5,i) = 0.0D0
         y(6,i) = 0.0D0
         y(7, i) = global%piclTemp
         y(8,i) = 0.0D0
         y(9,i) = 0.0D0
         y(10,i) = 0.0D0
  
         ! search for material
         matName = ADJUSTL(TRIM(matName))
         foundMat = .FALSE.
         DO iMat=1,global%nMaterials
            material => global%materials(iMat)
            IF (matName == material%name) THEN
               ppiclf_matname = matName
               IF (material%spht .GE. 1.0_RFREAL) THEN
                  ppiclf_rcp_part = material%spht
               ELSE
                 PRINT*, 'Material Specific Heat not found in input file' 
                 CALL ErrorStop(global,ERR_INRT_MISSPLAGMAT,374,matName)
               ENDIF
               rhop = material%dens
               foundMat = .TRUE.
               EXIT
            END IF
         END DO

         IF(.NOT. foundMat) THEN
            print*,global%myProcid,'stopping foundMat = False'
            CALL ErrorStop(global,ERR_INRT_MISSPLAGMAT,384,matName)
         END IF

         IF ( global%myProcid == MASTERPROC) then
            IF (i==1) THEN
               print*
               print*,'PPICLF MAT: ',trim(ppiclf_matname)
               print*,'   Density: ',rhop
               print*,'   C_p:     ',ppiclf_rcp_part
               print*,'   T_p:     ',global%piclTemp
               print*
            END IF
         END IF
    
         ! now set properties that are not interpolated from Rocflu onto the particles
         rprop(1,i) = rhop   ! particle density
         rprop(3,i)   = dp ! particle diameter
         rprop(4,i) = (4.0_RFREAL/3.0_RFREAL)*global%pi*&
                                   (0.5_RFREAL*dp)**3 ! particle volume
         ! Super Particle Loading (Real Number of particles = SPL * number of compuational particles)
         rprop(22,i) = 1.0_RFREAL

         ! Davin - added for burn rate model 02/22/2025
         rmass = rprop(4,i)*rhop
         rprop(34,i) = dp           ! Initial diameter
         rprop(35,i) = 0.0_RFREAL  ! Initial burntime
         y(11,i) = rmass            ! Initial AL mass
         y(12,i) = 0.0_RFREAL       ! Initial OX mass

         i = i + 1
      END IF
   END DO
   npart_local = i - 1

   CALL MPI_ALLREDUCE(xp_min,xp_min,1,MPI_RFREAL,MPI_MIN, &
        global%mpiComm,global%mpierr)

   CALL MPI_ALLREDUCE(xp_max,xp_max,1,MPI_RFREAL,MPI_MAX, &
        global%mpiComm,global%mpierr)

   CALL MPI_ALLREDUCE(yp_min,yp_min,1,MPI_RFREAL,MPI_MIN, &
        global%mpiComm,global%mpierr)

   CALL MPI_ALLREDUCE(yp_max,yp_max,1,MPI_RFREAL,MPI_MAX, &
        global%mpiComm,global%mpierr)

   CALL MPI_ALLREDUCE(zp_min,zp_min,1,MPI_RFREAL,MPI_MIN, &
        global%mpiComm,global%mpierr)

   CALL MPI_ALLREDUCE(zp_max,zp_max,1,MPI_RFREAL,MPI_MAX, &
        global%mpiComm,global%mpierr)

   IF(global%myProcid == MASTERPROC) THEN
      print*
      print*,'Particle domain boundaries at t=0'
      print*,'x - min, max, dx', xp_min, xp_max, xp_max - xp_min
      print*,'y - min, max, dx', yp_min, yp_max, yp_max - yp_min
      print*,'z - min, max, dx', zp_min, zp_max, zp_max - zp_min
      print*
   END IF


   IF(xp_min < x_per_min .OR. xp_max > x_per_max .OR. &
      yp_min < y_per_min .OR. yp_max > y_per_max .OR. &
      zp_min < z_per_min .OR. zp_max > z_per_max) THEN
     IF(global%myProcid == MASTERPROC) THEN
       WRITE(*,*) 'WARNING - Particles initalized outside of fluid domain'
       WRITE(*,*) 'Particle domain boundaries at t=0:'
       WRITE(*,*) 'x - min, max', xp_min, xp_max, xp_max - xp_min
       WRITE(*,*) 'y - min, max', yp_min, yp_max, yp_max - yp_min
       WRITE(*,*) 'z - min, max', zp_min, zp_max, zp_max - zp_min
       WRITE(*,*) 'Fluid domain boundaries:'
       WRITE(*,*) 'x fluid min/max', x_per_min, x_per_max
       WRITE(*,*) 'y fluid min/max', y_per_min, y_per_max
       WRITE(*,*) 'z fluid min/max', z_per_min, z_per_max
     END IF
     CALL MPI_Barrier(global%mpiComm,errorFlag)
     CALL ErrorStop(global,0,459,"rocpicl Init: Particles outside fluid domain")
   END IF

   ! Close points.dat file
   CLOSE(iFile, IOSTAT=errorFlag)
   global%error = errorFlag   
   IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,468,iFileName)
   END IF ! global%error  
ELSE
!  This is for a restart
   fexists = .true.
   ii=0 !last exists par.vtu file
   DO WHILE (fexists)
      ii = ii + 1
      vtuFile = ''
      WRITE(vtuFile,'(A3,I5.5,A4)') 'par',ii,'.vtu'
      INQUIRE(FILE=trim(vtuFile), EXIST=fexists)
      IF( global%myProcid == MASTERPROC) THEN
         WRITE(*,*) 'PPICLF par file: ',TRIM(vtuFile),ii,'  ',fexists
      END IF
   END DO
   ii = ii - 1
   vtuFile = ''
   WRITE(vtuFile,'(A3,I5.5,A4)') 'par',ii,'.vtu'
   IF( global%myProcid == MASTERPROC) THEN
      WRITE(*,*) 'Reading ', vtuFile, len(vtuFile)
   ENDIF

   IF(ii .lt. 0) THEN
      CALL ErrorStop(global,ERR_FILE_EXIST,491,vtuFile)
   END IF

   npart = -1
   dp_max = -1.0
   CALL ppiclf_io_ReadParticleVTU(trim(vtuFile), ii, npart, dp_max)
   CALL MPI_Allreduce(dp_max,dp_max,1,MPI_RFREAL,MPI_MAX, &
      global%mpiComm,global%mpierr )
   print*,global%myProcid,npart,dp_max,dp_max
   IF( global%myProcid == MASTERPROC) THEN
      print*, " "
      WRITE(*,*) 'Finished PPICLF Restart'
      print*, " "
      print*, " ======================================="
      print*, " "
   END IF
END IF ! global%restartFromScratch

! User sets up overlap grid:
nCells = pRegion%grid%nCells
 
ALLOCATE(xGrid(8,nCells),STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_ALLOCATE,515,'PPICLF:xGrid')
END IF ! global%error

ALLOCATE(yGrid(8,nCells),STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_ALLOCATE,521,'PPICLF:yGrid')
END IF ! global%error

ALLOCATE(zGrid(8,nCells),STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_ALLOCATE,527,'PPICLF:zGrid')
END IF ! global%error

ALLOCATE(rocGrid(7,nCells),STAT=errorFlag)
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_ALLOCATE,532,'PPICLF:rocGrid')
END IF ! global%error



DO i = 1, nCells 
  IF(pGrid%cellGlob2Loc(1,i) == 1) THEN
    ! Tetrahedral Cell
    CellVertices = 4
    DO k = 1, CellVertices
      vi = pRegion%grid%tet2v(k,i) 
             xGrid(k,i) = pRegion%grid%xyz(XCOORD,vi) 
             yGrid(k,i) = pRegion%grid%xyz(YCOORD,vi) 
             zGrid(k,i) = pRegion%grid%xyz(ZCOORD,vi) 
    END DO
  ELSE IF(pGrid%cellGlob2Loc(1,i) == 2) THEN
    ! Hexahedral Cell
    CellVertices = 8
    DO k = 1, CellVertices
      vi = pRegion%grid%hex2v(k,i) 
             xGrid(k,i) = pRegion%grid%xyz(XCOORD,vi) 
             yGrid(k,i) = pRegion%grid%xyz(YCOORD,vi) 
             zGrid(k,i) = pRegion%grid%xyz(ZCOORD,vi) 
    END DO !nCells
  ELSE
    CellVertices = 0
    WRITE(*,*) 'ERROR: Rocflupicl only support tetrahedral and hexahedral cell types.'
    CALL ErrorStop(global,ERR_ALLOCATE,559,'PPICLF:CellLen')
  END IF
END DO
! Find cell lengths
DO i = 1,nCells
  IF(pGrid%cellGlob2Loc(1,i) == 1) THEN
    ! Tetrahedral Cell
    CellVertices = 4
  ELSE IF(pGrid%cellGlob2Loc(1,i) == 2) THEN
    ! Hexahedral Cell
    CellVertices = 8
  ELSE
    CellVertices = 0
    WRITE(*,*) 'ERROR: Rocflupicl only support tetrahedral and hexahedral cell types.'
    CALL ErrorStop(global,ERR_ALLOCATE,573,'PPICLF:CellLen')
  END IF
  ! Initialize as zero for each element
  DO l = 1,3
    MaxPoint(l) = -1.0D10 
    MinPoint(l) =  1.0D10 
    CellLen(l)   =  0.0D0   
    Max_CellLen(l)  = 0.0D0
  END DO !l
  ! Add all x,y,z cell corners for centroid and find extremes
  DO k = 1,CellVertices
    IF (xGrid(k,i) > MaxPoint(1)) &
      MaxPoint(1) = xGrid(k,i)
    IF (xGrid(k,i) < MinPoint(1)) &
      MinPoint(1) = xGrid(k,i)
    IF (yGrid(k,i) > MaxPoint(2)) &
      MaxPoint(2) = yGrid(k,i)  
    IF (yGrid(k,i) < MinPoint(2)) &
      MinPoint(2) = yGrid(k,i)
    IF (zGrid(k,i) > MaxPoint(3)) &
      MaxPoint(3) = zGrid(k,i)  
    IF (zGrid(k,i) < MinPoint(3)) &
      MinPoint(3) = zGrid(k,i)
  END DO !k
  DO l = 1,3
    ! Find element length in all dimensions
    CellLen(l) = ABS(MaxPoint(l)-MinPoint(l))
    ! Find max lengths for grid cells in particle domain on this processor
    IF( (pRegion%grid%cofg(1,i) >= xp_min) .AND. (pRegion%grid%cofg(1,i) <= xp_max) .AND. & 
        (pRegion%grid%cofg(2,i) >= yp_min) .AND. (pRegion%grid%cofg(2,i) <= yp_max) .AND. & 
        (pRegion%grid%cofg(3,i) >= zp_min) .AND. (pRegion%grid%cofg(3,i) <= zp_max)) THEN 
      IF(CellLen(l) .GT. Max_CellLen(l)) Max_CellLen(l) = CellLen(l)
      ! For tets, make a rectangular box around it with maxCellLength equal in all dimensions
      IF(pGrid%cellGlob2Loc(1,i) == 1) Max_CellLen(l) = MAXVAL(Max_CellLen(:))
    END IF
  END DO !l
  IF(pGrid%cellGlob2Loc(1,i) == 1) THEN
    ! Make cell length same in all dimensions for tetrahedral
    ! This is due to their shape not being rectangular
    DO l = 1,3
      CellLen(l) = MAX(CellLen(1),CellLen(2),CellLen(3))
    END DO !l
  END IF

  rocGrid(1,i) = pRegion%grid%cofg(1,i) !Cell Centroid x position
  rocGrid(2,i) = pRegion%grid%cofg(2,i) !Cell Centroid y position
  rocGrid(3,i) = pRegion%grid%cofg(3,i) !Cell Centroid z position
  rocGrid(4,i) = CellLen(1)             !Cell largest dx
  rocGrid(5,i) = CellLen(2)             !Cell largest dy
  rocGrid(6,i) = CellLen(3)             !Cell largest dz
  rocGrid(7,i) = pRegion%grid%vol(i)    !Cell volume
END DO !i

MinFluidCells = 2.0D0 !Number of fluid cells for Minimum Bin Size and ppiclf_filter(1:3)
DO l = 1,3
  filter_local(l) = MinFluidCells*Max_CellLen(l) 
  ! Find max x,y,z cell lengths across all MPI ranks (entire rocflu domain)
  CALL MPI_Allreduce(filter_local(l),filter(l),1,MPI_RFREAL,MPI_MAX, &
      global%mpiComm,global%mpierr )
END DO

neighborWidth = 4.0_RFREAL*dp_max ! Minimum value. User can entire larger if desired
IF((neighborWidth .GT. global%piclNeighborWidth) &
    .AND. (global%myProcid == MASTERPROC)) THEN
    WRITE(STDOUT, '(A)') &
        '*** WARNING *** PICL NEIGHBORWIDTH too small, defaulting to 4*dp_max'
END IF
neighborWidth = MAX(neighborWidth, global%piclNeighborWidth)

IF(global%myProcid == MASTERPROC) THEN
   PRINT*,' '
   PRINT*,'PPICLF: '
   PRINT*,'  Inputs NEIGHBORWIDTH  : ',global%piclNeighborWidth
   PRINT*,'  dp_max (points.dat)   : ',dp_max
   PRINT*,'  ppiclf_filter(1)      : ',filter(1)
   PRINT*,'  ppiclf_filter(2)      : ',filter(2)
   PRINT*,'  ppiclf_filter(3)      : ',filter(3)
   PRINT*,'  ppiclf_nndist         : ',neighborWidth
   PRINT*,' '
END IF

! Sets ppiclf_nndist, ppiclf_filter(1:3), ppiclf_y(), ppiclf_rprop()
!(imethod (RK pick), nDimensions, iendian (IO format), .....)
CALL ppiclf_solve_InitParticle(2,3,0,npart_local,y,rprop,filter,neighborWidth) 
IF(global%myProcid == MASTERPROC) THEN
  PRINT*, 'x fluid min/max', x_per_min, x_per_max
  PRINT*, 'y fluid min/max', y_per_min, y_per_max
  PRINT*, 'z fluid min/max', z_per_min, z_per_max
END IF

CALL ppiclf_solve_Initialize( & 
           x_per_flag, x_per_min, x_per_max, &
           y_per_flag, y_per_min, y_per_max, &
           z_per_flag, z_per_min, z_per_max, &
           ang_per_flag, ang_per_angle, ang_per_xangle, ang_per_rin, ang_per_rout)

! *** AVERY - Update when general ang per is done ***
IF(((ang_per_flag.eq.1) .and. (x_per_flag.eq.1 .or. y_per_flag.eq.1)) .or. & 
    (ang_per_flag .gt. 1)) THEN
    CALL ErrorStop(global,ERR_PICL_INVALID_PERIODICITY,672,&
      'Wrong periodicity choices for ppiclF')
END IF 

! Creates OverlapGrid and Calls Init Solve
CALL ppiclf_comm_InitOverlapGrid(nCells,rocGrid)
CALL ppiclf_solve_InitSolve

INQUIRE(FILE='filein.vtk', EXIST=wall_exists)
IF(wall_exists) THEN
  CALL ppiclf_io_ReadWallVTK('filein.vtk')
  WRITE(*,*) 'Particle boundary wall read-in.'
ELSE IF (global%myProcid == MASTERPROC) THEN
  WRITE(*,*) '***WARNING: Could not find filein.vtk***'
END IF

! 03/24/2025 - Thierry - store the RocfluMP Flow Model chosen (Euler or NS)
!                        this is used in ppiclF for calculating the pressure gradient
!                        whether with or without the viscous part
!flow_model = 0
!flow_model = int(pRegion%mixtInput%flowModel)

DEALLOCATE(rocGrid,STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_DEALLOCATE,697,'PPICLF:rocGrid')
END IF ! global%error


DEALLOCATE(xGrid,STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_DEALLOCATE,704,'PPICLF:xGrid')
END IF ! global%error

DEALLOCATE(yGrid,STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_DEALLOCATE,710,'PPICLF:yGrid')
END IF ! global%error

DEALLOCATE(zGrid,STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_DEALLOCATE,716,'PPICLF:zGrid')
END IF ! global%error

ALLOCATE(volp(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,722,'PPICLF:volp')
    END IF ! global%error

DO i=1,pGrid%nCellsTot
        pRegion%mixt%piclVF(i) = 0.0_RFREAL
END DO

! Calling inisolve to gather particle volume fraction
CALL ppiclf_solve_initsolve

IF ( global%myProcid == MASTERPROC) WRITE(*,*) "PFINIT: Calc Init VolP"
DO i = 1, nCells
       CALL ppiclf_solve_GetProFld(i, 1, volp(i))
       IF (pRegion%mixtInput%axiFlag) THEN
           WRITE(*,*) "Need to properly implement axi-sym for phip init."
           CALL ErrorStop(global,ERR_OPTION_TYPE,737,'PPICLF:axi')
       END IF
       volp(i) = volp(i)/pRegion%grid%vol(i)
!*** VOL FRAC CAP
!*** I THINK WE SHOULD DELETE THIS CAP - AVERY ***
       IF(volp(i) .GT. 0.62) THEN
           volp(i) = 0.62
       END IF
       pRegion%mixt%piclVF(i) = volp(i) 
END DO

! TLJ:
! This section takes as input from utilities/init/RFLU_InitFlowHardCode.F90
!    (r,r*u,r*v,r*w,r*E) and changes to RocfluMP conserved variables
!    (phig*r,phig*r*u,phig*r*v,phig*r*w,phig*r*E), where phig is the gas
!    phase volume fraction and can only be computed after the particles
!    are read in.

DO icg = 1,pGrid%nCellsTot
    vFrac = 1.0_RFREAL - pRegion%mixt%piclVF(icg)
    pRegion%mixt%cv(CV_MIXT_DENS,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_DENS,icg)
    pRegion%mixt%cv(CV_MIXT_XMOM,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
    pRegion%mixt%cv(CV_MIXT_YMOM,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_YMOM,icg)
    pRegion%mixt%cv(CV_MIXT_ZMOM,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_ZMOM,icg)
    pRegion%mixt%cv(CV_MIXT_ENER,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_ENER,icg)
    IF(pRegion%mixt%cv(CV_MIXT_DENS,icg) .le. 0.0) THEN
         WRITE(*,*) "Error: negative density: ",pRegion%mixt%cv(CV_MIXT_DENS,icg)     
         PRINT*, 'From rocpicl/PICL_TEMP_InitSolver.F90' 
         CALL ErrorStop(global,ERR_INVALID_VALUE,765,'PPICLF:init')
    END IF    
END DO ! icg

DEALLOCATE(volp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,772,'PPICLF:zGrid')
    END IF ! global%error


IF ( global%myProcid == MASTERPROC) then
   print*, ' '
   print*, '***********************************************'
   print*, 'TLJ'
   print*, 'Starting PICL_TEMP_InitSolver.F90'
   print*, ' '
   print*, 'stationary           = ',global%piclStationaryFlag
   print*, 'qs_flag              = ',global%piclQsFlag
   print*, 'am_flag              = ',global%piclAmFlag
   print*, 'pg_flag              = ',global%piclPgFlag
   print*, 'collisional_flag     = ',global%piclCollisionFlag
   print*, 'ViscousUnsteady_flag = ',global%piclViscousUnsteady
   print*, 'heattransfer_flag    = ',global%piclHeatTransferFlag
   print*, 'feedback_flag        = ',global%piclFeedbackFlag
   print*, 'qs_fluct_flag        = ',global%piclQsFluctFlag
   print*, 'ppiclf_debug         = ',global%piclDebug
   print*, 'ppiclf_nUnsteadyData = ',ppiclf_nUnsteadyData
   print*, 'ppiclf_VU            = ',0
   print*, 'sbNearest_flag       = ',global%piclSBNearFlag
   print*, 'burnrate_flag        = ',global%piclBurnRateFlag

   IF (global%piclViscousUnsteady >=1) THEN
      print*,'  Using Viscous unsteady history term'
      print*,'    ppiclf_nTimeBH       = ',ppiclf_nTimeBH
      print*,'    ppiclf_nUnsteadyData = ',ppiclf_nUnsteadyData,0
   ENDIF

   print*, ' '
   print*, "XLOC = ", y(1 ,1) 
   print*, "YLOC = ", y(2 ,1) 
   print*, "ZLOC = ", y(3 ,1)   

   print*, ' '
   print*, 'Reading points.dat file...'
   print*, "part material  = ", TRIM(ppiclf_matname)
   print*, "npart          = ", npart
   print*, 'dp_max         = ', dp_max
   print*, 'rho_dens       = ', rhop
   print*, 'cv_particle    = ', ppiclf_rcp_part
   
   if(x_per_flag.eq.1) then
     print*, "ppiclF X-periodic (min, max) ", x_per_min, x_per_max
   endif  
   if(y_per_flag.eq.1) then
     print*, "ppiclF Y-periodic (min, max) ", y_per_min, y_per_max
   endif  
   if(z_per_flag.eq.1) then
     print*, "ppiclF Z-periodic (min, max) ", z_per_min, z_per_max
   endif  
   if(ang_per_flag.eq.1) then
     print*, "ppiclF Angular-periodic (axis, rin, rout, angle, x-angle) ", &
     ang_per_flag, ang_per_rin, ang_per_rout, &
     ang_per_angle, ang_per_xangle
   endif  

   print*, ' '
   print*, 'Ending PICL_TEMP_InitSolver.F90'
   print*, '***********************************************'
ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PICL_TEMP_InitSolver

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PICL_.F90,v $
!
!
!******************************************************************************


