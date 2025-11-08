










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
  USE ModRandom, ONLY: Rand1Uniform,Rand1Normal
  USE RFLU_ModInCellTest

  IMPLICIT NONE
!DEC$ NOFREEFORM

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
   LOGICAL :: notfoundFlag, pf_restart,pf_rpInit,pf_settle,pf_fluidInit,&
              wall_exists, fexists, foundMat
   INTEGER :: i,npart,nCells,lx,ly,lz,vi,vii,ii,jj,kk,loopCounter,ipart,icl
   INTEGER :: PPC,numPclCells,npart_local,i_global,i_global_min,i_global_max,&
              iFile,iMat, k, j, l, m
   REAL(RFREAL) :: dp_min,dp_max,rhop,tester,ratio,total_vol,filter,xMinCurt,&
                   xMaxCurt,yMinCurt,yMaxCurt,xMinCell,xMaxCell,yMinCell,&
                   yMaxCell,zMinCell,zMaxCell,x,vFrac,volpclsum,xLoc,yLoc,zLoc,yL, &
                   zpf_factor,xpf_factor,dp,neighborWidth,dp_max_l,xp_min,xp_max, &
                   xp_min_l,xp_max_l
   REAL(KIND=8) :: y(12, 20000), &
                   rprop(53, 20000)
   REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: xGrid, yGrid, zGrid,vfP
   REAL(RFREAL),ALLOCATABLE,DIMENSION(:) :: xData,yData,zData,rData,dumData     
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: volp,SPL 
   REAL(KIND=8), DIMENSION(3) :: tpw1,tpw2,tpw3         
   REAL(RFREAL) :: xin, wout, pi
   REAL(RFREAL) :: rmass

   INTEGER :: seed(33), isize

   integer*4 :: stationary, qs_flag, am_flag, pg_flag, &
        collisional_flag, heattransfer_flag, feedback_flag, &
        qs_fluct_flag, ppiclf_debug, rmu_flag, &
        rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag, &
        qs_fluct_filter_adapt_flag, &
        ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH, &
        sbNearest_flag, burnrate_flag, flow_model, pseudoTurb_flag, &
        HTUnsteady_flag
   real*8 :: rmu_ref, tref, suth, ksp, erest
   common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag, &
        collisional_flag, heattransfer_flag, feedback_flag, &
        qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth, &
        rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag, &
        qs_fluct_filter_adapt_flag, ksp, erest, &
        ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH, &
        sbNearest_flag, burnrate_flag, flow_model, pseudoTurb_flag, &
        HTUnsteady_flag
   real*8 :: ppiclf_rcp_part
   CHARACTER(12) :: ppiclf_matname
   common /RFLU_ppiclf_misc01/ ppiclf_rcp_part
   common /RFLU_ppiclf_misc02/ ppiclf_matname

   ! 08/19/24 - Thierry - added for Periodicity - begins here
   integer*4 x_per_flag, y_per_flag, z_per_flag, ang_per_flag 
   real*8 x_per_min, x_per_max, & 
          y_per_min, y_per_max, & 
          z_per_min, z_per_max, & 
          ang_per_angle, ang_per_xangle, &
          ang_per_rin, ang_per_rout
   ! 08/19/24 - Thierry - added for Periodicity - ends here

   REAL :: MaxPoint(3), MinPoint(3), EleLen(3), Max_EleLen(3)
   
   ! 04/04/2025 - TLJ - added min/max grid for periodicity
   integer :: errorFrag
   real*8 gridmin,gridmax

   
!******************************************************************************

  RCSIdentString = '$RCSfile: PICL_TEMP_InitSolver.F90,v $ $Revision: 1.1.1.1 $'
  
  global => pRegion%global !pRegion%global

!write(*,*) "Step 0:",global%myProcid  
  CALL RegisterFunction( global,'PICL_TEMP_InitSolver',"../rocpicl/PICL_TEMP_InitSolver.F90" )


! Set pointers ----------------------------------------------------------------

    !levels(0)=>pLevel    
  !  pRegion => regions(1)         !pLevel%regions(iReg)
    pGrid   => pRegion%grid!pRegion%grid

  CALL MPI_Barrier(global%mpiComm,errorFlag)


!MOVING 1 HERE TO AVOID SEG FAULT OF ROCFLU STORED VF

IF (global%rkscheme /= RK_SCHEME_3_WRAY) THEN
  CALL ErrorStop(global,ERR_PICL_WRONG_RK,199,'Wrong RK Scheme for ppiclf. Needs RK3')
END IF

stationary = global%piclStationaryFlag
qs_flag = global%piclQsFlag
am_flag = global%piclAmFlag
pg_flag = global%piclPgFlag
collisional_flag = global%piclCollisionFlag
ViscousUnsteady_flag = global%piclViscousUnsteady
heattransfer_flag = global%piclHeatTransferFlag
HTUnsteady_flag = global%piclHTUnsteadyFlag
feedback_flag = global%piclFeedbackFlag
qs_fluct_flag = global%piclQsFluctFlag
pseudoTurb_flag = global%piclPseudoTurbFlag
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

! 08/13/24 - Thierry - added for Periodicity - begins here
! 04/04/2025 - TLJ - modified to detemine min/max from Rocflu grid

x_per_min=-10000.0; x_per_max= 10000.0; ! set crazy value if not used
y_per_min=-10000.0; y_per_max= 10000.0; ! set crazy value if not used
z_per_min=-10000.0; z_per_max= 10000.0; ! set crazy value if not used

x_per_flag = global%piclPeriodicXFlag 
if (x_per_flag == 1) then
    gridmin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
    gridmax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
    CALL MPI_AllReduce(gridmin,x_per_min,1,MPI_RFREAL,MPI_MIN, &
            global%mpiComm,errorFlag)
    CALL MPI_AllReduce(gridmax,x_per_max,1,MPI_RFREAL,MPI_MAX, &
            global%mpiComm,errorFlag)
endif

y_per_flag = global%piclPeriodicYFlag    
if (y_per_flag == 1) then
    gridmin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
    gridmax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
    CALL MPI_AllReduce(gridmin,y_per_min,1,MPI_RFREAL,MPI_MIN, &
            global%mpiComm,errorFlag)
    CALL MPI_AllReduce(gridmax,y_per_max,1,MPI_RFREAL,MPI_MAX, &
            global%mpiComm,errorFlag)
endif

z_per_flag = global%piclPeriodicZFlag 
if (z_per_flag == 1) then
    gridmin = MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))
    gridmax = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))
    CALL MPI_AllReduce(gridmin,z_per_min,1,MPI_RFREAL,MPI_MIN, &
            global%mpiComm,errorFlag)
    CALL MPI_AllReduce(gridmax,z_per_max,1,MPI_RFREAL,MPI_MAX, &
            global%mpiComm,errorFlag)
endif

ang_per_flag   = global%piclAngularPeriodicFlag
ang_per_angle  = global%piclAngularPeriodicAngle
ang_per_xangle = global%piclAngularPeriodicXAngle
ang_per_rin    = global%piclAngularPeriodicRin
ang_per_rout   = global%piclAngularPeriodicRout

pi = acos(-1.0d0)
ang_per_angle  = global%piclAngularPeriodicAngle * pi / 180.0d0
ang_per_xangle = global%piclAngularPeriodicXAngle * pi / 180.0d0

call ppiclf_solve_Initialize( &
   x_per_flag, x_per_min, x_per_max, &
   y_per_flag, y_per_min, y_per_max, &
   z_per_flag, z_per_min, z_per_max, &
   ang_per_flag, ang_per_angle, ang_per_xangle, ang_per_rin, ang_per_rout)

! 08/13/24 - Thierry - added for Periodicity - ends here

! Sanity check for viscosity
if (rmu_ref .lt. 0.0d0) then
    CALL ErrorStop(global,ERR_PICL_INVALID_VISC,288,&
        'Negative viscosity for ppiclF')
end if

 ! Initialization for viscous unsteady term
 ppiclf_nTimeBH = 1
 ppiclf_nUnsteadyData = 50

 CALL MPI_Barrier(global%mpiComm,errorFlag)

seed = 1
call RANDOM_SEED(put=seed)
call RANDOM_SEED(size=isize)


! Sam - Making this as simple as possible for the initial commit. 
! We can later add back in some of these routines if we decide we
! need them. Cleaning house like Calvin, not Luther. 

pf_fluidInit = .true.

! Josh Gillis - Fixed restart probelm
! TLJ - we should probably use the .rin file instead
global%restartFromScratch = .true.
vtuFile1 = 'par00002.vtu'
INQUIRE(FILE=trim(vtuFile1), EXIST=fexists)
IF ( global%myProcid == MASTERPROC) then
   write(*,*) fexists, "fexists"
ENDIF
IF (fexists) THEN
   global%restartFromScratch = .false.
   IF ( global%myProcid == MASTERPROC) THEN
      print*, " "
      print*, " ======================================="
      print*, " "
      WRITE(*,*) 'Starting PPICLF Restart'
   ENDIF
ENDIF

IF (global%restartFromScratch) THEN
   ! This variable is stupid
   pRegion%mixt%piclGeom = 1.00_RFREAL
  
   ! Sam - initialization reading .dat file with points set all other
   ! properties manually except location
   IF ( global%myProcid == MASTERPROC) then
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading points.dat file...'
   END IF
  
   CALL MPI_Barrier(global%mpiComm,errorFlag)
   iFileName = 'points.dat'
   iFile = 0

   ! open data file
   OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
   global%error = errorFlag   
   IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,345,iFileName)
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
   IF (npart .gt. 20000*global%nProcs) THEN
      CALL ErrorStop(global,ERR_ILLEGAL_VALUE,360,'PPICLF:too &
        many particles to initialize')
   END IF
  
   npart_local = npart/global%nProcs+1
   i = 1
   i_global = 1
   i_global_min = npart_local*global%myProcid
   i_global_max = npart_local*(global%myProcid+1)
   IF(i_global_max > npart) i_global_max = npart
   !print*,global%myProcid,npart,i_global_min,i_global_max

   !IF ( global%myProcid == MASTERPROC) then
    !  print*,global%myProcid,npart,i_global_min,i_global_max
   !ENDIF

   rprop(1:53,1:20000) = 0.0d0
  
   dp_max = 0.0d0
   xp_min_l = +17400000.0
   xp_max_l = -17400000.0
   do i_global=1,npart
      READ(iFile,*) matName, (y(ii,i),ii=1,3),dp !points.dat not formated

      dp_max = max(dp_max,dp)
      xp_min_l = min(xp_min_l,y(1,i)-dp/2.0)!sqrt(y(1,i)*y(1,i)+y(2,i)*y(2,i)))
      xp_max_l = max(xp_max_l,y(1,i)+dp/2.0)!sqrt(y(1,i)*y(1,i)+y(2,i)*y(2,i)))
  
      ! if in range for this processor set all the other properties and increment i
      if ((i_global .gt. i_global_min) .and. (i_global .le. i_global_max)) then
         y(4,i) = 0.0d0
         y(5,i) = 0.0d0
         y(6,i) = 0.0d0
         y(7, i) = global%piclTemp
         y(8,i) = 0.0d0
         y(9,i) = 0.0d0
         y(10,i) = 0.0d0
  
         ! initially zero out all properties
         !do ii=1,53
         !  rprop(ii, i) = 0.0d0
         !end do

         ! search for material
         matName = ADJUSTL(TRIM(matName))
         foundMat = .FALSE.
         DO iMat=1,global%nMaterials
            material => global%materials(iMat)
            IF (matName == material%name) THEN
               ! TLJ - Set heat capacity of particle [J/kg-K]
               !     - If not set in *.inp file, set default
               ! Soda lime: Cv = 840; rhop = 2520.
               ppiclf_matname = matName
               IF (material%spht .GE. 10.0_RFREAL) THEN
                  ppiclf_rcp_part = material%spht
               ELSE
                  ppiclf_rcp_part = 840.0
               ENDIF
               rhop = material%dens
               foundMat = .TRUE.
               EXIT
            END IF
         END DO

         IF (.NOT. foundMat) THEN
            print*,global%myProcid,'stopping foundMat = False'
            CALL ErrorStop(global,ERR_INRT_MISSPLAGMAT,426,matName)
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
    
         rprop(22,i) = 1.0_RFREAL

         ! Davin - added for burn rate model 02/22/2025
         rmass = rprop(4,i)*rhop
         rprop(34,i) = dp           ! Initial diameter
         rprop(35,i) = 0.0_RFREAL  ! Initial burntime
         y(11,i) = rmass            ! Initial AL mass
         y(12,i) = 0.0_RFREAL       ! Initial OX mass

         i = i + 1
      endif
   enddo
   npart_local = i - 1
   CALL MPI_Allreduce(xp_min_l,xp_min,1,MPI_RFREAL,MPI_MIN, &
      global%mpiComm,global%mpierr )
   CALL MPI_Allreduce(xp_max_l,xp_max,1,MPI_RFREAL,MPI_MAX, &
      global%mpiComm,global%mpierr )

   IF ( global%myProcid == MASTERPROC) THEN
      print*
      print*,'TLJ starting location of particle bed (xp_min) = ',xp_min
      print*,'TLJ ending   location of particle bed (xp_max) = ',xp_max
      print*,'TLJ particle bed thickness (xp_max-xp_min)*1e3 = ',(xp_max-xp_min)*1e3
      print*
   ENDIF

   ! Close points.dat file
   CLOSE(iFile, IOSTAT=errorFlag)
   CALL MPI_Barrier(global%mpiComm,errorFlag)
   global%error = errorFlag   
   IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,477,iFileName)
   END IF ! global%error  

ELSE
!  This is for a restart
   pf_fluidInit = .false.
   fexists = .true.
   ii=0
   DO WHILE (fexists)
      ii = ii + 1
      vtuFile = ''
      write(vtuFile,'(A3,I5.5,A4)') 'par',ii,'.vtu'
      INQUIRE(FILE=trim(vtuFile), EXIST=fexists)
      IF ( global%myProcid == MASTERPROC) THEN
         WRITE(*,*) 'PPICLF par file: ',TRIM(vtuFile),ii,'  ',fexists
      ENDIF
   END DO
   ii = ii - 1
   vtuFile = ''
   write(vtuFile,'(A3,I5.5,A4)') 'par',ii,'.vtu'
   IF ( global%myProcid == MASTERPROC) THEN
      WRITE(*,*) 'Reading ', vtuFile, len(vtuFile)
   ENDIF

   IF (ii .lt. 0) THEN
      CALL ErrorStop(global,ERR_FILE_EXIST,502,vtuFile)
   END IF

   ! TLJ - 11/23/2024
   !  Fixed error in call to ReadParticleVTU:
   !    CALL ppiclf_io_ReadParticleVTU(trim(vtuFile), ii-1)
   !  Now passing back npart and dp_max
   !  Note that dp_max is needed for proper setting of neighborWidth
   npart = -1
   dp_max = -1.0
   CALL ppiclf_io_ReadParticleVTU(trim(vtuFile), ii, npart, dp_max_l)
   CALL MPI_Allreduce(dp_max_l,dp_max,1,MPI_RFREAL,MPI_MAX, &
      global%mpiComm,global%mpierr )
   print*,global%myProcid,npart,dp_max_l,dp_max

   i = 1 ! particles/rank to distribute??

   IF ( global%myProcid == MASTERPROC) THEN
      print*, " "
      WRITE(*,*) 'Finished PPICLF Restart'
      print*, " "
      print*, " ======================================="
      print*, " "
   ENDIF
END IF ! global%restartFromScratch


CALL MPI_Barrier(global%mpiComm,errorFlag)

!BRAD STARTS HERE
!Taking what was already done and building 1-way coupled first
!NEED TO BUILD overlap mesh 
!Loop through cells
 !Pull face id's
  !grab face center coords 
   !dump into approperate array
!How cooord is grab needs to match how props are grabbed

! User sets up overlap mesh:
nCells = pRegion%grid%nCells
lx = 2
ly = 2
lz = 2
 
ALLOCATE(xGrid(lx,ly,lz,nCells),STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_ALLOCATE,549,'PPICLF:xGrid')
END IF ! global%error

ALLOCATE(yGrid(lx,ly,lz,nCells),STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_ALLOCATE,555,'PPICLF:yGrid')
END IF ! global%error

ALLOCATE(zGrid(lx,ly,lz,nCells),STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_ALLOCATE,561,'PPICLF:zGrid')
END IF ! global%error

!Loop cells
DO i = 1, nCells 

   vi = pRegion%grid%hex2v(1,i) 
            xGrid(1,1,1,i) = pRegion%grid%xyz(XCOORD,vi) 
            yGrid(1,1,1,i) = pRegion%grid%xyz(YCOORD,vi) 
            zGrid(1,1,1,i) = pRegion%grid%xyz(ZCOORD,vi) 
   vi = pRegion%grid%hex2v(4,i) 
            xGrid(2,1,1,i) = pRegion%grid%xyz(XCOORD,vi)  
            yGrid(2,1,1,i) = pRegion%grid%xyz(YCOORD,vi) 
            zGrid(2,1,1,i) = pRegion%grid%xyz(ZCOORD,vi) 
   vi = pRegion%grid%hex2v(5,i) 
            xGrid(1,2,1,i) = pRegion%grid%xyz(XCOORD,vi) 
            yGrid(1,2,1,i) = pRegion%grid%xyz(YCOORD,vi) 
            zGrid(1,2,1,i) = pRegion%grid%xyz(ZCOORD,vi) 
   vi = pRegion%grid%hex2v(8,i) 
            xGrid(2,2,1,i) = pRegion%grid%xyz(XCOORD,vi) 
            yGrid(2,2,1,i) = pRegion%grid%xyz(YCOORD,vi) 
            zGrid(2,2,1,i) = pRegion%grid%xyz(ZCOORD,vi) 
   vi = pRegion%grid%hex2v(2,i) 
            xGrid(1,1,2,i) = pRegion%grid%xyz(XCOORD,vi) 
            yGrid(1,1,2,i) = pRegion%grid%xyz(YCOORD,vi)
            zGrid(1,1,2,i) = pRegion%grid%xyz(ZCOORD,vi)
   vi = pRegion%grid%hex2v(3,i) 
            xGrid(2,1,2,i) = pRegion%grid%xyz(XCOORD,vi) 
            yGrid(2,1,2,i) = pRegion%grid%xyz(YCOORD,vi) 
            zGrid(2,1,2,i) = pRegion%grid%xyz(ZCOORD,vi) 
   vi = pRegion%grid%hex2v(6,i) 
            xGrid(1,2,2,i) = pRegion%grid%xyz(XCOORD,vi) 
            yGrid(1,2,2,i) = pRegion%grid%xyz(YCOORD,vi) 
            zGrid(1,2,2,i) = pRegion%grid%xyz(ZCOORD,vi) 
   vi = pRegion%grid%hex2v(7,i) 
            xGrid(2,2,2,i) = pRegion%grid%xyz(XCOORD,vi) 
            yGrid(2,2,2,i) = pRegion%grid%xyz(YCOORD,vi)
            zGrid(2,2,2,i) = pRegion%grid%xyz(ZCOORD,vi) 

END DO !nCells
! Max cell lengths
DO i = 1,nCells
  ! Initialize as zero for each element
  DO l = 1,3
    MaxPoint(l) = -1.0E10 
    MinPoint(l) =  1.0E10 
    EleLen(l)   =  0.0   
    Max_EleLen(l)  =  1.23456789D-10
  ENDDO !l
  ! Add all x,y,z mesh points for centroid and find extremes
  ! AVERY FUTURE WORK - update filter to be 3 dimensional.
    DO k = 1,2!2
      DO j = 1,2!2
        DO m = 1,2!2
          IF (xGrid(k,j,m,i) > MaxPoint(1)) &
            MaxPoint(1) = xGrid(k,j,m,i)
          IF (xGrid(k,j,m,i) < MinPoint(1)) &
            MinPoint(1) = xGrid(k,j,m,i)
          IF (yGrid(k,j,m,i) > MaxPoint(2)) &
            MaxPoint(2) = yGrid(k,j,m,i)  
          IF (yGrid(k,j,m,i) < MinPoint(2)) &
            MinPoint(2) = yGrid(k,j,m,i)
          IF (zGrid(k,j,m,i) > MaxPoint(3)) &
            MaxPoint(3) = zGrid(k,j,m,i)  
          IF (zGrid(k,j,m,i) < MinPoint(3)) &
            MinPoint(3) = zGrid(k,j,m,i)
        ENDDO !i
      ENDDO !j
    ENDDO !k
  DO l = 1,3
    ! Find max element length in all dimensions
    EleLen(l) = ABS(MaxPoint(l)-MinPoint(l))
    ! Find max lengths for all mesh elements in all directions
    IF (EleLen(l) > 1D-2 .OR. EleLen(l) < 1D-7) THEN
      !WRITE(*,*) 'AVERY - Extreme points:', MaxPoint(l), MinPoint(l), 'Dimension:',l
      CYCLE
    END IF
    IF (EleLen(l) .GT. Max_EleLen(l)) Max_EleLen(l) = EleLen(l)
  ENDDO !l
ENDDO !i

! Old filter
filter = global%piclFilterWidth
filter=filter/2.0d0
if(global%myProcId .eq. MASTERPROC) print*, "Old filter =", filter

! Avery's modified filter
filter = 4*MAXVAL(Max_EleLen) ! Minimum two cells per ppiclf_bin
filter=filter/2.0d0
if(global%myProcId .eq. MASTERPROC) print*, "Avery filter =", filter

! TLJ compute ppcilf_d2chk here in rocpicl
!     this is needed to have the bin at t=0 be
!     the correct size
!
! TLJ: Here we need dp_min to be defined
! Sam - switching to dp_max which is the worst case
neighborWidth = 4.0_RFREAL*dp_max
if ((neighborWidth .gt. global%piclNeighborWidth) &
    .and. (global%myProcid == MASTERPROC)) then
    WRITE(STDOUT, '(A)') &
        '*** WARNING *** PICL NEIGHBORWIDTH too small, defaulting to 4*dp_max'
end if
neighborWidth = MAX(neighborWidth, global%piclNeighborWidth)

! Thierry - checking is this is the issue for the quarter cylinder shock tube not working


if (global%myProcid == MASTERPROC)  then
   print*,' '
   print*,'PPICLF: '
   print*,'  Inputs FILTERWIDTH    : ',global%piclFilterWidth
   print*,'  Inputs NEIGHBORWIDTH  : ',global%piclNeighborWidth
   print*,'  dp_max (points.dat)               : ',dp_max
   print*,'  d2chk(2) = FILTERWIDTH/2          : ',filter
   print*,'  d2chk(3) = max(4dp,NEIGHBORWIDTH) : ',neighborWidth
   print*,'  d2chk(1) = max(d2chk(2),d2chk(3)) : ',max(filter,neighborWidth)
   print*,'    d2chk(1) = max(d2chk(1),d2chk(2); used in CreateBin'
   print*,'    d2chk(2) = filter; used in filters'
   print*,'    d2chk(3) = neighborWidth; used in nearestNeighbor'
   print*,' '
endif

! TLJ after computing d2chk, we can initialize bins, etc.
call ppiclf_solve_InitParticle(2,3,0,npart_local,y,rprop,filter,neighborWidth) 

! TLJ: CAUTION - Gaussian filter needs to be fixed
! TLJ: Initialize Box Filter
!      This sets ppiclf_d2chk(2) used in Nearest Neighbors
! subroutine ppiclf_solve_InitBoxFilter(filt,iwallm,sngl_elem)
call ppiclf_solve_InitBoxFilter(filter,0,1)!global%piclFilterWidth,0,1)


! Sam - must call initneighborbin before initoverlap for new interpolation scheme
call ppiclf_solve_InitNeighborBin(neighborWidth)

call ppiclf_comm_InitOverlapMesh(nCells,lx,ly,lz,xGrid,yGrid,zGrid)

! 08/13/24 - Thierry - added for Periodicity - begins here

! Thierry - user cannot invoke linear periodicity in x or y when invoking 
!           angular periodicity around z-axis
IF(((ang_per_flag.eq.1) .and. (x_per_flag.eq.1 .or. y_per_flag.eq.1)) .or. & 
    (ang_per_flag .gt. 1)) THEN
    CALL ErrorStop(global,ERR_PICL_INVALID_PERIODICITY,705,&
      'Wrong periodicity choices for ppiclF')
END IF 

! Angular-Periodic
IF(ang_per_flag .eq. 1) then
   IF(global%myProcid == MASTERPROC) print*, "PPICLF Angular Periodic Invoked"
     call ppiclf_solve_InitAngularPeriodic(ang_per_flag , &
                                           ang_per_rin  , ang_per_rout, &
                                           ang_per_angle, ang_per_xangle)
   IF(global%myProcid == MASTERPROC) print*, "PPICLF Angular Periodic Done"
END IF

! Linear X-Periodic
IF(x_per_flag .eq. 1) then  
   IF(global%myProcid == MASTERPROC) print*, "PPICLF PeriodicX Invoked"
   call ppiclf_solve_InitPeriodicX(x_per_min, x_per_max) 
   IF(global%myProcid == MASTERPROC) print*, "PPICLF PeriodicX Done" 
END IF 

! Linear Y-Periodic
IF(y_per_flag .eq. 1) then
   IF(global%myProcid == MASTERPROC) print*, "PPICLF PeriodicY Invoked"
  call ppiclf_solve_InitPeriodicY(y_per_min, y_per_max)
   IF(global%myProcid == MASTERPROC) print*, "PPICLF PeriodicY Done" 
END IF

! Linear Z-Periodic
IF(z_per_flag .eq. 1) then 
   IF(global%myProcid == MASTERPROC) print*, "PPICLF PeriodicZ Invoked"
 call ppiclf_solve_InitPeriodicZ(z_per_min, z_per_max)
   IF(global%myProcid == MASTERPROC) print*, "PPICLF PeriodicZ Done" 
END IF 

! 08/13/24 - Thierry - added for Periodicity - ends here

INQUIRE(FILE='filein.vtk', EXIST=wall_exists)
if (wall_exists) then
  call ppiclf_io_ReadWallVTK('filein.vtk')
else if (global%myProcid == MASTERPROC) then
  WRITE(*,*) 'Could not find filein.vtk'
end if


! 03/24/2025 - Thierry - store the RocfluMP Flow Model chosen (Euler or NS)
!                        this is used in ppiclF for calculating the pressure gradient
!                        whether with or without the viscous part
!flow_model = 0
!flow_model = int(pRegion%mixtInput%flowModel)

! Important note from BRAD:
!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!
!Make sure this call is done specifically when starting the simulation
!for the First Time. As in t = 0.000 else the collision factor will be messed up 
!TLJ commented this out because not needed
!if (pf_fluidInit .eqv. .true. ) call ppiclf_solve_InterpFieldUser(20,xGrid)    

DEALLOCATE(xGrid,STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_DEALLOCATE,765,'PPICLF:xGrid')
END IF ! global%error

DEALLOCATE(yGrid,STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_DEALLOCATE,771,'PPICLF:yGrid')
END IF ! global%error

DEALLOCATE(zGrid,STAT=errorFlag)
global%error = errorFlag
IF ( global%error /= ERR_NONE ) THEN
  CALL ErrorStop(global,ERR_DEALLOCATE,777,'PPICLF:zGrid')
END IF ! global%error

ALLOCATE(vfP(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,783,'PPICLF:xGrid')
    END IF ! global%error

ALLOCATE(volp(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,789,'PPICLF:xGrid')
    END IF ! global%error

do i=1,pGrid%nCellsTot
        pRegion%mixt%piclVF(i) = 0.0_RFREAL
end do

IF ( global%myProcid == MASTERPROC) write(*,*) "PFINIT: Calc Init VolP"
DO i = 1,pRegion%grid%nCells
       volp(i) = 0 
       tester = 0
       do lz=1,2
       do ly=1,2
       do lx=1,2
       call ppiclf_solve_GetProFldIJKEF(lx, ly, lz, i, 1,&
                        vfP(lx,ly,lz,i))

       IF (pRegion%mixtInput%axiFlag) THEN
           WRITE(*,*) "Need to properly implement axi-sym for phip init."
           CALL ErrorStop(global,ERR_OPTION_TYPE,808,'PPICLF:axi')
       END IF

       tester = tester + (0.125*vfP(lx,ly,lz,i))*pRegion%grid%vol(i)

       end do 
       end do 
       end do
       volp(i) = (tester/(pRegion%grid%vol(i)))
!*** VOL FRAC CAP
! TLJ: Increased from 0.6 to 0.62
       if (volp(i) .gt. 0.62) then
           volp(i) = 0.62 
       endif
       pRegion%mixt%piclVF(i) = volp(i) 
end do

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
    if (pRegion%mixt%cv(CV_MIXT_DENS,icg) .le. 0.0) then
         WRITE(*,*) "Error: negative density: ",pRegion%mixt%cv(CV_MIXT_DENS,icg)      
         CALL ErrorStop(global,ERR_INVALID_VALUE,840,'PPICLF:init')
    end if    

END DO ! icg

DEALLOCATE(vfP,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,848,'PPICLF:zGrid')
    END IF ! global%error
DEALLOCATE(volp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,853,'PPICLF:zGrid')
    END IF ! global%error

!!Josh - Removed Brad Comments and Restart Section
!!Can now restart using last par output and converting to .dat format
!!Note that par naming will restart from 00001 (this can be changed)


!1 MOVE HERE END

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
     print*, 'HTUnsteady_flag      = ',global%piclHTUnsteadyFlag
     print*, 'feedback_flag        = ',global%piclFeedbackFlag
     print*, 'qs_fluct_flag        = ',global%piclQsFluctFlag
     print*, 'ppiclf_debug         = ',global%piclDebug
     print*, 'ppiclf_nUnsteadyData = ',ppiclf_nUnsteadyData
     print*, 'ppiclf_VU            = ',50
     print*, 'sbNearest_flag       = ',global%piclSBNearFlag
     print*, 'burnrate_flag        = ',global%piclBurnRateFlag

     IF (global%piclViscousUnsteady >=1) THEN
        print*,'  Using Viscous unsteady history term'
        print*,'    ppiclf_nTimeBH       = ',ppiclf_nTimeBH
        print*,'    ppiclf_nUnsteadyData = ',ppiclf_nUnsteadyData,50
     ENDIF

     print*, ' '
     print*, "BOXF = ", filter 
     print*, "PHIPF = ",  global%pi/6.0*(filter)**3
     print*, "XLOC = ", y(1 ,1) 
     print*, "YLOC = ", y(2 ,1) 
     print*, "ZLOC = ", y(3 ,1)   
     print*, 'Box filter rwidth (in meters) = ',global%piclFilterWidth

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

