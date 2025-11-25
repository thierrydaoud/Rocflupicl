#include "../../ppiclF/source/PPICLF_USER.h"
#include "../../ppiclF/source/PPICLF_STD.h"
!----------------------------------------------------------------------
      PROGRAM main

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF_UNIT_TEST'

      INTEGER*4 i, j, k, l, ie, BinMaxNum(3), loop, wf
      REAL*8    numBins, xmin, xmax, ymin, ymax, zmin, zmax,
     >          BinMinLen(3)
      LOGICAL   BinsCreated
      CHARACTER*50 filename, testNum

      rootProc = 0
      PI       = 4.0D0*ATAN(1.0) ! pi
      BinsCreated   = .TRUE.

! MPI Setup
!**********************************************************************
      CALL MPI_INIT(ierr)
      icomm = MPI_COMM_WORLD
      CALL MPI_COMM_SIZE(icomm,nproc,ierr)
      CALL MPI_COMM_RANK(icomm,nid,ierr)
      CALL ppiclf_comm_InitMPI(icomm, nid, nproc)
! Write Output Setup 
!**********************************************************************
      IF(nid .EQ. rootProc) THEN
        PRINT*,' CreateBin Test - STARTED'
      END IF
      xmin = -1.0D0
      ymin = -1.0D0
      zmin = -1.0D0

      xmax = 1.0D0
      ymax = 1.0D0
      zmax = 1.0D0

      DO loop = 1,4
        IF(loop .EQ. 2) THEN
          ymax = 40.0D0
        ELSEIF(loop .EQ. 3) THEN
          zmax = 20.0D0
        ELSEIF(loop .EQ. 4) THEN
          xmax = 40.0D0
        END IF
        ! UT_setup has default grid bounds
        ! unless this logical is T
        gridBoundsDefined = .TRUE.
        gridDomain(1,1) = xmin
        gridDomain(2,1) = xmax
        gridDomain(1,2) = ymin
        gridDomain(2,2) = ymax
        gridDomain(1,3) = zmin
        gridDomain(2,3) = zmax
        CALL UT_setup

! ppiclF Inputs and test case setup
!**********************************************************************
        ang_per_flag   = 0
        ang_per_angle  = 0.0D0
        ang_per_xangle = 0.0D0
        ang_per_rin    = 0.0D0
        ang_per_rout   = 0.0D0  
        x_per_flag     = 0  
        y_per_flag     = 0
        z_per_flag     = 0

! Start ppiclF Calls
!**********************************************************************
        PPICLF_OVERLAP = .FALSE.
        CALL ppiclf_solve_InitParticle(2,3,0,npart_local,
     >                                 p_part_y,p_part_r,filter,nndist)
        CALL ppiclf_solve_Initialize(x_per_flag, x_per_min, x_per_max,
     >                               y_per_flag, y_per_min, y_per_max, 
     >                               z_per_flag, z_per_min, z_per_max, 
     >                               ang_per_flag, ang_per_angle, 
     >                               ang_per_xangle, ang_per_rin,
     >                                                    ang_per_rout)
        CALL ppiclf_comm_InitOverlapGrid(proc_ncells,p_grid)
        CALL MPI_BARRIER(icomm,ierr)

! Test CreateBin
!********************************************************************** 
        DO i = 1,3
          BinMinLen(i) = MAX(ppiclf_filter(i),ppiclf_nndist)
          BinMaxNum(i) = INT((gridDomain(2,i) - gridDomain(1,i))
     >                   / BinMinLen(i))
        END DO
        WRITE(testNum, '(I0)') loop
        wf = 100 + loop
        filename = 'Bin_Results_' // TRIM(testNum) //'.txt'
        OPEN(UNIT=wf,FILE=filename, STATUS='REPLACE',ACTION='WRITE')
        WRITE(wf,*) 'New domain setup'
        WRITE(wf,*) 'x min:', xmin,'x max:', xmax
        WRITE(wf,*) 'y min:', xmin,'y max:', ymax
        WRITE(wf,*) 'z min:', xmin,'z max:', zmax
        WRITE(wf,*) 'Max Numbers of bins',
     >               ' based on filter size:',
     >               BinMaxNum(1), BinMaxNum(2), BinMaxNum(3)
        WRITE(wf,*) 'Number of Processors, x bins, y bins, z bins,',
     >                ', total bins, Percent of Processors In Use:' 

        DO i = 1,100000
          ppiclf_np = i
          DO j = 1,3
            ppiclf_n_bins(j) = 0
          END DO
          CALL ppiclf_comm_CreateBin
          ! Catches errors
          numBins = ppiclf_n_bins(1)*
     >               ppiclf_n_bins(2)*ppiclf_n_bins(3)
          IF(numBins .LT. 1 .OR. numBins .GT. ppiclf_np) 
     >       BinsCreated = .FALSE.
          WRITE(wf,*) ppiclf_np, 
     >       ppiclf_n_bins(1), ppiclf_n_bins(2), ppiclf_n_bins(3),
     >       numBins , numBins/ppiclf_np*100
          IF(BinMaxNum(1) .EQ. ppiclf_n_bins(1) .AND.
     >       BinMaxNum(2) .EQ. ppiclf_n_bins(2) .AND.
     >       BinMaxNum(3) .EQ. ppiclf_n_bins(3)) THEN
            WRITE(wf,*) 'Max Number of bins reached based',
     >                   ' on minimum bin lengths!'
            EXIT
          END IF
        END DO
        CLOSE(UNIT=wf)
      END DO
      ! Only passing criteria is the program didn't crash.
      IF(BinsCreated .AND. nid .EQ. rootProc) THEN
        PRINT*,' CreateBin Test - PASSED'
      ELSE
        PRINT*,' CreateBin Test - FAILED'
      END IF
      CALL MPI_FINALIZE(ierr)
      END PROGRAM

