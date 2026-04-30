#include "../../ppiclF/source/PPICLF_USER.h"
#include "../../ppiclF/source/PPICLF_STD.h"
!----------------------------------------------------------------------
      PROGRAM main

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF_UNIT_TEST'
      INTEGER*4 i, j, k, l, m, ie

      ! Timer variables
      REAL*8    Tstart, Tend, SBtime, Regtime

      ! Grid variables

      REAL*8    tpF(PPICLF_LEE), numBins,
     >          feedback1(PPICLF_LEE), feedback2(PPICLF_LEE),
     >          x_norm, y_norm, z_norm

      ! Particle variables
      REAL*8    xp(3), totErr, InteriorErr, xFaceErr,
     >          yFaceErr, zFaceErr, xyEdgeErr, xzEdgeErr, yzEdgeErr,
     >          xyzCornerErr, totCnt, InteriorCnt,
     >          xFaceCnt, yFaceCnt, zFaceCnt, xyEdgeCnt, xzEdgeCnt,
     >          yzEdgeCnt, xyzCornerCnt, wInterp(27),
     >          wInterpTot, xcell(3), T_calc,
     >          calcErr, eps, T_analytic, NNDistSQ, nndiff, mindiff,
     >          binblen(3)

      INTEGER*4 ip, np(3), part_cell(3), xstart, xend,
     >          ystart, yend, zstart, zend, icount, loopcount, 
     >          NNCount, id1, id2, projCells, 
     >          cellNumber, icellNumber, ix, iy, iz, 
     >          iParticle(PPICLF_LPART), iteration, series

      LOGICAL   xFace(2), yFace(2), zFace(2), farAway, interpolation(8),
     >          projection(8), binGen, nnpart(8), interp_logical(8),
     >          proj_logical(8), nn_logical(8), PartOnProc, PPInteract

      ! Projection variables
      REAL*8    wProjTot, dSQl, dSQi, dist, CellVol, GaussianConst,
     >          wProj(27),part_feedbk1(PPICLF_LPART),
     >          TrueFeedback1(PPICLF_LEE), TrueFeedback2(PPICLF_LEE),
     >          part_feedbk2(PPICLF_LPART), error1, error2, e1avg, e2avg

      CHARACTER*50 filename, testcase, procString, npstring, ARG1, ARG2

      CHARACTER*2 par 

      CALL GETARG(1, ARG1)
      READ(ARG1,*) iteration

      CALL GETARG(2, ARG2)
      READ(ARG2,*) series

      rootProc = 0
      PI = 4.0D0*ATAN(1.0) ! pi
      ! Set to TRUE and define gridDomain(1:2,1:3)
      ! if you want a particular size for this unit test
      gridBoundsDefined  = .FALSE.

! MPI Setup
!**********************************************************************
      CALL MPI_Init(ierr)
      icomm = MPI_COMM_WORLD
      CALL MPI_Comm_size(icomm,nproc,ierr)
      CALL MPI_Comm_rank(icomm,nid,ierr)
      WRITE(npstring, '(I0)') nproc
      WRITE(procString, '(I0)') nid
      IF(nid .LT. 10) THEN
        procString = '0' // TRIM(procString) //'_' // TRIM(npstring) 
      ELSE
        procString = TRIM(procString) // '_' // TRIM(npstring)
      END IF

      CALL UT_setup(series)
!      IF(nid .EQ. 0) THEN
!        PRINT*, 'Number of particles:',totalParticles
!        PRINT*, 'Number of cells; x: ',nCells(1),', y: ',nCells(2),
!     >    ', z: ',nCells(3), ', total: ', nCells(1)*nCells(2)*nCells(3)
!      END IF
      CALL ppiclf_comm_InitMPI(icomm, nid, nproc)

!      DO l = 1,2
!        IF(l .EQ. 1) test = 1
!        IF(l .EQ. 2) test = 8
        test = 1 
       ! Periodicity flag Setup
        CALL test_setperiodic(x_per_flag,y_per_flag,z_per_flag,
     >                                            test,testcase)
        IF(nid .EQ. rootProc) THEN
          !PRINT*, 'test case: ', TRIM(testcase)
        END IF
        ! Will handle angular periodicity separately
        ang_per_flag   = 0
        ang_per_angle  = 0.0D0
        ang_per_xangle = 0.0D0
        ang_per_rin    = 0.0D0
        ang_per_rout   = 0.0D0  

        CALL MPI_BARRIER(icomm,ierr)
 
! Start ppiclF Calls
!**********************************************************************
        PPICLF_OVERLAP = .FALSE.
        PPInteract = .TRUE.
        CALL ppiclf_solve_InitParticle(2,3,0,npart_local,
     >                                 p_part_y,p_part_r,filter,nndist)
        CALL ppiclf_solve_Initialize(PPInteract,
     >                               x_per_flag, x_per_min, x_per_max,
     >                               y_per_flag, y_per_min, y_per_max, 
     >                               z_per_flag, z_per_min, z_per_max, 
     >                               ang_per_flag, ang_per_angle, 
     >                               ang_per_xangle, ang_per_rin,
     >                                                    ang_per_rout)
        CALL ppiclf_comm_InitOverlapGrid(proc_ncells,p_grid)
        ! Setup fluid temperature field for ppiclf input
        ! This is temperature for cells in processor's grid domain
!        DO i = 1,proc_ncells
!          x_norm = (p_grid(1,i) - gridDomain(1,1))
!     >             /(gridDomain(2,1)-gridDomain(1,1))
!          y_norm = (p_grid(2,i) - gridDomain(1,2))
!     >             /(gridDomain(2,2)-gridDomain(1,2))
!          z_norm = (p_grid(3,i) - gridDomain(1,3))
!     >             /(gridDomain(2,3)-gridDomain(1,3))
!          tpF(i) =   COS(2*PI*x_norm) +
!     >               COS(2*PI*y_norm) +
!     >               COS(2*PI*z_norm)
!        END DO
!        DO i = 1,2
!          CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JT,tpF)
!          CALL ppiclf_comm_CreateBin
!          CALL MPI_BARRIER(ppiclf_comm,ierr)
!          CALL ppiclf_comm_FindParticle
!          CALL ppiclf_comm_MoveParticle
!          CALL ppiclf_comm_MapOverlapGrid
!          CALL MPI_BARRIER(ppiclf_comm,ierr)
!          CALL ppiclf_solve_InitInterp
!          DO j = 1,PPICLF_INT_ICNT
!            CALL ppiclf_solve_InterpField(j)
!          END DO
!          CALL ppiclf_solve_InterpTupleTransfer
!          CALL MPI_BARRIER(ppiclf_comm,ierr)
!          IF(nid .EQ. rootProc) CALL CPU_TIME(Tstart)
!          IF(i .EQ. 1) THEN 
!          ! SB
!            CALL ppiclf_solve_SBParticleToCellMap
!          ELSE
!          !brute force particle to cell mapping
!            CALL ppiclf_solve_ParticleToCellMap
!          END IF
!          CALL MPI_BARRIER(ppiclf_comm,ierr)
!          IF(nid .EQ. rootProc) CALL CPU_TIME(Tend)
!          IF(nid .EQ. rootProc) THEN
!            IF(i .EQ. 1) THEN
!              SBtime = Tend-Tstart
!            ELSE
!              Regtime = Tend-Tstart
!              PRINT*, 'Reg time: ',Regtime, 
!     >                'SB time: ',SBtime,
!     >                'Reg/SB ratio: ', Regtime/SBtime
!            END IF
!          END IF
!        END DO !i
!      END DO !l
!
!      CALL MPI_BARRIER(ppiclf_comm,ierr)
!
!      DO i = 1,10000
!        Regtime = SQRT(Regtime)
!        Regtime = Regtime**2
!      END DO
!      CALL MPI_BARRIER(ppiclf_comm,ierr)

      DO i = 1,2000
        !PRINT*, 'start loop'
        ! Manually set number of MPI ranks
        ppiclf_np = i
        ! ppiclf_binchanged set in CreateBin
        ! ppiclf_binchanged .TRUE. means
        ! bin coordinates changed
        CALL ppiclf_comm_CreateBinPartLB

        ! ppiclf_particleMoved set in FindParticle
        ! ppiclf_particleMoved FALSE means all particles
        ! stayed in same rank as previous RK Stage.
        CALL ppiclf_comm_FindParticlePartLB
        CALL MPI_BARRIER(ppiclf_comm,ierr)

        CALL ppiclf_comm_PartLoadBalance

        CALL ppiclf_comm_setRankBoundaries

        CALL ppiclf_comm_setEmptyIndicator

        CALL ppiclf_comm_setInterfaceIndicator

        ! Make-shift wait timer
        CALL MPI_BARRIER(icomm,ierr)
        DO j = 1,10000
          Regtime = SQRT(Regtime)
          Regtime = Regtime**2
        END DO
        CALL ppiclf_comm_LBCheck(series,iteration)
        ! Make-shift wait timer
        DO j = 1,10000
          Regtime = SQRT(Regtime)
          Regtime = Regtime**2
        END DO
        CALL MPI_BARRIER(icomm,ierr)
      END DO
      !CALL ppiclf_io_WriteParticleVTU('1') 
      CALL MPI_FINALIZE(ierr)

      END PROGRAM

!----------------------------------------------------------------------
      SUBROUTINE test_setperiodic(xf,yf,zf,i,tc)

      ! Input/Output
      ! xf - x periodic flag
      ! yf - y periodic flag
      ! zf - z periodic flag
      ! i  - test iteration
      ! tc - test case name
      IMPLICIT NONE

      INTEGER*4 xf, yf, zf, i
      CHARACTER*50 tc

      xf = 0
      yf = 0
      zf = 0

      IF(i .EQ. 1) THEN
        tc = 'NonPeriodic' 
      ELSE IF(i .EQ. 2) THEN
        tc = 'Periodic_x' 
        xf = 1
      ELSE IF(i .EQ. 3) THEN
        tc = 'Periodic_y' 
        yf = 1
      ELSE IF(i .EQ. 4) THEN
        tc = 'Periodic_z' 
        zf = 1
      ELSE IF(i .EQ. 5) THEN
        tc = 'Periodic_xy' 
        xf = 1
        yf = 1
      ELSE IF(i .EQ. 6) THEN
        tc = 'Periodic_xz' 
        xf = 1
        zf = 1
      ELSE IF(i .EQ. 7) THEN
        tc = 'Periodic_yz'
        yf = 1
        zf = 1 
      ELSE IF(i .EQ. 8) THEN
        tc = 'Periodic_xyz' 
        xf = 1
        yf = 1
        zf = 1           
      END IF

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------

 
