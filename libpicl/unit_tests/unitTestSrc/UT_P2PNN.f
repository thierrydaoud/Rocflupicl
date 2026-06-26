#include "../../ppiclF/source/PPICLF_USER.h"
#include "../../ppiclF/source/PPICLF_STD.h"
!----------------------------------------------------------------------
      PROGRAM main

      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF_UNIT_TEST'
      INTEGER*4 i, j, k, l, ie, iii

      ! Grid variables

      REAL*8    tpF(PPICLF_LEE), numBins,
     >          feedback1(PPICLF_LEE), feedback2(PPICLF_LEE),
     >          x_norm, y_norm, z_norm

      ! Particle variables
      REAL*8    xp, yp, zp, totErr, InteriorErr, xFaceErr,
     >          yFaceErr, zFaceErr, xyEdgeErr, xzEdgeErr, yzEdgeErr,
     >          xyzCornerErr, totCnt, InteriorCnt,
     >          xFaceCnt, yFaceCnt, zFaceCnt, xyEdgeCnt, xzEdgeCnt,
     >          yzEdgeCnt, xyzCornerCnt, weightCell(27),
     >          weightTot, xcell, ycell, zcell, T_calc,
     >          calcErr, eps, T_analytic, NNDistSQ

      INTEGER*4 ip, np(3), part_cell(3), xstart, xend,
     >          ystart, yend, zstart, zend, icount, loopcount, 
     >          NNCount, id1, id2, projCells

      LOGICAL   xFace, yFace, zFace, farAway, interpolation(8), 
     >          projection(8), binGen, nnpart(8), interp_logical(8),
     >          proj_logical(8), nn_logical(8), fineGrid

      ! Projection variables
      REAL*8    wsum, dSQl, dSQi, dist, CellVol, GaussianConst,
     >          w(PPICLF_LEE),part_feedbk1(PPICLF_LPART),
     >          TrueFeedback1(PPICLF_LEE), TrueFeedback2(PPICLF_LEE),
     >          part_feedbk2(PPICLF_LPART), error1, error2, e1avg, e2avg

      CHARACTER*50 filename, testcase, procString, npstring, par

      REAL*8    SBtime, Regtime, Tend, Tstart

      IF(nproc .GT. 1) THEN
        par = ' PARALLEL'
      ELSE
        par = ''
      END IF
      rootProc = 0
      PI = 4.0D0*ATAN(1.0) ! pi
      DO i = 1,8
        nnpart(i)        = .TRUE.
      END DO
      gridBoundsDefined  = .FALSE.
! MPI Setup
!**********************************************************************
      CALL MPI_Init(ierr)
      icomm = MPI_COMM_WORLD
      CALL MPI_Comm_size(icomm,nproc,ierr)
      CALL MPI_Comm_rank(icomm,nid,ierr)
      CALL ppiclf_comm_InitMPI(icomm, nid, nproc)
      WRITE(npstring, '(I0)') nproc
      WRITE(procString, '(I0)') nid
      IF(nid .LT. 10) THEN
        procString = '0' // TRIM(procString) //'_' // TRIM(npstring) 
      ELSE
        procString = TRIM(procString) // '_' // TRIM(npstring)
      END IF

      CALL UT_setup

! Particle to Particle Nearest Neighbor Test
!**********************************************************************
! Start loop for varying periodicity cases
      DO test = 1,8
        ! Periodicity flag Setup
        CALL test_setperiodic(x_per_flag,y_per_flag,z_per_flag,
     >                                            test,testcase)
!        IF(nid .EQ. rootProc) PRINT*, 'Case:',testcase
        ! Will handle angular periodicity separately
        ang_per_flag   = 0
        ang_per_angle  = 0.0D0
        ang_per_xangle = 0.0D0
        ang_per_rin    = 0.0D0
        ang_per_rout   = 0.0D0  

        CALL MPI_BARRIER(icomm,ierr)
! ppiclF Calls
!**********************************************************************
        PPICLF_OVERLAP = .FALSE.
        CALL ppiclf_solve_InitParticle(2,3,0,npart_local,
     >                                 p_part_y,p_part_r,filter,nndist)
        CALL ppiclf_solve_Initialize(.TRUE.,
     >                               x_per_flag, x_per_min, x_per_max,
     >                               y_per_flag, y_per_min, y_per_max, 
     >                               z_per_flag, z_per_min, z_per_max, 
     >                               ang_per_flag, ang_per_angle, 
     >                               ang_per_xangle, ang_per_rin,
     >                                                    ang_per_rout)
        CALL ppiclf_comm_InitOverlapGrid(proc_ncells,p_grid)
        CALL ppiclf_comm_CreateBinPartLB
        CALL ppiclf_comm_FindParticlePartLB
        CALL ppiclf_comm_PartLoadBalance
        CALL ppiclf_comm_MoveParticlePartLB
        CALL ppiclf_comm_subbinRealParticleMap
        CALL ppiclf_comm_CreateGhostPartLB
        CALL ppiclf_comm_MoveGhostPartLB

        IF(ppiclf_npart .GT. 0) THEN 
          ! SetYdot only runs with nearest neighbor subroutine
          ! due to preprocess compile flag.  PARTICLE_NN(i) saves
          ! the number of nearest neighbors per particle in ppiclF.

!          !warming paths once
!          CALL ppiclf_comm_subbinFineParticleMap
!          DO iii = 1,ppiclf_npart
!            CALL ppiclf_solve_NearestNeighborSB(iii,.TRUE.)
!          END DO
!          CALL ppiclf_comm_subbinRealParticleMap
!          CALL ppiclf_comm_subbinGhostParticleMap
!          DO iii = 1,ppiclf_npart
!            CALL ppiclf_solve_NearestNeighborSB(iii,.FALSE.)
!          END DO

          DO i = 1,1
            CALL MPI_BARRIER(icomm,ierr)
            Tstart = MPI_WTIME()
            IF(i .EQ. 1) THEN 
              ! SB
              ! ppiclf_useFineGrid set false when compiled with TEST=1
              CALL ppiclf_comm_subbinFineParticleMap
              fineGrid = .TRUE.
              DO iii = 1,ppiclf_npart
                CALL ppiclf_solve_NearestNeighborSB(iii,fineGrid)
              END DO
            ELSE
              !CALL ppiclf_comm_subbinRealParticleMap
              CALL ppiclf_comm_subbinGhostParticleMap
              fineGrid = .FALSE.
              DO iii = 1,ppiclf_npart
                CALL ppiclf_solve_NearestNeighborSB(iii,fineGrid)
              END DO
            END IF
            CALL MPI_BARRIER(ppiclf_comm,ierr)
            Tend = MPI_WTIME()
            IF(nid .EQ. rootProc) THEN
              IF(i .EQ. 1) THEN
                SBtime = Tend-Tstart
              ELSE
                Regtime = Tend-Tstart
                PRINT*, 'Num Proc', ppiclf_np
                PRINT*, 'Reg time: ',Regtime 
                PRINT*, 'SB time: ',SBtime
                PRINT*, 'Reg/SB ratio: ', Regtime/SBtime
              END IF
              Tend = 0.0D0
            END IF
 
          END DO !i

          CALL MPI_BARRIER(icomm,ierr)
          ! Calculate number of nearest neighbors
          ! for all particles in domain, and compare
          ! with results from ppiclF
          DO i = 1,ppiclf_npart 
            id1 = ppiclf_iprop(1,i)
            id2 = ppiclf_iprop(2,i)
            NNCount   = 0
            NNDistSQ  = 0.0D0
            ! Loop through all particles in domain
            DO j = 1, totalparticles
              IF(ppiclf_y(1,i) .EQ. part_y(1,j) .AND.
     >           ppiclf_y(2,i) .EQ. part_y(2,j) .AND.
     >           ppiclf_y(3,i) .EQ. part_y(3,j)) CYCLE

!              IF(j .EQ. (id1+id2*particlesPerProc)) CYCLE 
              dSQi = 0.0D0
              dSQl = 0.0D0
              DO l = 1,3
                IF(ppiclf_linperiodic(l) .AND. 
     >                               ppiclf_EqualDomain(l)) THEN
                  dSQl = MIN( (ppiclf_y(l,i) - part_y(l,j))**2, 
     >             ( (gridDomain(2,l) - gridDomain(1,l)) -
     >                   ABS(ppiclf_y(l,i) - part_y(l,j)) )**2 )
                ELSE
                  dSQl = (ppiclf_y(l,i) - part_y(l,j))**2
                END IF
                dSQi = dSQi + dSQl
              END DO !l
              IF(dSQi .GT. ppiclf_nndist**2) CYCLE
              NNCount = NNCount + 1
              NNDistSQ = NNDistSQ + dSQi
            END DO
            IF(NNCount .NE. PARTICLE_NN(i)) THEN
              IF((ABS((ABS(NNDistSQ - PPICLF_TOTNNDIST(i))
     >             /ABS(NNCount-PARTICLE_NN(i)))-ppiclf_nndist**2)
     >                 /ppiclf_nndist**2) .GT. 0.1) THEN
                nnpart(test) = .FALSE.
                ! Check to see if periodicity causing issue
                IF(ABS(ppiclf_y(1,i)-gridDomain(2,1)) .LT. 
     >             ppiclf_nndist .OR.
     >             ABS(ppiclf_y(1,i)-gridDomain(1,1)) .LT. 
     >             ppiclf_nndist .OR.
     >             ABS(ppiclf_y(2,i)-gridDomain(2,2)) .LT. 
     >             ppiclf_nndist .OR.
     >             ABS(ppiclf_y(2,i)-gridDomain(1,2)) .LT. 
     >             ppiclf_nndist .OR.
     >             ABS(ppiclf_y(3,i)-gridDomain(2,3)) .LT. 
     >             ppiclf_nndist .OR.
     >             ABS(ppiclf_y(3,i)-gridDomain(1,3)) .LT. 
     >             ppiclf_nndist) THEN
!                    PRINT*, 'PERIODIC ISSUE!'
                ELSE
!                     PRINT*, 'NonPeriodic Interface ISSUE!'
                END IF
!                PRINT*, 'Count is different! ppiclf , truth:', 
!     >                  PARTICLE_NN(i),NNCount
!                PRINT*, 'Particle Location:', ppiclf_y(1,i),
!     >                  ppiclf_y(2,i), ppiclf_y(3,i)
!                PRINT*, '((RefDist^2-ppiclf_dist^2)/#NNdiff)/nndist**2'
!     >            ,(ABS((ABS(NNDistSQ - PPICLF_TOTNNDIST(i))
!     >             /ABS((NNCount-PARTICLE_NN(i)))))
!     >                 /ppiclf_nndist**2)
!                PRINT*, ''
!                PRINT*, ''
!                IF(PARTICLE_NN(i) .GT. NNCount) THEN
!                  PRINT*, 'More on ppiclf, distdiff**2/CountDiff:',
!     >                    (PPICLF_TOTNNDIST(i)-NNDistSQ)
!     >                    /(PARTICLE_NN(i)-NNCount) ,
!     >                    'nndist^2:',ppiclf_nndist**2
!                  PRINT*, PARTICLE_NN(i),NNCount
!                ELSE
!                  PRINT*, 'More on ref, distdiff**2/CountDiff:',
!     >                    (- PPICLF_TOTNNDIST(i)+NNDistSQ)
!     >                    /(-PARTICLE_NN(i)+NNCount) ,
!     >                    'nndist:',ppiclf_nndist**2
!                  PRINT*, NNCount, PARTICLE_NN(i)
!                END IF
!                EXIT
              END IF
            END IF
          END DO
        END IF
        ! Ensure nnpart(test) is true across all processors
        CALL MPI_ALLREDUCE(nnpart(test), nn_logical(test), 1,
     >                     MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(icomm,ierr)
        IF(test .EQ. 1) THEN
          IF(nid .EQ. rootProc) THEN
            PRINT*, 'Non-Periodic Case:'
          END IF
          CALL ppiclf_solve_PrintQuantities
        ELSE IF(test .EQ. 8) THEN
          IF(nid .EQ. rootProc) THEN
            PRINT*, 'Triply-Perioidc Case:'
          END IF
          CALL ppiclf_solve_PrintQuantities
        END IF
        ! Print Pas/Fail
        IF(nid .EQ. rootProc) THEN
          IF(test .EQ. 1) THEN
            PRINT*, ''
            PRINT*, 'NonPeriodic Results:'
            IF(nproc .EQ. 1) THEN
              IF(nn_logical(1)) THEN
                PRINT*, ' PARTICLE NEAREST NEIGHBOR SEARCH - PASSED'
              ELSE
                PRINT*, ' PARTICLE NEAREST NEIGHBOR SEARCH - FAILED'
              END IF
            ELSE
              IF(nn_logical(1)) THEN
                PRINT*, TRIM(par) // ' PARTILCE NN SEARCH ',
     >                               'WITH GHOST PARTICLES - PASSED'
              ELSE
                PRINT*, TRIM(par) // ' PARTILCE NN SEARCH ',
     >                               'WITH GHOST PARTICLES - FAILED'
              END IF
            END IF
          ELSE IF(test .EQ. 8) THEN
            PRINT*, 'Linear Periodicity Results:'
            IF(nn_logical(2) .AND.
     >         nn_logical(3) .AND.
     >         nn_logical(4) .AND.
     >         nn_logical(5) .AND.
     >         nn_logical(6) .AND.
     >         nn_logical(7) .AND.
     >         nn_logical(8)      ) THEN
              PRINT*, TRIM(par) // ' PARTILCE NN SEARCH WITH '
     >                            ,'GHOST PARTICLES - PASSED'
            ELSE
              PRINT*, TRIM(par) // ' PARTILCE NN SEARCH WITH '
     >                            ,'GHOST PARTICLES - FAILED'
            END IF
          END IF
        END IF
      END DO !test
!********************************************************************** 
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

 
