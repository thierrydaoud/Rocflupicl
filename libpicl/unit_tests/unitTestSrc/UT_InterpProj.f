#include "../../ppiclF/source/PPICLF_USER.h"
#include "../../ppiclF/source/PPICLF_STD.h"
!----------------------------------------------------------------------
      PROGRAM main

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF_UNIT_TEST'
      INTEGER*4 i, j, k, l, m, ie

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
     >          iParticle(PPICLF_LPART)

      LOGICAL   xFace(2), yFace(2), zFace(2), farAway, interpolation(8),
     >          projection(8), binGen, nnpart(8), interp_logical(8),
     >          proj_logical(8), nn_logical(8), PartOnProc, PPInteract

      ! Projection variables
      REAL*8    wProjTot, dSQl, dSQi, dist, CellVol, GaussianConst,
     >          wProj(27),part_feedbk1(PPICLF_LPART),
     >          TrueFeedback1(PPICLF_LEE), TrueFeedback2(PPICLF_LEE),
     >          part_feedbk2(PPICLF_LPART), error1, error2, e1avg, e2avg

      CHARACTER*50 filename, testcase, procString, npstring, par

      rootProc = 0
      PI = 4.0D0*ATAN(1.0) ! pi
      DO i = 1,8
        interpolation(i) = .FALSE.
        projection(i)    = .FALSE.
      END DO
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

      CALL UT_setup
      CALL ppiclf_comm_InitMPI(icomm, nid, nproc)

      DO test = 1,8
       
        ! Periodicity flag Setup
        CALL test_setperiodic(x_per_flag,y_per_flag,z_per_flag,
     >                                            test,testcase)
        ! Will handle angular periodicity separately
        ang_per_flag   = 0
        ang_per_angle  = 0.0D0
        ang_per_xangle = 0.0D0
        ang_per_rin    = 0.0D0
        ang_per_rout   = 0.0D0  

 
! Start ppiclF Calls
!**********************************************************************
        PPICLF_OVERLAP = .FALSE.
        PPInteract = .TRUE.
        CALL ppiclf_solve_InitParticle(2,3,0,npart_local,
     >                                 p_part_y,p_part_r,filter,nndist)
        CALL MPI_BARRIER(icomm,ierr)
        CALL ppiclf_solve_Initialize(PPInteract,
     >                               x_per_flag, x_per_min, x_per_max,
     >                               y_per_flag, y_per_min, y_per_max, 
     >                               z_per_flag, z_per_min, z_per_max, 
     >                               ang_per_flag, ang_per_angle, 
     >                               ang_per_xangle, ang_per_rin,
     >                                                    ang_per_rout)
        CALL ppiclf_comm_InitOverlapGrid(proc_ncells,p_grid)
        CALL MPI_BARRIER(icomm,ierr)
        ! Setup fluid temperature field for ppiclf input
        ! This is temperature for cells in processor's grid domain
        DO i = 1,proc_ncells
          x_norm = (p_grid(1,i) - gridDomain(1,1))
     >             /(gridDomain(2,1)-gridDomain(1,1))
          y_norm = (p_grid(2,i) - gridDomain(1,2))
     >             /(gridDomain(2,2)-gridDomain(1,2))
          z_norm = (p_grid(3,i) - gridDomain(1,3))
     >             /(gridDomain(2,3)-gridDomain(1,3))
          tpF(i) =   COS(2*PI*x_norm) +
     >               COS(2*PI*y_norm) +
     >               COS(2*PI*z_norm)
        END DO
        CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JT,tpF)
        CALL MPI_BARRIER(icomm,ierr)
        CALL ppiclf_solve_InitSolve
        CALL MPI_BARRIER(icomm,ierr)
        DO ie = 1,proc_ncells
          CALL ppiclf_solve_GetProFld(ie,1,feedback1(ie))
          CALL ppiclf_solve_GetProFld(ie,2,feedback2(ie))
        END DO
        CALL MPI_BARRIER(icomm,ierr)

        ! Store bin boundary length for periodicity purposes
        DO i = 1,3
          binblen(i) = ppiclf_binb(2*i) - ppiclf_binb(2*i-1)
        END DO

        ! Setup fluid temperature field for error checking
        ! This is temperature for every cell on entire grid domain
        tpF = 0.0
        DO i = 1,numCells
          x_norm = (grid(1,i) - gridDomain(1,1))
     >             /(gridDomain(2,1)-gridDomain(1,1))
          y_norm = (grid(2,i) - gridDomain(1,2))
     >             /(gridDomain(2,2)-gridDomain(1,2))
          z_norm = (grid(3,i) - gridDomain(1,3))
     >             /(gridDomain(2,3)-gridDomain(1,3))
          tpF(i) =   COS(2*PI*x_norm) +
     >               COS(2*PI*y_norm) +
     >               COS(2*PI*z_norm)
        END DO

! Calculate Interpolation Reference and Errors 
!**********************************************************************
        ! Zero out errors and counters
        calcErr      = 0.0D0
        totErr       = 0.0D0
        InteriorErr  = 0.0D0
        xFaceErr     = 0.0D0  
        yFaceErr     = 0.0D0    
        zFaceErr     = 0.0D0    
        xyEdgeErr    = 0.0D0    
        xzEdgeErr    = 0.0D0  
        yzEdgeErr    = 0.0D0   
        xyzCornerErr = 0.0D0    
        totCnt       = 0.0D0          
        InteriorCnt  = 0.0D0          
        xFaceCnt     = 0.0D0          
        yFaceCnt     = 0.0D0          
        zFaceCnt     = 0.0D0          
        xyEdgeCnt    = 0.0D0          
        xzEdgeCnt    = 0.0D0          
        yzEdgeCnt    = 0.0D0          
        xyzCornerCnt = 0.0D0   

        eps = 1.0D-60 !same as ppiclF
        IF(ppiclf_npart .GT. 0) THEN
          T_analytic = 0.0D0
          calcErr    = 0.0D0
          DO i = 1, ppiclf_npart
            xFace(1)     = .FALSE.
            xFace(2)     = .FALSE.
            yFace(1)     = .FALSE.
            yFace(2)     = .FALSE.
            zFace(1)     = .FALSE.
            zFace(2)     = .FALSE.

            xp(1) = ppiclf_y(1,i)
            xp(2) = ppiclf_y(2,i)
            xp(3) = ppiclf_y(3,i)
            
            !adjust gridDX to first and last gridDX if varying griddx
            IF(ABS(xp(1) - gridDomain(1,1)) .LT. gridDX(1)) 
     >           xFace(1)=.TRUE.
            IF(ABS(gridDomain(2,1) - xp(1)) .LT. gridDX(1))  
     >           xFace(2)=.TRUE.

            IF(ABS(xp(2) - gridDomain(1,2)) .LT. gridDX(2))  
     >           yFace(1)=.TRUE.
            IF(ABS(gridDomain(2,2) - xp(2)) .LT. gridDX(2))  
     >           yFace(2)=.TRUE.

            IF(ABS(xp(3) - gridDomain(1,3)) .LT. gridDX(3))  
     >           zFace(1)=.TRUE.
            IF(ABS(gridDomain(2,3) - xp(3)) .LT. gridDX(3))  
     >           zFace(2)=.TRUE.

            ! Find cell that particle resides in
            mindiff = 1.0D10
            DO j = 1,numCells
              nndiff = SQRT((xp(1)-grid(1,j))**2 + (xp(2)-grid(2,j))**2 
     >                                           + (xp(3)-grid(3,j))**2)
              IF(nndiff .LT. mindiff) THEN
                mindiff = nndiff
                cellNumber = j
              END IF
            END DO

            ! Nominally look at cells +/- 1, for 27 total
            xstart = -1
            xend   =  1
            ystart = -1
            yend   =  1
            zstart = -1
            zend   =  1

            ! Don't go +/- 1 for nonperiodic and on face
            IF(.NOT.ppiclf_linperiodic(1) .AND. xFace(1)) xstart = 0
            IF(.NOT.ppiclf_linperiodic(1) .AND. xFace(2)) xend   = 0
            IF(.NOT.ppiclf_linperiodic(2) .AND. yFace(1)) ystart = 0
            IF(.NOT.ppiclf_linperiodic(2) .AND. yFace(2)) yend   = 0 
            IF(.NOT.ppiclf_linperiodic(3) .AND. zFace(1)) zstart = 0
            IF(.NOT.ppiclf_linperiodic(3) .AND. zFace(2)) zend   = 0 

            wInterp    = 0.0D0
            wInterpTot = 0.0D0
            T_calc     = 0.0D0
            DO loopcount = 1,2
              icount = 0
              DO j = xstart,xend
                ix = j
                ! Periodic shift.  If non-periodic, won't be -/+1.
                IF(ix .EQ. -1 .AND. xFace(1)) ix =  nCells(1) - 1 
                IF(ix .EQ.  1 .AND. xFace(2)) ix = -nCells(1) + 1
                DO k = ystart,yend
                  iy = k
                  ! Periodic shift.  If non-periodic, won't be -/+1.
                  IF(iy .EQ. -1 .AND. yFace(1)) iy =  nCells(2) - 1
                  IF(iy .EQ.  1 .AND. yFace(2)) iy = -nCells(2) + 1
                  DO l = zstart,zend
                    iz = l
                    ! Periodic shift.  If non-periodic, won't be -/+1.
                    IF(iz .EQ. -1 .AND. zFace(1)) iz =  nCells(3)-1
                    IF(iz .EQ.  1 .AND. zFace(2)) iz = -nCells(3)+1
                    icount = icount + 1
                    icellNumber = cellNumber + iz + iy*nCells(3) +
     >                            ix*nCells(2)*nCells(3)
                    xcell(1)  = grid(1,icellNumber) 
                    xcell(2)  = grid(2,icellNumber)
                    xcell(3)  = grid(3,icellNumber)
                    IF(loopcount .EQ. 1) THEN
                      ! Establish weighting per cell
                      dist = 0.0D0
                      DO m = 1,3
                        IF(ppiclf_linperiodic(m)) THEN
                          dist =  dist +
     >                      MIN((xp(m)-xcell(m))**2,  
     >                          (binblen(m) - ABS(xp(m)-xcell(m)))**2)
                        ELSE
                          dist =  dist + (xp(m)-xcell(m))**2  
                        END IF
                      END DO
                      dist = SQRT(dist)
                      wInterp(icount) = 1.0D0/((dist + eps))**3
                      wInterpTot  = wInterp(icount) + wInterpTot 

                    ELSE IF(loopcount .EQ. 2) THEN
                        T_analytic = tpF(icellNumber)
                        T_calc = T_calc + 
     >                        (wInterp(icount)/wInterpTot)*T_analytic
                    END IF
                  END DO !l
                END DO !k
              END DO !j
            END DO !loopcount

            ! All particles
            calcErr = calcErr + 
     >          ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) -
     >          ABS(T_calc))
            totCnt = totCnt + 1.0D0

            ! Particles in Interior cells
            IF(.NOT. xFace(1) .AND. .NOT. xFace(2) .AND.
     >         .NOT. yFace(1) .AND. .NOT. yFace(2) .AND.
     >         .NOT. zFace(1) .AND. .NOT. zFace(2)) THEN
              InteriorErr = InteriorErr +
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) -
     >             ABS(T_calc))
              InteriorCnt = InteriorCnt + 1.0D0

            ! Particles in single Face cells
            ELSE IF((xFace(1) .OR. xFace(2)) .AND.
     >              .NOT. yFace(1) .AND. .NOT. yFace(2) .AND.
     >              .NOT. zFace(1) .AND. .NOT. zFace(2)) THEN
              xFaceErr = xFaceErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - 
     >             ABS(T_calc))
              xFaceCnt = xFaceCnt + 1.0D0
            ELSE IF(.NOT. xFace(1) .AND. .NOT. xFace(2) .AND.
     >              (yFace(1) .OR. yFace(2)) .AND.
     >              .NOT. zFace(1) .AND. .NOT. zFace(2)) THEN
              yFaceErr = yFaceErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - 
     >             ABS(T_calc))
              yFaceCnt = yFaceCnt + 1.0D0
            ELSE IF(.NOT. xFace(1) .AND. .NOT. xFace(2) .AND.
     >              .NOT. yFace(1) .AND. .NOT. yFace(2) .AND.
     >              (zFace(1) .OR. zFace(2))) THEN
              zFaceErr = zFaceErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - 
     >             ABS(T_calc))
              zFaceCnt = zFaceCnt + 1.0D0

            ! Particles in Two Face (Edge) cells
            ELSE IF((xFace(1) .OR. xFace(2)) .AND.
     >              (yFace(1) .OR. yFace(2)) .AND.
     >              .NOT. zFace(1) .AND. .NOT. zFace(2)) THEN
              xyEdgeErr = xyEdgeErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - 
     >             ABS(T_calc))
              xyEdgeCnt = xyEdgeCnt + 1.0D0
            ELSE IF((xFace(1) .OR. xFace(2)) .AND.
     >              .NOT. yFace(1) .AND. .NOT. yFace(2) .AND.
     >              (zFace(1) .OR. zFace(2))) THEN
              xzEdgeErr = xzEdgeErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - 
     >             ABS(T_calc))
              xzEdgeCnt = xzEdgeCnt + 1.0D0
            ELSE IF(.NOT. xFace(1) .AND. .NOT. xFace(2) .AND.
     >              (yFace(1) .OR. yFace(2)) .AND.
     >              (zFace(1) .OR. zFace(2))) THEN
              yzEdgeErr = yzEdgeErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - 
     >             ABS(T_calc))
              yzEdgeCnt = yzEdgeCnt + 1.0D0

            ! Particles in Three Face (Corner) cells
            ELSE IF((xFace(1) .OR. xFace(2)) .AND.
     >              (yFace(1) .OR. yFace(2)) .AND.
     >              (zFace(1) .OR. zFace(2))) THEN
              xyzCornerErr = xyzCornerErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - 
     >             ABS(T_calc))
              xyzCornerCnt = xyzCornerCnt + 1.0D0
            ELSE
              PRINT*, 'ERROR!!!! did not classify particles correctly'
            END IF
          END DO !i
        END IF
        CALL MPI_BARRIER(icomm,ierr)
        ! Percent Error Sums
        CALL MPI_Allreduce(calcErr,calcErr,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(InteriorErr,InteriorErr,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(xFaceErr,xFaceErr,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(yFaceErr,yFaceErr,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(zFaceErr,zFaceErr,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(xyEdgeErr,xyEdgeErr,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(xzEdgeErr,xzEdgeErr,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(yzEdgeErr,yzEdgeErr,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(xyzCornerErr,xyzCornerErr,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        ! Counter MPI Sums
        CALL MPI_Allreduce(totCnt,totCnt,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(InteriorCnt,InteriorCnt,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(xFaceCnt,xFaceCnt,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(yFaceCnt,yFaceCnt,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(zFaceCnt,zFaceCnt,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(xyEdgeCnt,xyEdgeCnt,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(xzEdgeCnt,xzEdgeCnt,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(yzEdgeCnt,yzEdgeCnt,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(xyzCornerCnt,xyzCornerCnt,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)

        CALL MPI_BARRIER(icomm,ierr)

        IF(calcErr/totCnt .LT. 1.0D-4) interpolation(test) = .TRUE.
        CALL MPI_ALLREDUCE(interpolation(test), interp_logical(test), 1,
     >                     MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(icomm,ierr)

        IF(nid .EQ. rootProc) THEN
          IF(.NOT. interpolation(test)) PRINT*, 
     >         'Interpolation failed for case: ', TRIM(testcase), 
     >         ' and number of processors: ', nproc
        END IF

        ! Only print out error results for x,y,z periodicitiy on.
        ! This shows all periodicity and interior particles.
        IF(test .EQ. 8) THEN
          IF(nid .EQ. rootProc .AND. nproc .EQ. 1) THEN
            PRINT*,''
            PRINT*,'Linear Periodicity Implemented in x, y, & z'
            PRINT*,'Particle average interpolation errors:'
            PRINT*,'All particles error:      ' , calcErr/totCnt
            PRINT*,'Number of particles:      ' , INT(totCnt)
            PRINT*,'Interior particles error: ' , InteriorErr/
     >                                            InteriorCnt
            PRINT*,'Number of internal:       ' , INT(InteriorCnt)
            PRINT*,'xFace particles error:    ' , xFaceErr/xFaceCnt
            PRINT*,'Number of xFace:          ' , INT(xFaceCnt)
            PRINT*,'yFace particles error:    ' , yFaceErr/yFaceCnt
            PRINT*,'Number of yFace:          ' , INT(yFaceCnt)
            PRINT*,'zFace particles error:    ' , zFaceErr/zFaceCnt
            PRINT*,'Number of zFace:          ' , INT(zFaceCnt)
            PRINT*,'xyEdge particles error:   ' , xyEdgeErr/xyEdgeCnt
            PRINT*,'Number of xyEdge:         ' , INT(xyEdgeCnt)
            PRINT*,'xzEdge particles error:   ' , xzEdgeErr/xzEdgeCnt
            PRINT*,'Number of xzEdge:         ' , INT(xzEdgeCnt)
            PRINT*,'yzEdge particles error:   ' , yzEdgeErr/yzEdgeCnt
            PRINT*,'Number of yzEdge:         ' , INT(yzEdgeCnt)
            PRINT*,'xyzCorner particles error:' , xyzCornerErr/
     >                                            xyzCornerCnt
            PRINT*,'Number of Corner:         ' , INT(xyzCornerCnt)
            PRINT*,''
          END IF
        END IF

! Calculate Projection Error
!**********************************************************************
        TrueFeedback1 = 0.0D0
        TrueFeedback2 = 0.0D0
        DO i = 1, totalParticles 
          xFace(1)     = .FALSE.
          xFace(2)     = .FALSE.
          yFace(1)     = .FALSE.
          yFace(2)     = .FALSE.
          zFace(1)     = .FALSE.
          zFace(2)     = .FALSE.

          xp(1) = part_y(1,i)
          xp(2) = part_y(2,i)
          xp(3) = part_y(3,i)
          
          !adjust gridDX to first and last gridDX if varying griddx
          IF(ABS(xp(1) - gridDomain(1,1)) .LT. gridDX(1)) 
     >         xFace(1)=.TRUE.
          IF(ABS(gridDomain(2,1) - xp(1)) .LT. gridDX(1))  
     >         xFace(2)=.TRUE.

          IF(ABS(xp(2) - gridDomain(1,2)) .LT. gridDX(2))  
     >         yFace(1)=.TRUE.
          IF(ABS(gridDomain(2,2) - xp(2)) .LT. gridDX(2))  
     >         yFace(2)=.TRUE.

          IF(ABS(xp(3) - gridDomain(1,3)) .LT. gridDX(3))  
     >         zFace(1)=.TRUE.
          IF(ABS(gridDomain(2,3) - xp(3)) .LT. gridDX(3))  
     >         zFace(2)=.TRUE.

          ! Find cell that particle resides in
          mindiff = 1.0D10
          DO j = 1,numCells
            nndiff = SQRT((xp(1)-grid(1,j))**2 + (xp(2)-grid(2,j))**2 
     >                                         + (xp(3)-grid(3,j))**2)
            IF(nndiff .LT. mindiff) THEN
              mindiff = nndiff
              cellNumber = j
            END IF
          END DO

          ! Nominally look at cells +/- 1, for 27 total
          xstart = -1
          xend   =  1
          ystart = -1
          yend   =  1
          zstart = -1
          zend   =  1

          ! Don't go +/- 1 for nonperiodic and on face
          IF(.NOT.ppiclf_linperiodic(1) .AND. xFace(1)) xstart = 0
          IF(.NOT.ppiclf_linperiodic(1) .AND. xFace(2)) xend   = 0
          IF(.NOT.ppiclf_linperiodic(2) .AND. yFace(1)) ystart = 0
          IF(.NOT.ppiclf_linperiodic(2) .AND. yFace(2)) yend   = 0 
          IF(.NOT.ppiclf_linperiodic(3) .AND. zFace(1)) zstart = 0
          IF(.NOT.ppiclf_linperiodic(3) .AND. zFace(2)) zend   = 0 

          wProj      = 0.0D0
          wProjTot   = 0.0D0
          DO loopcount = 1,2
            icount = 0
            DO j = xstart,xend
              ix = j
              ! Periodic shift.  If non-periodic, won't be -/+1.
              IF(ix .EQ. -1 .AND. xFace(1)) ix =  nCells(1) - 1 
              IF(ix .EQ.  1 .AND. xFace(2)) ix = -nCells(1) + 1
              DO k = ystart,yend
                iy = k
                ! Periodic shift.  If non-periodic, won't be -/+1.
                IF(iy .EQ. -1 .AND. yFace(1)) iy =  nCells(2) - 1
                IF(iy .EQ.  1 .AND. yFace(2)) iy = -nCells(2) + 1
                DO l = zstart,zend
                  iz = l
                  ! Periodic shift.  If non-periodic, won't be -/+1.
                  IF(iz .EQ. -1 .AND. zFace(1)) iz =  nCells(3)-1
                  IF(iz .EQ.  1 .AND. zFace(2)) iz = -nCells(3)+1
                  icount = icount + 1
                  icellNumber = cellNumber + iz + iy*nCells(3) +
     >                          ix*nCells(2)*nCells(3)
     
                  xcell(1)  = grid(1,icellNumber) 
                  xcell(2)  = grid(2,icellNumber)
                  xcell(3)  = grid(3,icellNumber)
                  IF(loopcount .EQ. 1) THEN
                    ! Establish weighting per cell
                    dist = 0.0D0
                    DO m = 1,3
                      IF(ppiclf_linperiodic(m)) THEN
                        dist =  dist +
     >                    MIN((xp(m)-xcell(m))**2,  
     >                        (binblen(m) - ABS(xp(m)-xcell(m)))**2)
                      ELSE
                        dist =  dist + (xp(m)-xcell(m))**2  
                      END IF
                    END DO
                    dist = SQRT(dist)
                    CellVol = grid(7,icellNumber)
                    GaussianConst = 2.305D0
                    wProj(icount) = ABS(CellVol*EXP(-GaussianConst*
     >                         (dist**2) / (CellVol**(2.0D0/3.0D0))))
                    wProjTot = wProjTot + wProj(icount)          

                  ELSE IF(loopcount .EQ. 2) THEN
                    ! calculate particle feedback to cell
                    x_norm  = (xp(1) - gridDomain(1,1))
     >                        /(gridDomain(2,1)-gridDomain(1,1))
                    y_norm  = (xp(2) - gridDomain(1,2))
     >                        /(gridDomain(2,2)-gridDomain(1,2))
                    z_norm  = (xp(3) - gridDomain(1,3))
     >                        /(gridDomain(2,3)-gridDomain(1,3))
                    part_feedbk1(i) = 1.0
                    part_feedbk2(i) = SIN(2*PI*x_norm) +
     >                                SIN(2*PI*y_norm) +
     >                                SIN(2*PI*z_norm)
                    TrueFeedback1(icellNumber) = 
     >                                   TrueFeedback1(icellNumber) + 
     >                                   (wProj(icount)/wProjTot) *
     >                                   part_feedbk1(i)
                    TrueFeedback2(icellNumber) = 
     >                                   TrueFeedback2(icellNumber) +
     >                                   (wProj(icount)/wProjTot) *
     >                                   part_feedbk2(i)
                   END IF
                END DO !l
              END DO !k
            END DO !j
          END DO !loopcount
        END DO
        e1avg = 0.0D0
        e2avg = 0.0D0
        IF(proc_ncells .GT. 0) THEN
          DO ie = 1, proc_ncells
            error1 = ABS(ABS(feedback1(ie)) -
     >                  ABS(TrueFeedback1(iCstart - 1 + ie)))
            e1avg = e1avg + error1
            error2 = ABS(ABS(feedback2(ie)) -
     >                  ABS(TrueFeedback2(iCstart - 1 + ie)))
            e2avg = e2avg + error2
          END DO
          e1avg = e1avg / REAL(proc_ncells)
          e2avg = e2avg / REAL(proc_ncells)
          IF(e1avg .LT. 1.0D-4 .AND. e2avg .LT. 1.0D-4)
     >      projection(test) = .TRUE.
        ELSE
          projection(test) = .TRUE.
        END IF
        CALL MPI_ALLREDUCE(projection(test), proj_logical(test), 1,
     >                     MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
        e1avg = e1avg * REAL(proc_ncells)
        e2avg = e2avg * REAL(proc_ncells)
        projCells = proc_ncells

        ! Print out average feedback error for x,y,z periodic case
        CALL MPI_Allreduce(e1avg,e1avg,1,MPI_DOUBLE,
     >                                      MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(e2avg,e2avg,1,MPI_DOUBLE,
     >                                      MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(projCells,projCells,1,MPI_INTEGER,
     >                                      MPI_SUM,iComm,ierr)
        e1avg = e1avg / REAL(projCells)
        e2avg = e2avg / REAL(projCells)
        IF(test .EQ. 8) THEN
          IF(nid .EQ. rootProc .AND. nproc .EQ. 1) THEN
            PRINT*,'Linear Periodicity Implemented in x, y, & z'
            PRINT*,'Cell average projection errors:'
            PRINT*,''          
            PRINT*, 'Unity projection:', e1avg
            PRINT*, 'Sine projeciton:', e2avg
            PRINT*, ''
          END IF
        END IF
        IF(nid .EQ. rootProc) THEN
          IF(.NOT. proj_logical(test)) PRINT*, 
     >         'Projection failed for case: ', TRIM(testcase), 
     >         ' and number of processors: ', nproc
        END IF

        CALL MPI_BARRIER(icomm,ierr)
!********************************************************************** 
! End periodicity loop testing
      END DO !test
! Print final results & close out program
!********************************************************************** 
      IF(nproc .GT. 1) THEN
        par = ' PARALLEL'
      ELSE
        par = ''
      END IF
      IF(nid .EQ. rootProc) THEN
        PRINT*, 'NonPeriodic Results:'

        IF(interp_logical(1)) THEN
          PRINT*, TRIM(par) // ' INTERPOLATION - PASSED'
        ELSE
          PRINT*, TRIM(par) // ' INTERPOLATION - FAILED'
        END IF

        IF(proj_logical(1)) THEN
          PRINT*, TRIM(par) // ' PROJECTION - PASSED'
        ELSE
          PRINT*, TRIM(par) // ' PROJECTION - FAILED'
        END IF


        PRINT*, 'Linear Periodicity Results:'

        IF(interp_logical(2) .AND.
     >     interp_logical(3) .AND.
     >     interp_logical(4) .AND.
     >     interp_logical(5) .AND.
     >     interp_logical(6) .AND.
     >     interp_logical(7) .AND.
     >     interp_logical(8)      ) THEN
          PRINT*, TRIM(par) // ' INTERPOLATION - PASSED'
        ELSE
          PRINT*, TRIM(par) // ' INTERPOLATION - FAILED'
        END IF

        IF(proj_logical(2) .AND.
     >     proj_logical(3) .AND.
     >     proj_logical(4) .AND.
     >     proj_logical(5) .AND.
     >     proj_logical(6) .AND.
     >     proj_logical(7) .AND.
     >     proj_logical(8)      ) THEN
          PRINT*, TRIM(par) // ' PROJECTION - PASSED'
        ELSE
          PRINT*, TRIM(par) // ' PROJECTION - FAILED'
        END IF
      END IF

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

 
