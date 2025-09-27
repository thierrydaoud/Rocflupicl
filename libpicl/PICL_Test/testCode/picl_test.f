#include "../../ppiclF/source/PPICLF_USER.h"
#include "../../ppiclF/source/PPICLF_STD.h"
!----------------------------------------------------------------------
      PROGRAM main

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF'

      INTEGER*4 i, j, k, l, nproc, nid, icomm, ierr, test

      REAL*8    nndistTemp, nndist, PI

      ! Grid variables
      INTEGER*4 nCells(3), proc_ncells, ie, numCells,
     >          MapGrid(PPICLF_LEE), iCend, iCstart,
     >          cellsPerProc, iPend, iPstart

      REAL*8    p_grid(7,PPICLF_LEE), grid(7,PPICLF_LEE),
     >          gridDomain(2,3), gridDX(3), filter(3), 
     >          nFilterCells, tpF(PPICLF_LEE), numBins,
     >          dx_min(3), feedback1(PPICLF_LEE), feedback2(PPICLF_LEE),
     >          x_norm, y_norm, z_norm

      ! Particle variables
      REAL*8    part_y(PPICLF_LRS,PPICLF_LPART), pdia, 
     >          p_part_y(PPICLF_LRS,PPICLF_LPART), part_dx(3),
     >          p_part_r(PPICLF_LRP,PPICLF_LPART),
     >          xp, yp, zp, totErr, InteriorErr, xFaceErr,
     >          yFaceErr, zFaceErr, xyEdgeErr, xzEdgeErr, yzEdgeErr,
     >          xyzCornerErr, totCnt, InteriorCnt,
     >          xFaceCnt, yFaceCnt, zFaceCnt, xyEdgeCnt, xzEdgeCnt,
     >          yzEdgeCnt, xyzCornerCnt, weightCell(27),
     >          weightTot, xcell, ycell, zcell, T_calc,
     >          calcErr, eps, T_analytic, NNDistSQ, randNum 

      INTEGER*4 particlesPerCell(3), particlesPerProc, npart_local,
     >          totalParticles, ip, np(3), part_cell(3), xstart, xend,
     >          ystart, yend, zstart, zend, icount, loopcount, 
     >          NNCount, id1, id2, projCells

      LOGICAL   xFace, yFace, zFace, farAway, interpolation(8), 
     >          projection(8), binGen, nnpart(8), interp_logical(8),
     >          proj_logical(8), nn_logical(8)

      ! Projection variables
      REAL*8    wsum, dSQl, dSQi, dist, CellVol, GaussianConst,
     >          w(PPICLF_LEE),part_feedbk1(PPICLF_LPART),
     >          TrueFeedback1(PPICLF_LEE), TrueFeedback2(PPICLF_LEE),
     >          part_feedbk2(PPICLF_LPART), error1, error2, e1avg, e2avg

      CHARACTER*50 filename, testcase, procString, npstring, par


      PI = 4.0D0*ATAN(1.0) ! pi
      DO i = 1,8
        interpolation(i) = .FALSE.
        projection(i)    = .FALSE.
        nnpart(i)        = .TRUE.
      END DO
      binGen             = .FALSE.

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
! Grid Setup
!**********************************************************************
      ! Create rectangular grid
      gridDomain(1,1) = 0.0D0 !x domain min
      gridDomain(2,1) = 10.0D0 !x domain max
      nCells(1)       = 20    !Number of x cells in domain
      gridDX(1) = (gridDomain(2,1) - gridDomain(1,1))/REAL(nCells(1))

      gridDomain(1,2) = 0.0D0 !y domain min
      gridDomain(2,2) = 1.5D0 !y domain max     
      nCells(2)       = 10    !Number of y cells in domain
      gridDX(2) = (gridDomain(2,2) - gridDomain(1,2))/REAL(nCells(2))

      gridDomain(1,3) = 0.0D0 !z domain min
      gridDomain(2,3) = 1.0D0 !z domain max     
      nCells(3)       = 10    !Number of z cells in domain
      gridDX(3) = (gridDomain(2,3) - gridDomain(1,3))/REAL(nCells(3))

      IF(nCells(1)*nCells(2)*nCells(3) .GT. PPICLF_LEE) THEN
        PRINT*, 'ERROR, PPICLF_LEE = ',PPICLF_LEE, 'Cells to test =',
     >          nCells(1)*nCells(2)*nCells(3)
        CALL MPI_FINALIZE(ierr)
        STOP
      END IF

      ! Build full grid on each processor 
      ! Build in order different from bin numbering
      ! to ensure GSLIB calls are fully tested.
      numCells = 0
      DO k = 1,nCells(3)
        DO j = 1,nCells(2)
          DO i = 1,nCells(1)
            numCells = numCells + 1
            grid(1,numCells) = gridDomain(1,1) + (i-0.5)*gridDX(1) !x centroid
            grid(2,numCells) = gridDomain(1,2) + (j-0.5)*gridDX(2) !y centroid
            grid(3,numCells) = gridDomain(1,3) + (k-0.5)*gridDX(3) !z centroid
            grid(4,numCells) = gridDX(1)
            grid(5,numCells) = gridDX(2)
            grid(6,numCells) = gridDX(3)
            grid(7,numCells) = gridDX(1)*gridDX(2)*gridDX(3)
          END DO !k
        END DO !j
      END DO !i

      ! Divide grid for ppiclf input by number of parallel processors
      cellsPerProc = INT(REAL(numCells)/REAL(nproc))
      IF(cellsPerProc .EQ. 0) cellsPerProc = 1
      iCstart = cellsPerProc*(nid  ) + 1
      iCend   = cellsPerProc*(nid+1)
      proc_ncells = 0
      IF(iCstart .LE. numCells) THEN
        IF(iCend .GT. numCells) iCend = numCells
        IF(nid .EQ. nproc-1) iCend = numCells
!        IF(nproc .EQ. 1) THEN
!          filename = 'GridPoints.txt'
!          OPEN(UNIT=1,FILE=filename, STATUS='REPLACE', ACTION='WRITE')
!        END IF
        DO i = iCstart,iCend  
          proc_ncells = proc_ncells + 1
          MapGrid(proc_ncells) = i
          DO j = 1,7
            p_grid(j,proc_ncells) = grid(j,i)
          END DO !j
!          IF(nproc .EQ. 1) 
!     >      WRITE(1,*) i, p_grid(1,i), p_grid(2,i), p_grid(3,i),
!     >                    p_grid(4,i), p_grid(5,i), p_grid(6,i),
!     >                    p_grid(7,i)
        END DO !i
!        IF(nproc .EQ. 1) CLOSE(UNIT=1)
      END IF

      CALL MPI_BARRIER(icomm, ierr)

      ! Find cell filter search distance
      filter       = 1.0D-9 !dummy
      dx_min       = 1.0D9  !dummy
      nFilterCells = 2.0
      DO j = 1,proc_ncells
        DO i = 4,6
          ! Find largest & smallest grid dx, dy, dz
          IF(p_grid(i,j) > filter(i-3))
     >      filter(i-3) = p_grid(i,j)
          IF(p_grid(i,j) < dx_min(i-3)) dx_min(i-3) = p_grid(i,j)
        END DO
      END DO

      DO i = 1,3
        filter(i) = nFilterCells*filter(i)
      END DO
      
      ! Setup fluid temperature field for ppiclf input
      DO i = 1,proc_ncells
        x_norm = (p_grid(1,i) - gridDomain(1,1))
     >           /(gridDomain(2,1)-gridDomain(1,1))
        y_norm = (p_grid(2,i) - gridDomain(1,2))
     >           /(gridDomain(2,2)-gridDomain(1,2))
        z_norm = (p_grid(3,i) - gridDomain(1,3))
     >           /(gridDomain(2,3)-gridDomain(1,3))
        tpF(i) =   COS(2*PI*x_norm) +
     >             COS(2*PI*y_norm) +
     >             COS(2*PI*z_norm)
      END DO

      ! Setup filter(1:3) and smallest cell dx across processors 
      DO i = 1,3
        CALL MPI_Allreduce(filter(i),filter(i),1,MPI_DOUBLE,
     >                                      MPI_MAX,iComm,ierr)
        CALL MPI_Allreduce(dx_min(i),dx_min(i),1,MPI_DOUBLE,
     >                                      MPI_MIN,iComm,ierr)
      END DO

      ! Fluid Domain Min/Max
      x_per_min = gridDomain(1,1)
      x_per_max = gridDomain(2,1)
      y_per_min = gridDomain(1,2)
      y_per_max = gridDomain(2,2)
      z_per_min = gridDomain(1,3)
      z_per_max = gridDomain(2,3)

      CALL MPI_BARRIER(icomm,ierr)

! Particle Setup   
!********************************************************************** 
      particlesPerCell(1) = 3.0 ! in x dimension
      particlesPerCell(2) = 3.0 ! in y dimension
      particlesPerCell(3) = 3.0 ! in z dimension

      ! set particle diameter to prevent overlapping
      pdia  = MIN(dx_min(1)/ParticlesPerCell(1),
     >            dx_min(2)/ParticlesPerCell(2),
     >            dx_min(3)/ParticlesPerCell(3))

      totalParticles = 1 ! global number of particles
      DO i = 1,3
        part_dx(i) = gridDX(i)/ParticlesPerCell(i)
        np(i) = particlesPerCell(i)*nCells(i) ! Num particles per dimension
        totalParticles = totalParticles*np(i)
      END DO

      IF(totalParticles .GT. PPICLF_LPART) THEN
        PRINT*, 'Tried to make too many particles!'
        PRINT*, 'PPICLF_LPART =',PPICLF_LPART,'TotalParticles=',
     >          totalParticles
        CALL MPI_FINALIZE(ierr)
        STOP
      END IF

      ! Build full particle dispersion on each processor 
      ! Build in order different from bin & grid numbering
      ! to ensure GSLIB calls are fully tested.
      totalParticles = 0
      DO j = 1,np(2)
        DO k = 1,np(3)
          DO i = 1,np(1)
            totalParticles = totalParticles + 1
            part_y(1,totalParticles) = gridDomain(1,1) 
     >                                 + part_dx(1)*(i-0.5)
            part_y(2,totalParticles) = gridDomain(1,2)
     >                                 + part_dx(2)*(j-0.5)
            part_y(3,totalParticles) = gridDomain(1,3) 
     >                                 + part_dx(3)*(k-0.5)
!This can test out the remove particle function
!            IF(i .EQ. INT(np(1)/2) .AND. j .EQ. 1 .AND. k .EQ. 1)
!     >            part_y(1,totalParticles) = 2.0*gridDomain(2,1) 
          END DO
        END DO
      END DO
      particlesPerProc = INT(REAL(totalParticles)/REAL(nproc))
      npart_local = 0
      iPstart     = particlesPerProc*(nid  ) + 1 
      iPend       = particlesPerProc*(nid+1)
      IF(iPstart .LE. totalParticles) THEN
        IF(iPend .GT. totalParticles) iPend = totalParticles
        IF(nid  .EQ. nproc-1)         iPend = totalParticles
        DO i = iPstart,iPend 
          npart_local = npart_local + 1
          DO j = 1,3
            p_part_y(j,npart_local) = part_y(j,i)
          END DO
        END DO
      END IF

      CALL MPI_BARRIER(icomm, ierr)
      p_part_r = 0.0D0
      rhop   = 7730.0D0 ! steel particles
      DO i = 1,npart_local
        p_part_y(PPICLF_JVX,i) = 0.0D0
        p_part_y(PPICLF_JVY,i) = 0.0D0
        p_part_y(PPICLF_JVZ,i) = 0.0D0
        p_part_y(PPICLF_JT, i) = 0.0D0 ! particle temp
        p_part_y(PPICLF_JOX,i) = 0.0D0
        p_part_y(PPICLF_JOY,i) = 0.0D0
        p_part_y(PPICLF_JOZ,i) = 0.0D0
        p_part_r(PPICLF_R_JRHOP,i) = rhop ! particle density
        p_part_r(PPICLF_R_JDP,i)   = pdia ! particle diameter
        p_part_r(PPICLF_R_JVOLP,i) = (4.0D0/3.0D0)*PI
     >                              *(0.5D0*pdia)**3 ! particle volume
        p_part_r(PPICLF_R_JSPL,i) = 1.0D0 ! Super Particle Loading 
      END DO

      nndistTemp  = 4.0D0*pdia
      ! If the above is an even number, then it may have an issue due to
      ! precision error. 
      CALL MPI_Allreduce(nndistTemp,nndist,1,MPI_DOUBLE,
     >                                      MPI_MAX,iComm,ierr)

      CALL MPI_BARRIER(icomm,ierr)
! ppiclF Inputs and test case setup
!**********************************************************************

      IF(nproc .EQ. 1) THEN
        PRINT*, 'Number of grid cells in x direction:        ',
     >          nCells(1)
        PRINT*, 'Number of grid cells in y direction:        ',
     >          nCells(2)
        PRINT*, 'Number of grid cells in z direction:        ',
     >          nCells(3)
        PRINT*, 'Number of particles per cell in x direction:',
     >          particlesPerCell(1)
        PRINT*, 'Number of particles per cell in y direction:',
     >          particlesPerCell(2)        
        PRINT*, 'Number of particles per cell in z direction:',
     >          particlesPerCell(3)      
        PRINT*, ''
        PRINT*, 'Total Cells:    ', numCells
        PRINT*, 'Total Particles:', totalParticles 
        PRINT*, ''  
        PRINT*, 'Interpolation and Projection passes if the '
        PRINT*, 'average absolute error is less than 1.0D-5.'
        PRINT*, ''
      END IF

      IF(nid .EQ. 0)
     >  PRINT*,'****************************************************'
      IF(nid .EQ. 0)  PRINT*, 'Number of Processors:',nproc

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

        CALL test_PrintBanner(1,nid,nproc,testcase)
        CALL MPI_BARRIER(icomm,ierr)
 
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
        CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JT,tpF)
        CALL ppiclf_solve_InitSolve
        CALL MPI_BARRIER(icomm,ierr)
        DO ie = 1,proc_ncells
          CALL ppiclf_solve_GetProFld(ie,1,feedback1(ie))
          CALL ppiclf_solve_GetProFld(ie,2,feedback2(ie))
        END DO
        CALL MPI_BARRIER(icomm,ierr)

! Calculate Interpolation results 
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
        eps = 1.0D-60 
        IF(ppiclf_npart .GT. 0) THEN
          T_analytic = 0.0D0
          calcErr    = 0.0D0
          DO i = 1,ppiclf_npart
            xFace = .FALSE.
            yFace = .FALSE.
            zFace = .FALSE.

            xp = ppiclf_y(PPICLF_JX,i)
            yp = ppiclf_y(PPICLF_JY,i)
            zp = ppiclf_y(PPICLF_JZ,i) 

            IF( (xp - gridDomain(1,1)) .LT. gridDX(1) .OR.
     >          (gridDomain(2,1) - xp) .LT. gridDX(1) ) xFace = .TRUE.
            IF( (yp - gridDomain(1,2)) .LT. gridDX(2) .OR.
     >          (gridDomain(2,2) - yp) .LT. gridDX(2) ) yFace = .TRUE.
            IF( (zp - gridDomain(1,3)) .LT. gridDX(3) .OR.
     >          (gridDomain(2,3) - zp) .LT. gridDX(3) ) zFace = .TRUE.

            ! First find "calculation error", based on a inverse distance 
            ! interpolation with surrounding cells (~27 cells)
            ! This assumes a constant grid dx per direction!!
            part_cell(1) = CEILING(xp/gridDX(1))
            part_cell(2) = CEILING(yp/gridDX(2))
            part_cell(3) = CEILING(zp/gridDX(3))
            xstart = -1
            xend   = 1
            IF(.NOT.ppiclf_linperiodic(1) .AND. xFace) THEN
              IF(part_cell(1) .LE. 1)         xstart = 0
              IF(part_cell(1) .GE. nCells(1)) xend   = 0 
            END IF
            ystart = -1
            yend   = 1
            IF(.NOT.ppiclf_linperiodic(2) .AND. yFace) THEN
              IF(part_cell(2) .LE. 1)         ystart = 0
              IF(part_cell(2) .GE. nCells(2)) yend   = 0 
            END IF
            zstart = -1
            zend   = 1
            IF(.NOT.ppiclf_linperiodic(3) .AND. zFace) THEN
              IF(part_cell(3) .LE. 1)         zstart = 0
              IF(part_cell(3) .GE. nCells(3)) zend   = 0 
            END IF
            weightCell = 0.0D0
            weightTot  = 0.0D0
            T_calc     = 0.0D0
            DO loopcount = 1,2
              icount = 0
              DO j = xstart,xend
                DO k = ystart,yend
                  DO l = zstart,zend
                    icount = icount + 1
                    xcell  = (part_cell(1) - 0.5D0 + j)*gridDX(1) 
                    ycell  = (part_cell(2) - 0.5D0 + k)*gridDX(2)
                    zcell  = (part_cell(3) - 0.5D0 + l)*gridDX(3)
                    IF(loopcount .EQ. 1) THEN
                      weightCell(icount) = 1.0D0/
     >                                    (SQRT((xp-xcell)**2 + 
     >                                          (yp-ycell)**2 +
     >                                          (zp-zcell)**2 ) + eps)
                      weightTot  = weightCell(icount) + weightTot 
                    END IF
                    IF(loopcount .EQ. 2) THEN
                      x_norm  = (xcell - gridDomain(1,1))
     >                          /(gridDomain(2,1)-gridDomain(1,1))
                      y_norm  = (ycell - gridDomain(1,2))
     >                          /(gridDomain(2,2)-gridDomain(1,2))
                      z_norm  = (zcell - gridDomain(1,3))
     >                          /(gridDomain(2,3)-gridDomain(1,3))
                      T_analytic = COS(2*PI*x_norm) +
     >                             COS(2*PI*y_norm) +
     >                             COS(2*PI*z_norm)
                      T_calc = T_calc + 
     >                        (weightCell(icount)/weightTot)*T_analytic
                    END IF
                  END DO !l
                END DO !k
              END DO !j
            END DO !loopcount

            ! All particles
            calcErr = calcErr + 
     >          ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - ABS(T_calc))
            totCnt = totCnt + 1.0D0

            ! Particles in Interior cells
            IF(.NOT. xFace .AND. .NOT. yFace .AND. .NOT. zFace) THEN
              InteriorErr = InteriorErr +
     >               ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - ABS(T_calc))
              InteriorCnt = InteriorCnt + 1.0D0

            ! Particles in single Face cells
            ELSE IF(xFace .AND. .NOT. yFace .AND. .NOT. zFace) THEN
              xFaceErr = xFaceErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - ABS(T_calc))
              xFaceCnt = xFaceCnt + 1.0D0
            ELSE IF(.NOT. xFace .AND. yFace .AND. .NOT. zFace) THEN
              yFaceErr = yFaceErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - ABS(T_calc))
              yFaceCnt = yFaceCnt + 1.0D0
            ELSE IF(.NOT. xFace .AND. .NOT. yFace .AND. zFace) THEN
              zFaceErr = zFaceErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - ABS(T_calc))
              zFaceCnt = zFaceCnt + 1.0D0

            ! Particles in Two Face (Edge) cells
            ELSE IF(xFace .AND. yFace .AND. .NOT. zFace) THEN
              xyEdgeErr = xyEdgeErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - ABS(T_calc))
              xyEdgeCnt = xyEdgeCnt + 1.0D0
            ELSE IF(xFace .AND. .NOT. yFace .AND. zFace) THEN
              xzEdgeErr = xzEdgeErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - ABS(T_calc))
              xzEdgeCnt = xzEdgeCnt + 1.0D0
            ELSE IF(.NOT. xFace .AND. yFace .AND. zFace) THEN
              yzEdgeErr = yzEdgeErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - ABS(T_calc))
              yzEdgeCnt = yzEdgeCnt + 1.0D0

            ! Particles in Three Face (Corner) cells
            ELSE IF(xFace .AND. yFace .AND. zFace) THEN
              xyzCornerErr = xyzCornerErr + 
     >             ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i)) - ABS(T_calc))
              xyzCornerCnt = xyzCornerCnt + 1.0D0
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

        ! Only print out error results for x,y,z periodicitiy on.
        ! This shows all periodicity and interior particles.
        IF(test .EQ. 8) THEN
          IF(nproc .EQ. 1) THEN
            PRINT*,''
            PRINT*,'Linear Periodicity Implemented in x, y, & z'
            PRINT*,'Particle average interpolation errors:'
            PRINT*,''
            PRINT*,'All particles:      ' , calcErr/totCnt
            PRINT*,'Interior particles: ' , InteriorErr/InteriorCnt
            PRINT*,'xFace particles:    ' , xFaceErr/xFaceCnt
            PRINT*,'yFace particles:    ' , yFaceErr/yFaceCnt
            PRINT*,'zFace particles:    ' , zFaceErr/zFaceCnt
            PRINT*,'xyEdge particles:   ' , xyEdgeErr/xyEdgeCnt
            PRINT*,'xzEdge particles:   ' , xzEdgeErr/xzEdgeCnt
            PRINT*,'yzEdge particles:   ' , yzEdgeErr/yzEdgeCnt
            PRINT*,'xyzCorner particles:' , xyzCornerErr/xyzCornerCnt
            PRINT*,''
          END IF
        END IF
        IF(calcErr/totCnt .LT. 1.0D-5) interpolation(test) = .TRUE.
        CALL MPI_ALLREDUCE(interpolation(test), interp_logical(test), 1,
     >                     MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
        IF(nid .EQ. 0) THEN
          IF(.NOT. interp_logical(test)) PRINT*, 
     >         'Interpolation failed for:',testcase
        END IF
        CALL MPI_BARRIER(icomm,ierr)

! Calculate Projection Results
!**********************************************************************
        ! Find true projection result for all cells on each processor
        TrueFeedback1 = 0.0D0
        TrueFeedback2 = 0.0D0
        IF(nproc .EQ. 1) THEN
!          filename = TRIM(testcase) // '_' //'Solution_Feedback.txt'
!          OPEN(UNIT=499,FILE=filename, STATUS='REPLACE',ACTION='WRITE')
!          WRITE(499,*) 'Cell ID      ',
!     >                 'x_centroid                ',
!     >                 'y_centroid                ',
!     >                 'z_centroid                ',' ',
!     >                 'Unity Absolute Error      ',
!     >                 'Sine Absolute Error       ',
!     >                 'Unity Solution            ',
!     >                 'Sine Solution'
        END IF

        ! This assumes uniform grid dx per dimension!!
        DO ip = 1,totalParticles
          w     = 0.0D0
          wsum  = 0.0D0
!          ! This takes care of particle that was intentionally removed
!          IF(part_y(1,ip) .GT. gridDomain(2,1)) CYCLE

          ! Loop to find individual cell weightings
          DO ie = 1,numCells
            dSQi = 0.0D0
            dSQl = 0.0D0
            farAway = .FALSE.
            DO l = 1,3
              IF(ppiclf_linperiodic(l) .AND. 
     >                             ppiclf_EqualDomain(l)) THEN
                dSQl = MIN( (grid(l,ie) - part_y(l,ip))**2, 
     >           ( (gridDomain(2,l) - gridDomain(1,l)) -
     >                 ABS(grid(l,ie) - part_y(l,ip)) )**2 )
              ELSE
                dSQl = (grid(l,ie) - part_y(l,ip))**2
              END IF
              IF(dSQl .GT. (gridDX(l)*1.50001D0)**2) farAway = .TRUE.
              dSQi = dSQi + dSQl
            END DO !l
            IF(farAway) CYCLE
            dist = SQRT(dSQi)
            CellVol = grid(7,ie)
            GaussianConst = 2.305D0
            w(ie) = ABS(CellVol*EXP(-GaussianConst*(dist**2)
     >                / (CellVol**(2.0D0/3.0D0))))
            wsum = wsum + w(ie)
          END DO !ie

          x_norm = (part_y(1,ip) - gridDomain(1,1))
     >          / (gridDomain(2,1) - gridDomain(1,1))
          y_norm = (part_y(2,ip) - gridDomain(1,2))
     >         /  (gridDomain(2,2) - gridDomain(1,2))
          z_norm = (part_y(3,ip) - gridDomain(1,3))
     >         /  (gridDomain(2,3) - gridDomain(1,3))

          part_feedbk1(ip) = 1.0
          part_feedbk2(ip) = SIN(2*PI*x_norm) + SIN(2*PI*y_norm)
     >                     + SIN(2*PI*z_norm)
          DO ie = 1, numCells
            TrueFeedback1(ie) = TrueFeedback1(ie) + 
     >                         (w(ie)/wsum) * part_feedbk1(ip)
            TrueFeedback2(ie) = TrueFeedback2(ie) +
     >                         (w(ie)/wsum) * part_feedbk2(ip)
          END DO !ie
        END DO !ip
        ! Now only loop through cells that this processor
        ! input to ppiclF
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
!            IF(nproc .EQ. 1) THEN
!              WRITE(499,*) ie, grid(1,ie), grid(2,ie), grid(3,ie),
!     >                     error1, error2, 
!     >                     TrueFeedback1(iCstart-1+ie),
!     >                     TrueFeedback2(iCstart-1+ie) 
!            END IF
          END DO
!          IF(nproc .EQ. 1) CLOSE(UNIT=499)

          e1avg = e1avg / proc_ncells
          e2avg = e2avg / proc_ncells
          IF(e1avg .LT. 1.0D-5 .AND. e2avg .LT. 1.0D-5)
     >      projection(test) = .TRUE.
        ELSE
          projection(test) = .TRUE.
        END IF
        CALL MPI_ALLREDUCE(projection(test), proj_logical(test), 1,
     >                     MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
        e1avg = e1avg * proc_ncells
        e2avg = e2avg * proc_ncells
        projCells = proc_ncells

        ! Print out average feedback error for x,y,z periodic case
        IF(test .EQ. 8) THEN
          IF(nproc .EQ. 1) THEN
            CALL MPI_Allreduce(e1avg,e1avg,1,MPI_DOUBLE,
     >                                          MPI_SUM,iComm,ierr)
            CALL MPI_Allreduce(e2avg,e2avg,1,MPI_DOUBLE,
     >                                          MPI_SUM,iComm,ierr)
            CALL MPI_Allreduce(projCells,projCells,1,MPI_INTEGER,
     >                                          MPI_SUM,iComm,ierr)
            e1avg = e1avg / projCells
            e2avg = e2avg / projCells
            PRINT*,'Linear Periodicity Implemented in x, y, & z'
            PRINT*,'Cell average projection errors:'
            PRINT*,''          
            PRINT*, 'Unity projection:', e1avg
            PRINT*, 'Sine projeciton:', e2avg
            PRINT*, ''
          END IF
        END IF
        IF(nid .EQ. 0) THEN
          IF(.NOT. proj_logical(test)) PRINT*, 
     >         'Projection failed for:',testcase
        END IF

        CALL MPI_BARRIER(icomm,ierr)
!********************************************************************** 
! End periodicity loop testing
      END DO !test

! Particle to Particle Nearest Neighbor Test
!**********************************************************************

! Add random particles to ensure results are robust
          DO i = 1,1000
            totalParticles = totalParticles + 1
            CALL RANDOM_NUMBER(randNum)
            part_y(1,totalParticles) = gridDomain(1,1) 
     >              + (gridDomain(2,1) - gridDomain(1,1))*randNum
            CALL RANDOM_NUMBER(randNum)
            part_y(2,totalParticles) = gridDomain(1,2)
     >              + (gridDomain(2,2) - gridDomain(1,2))*randNum
            CALL RANDOM_NUMBER(randNum)
            part_y(3,totalParticles) = gridDomain(1,3) 
     >              + (gridDomain(2,3) - gridDomain(1,3))*randNum
          END DO
        IF(nid .EQ. nproc-1) THEN
          DO i = totalParticles-999,totalParticles 
              npart_local = npart_local + 1
              DO j = 1,3
                p_part_y(j,npart_local) = part_y(j,i)
              END DO
          END DO

          p_part_r = 0.0D0
          rhop   = 7730.0D0 ! steel particles
          DO i = 1,npart_local
            p_part_y(PPICLF_JVX,i) = 0.0D0
            p_part_y(PPICLF_JVY,i) = 0.0D0
            p_part_y(PPICLF_JVZ,i) = 0.0D0
            p_part_y(PPICLF_JT, i) = 0.0D0 
            p_part_y(PPICLF_JOX,i) = 0.0D0
            p_part_y(PPICLF_JOY,i) = 0.0D0
            p_part_y(PPICLF_JOZ,i) = 0.0D0
            p_part_r(PPICLF_R_JRHOP,i) = rhop ! particle density
            p_part_r(PPICLF_R_JDP,i)   = pdia ! particle diameter
            p_part_r(PPICLF_R_JVOLP,i) = (4.0D0/3.0D0)*PI
     >                                  *(0.5D0*pdia)**3 ! particle volume
            p_part_r(PPICLF_R_JSPL,i) = 1.0D0 ! Super Particle Loading 
          END DO
        END IF

! Start loop for varying periodicity cases
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

        CALL MPI_BARRIER(icomm,ierr)
! ppiclF Calls
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
        CALL ppiclf_comm_CreateBin
        CALL ppiclf_comm_FindParticle
        CALL ppiclf_comm_MoveParticle
        CALL ppiclf_comm_CreateGhost
        CALL ppiclf_comm_MoveGhost
        CALL MPI_BARRIER(icomm,ierr)

        IF(ppiclf_npart .GT. 0) THEN 
          ! SetYdot only runs with nearest neighbor subroutine
          ! due to preprocess compile flag.  PARTICLE_NN(i) saves
          ! the number of nearest neighbors per particle in ppiclF.
          CALL ppiclf_user_SetYdot
 
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
              IF(j .EQ. (id1+id2*particlesPerProc)) CYCLE 
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
            IF(NNCount .NE. PARTICLE_NN(i) .OR.
     >         (ABS(NNDistSQ - PPICLF_TOTNNDIST(i))
     >                            /  ABS(NNDistSQ)) .GT. 1D0) THEN
              nnpart(test) = .FALSE.
!              PRINT*, 'Count Diff (ppiclf,ref):', PARTICLE_NN(i),NNCount
!              PRINT*, 'Dist SQ percent Error:', (ABS(NNDistSQ - 
!     >                 PPICLF_TOTNNDIST(i)) /  ABS(NNDistSQ))
              IF(PARTICLE_NN(i) .GT. NNCount) THEN
                PRINT*, 'diff/CountDiff:',
     >                  (PPICLF_TOTNNDIST(i)-NNDistSQ)
     >                  /(PARTICLE_NN(i)-NNCount) ,
     >                  'nndist^2:',ppiclf_nndist**2
                PRINT*, PARTICLE_NN(i),NNCount
              ELSE
                PRINT*, 'diff/CountDiff:',
     >                  (- PPICLF_TOTNNDIST(i)+NNDistSQ)
     >                  /(-PARTICLE_NN(i)+NNCount) ,
     >                  'nndist^2:',ppiclf_nndist**2
                PRINT*, NNCount, PARTICLE_NN(i)
              END IF
              EXIT
            END IF
          END DO
        END IF
        ! Ensure nnpart(test) is true across all processors
        CALL MPI_ALLREDUCE(nnpart(test), nn_logical(test), 1,
     >                     MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(icomm,ierr)
      END DO !test
! Test CreateBin variations
!********************************************************************** 
      IF(nproc .EQ. 1) THEN
        filename = 'CreateBin_Results.txt'
        OPEN(UNIT=100,FILE=filename, STATUS='REPLACE',ACTION='WRITE')
        WRITE(100,*) 'Number of Processors, x bins, y bins, z bins,',
     >                ', total bins, Percent of Processors In Use:' 
        DO i = 1,125
          ppiclf_np = i
          CALL ppiclf_comm_CreateBin
          numBins = ppiclf_n_bins(1)*
     >               ppiclf_n_bins(2)*ppiclf_n_bins(3)
          WRITE(100,*) ppiclf_np, 
     >       ppiclf_n_bins(1), ppiclf_n_bins(2), ppiclf_n_bins(3),
     >       numBins , numBins/ppiclf_np*100
        END DO
        CLOSE(UNIT=100)
        binGen = .TRUE. ! only criteria is the program didn't crash.
        PRINT*,'Bin Stress test: (only ran on single processor case)'
        IF(binGen) THEN
          PRINT*,' CreateBin - PASSED'
        ELSE
          PRINT*,' CreateBin - FAILED'
        END IF
      END IF
! Print final results & close out program
!********************************************************************** 
      IF(nproc .GT. 1) THEN
        par = ' PARALLEL'
      ELSE
        par = ''
      END IF
      IF(nid .EQ. 0) THEN
        IF(nproc .EQ. 1) THEN
          PRINT*, ''
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

          IF(nproc .EQ. 1) THEN
            IF(nn_logical(1)) THEN
              PRINT*, ' PARTICLE NEAREST NEIGHBOR SEARCH - PASSED'
            ELSE
              PRINT*, ' PARTICLE NEAREST NEIGHBOR SEARCH - FAILED'
            END IF
          ELSE
            IF(nn_logical(1)) THEN
              PRINT*, TRIM(par) // ' PARTILCE NN SEARCH ',
     >                             'WITH GHOST PARTICLES - PASSED'
            ELSE
              PRINT*, TRIM(par) // ' PARTILCE NN SEARCH ',
     >                             'WITH GHOST PARTICLES - FAILED'
            END IF

          END IF

          PRINT*,''
          PRINT*, 'Linearly Periodic Results:'

          IF(interp_logical(2) .AND.
     >       interp_logical(3) .AND.
     >       interp_logical(4) .AND.
     >       interp_logical(5) .AND.
     >       interp_logical(6) .AND.
     >       interp_logical(7) .AND.
     >       interp_logical(8)      ) THEN
            PRINT*, TRIM(par) // ' INTERPOLATION - PASSED'
          ELSE
            PRINT*, TRIM(par) // ' INTERPOLATION - FAILED'
          END IF

          IF(proj_logical(2) .AND.
     >       proj_logical(3) .AND.
     >       proj_logical(4) .AND.
     >       proj_logical(5) .AND.
     >       proj_logical(6) .AND.
     >       proj_logical(7) .AND.
     >       proj_logical(8)      ) THEN
            PRINT*, TRIM(par) // ' PROJECTION - PASSED'
          ELSE
            PRINT*, TRIM(par) // ' PROJECTION - FAILED'
          END IF

          IF(nn_logical(2) .AND.
     >       nn_logical(3) .AND.
     >       nn_logical(4) .AND.
     >       nn_logical(5) .AND.
     >       nn_logical(6) .AND.
     >       nn_logical(7) .AND.
     >       nn_logical(8)      ) THEN
            PRINT*, TRIM(par) // ' PARTILCE NN SEARCH WITH '
     >                          ,'GHOST PARTICLES - PASSED'
          ELSE
            PRINT*, TRIM(par) // ' PARTILCE NN SEARCH WITH '
     >                          ,'GHOST PARTICLES - FAILED'
          END IF
          PRINT*,''
        END IF
      
        IF(interp_logical(1) .AND.
     >     interp_logical(2) .AND.
     >     interp_logical(3) .AND.
     >     interp_logical(4) .AND.
     >     interp_logical(5) .AND.
     >     interp_logical(6) .AND.
     >     interp_logical(7) .AND.
     >     interp_logical(8) .AND.
     >     proj_logical(1)  .AND.
     >     proj_logical(2)  .AND.
     >     proj_logical(3)  .AND.
     >     proj_logical(4)  .AND.
     >     proj_logical(5)  .AND.
     >     proj_logical(6)  .AND.
     >     proj_logical(7)  .AND.
     >     proj_logical(8)  .AND.
     >     nn_logical(1)    .AND.
     >     nn_logical(2)    .AND.
     >     nn_logical(3)    .AND.
     >     nn_logical(4)    .AND.
     >     nn_logical(5)    .AND.
     >     nn_logical(6)    .AND.
     >     nn_logical(7)    .AND.
     >     nn_logical(8)         ) THEN
           PRINT*, 
     >   'ALL TESTING PASSED - number of processors:'
     >   ,nproc
        ELSE
          PRINT*, 'TEST FAILED - number of processors:', nproc
          PRINT*,''
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

          IF(nproc .EQ. 1) THEN
            IF(nn_logical(1)) THEN
              PRINT*, ' PARTICLE NEAREST NEIGHBOR SEARCH - PASSED'
            ELSE
              PRINT*, ' PARTICLE NEAREST NEIGHBOR SEARCH - FAILED'
            END IF
          ELSE
            IF(nn_logical(1)) THEN
              PRINT*, TRIM(par) // ' PARTILCE NN SEARCH ',
     >                             'WITH GHOST PARTICLES - PASSED'
            ELSE
              PRINT*, TRIM(par) // ' PARTILCE NN SEARCH ',
     >                             'WITH GHOST PARTICLES - FAILED'
            END IF

          END IF

          PRINT*,''
          PRINT*, 'Linearly Periodic Results:'

          IF(interp_logical(2) .AND.
     >       interp_logical(3) .AND.
     >       interp_logical(4) .AND.
     >       interp_logical(5) .AND.
     >       interp_logical(6) .AND.
     >       interp_logical(7) .AND.
     >       interp_logical(8)      ) THEN
            PRINT*, TRIM(par) // ' INTERPOLATION - PASSED'
          ELSE
            PRINT*, TRIM(par) // ' INTERPOLATION - FAILED'
          END IF

          IF(proj_logical(2) .AND.
     >       proj_logical(3) .AND.
     >       proj_logical(4) .AND.
     >       proj_logical(5) .AND.
     >       proj_logical(6) .AND.
     >       proj_logical(7) .AND.
     >       proj_logical(8)      ) THEN
            PRINT*, TRIM(par) // ' PROJECTION - PASSED'
          ELSE
            PRINT*, TRIM(par) // ' PROJECTION - FAILED'
          END IF

          IF(nn_logical(2) .AND.
     >       nn_logical(3) .AND.
     >       nn_logical(4) .AND.
     >       nn_logical(5) .AND.
     >       nn_logical(6) .AND.
     >       nn_logical(7) .AND.
     >       nn_logical(8)      ) THEN
            PRINT*, TRIM(par) // ' PARTILCE NN SEARCH WITH '
     >                          ,'GHOST PARTICLES - PASSED'
          ELSE
            PRINT*, TRIM(par) // ' PARTILCE NN SEARCH WITH '
     >                          ,'GHOST PARTICLES - FAILED'
          END IF
          PRINT*,''
        END IF
      END IF

      CALL MPI_FINALIZE(ierr)
      CALL test_PrintBanner(0,nid,nproc,testcase)
      END PROGRAM

!----------------------------------------------------------------------

      SUBROUTINE test_PrintBanner(i,nid,nproc,tc)
      ! Input/Output
      ! i  - 1:Start of test, 2:End of test
      ! nid - Processor ID
      ! nproc - Number of Processors
      ! tc - periodicity test case
  
      IMPLICIT NONE

      INTEGER*4 i, nid, nproc
      CHARACTER*50 tc

      IF(i .EQ. 1) THEN
        IF(nid .EQ. 0) THEN
          !PRINT*, '*** Test Case:',TRIM(tc), ' ***'
        END IF
      ELSE IF(i .EQ. 0) THEN
        IF(nid .EQ. 0) THEN
          PRINT*, 'ppiclF test run completed for Processors:',nproc
          PRINT*, '****************************************************'
          PRINT*, ''
        END IF
      ELSE
        PRINT*, 'PrintBanner error'
      END IF

      RETURN
      END SUBROUTINE
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

 
