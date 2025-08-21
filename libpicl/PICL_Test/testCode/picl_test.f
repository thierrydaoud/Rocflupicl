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
     >          MapGrid(PPICLF_LEE), iend, istart, cellsPerProc

      REAL*8    p_grid(7,PPICLF_LEE), grid(7,PPICLF_LEE),
     >          gridDomain(2,3), gridDX(3), filter(3), 
     >          nFilterCells, tpF(PPICLF_LEE), numBins,
     >          dx_min(3), feedback1(PPICLF_LEE), feedback2(PPICLF_LEE),
     >          x_norm, y_norm, z_norm

      ! Particle variables
      REAL*8    part_y(PPICLF_LRS,PPICLF_LPART), pdia, 
     >          p_part_y(PPICLF_LRS,PPICLF_LPART), part_dx(3),
     >          p_part_r(PPICLF_LRP,PPICLF_LPART),T_truth(PPICLF_LPART),
     >          xp, yp, zp, i_err, p_err, totErr, InteriorErr, xFaceErr,
     >          yFaceErr, zFaceErr, xyEdgeErr, xzEdgeErr, yzEdgeErr,
     >          xyzCornerErr, totCnt, InteriorCnt,
     >          xFaceCnt, yFaceCnt, zFaceCnt, xyEdgeCnt, xzEdgeCnt,
     >          yzEdgeCnt, xyzCornerCnt, CalcErr, 

      INTEGER*4 particlesPerCell(3), particlesPerProc, npart_local,
     >          totalParticles, ip, np(3), part_cell(3)

      LOGICAL   xFace, yFace, zFace

      ! Projection variables
      REAL*8    wsum, dSQl, dSQi, dist, CellVol, GaussianConst,
     >          w(PPICLF_LEE),part_feedbk1(PPICLF_LPART),
     >          TrueFeedback1(PPICLF_LEE), TrueFeedback2(PPICLF_LEE),
     >          part_feedbk2(PPICLF_LPART), error1, error2

      CHARACTER*50 filename, testcase, procString


      PI = 4.0D0*ATAN(1.0) ! pi
      k  = 3.0D0 ! wave number

! MPI Setup
!**********************************************************************
      CALL MPI_Init(ierr)
      icomm = MPI_COMM_WORLD
      CALL MPI_Comm_size(icomm,nproc,ierr)
      CALL MPI_Comm_rank(icomm,nid,ierr)
      WRITE(procString, '(I0)') nid
      IF(nid .LT. 10) procString = '0' // TRIM(procString)

! Grid Setup
!**********************************************************************
      ! Create rectangular grid
      gridDomain(1,1) = 0.0D0 !x domain min
      gridDomain(2,1) = 1.0D0 !x domain max
      nCells(1)       = 5    !Number of x cells in domain
      gridDX(1) = (gridDomain(2,1) - gridDomain(1,1))/REAL(nCells(1))

      gridDomain(1,2) = 0.0D0 !y domain min
      gridDomain(2,2) = 1.0D0 !y domain max     
      nCells(2)       = 5    !Number of y cells in domain
      gridDX(2) = (gridDomain(2,2) - gridDomain(1,2))/REAL(nCells(2))

      gridDomain(1,3) = 0.0D0 !z domain min
      gridDomain(2,3) = 1.0D0 !z domain max     
      nCells(3)       = 5    !Number of z cells in domain
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
      istart = cellsPerProc*(nid  ) + 1
      iend   = cellsPerProc*(nid+1)
      proc_ncells = 0
      IF(istart .LE. numCells) THEN
        IF(iend .GT. numCells) iend = numCells
        IF(nid .EQ. nproc-1) iend = numCells
        filename = 'Grid_Proc_' // TRIM(procString) // '.txt'
        OPEN(UNIT=1,FILE=filename, STATUS='REPLACE', ACTION='WRITE')
        DO i = istart,iend  
          proc_ncells = proc_ncells + 1
          MapGrid(proc_ncells) = i
          DO j = 1,7
            p_grid(j,proc_ncells) = grid(j,i)
          END DO !j
          WRITE(1,*) i, p_grid(1,i), p_grid(2,i), p_grid(3,i),
     >                  p_grid(4,i), p_grid(5,i), p_grid(6,i),
     >                  p_grid(7,i)
        END DO !i
        CLOSE(UNIT=1)
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
      particlesPerCell(1) = 2.0 ! in x dimension
      particlesPerCell(2) = 2.0 ! in y dimension
      particlesPerCell(3) = 2.0 ! in z dimension

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
          END DO
        END DO
      END DO
     
      particlesPerProc = INT(REAL(totalParticles)/REAL(nproc)) 
      npart_local = 0
      istart     = particlesPerProc*(nid  ) + 1 
      iend       = particlesPerProc*(nid+1)
      IF(istart .LE. totalParticles) THEN
        IF(iend .GT. totalParticles) iend = totalParticles
        IF(nid  .EQ. nproc-1)        iend = totalParticles
        DO i = istart,iend 
          npart_local = npart_local + 1
          DO j = 1,3
            p_part_y(j,npart_local) = part_y(j,i)
          END DO
        END DO
      END IF

      CALL MPI_BARRIER(icomm, ierr)

      IF(nid .EQ. 0) PRINT*,'Total Particles:',totalParticles 

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

      nndistTemp  = 1.3D0*pdia
      CALL MPI_Allreduce(nndistTemp,nndist,1,MPI_DOUBLE,
     >                                      MPI_MAX,iComm,ierr)

      CALL MPI_BARRIER(icomm,ierr)
! ppiclF Inputs and test case setup
!**********************************************************************

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
        PPICLF_TEST    = .TRUE.
        PPICLF_PERTEST = .TRUE.
        PPICLF_OVERLAP = .FALSE.
        CALL ppiclf_solve_InitParticle(2,3,0,npart_local,
     >                                 p_part_y,p_part_r,filter,nndist)
        PPICLF_TEST = .TRUE.
        PPICLF_PERTEST = .TRUE.
        CALL ppiclf_solve_Initialize(x_per_flag, x_per_min, x_per_max,
     >                               y_per_flag, y_per_min, y_per_max, 
     >                               z_per_flag, z_per_min, z_per_max, 
     >                               ang_per_flag, ang_per_angle, 
     >                               ang_per_xangle, ang_per_rin,
     >                                                    ang_per_rout)
        PPICLF_TEST = .TRUE.
        PPICLF_PERTEST = .TRUE.
        CALL ppiclf_comm_InitOverlapGrid(proc_ncells,p_grid)
        CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JT,tpF)
        CALL ppiclf_solve_InitSolve
        CALL MPI_BARRIER(icomm,ierr)
        DO ie = 1,proc_ncells
          CALL ppiclf_solve_GetProFld(ie,1,feedback1(ie))
          CALL ppiclf_solve_GetProFld(ie,2,feedback2(ie))
        END DO
        CALL MPI_BARRIER(icomm,ierr)
        PRINT*, 'Number of ghost particles:',ppiclf_npart_gp
! Particle to Particle Nearest Neighbor Test
!**********************************************************************
        ! This only prints out particle to particle
        ! nearest neighbor results since PPICLF_TEST = .TRUE.
        CALL ppiclf_user_SetYdot
        ! This finds true solution
        IF(nproc .EQ. 1) THEN
          DO i = 1,ppiclf_npart
            DO j = 1,ppiclf_npart
              IF(i .EQ. j) CYCLE
              dSQi = 0.0D0
              dSQl = 0.0D0
              DO l = 1,3
                IF(ppiclf_linperiodic(l) .AND. 
     >                               ppiclf_EqualDomain(l)) THEN
                  dSQl = MIN( (ppiclf_y(l,i) - ppiclf_y(l,j))**2, 
     >             ( (gridDomain(2,l) - gridDomain(1,l)) -
     >                   ABS(ppiclf_y(l,i) - ppiclf_y(l,j)) )**2 )
                ELSE
                  dSQl = (ppiclf_y(l,i) - ppiclf_y(l,j))**2
                END IF
                dSQi = dSQi + dSQl
              END DO !l
              IF(dSQi .GT. ppiclf_nndist**2) CYCLE
              PRINT*,'SOL Particles close',i,j
              PRINT*, 'x:', ppiclf_y(1,i), ppiclf_y(1,j)
              PRINT*, 'y:', ppiclf_y(2,i), ppiclf_y(2,j)
              PRINT*, 'z:', ppiclf_y(3,i), ppiclf_y(3,j)
            END DO
          END DO
        END IF
        CALL MPI_BARRIER(icomm,ierr)

! Print Interpolation results 
!**********************************************************************
        ! Zero out errors and counters
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
  
        IF(ppiclf_npart .GT. 0) THEN

          filename = TRIM(testcase) // '_' // 'Interpolation_Proc_'
     >                  // TRIM(procString) // '.txt'
          OPEN(UNIT=300,FILE=TRIM(filename), STATUS='REPLACE',
     >                                      ACTION='WRITE')
          WRITE(300,*) 'Particle ID, x, y,',
     >                        'z, Interpolation Error (%).'
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
            IF(.NOT.ppiclf_linperiodic(1) .AND. .NOT. 
     >                               ppiclf_EqualDomain(1)) THEN

              IF(part_cell(1) .LE. 1)         xstart = 0
              IF(part_cell(1) .GE. nCells(1)) xend   = 0 
            END IF
            ystart = -1
            yend   = 1
            IF(.NOT.ppiclf_linperiodic(2) .AND. .NOT. 
     >                               ppiclf_EqualDomain(2)) THEN

              IF(part_cell(2) .LE. 1)         ystart = 0
              IF(part_cell(2) .GE. nCells(2)) yend   = 0 
            END IF
            zstart = -1
            zend   = 1
            IF(.NOT.ppiclf_linperiodic(3) .AND. .NOT. 
     >                               ppiclf_EqualDomain(3)) THEN

              IF(part_cell(3) .LE. 1)         zstart = 0
              IF(part_cell(3) .GE. nCells(3)) zend   = 0 
            END IF
            weightCell = 0.0D0
            weightTot  = 0.0D0
            DO loopcount = 1,2
              icount = 0
              DO j = xstart,xend
                DO k = ystart,yend
                  DO l = zstart,zend
                    icount = icount + 1
                    IF(loopcount .EQ. 1) THEN
                      xcell  = (part_cell(1) - 0.5D0 + j)*gridDX(1) 
                      ycell  = (part_cell(2) - 0.5D0 + k)*gridDX(2)
                      zcell  = (part_cell(3) - 0.5D0 + l)*gridDX(3)
                      IF(xcell .LT. 
                      weightCell(icount) = 1.0D0/(SQRT((xp-xcell)**2 + 
     >                                                 (yp-ycell)**2 +
     >                                                 (zp-zcell)**2))
                      weightTot  = weightCell(icount) + weightTot 
                    END IF
                    IF(loopcount .EQ. 2) THEN
                      x_norm  = (xcell - gridDomain(1,1))
     >                          /(gridDomain(2,1)-gridDomain(1,1))
                      y_norm  = (ycell - gridDomain(1,2))
     >                          /(gridDomain(2,2)-gridDomain(1,2))
                      z_norm  = (zcell - gridDomain(1,3))
     >                          /(gridDomain(2,3)-gridDomain(1,3))
                      T_truth = COS(2*PI*x_norm) +
     >                             COS(2*PI*y_norm) +
     >                             COS(2*PI*z_norm)
                      Tcalc = (weightCell(icount)/weightTot)*T_truth 
                    END IF
                  END DO
                END DO
              END DO
            END DO
            calcErr = ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i))-ABS(Tcalc))
            ! Find ppiclF interpolation error vs analytical funciton
            x_norm = (xp - gridDomain(1,1))
     >              /(gridDomain(2,1)-gridDomain(1,1))
            y_norm = (yp - gridDomain(1,2))
     >              /(gridDomain(2,2)-gridDomain(1,2))
            z_norm = (zp - gridDomain(1,3))
     >              /(gridDomain(2,3)-gridDomain(1,3))

            T_truth(i) = COS(2*PI*x_norm) +
     >                   COS(2*PI*y_norm) +
     >                   COS(2*PI*z_norm)

            IF(T_truth(i) .LT. 0.1) CYCLE

            i_err = ABS(ABS(ppiclf_rprop(PPICLF_R_JT,i))
     >                     - ABS(T_truth(i)))
            p_err = i_err!/ABS(T_truth(i))*100
            totErr = totErr + p_err 
            totCnt = totCnt + 1.0D0
            ! Interior cell particles
            IF(.NOT. xFace .AND. .NOT. yFace .AND. .NOT. zFace) THEN
              InteriorErr = InteriorErr + p_err
              InteriorCnt = InteriorCnt + 1.0D0

            ! Single Face cell particles
            ELSE IF(xFace .AND. .NOT. yFace .AND. .NOT. zFace) THEN
              xFaceErr = xFaceErr + p_err
              xFaceCnt = xFaceCnt + 1.0D0
            ELSE IF(.NOT. xFace .AND. yFace .AND. .NOT. zFace) THEN
              yFaceErr = yFaceErr + p_err
              yFaceCnt = yFaceCnt + 1.0D0
            ELSE IF(.NOT. xFace .AND. .NOT. yFace .AND. zFace) THEN
              zFaceErr = zFaceErr + p_err
              zFaceCnt = zFaceCnt + 1.0D0

            ! Two Face (Edge) cell particles
            ELSE IF(xFace .AND. yFace .AND. .NOT. zFace) THEN
              xyEdgeErr = xyEdgeErr + p_err
              xyEdgeCnt = xyEdgeCnt + 1.0D0
            ELSE IF(xFace .AND. .NOT. yFace .AND. zFace) THEN
              xzEdgeErr = xzEdgeErr + p_err
              xzEdgeCnt = xzEdgeCnt + 1.0D0
            ELSE IF(.NOT. xFace .AND. yFace .AND. zFace) THEN
              yzEdgeErr = yzEdgeErr + p_err
              yzEdgeCnt = yzEdgeCnt + 1.0D0

            ! Three Face (Corner) cell particles
            ELSE IF(xFace .AND. yFace .AND. zFace) THEN
              xyzCornerErr = xyzCornerErr + p_err
              xyzCornerCnt = xyzCornerCnt + 1.0D0
            ! Something went wrong if not caught yet
            ELSE
              PRINT*,'Interp boolean error. Particle not',
     >               'counted right. Postion:', xp, yp, zp
            END IF

            WRITE(300,*) i, xp, yp, zp, i_err, p_err,'%'
          END DO
          CLOSE(UNIT=300)
        END IF
        CALL MPI_BARRIER(icomm,ierr)
        ! Percent Error Sums
        CALL MPI_Allreduce(totErr,totErr,1,MPI_DOUBLE,
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
        ! Print Bin Data
        IF(nid .EQ. 0) THEN
          PRINT*,'Average Interpolation Error for all particles:'
     >           ,totErr/totCnt
          PRINT*,'Average Interpolation Error Interior particles:'
     >           ,InteriorErr/InteriorCnt
          PRINT*,'Average Interpolation Error xFace particles:'
     >           ,xFaceErr/xFaceCnt
          PRINT*,'Average Interpolation Error yFace particles:'
     >           ,yFaceErr/yFaceCnt
          PRINT*,'Average Interpolation Error zFace particles:'
     >           ,zFaceErr/zFaceCnt
          PRINT*,'Average Interpolation Error xyEdge particles:'
     >           ,xyEdgeErr/xyEdgeCnt
          PRINT*,'Average Interpolation Error xzEdge particles:'
     >           ,xzEdgeErr/xzEdgeCnt
          PRINT*,'Average Interpolation Error yzEdge particles:'
     >           ,yzEdgeErr/yzEdgeCnt
          PRINT*,'Average Interpolation Error xyzCorner particles:'
     >           ,xyzCornerErr/xyzCornerCnt
         END IF
        CALL MPI_BARRIER(icomm,ierr)
! Print Projection Results
!**********************************************************************
        ! Find true projection result
        IF(nproc .EQ. 1) THEN
          TrueFeedback1 = 0.0D0
          TrueFeedback2 = 0.0D0
          DO ip = 1,npart_local
            w     = 0.0D0
            wsum  = 0.0D0
            ! Loop to find individual cell weightings
            DO ie = 1,proc_ncells
              dSQi = 0.0D0
              dSQl = 0.0D0
              DO l = 1,3
                IF(ppiclf_linperiodic(l) .AND. 
     >                               ppiclf_EqualDomain(l)) THEN
                  dSQl = MIN( (p_grid(l,ie) - part_y(l,ip))**2, 
     >             ( (gridDomain(2,l) - gridDomain(1,l)) -
     >                   ABS(p_grid(l,ie) - part_y(l,ip)) )**2 )
                ELSE
                  dSQl = (p_grid(l,ie) - part_y(l,ip))**2
                END IF
                dSQi = dSQi + dSQl
              END DO !l

              dist = SQRT(dSQi)
              CellVol = p_grid(7,ie)
              GaussianConst = 2.305D0
              w(ie) = ABS(CellVol*EXP(-GaussianConst*(dist**2)
     >                  / (CellVol**(2.0D0/3.0D0))))
              wsum = wsum + w(ie)
            END DO !ie

            x_norm = (part_y(1,ip) - gridDomain(1,1))
     >            / (gridDomain(2,1) - gridDomain(1,1))
            y_norm = (part_y(2,ip) - gridDomain(1,2))
     >           /  (gridDomain(2,2) - gridDomain(1,2))
            z_norm = (part_y(3,ip) - gridDomain(1,3))
     >           /  (gridDomain(2,3) - gridDomain(1,3))

            part_feedbk1(ip) = 1.0
            part_feedbk2(ip) = SIN(2*PI*x_norm) + SIN(2*PI*y_norm)
     >                       + SIN(2*PI*z_norm)
            DO ie = 1,proc_ncells
              TrueFeedback1(ie) = TrueFeedback1(ie) + 
     >                           (w(ie)/wsum) * part_feedbk1(ip)
              TrueFeedback2(ie) = TrueFeedback2(ie) +
     >                           (w(ie)/wsum) * part_feedbk2(ip)
            END DO !ie
          END DO !ip

          filename = TRIM(testcase) // '_' //'FeedbackSolution.txt'
          OPEN(UNIT=499,FILE=filename, STATUS='REPLACE',ACTION='WRITE')
          WRITE(499,*) 'Cell ID, x_centroid, y_centroid, ',
     >                 'z_centroid, Unity Percent Error, ',
     >                 'SIN Percent Error, Unity Solution, ',
     >                 'SIN Solution'
          DO ie = 1,proc_ncells
            error1 = ABS(ABS(feedback1(ie)) - ABS(TrueFeedback1(ie)))
     >              / ABS(TrueFeedback1(ie)) * 100.0D0
            error2 = ABS(ABS(feedback2(ie)) - ABS(TrueFeedback2(ie)))
     >              / ABS(TrueFeedback2(ie)) * 100.0D0
            WRITE(499,*) ie, p_grid(1,ie), p_grid(2,ie),
     >                   p_grid(3,ie), error1, '%', error2, '%',
     >                   TrueFeedback1(ie), TrueFeedback2(ie) 
          END DO
          CLOSE(UNIT=499)
        END IF


        CALL MPI_BARRIER(icomm,ierr)
        IF(proc_ncells .GT. 0) THEN
          filename = TRIM(testcase) // '_' // 'Feedback_Proc_' //
     >                              TRIM(procString) // '.txt'
          OPEN(UNIT=400,FILE=filename, STATUS='REPLACE',ACTION='WRITE')
          WRITE(400,*) 'Cell ID, x_centroid, y_centroid,',
     >                 'z_centroid, Unity Feedback, SIN Feedback'
          DO ie = 1,proc_ncells
            WRITE(400,*) ie, p_grid(1,ie), p_grid(2,ie),
     >                   p_grid(3,ie), feedback1(ie), feedback2(ie)
          END DO
          CLOSE(UNIT=400)
        END IF
      END DO !test
! Test CreateBin variations
!********************************************************************** 
      IF(nproc .EQ. 1) THEN
        PRINT*,'******************************************************'
        PRINT*,'CreateBin Testing:'
!        DO j = 0,3
!          ! Tests 4 cases
!          ! a) bin x = 1.00; bin y = 1.00; bin z = 1.00
!          !    Number of bins equal in all direcitons.
!          ! b) bin x = 2.01; bin y = 1.00; bin z = 1.00
!          !    More bins in x.  Bins and y and z equal.
!          ! c) bin x = 2.01; bin y = 2.02; bin z = 1.00
!          !    More bins in x, then y, then z
!          ! d) bin x = 2.03; bin y = 2.01, bin z = 2.03
!          !    More ins in x or z.  less Bins in y
!          IF(j .NE. 0)
!     >     ppiclf_xdrange(2,j) = 2.0*gridDomain(2,j)*(j/100) 
!          IF(j .EQ. 3) ppiclf_xdrange(2,1) = ppiclf_xdrange(2,j)
!
!          WRITE(100+j,*) 'x domain:', ppiclf_xdrange(2,1)
!          WRITE(100+j,*) 'y domain:', ppiclf_xdrange(2,2)
!          WRITE(100+j,*) 'z domain:', ppiclf_xdrange(2,3)
!
          filename = 'CreateBin_Results.txt'
          OPEN(UNIT=100,FILE=filename, STATUS='REPLACE',ACTION='WRITE')
          WRITE(100,*) 'Number of Processors, x bins, y bins, z bins,',
     >                  ', total bins, Percent of Processors In Use:' 
          DO i = 1,512
            ppiclf_np = i
            numBins = ppiclf_n_bins(1)*
     >                 ppiclf_n_bins(2)*ppiclf_n_bins(3)
            CALL ppiclf_comm_CreateBin
            WRITE(100,*) ppiclf_np, 
     >         ppiclf_n_bins(1), ppiclf_n_bins(2), ppiclf_n_bins(3),
     >         numBins , numBins/ppiclf_np*100
          END DO
          CLOSE(UNIT=100)
!        END DO
      END IF


! Close out program
!********************************************************************** 
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
          PRINT*, ''
          PRINT*, '****************************************************'
          PRINT*, 'ppiclF test run starting'
          PRINT*, 'Number of Processors:',nproc
          PRINT*, 'Test Case:',tc
          PRINT*, ''
        END IF
      ELSE IF(i .EQ. 0) THEN
        IF(nid .EQ. 0) THEN
          PRINT*, ''
          PRINT*, 'ppiclF test run completed'
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
