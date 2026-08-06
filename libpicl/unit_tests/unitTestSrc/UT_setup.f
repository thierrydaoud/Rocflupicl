#include"../../ppiclF/source/PPICLF_USER.h"
#include"../../ppiclF/source/PPICLF_STD.h"

!**********************************************************************
! This setup_Unit Test should be called at the beginning of each unit
! test.  It creates a rectangular grid and random particle domain. 
! It ensures particles are located at each domain face for linear
! periodicity.  It also applies a temperature field to the cells 
! for the interpolation unit test.
 
      SUBROUTINE UT_setup

      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF_UNIT_TEST'

      ! local variables
      INTEGER*4 i, j, k, NumXCells, NumYCells, NumZCells, faceParticles 
      REAL*8    xmin, xmax, ymin, ymax, zmin, zmax, nndistTemp

! Setup Inputs
!**********************************************************************
      ! DEFAULT VALUES
      ! Rectangular Grid Input
      NumXCells = 15  
      NumYCells = 45 
      NumZCells = 15 

      xmin = -1.0D0
      ymin = -1.0D0
      zmin = -1.0D0

      xmax =  1.0D0
      ymax =  1.0D0
      zmax =  1.0D0

      ! Particle Input
      totalParticles = 25000
      ! Particles to calculate periodicity error
      faceParticles = 0!4200
      
      ! Ensuring PPICLF array limits aren't violated
      IF(NumXCells*NumYCells*NumZCells .GT. PPICLF_LEE) THEN
        PRINT*, 'ERROR, PPICLF_LEE = ',PPICLF_LEE, 'Cells to test =',
     >          NumXCells*NumYCells*NumZCells
        CALL MPI_FINALIZE(ierr)
        STOP
      END IF

      IF(totalParticles .GT. PPICLF_LPART) THEN
        PRINT*, 'Tried to make too many particles!'
        PRINT*, 'PPICLF_LPART =',PPICLF_LPART,'TotalParticles=',
     >          totalParticles
        CALL MPI_FINALIZE(ierr)
        STOP
      END IF
      ppiclf_glnpart = totalParticles
! Grid Setup
!**********************************************************************
      ! Create rectangular grid
      IF(.NOT.gridBoundsDefined) THEN
        gridDomain(1,1) = xmin
        gridDomain(2,1) = xmax
      END IF
      nCells(1) = NumXCells 
      gridDX(1) = (gridDomain(2,1) - gridDomain(1,1))/REAL(nCells(1))

      IF(.NOT.gridBoundsDefined) THEN
        gridDomain(1,2) = ymin 
        gridDomain(2,2) = ymax 
      END IF
      nCells(2) = NumYCells
      gridDX(2) = (gridDomain(2,2) - gridDomain(1,2))/REAL(nCells(2))

      IF(.NOT.gridBoundsDefined) THEN
        gridDomain(1,3) = zmin
        gridDomain(2,3) = zmax
      END IF
      nCells(3) = NumZCells
      gridDX(3) = (gridDomain(2,3) - gridDomain(1,3))/REAL(nCells(3))

      ! Build full grid on each processor 
      ! Build in order different from bin numbering
      ! to ensure GSLIB calls are fully tested.
      numCells = 0
      DO i = 1,nCells(1)
        DO j = 1,nCells(2)
          DO k = 1,nCells(3)
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
        DO i = iCstart,iCend  
          proc_ncells = proc_ncells + 1
          DO j = 1,7
            p_grid(j,proc_ncells) = grid(j,i)
          END DO !j
        END DO !i
      END IF
      CALL MPI_BARRIER(icomm, ierr)

      ! Find cell filter search distance
      filter       = 1.0D-9 !dummy
      dx_min       = 1.0D9  !dummy
      nFilterCells = 1.5D0
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
      ! A good rule of thumb is particle diameter <= dx/5
      pdia = MIN(dx_min(1)/10.0D0, dx_min(2)/10.0D0, dx_min(3)/10.0D0)

      ! Build full random particle dispersion on root processor
      IF(nid .EQ. rootProc) THEN 
        DO i = 1,totalParticles
          DO j = 1,3 ! loop through x, y, & z coordinates
            CALL RANDOM_NUMBER(randNum)
            part_y(j,i) = (gridDomain(2,j)-gridDomain(1,j))
     >                    * randNum
     >                    + gridDomain(1,j)
            IF(part_y(j,i) .GT. gridDomain(2,j) .OR. 
     >         part_y(j,i) .LT. gridDomain(1,j)     ) THEN
              PRINT*, 'Particle printed outside domain:',j, part_y(j,i),
     >                'Max:', gridDomain(2,j), 'Min:', gridDomain(1,j)
            END IF
            ! IF's below put first faceParticles within 1/3 a particle
            ! diameter from a face for linear periodicity testing.
            IF(i .LE. faceParticles) THEN
              ! Choose face (k)
              k = 3 !zface
              IF(i .LT. INT(REAL(faceParticles*2.0/3.0)+0.49D0) + 1) 
     >          k = 2 ! yface
              IF(i .LT. INT(REAL(faceParticles*1.0/3.0)+0.49D0) + 1) 
     >          k = 1 ! xface
              ! only modify one coordinate to put on face 
              IF(j .EQ. k) THEN
                IF(i - (k-1)*INT(REAL(faceParticles/3.0)+0.49D0) .LT.
     >                       INT(REAL(faceParticles/6.0)+1.49D0)) THEN
                  ! 1/6 faceParticles on min face
                  part_y(j,i) = gridDomain(1,j) + (pdia/2.0)*randNum
                ELSE
                  ! 1/6 faceParticles on max face
                  part_y(j,i) = gridDomain(2,j) - (pdia/2.0)*randNum
                END IF
              END IF ! j == k
            END IF ! i < faceParticles
            ! Make sure that there are some corner particles!
            IF(i .GT. faceParticles .AND.
     >         i .LE. INT(REAL(faceParticles)*(1.0D0+1.0D0/6.0D0))) 
     >        part_y(j,i) = gridDomain(2,j) - pdia/2.0D0*randNum
          END DO ! j
        END DO ! i 
      END IF 


      ! Share part_y from root processor to all processors
      CALL MPI_BCAST(part_y,totalParticles*SIZE(part_y,DIM=1),
     >               MPI_DOUBLE,rootProc,MPI_COMM_WORLD,ierr)

      ! Distribute particles across processors for ppiclf input
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
      rhoP   = 7730.0D0 ! steel particles
      DO i = 1,npart_local
        p_part_y(PPICLF_JVX,i) = 0.0D0
        p_part_y(PPICLF_JVY,i) = 0.0D0
        p_part_y(PPICLF_JVZ,i) = 0.0D0
        p_part_y(PPICLF_JT, i) = 0.0D0 
        p_part_y(PPICLF_JOX,i) = 0.0D0
        p_part_y(PPICLF_JOY,i) = 0.0D0
        p_part_y(PPICLF_JOZ,i) = 0.0D0
        p_part_r(PPICLF_R_JRHOP,i) = rhoP ! particle density
        p_part_r(PPICLF_R_JDP,i)   = pdia ! particle diameter
        p_part_r(PPICLF_R_JVOLP,i) = (4.0D0/3.0D0)*PI
     >                              *(0.5D0*pdia)**3 ! particle volume
        p_part_r(PPICLF_R_JSPL,i) = 1.0D0 ! Super Particle Loading 
      END DO

      nndistTemp  = 2.0D0*pdia
      CALL MPI_Allreduce(nndistTemp,nndist,1,MPI_DOUBLE,
     >                                      MPI_MAX,iComm,ierr)

      CALL MPI_BARRIER(icomm,ierr)

      RETURN
      END SUBROUTINE
!**********************************************************************
