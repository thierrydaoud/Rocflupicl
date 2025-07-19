      PROGRAM main

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF'

      INTEGER*4 i, j, k, nproc, nid, icomm, ierr

      ! Grid variables
      INTEGER*4 nCells
      REAL*8    grid(7,8), gridDomain(2,3), xMin(3),xMax(3)

      ! Particle variables
      REAL*8    part_y(3,30000),pdia,dx_part,dx_ratio, x0(3)
      INTEGER*4 nPcells

      CALL MPI_Init(ierr)
      icomm = MPI_COMM_WORLD
      CALL MPI_Comm_size(icomm,nproc,ierr)
      CALL MPI_Comm_rank(icomm,nid,ierr)  

      IF(nid .EQ. 0) THEN
        PRINT*, ''
        PRINT*, '*****************************************************'
        PRINT*, 'ppiclF test run starting'
        PRINT*, 'Number of Processors:',nproc
        PRINT*, ''
      END IF

      ! Create uniform rectangular grid
      DO i = 1,3 !3D uniform grid
        gridDomain(1,i) = 0.0D-2 !Min
        gridDomain(2,i) = 1.0D-2 !Max
      END DO
      nCells = 8 !Per processor & must be equal to grid() Num columns
      CALL test_CreateGrid(gridDomain,nCells,nid,nproc,grid)
      CALL MPI_BARRIER(icomm,ierr)

      ! Create Particle Pack positions
      dx_ratio = 10.0D0 !cell_dx/part_dx = 10
      pdia = grid(4,1)/dx_ratio ! 10**3=1000 particles per cell (grid(4,1) = cell dx)
      dx_part = pdia ! particles are not overlapping
      !Numbers of cells to fill with particles per dimension = 3
      nPcells = 3
      !27,007 particles total
      CALL test_CreateParticles(gridDomain,dx_ratio,pdia,dx_part
     >                          ,nPcells,part_y)


      CALL ppiclf_comm_InitMPI(icomm, nid, nproc)

      IF(nid .EQ. 0) THEN
        PRINT*, ''
        PRINT*, 'ppiclF test run completed'
        PRINT*, '*****************************************************'
        PRINT*, ''
      END IF
      CALL MPI_BARRIER(icomm,ierr)
      CALL MPI_FINALIZE(ierr)

      END PROGRAM


!!*** START grid test ***
!      DO i = 1,3
!        xMin(i) = 1.0D9
!        xMax(i) = 1.0D-9
!      END DO
!      DO j = 1,nCells
!        DO i = 1,3
!          IF(grid(i,j) .LT. xMin(i)) xMin(i) = grid(i,j)
!          IF(grid(i,j) .GT. xMax(i)) xMax(i) = grid(i,j)
!        END DO
!      END DO
!      DO i = 1,3 !shift from min/max center to min/max face     
!        xMin(i) = xMin(i) - grid(3+i,1)/2
!        xMax(i) = xMax(i) + grid(3+i,1)/2
!      END DO
!      IF(nid .EQ. 0) THEN
!        PRINT*, 'dx, dy, dz, vol:', grid(4,1), grid(5,1),
!     >            grid(6,1), grid(7,1)
!        PRINT*, 'Input Domain:', gridDomain(1,1), gridDomain(2,1)
!      END IF
!      CALL MPI_BARRIER(icomm,ierr)
!      DO i = 1,1
!        PRINT*, 'Proc, xMin, xMax:', nid, i, xMin(i), xMax(i)
!      END DO 
!      CALL MPI_BARRIER(icomm,ierr)
!!*** END grid test ***


