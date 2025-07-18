      PROGRAM main

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF'

      INTEGER*4 i, j, k, nproc, nid, icomm, ierr

      ! Grid variables
      INTEGER*4 nCells
      REAL*8    grid(7,100), gridDomain(2,3), xMin(3),xMax(3)

      ! Particle variables
      !REAL*8 

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

      ! This creates a uniform rectangular grid
      DO i = 1,3 !3D uniform grid
        gridDomain(1,i) = 0.0D0 !Min
        gridDomain(2,i) = 10.0D0 !Max
      END DO
      nCells = 1000 !Per processor & must be equal to grid() Num columns
      CALL test_CreateGrid(gridDomain,nCells,nid,nproc,grid)
 
      ! START grid test
      DO i = 1,3
        xMin(i) = 1.0D9
        xMax(i) = 1.0D-9
      END DO
      DO j = 1,nCells
        DO i = 1,3
          IF(grid(i,j) .LT. xMin(i)) xMin(i) = grid(i,j)
          IF(grid(i,j) .GT. xMax(i)) xMax(i) = grid(i,j)
        END DO
      END DO     
      CALL MPI_BARRIER(icomm,ierr)
      DO i = 1,3
        PRINT*, 'Proc, Dim, xMin, xMax:', nid, i, xMin(i), xMax(i)
      END DO 
      ! END grid test

      CALL ppiclf_comm_InitMPI(icomm, nid, nproc)
      ! Create particle array

      IF(nid .EQ. 0) THEN
        PRINT*, ''
        PRINT*, 'ppiclF test run completed'
        PRINT*, '*****************************************************'
        PRINT*, ''
      END IF
      CALL MPI_FINALIZE(ierr)

      END PROGRAM
