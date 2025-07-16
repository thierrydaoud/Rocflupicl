      PROGRAM main

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF'

      INTEGER*4 test, np, nid, icomm, ierr
      REAL*8 dia_p, x(2), y(2), z(2), dx, dy, dz


      CALL MPI_Init(ierr)
      icomm = MPI_COMM_WORLD
      CALL MPI_Comm_size(icomm,np,ierr)
      CALL MPI_Comm_rank(icomm,nid,ierr)  

      IF(nid .EQ. 0) THEN
        PRINT*, ''
        PRINT*, '*****************************************************'
        PRINT*, 'ppiclF test run starting'
        PRINT*, 'Number of Processors:',np
        PRINT*, ''
       END IF
      CALL ppiclf_comm_InitMPI(icomm, nid, np)

      IF(nid .EQ. 0) THEN
        dia_p = 0.1
        dx = 1.0
        dy = 1.0
        dz = 1.0
        x(1) = 0.0
        x(2) = 1.0
        y(1) = 0.0
        y(2) = 1.0
        z(1) = 0.0
        z(2) = 1.0 
        CALL test_particles_CreateParticles(dia_p)
        CALL test_grid_CreateGrid(x,y,z,dx,dy,dz)
      END IF

      IF(nid .EQ. 0) THEN
        PRINT*, ''
        PRINT*, 'ppiclF test run completed'
        PRINT*, '*****************************************************'
        PRINT*, ''
      END IF
      CALL MPI_FINALIZE(ierr)

      END PROGRAM
