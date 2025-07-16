      PROGRAM main

      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF'

      INTEGER*4 test, np, nid, icomm, ierr


      CALL MPI_Init(ierr)
      icomm = MPI_COMM_WORLD
      CALL MPI_Comm_size(icomm,np,ierr)
      CALL MPI_Comm_rank(icomm,nid,ierr)  

      PRINT*, 'Number of Processors:',np
      PRINT*, 'Processor ID:',nid

      CALL ppiclf_comm_InitMPI(icomm, nid, np)

      CALL MPI_FINALIZE(ierr)

      END PROGRAM
