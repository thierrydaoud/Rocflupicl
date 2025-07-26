#include "../../ppiclF/source/PPICLF_USER.h"
#include "../../ppiclF/source/PPICLF_STD.h"

!----------------------------------------------------------------------

      PROGRAM main

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF'

      INTEGER*4 i, j, nproc, nid, icomm, ierr
      REAL*8    nndistTemp, nndist, PI, k
      ! Grid variables
      INTEGER*4 nCells(3), proc_ncells
      REAL*8    grid(7,PPICLF_LEE), gridDomain(2,3), filter(3), 
     >          filterTemp(3), nFilterCells, tpF(PPICLF_LEE) 
      ! Particle variables
      REAL*8    part_y(PPICLF_LRS,PPICLF_LPART), pdia, dx_part(3), 
     >          dx_ratio, xp, yp, zp,
     >          part_r(PPICLF_LRP,PPICLF_LPART),T_truth(PPICLF_LPART),
     >          totErr, numErr
      INTEGER*4 nPcells(3), npart_local

      PI = 4.0D0*ATAN(1.0) ! pi
      k  = 1.0D0 ! wave number
! MPI Setup
!**********************************************************************
      CALL MPI_Init(ierr)
      icomm = MPI_COMM_WORLD
      CALL MPI_Comm_size(icomm,nproc,ierr)
      CALL MPI_Comm_rank(icomm,nid,ierr) 
      CALL ppiclf_comm_InitMPI(icomm, nid, nproc)
      CALL test_PrintBanner(1,nid,nproc)
      CALL MPI_BARRIER(icomm,ierr)

! Grid Setup
!**********************************************************************
      ! Create rectangular grid
      DO i = 1,3
        gridDomain(1,i) = 0.0D0  !x Min
        gridDomain(2,i) = 1.0D0 !x Max
      END DO
      nCells(1) = 3 !Number of x cells per processor
      nCells(2) = 3 !Number of y cells per processor
      nCells(3) = 3 !Number of z cells per processor
      CALL test_CreateGrid(gridDomain,nCells,nid,nproc,grid,proc_ncells)
      CALL MPI_BARRIER(icomm, ierr)

! Particle Setup   
!**********************************************************************
      dx_ratio = 5.0D0 !cell_dx/part_dx
      !Assumes dx, dy, dz are constant
      dx_part(1) = grid(4,1)/dx_ratio 
      dx_part(2) = grid(5,1)/dx_ratio 
      dx_part(3) = grid(6,1)/dx_ratio 
      pdia = MIN(dx_part(1),dx_part(2),dx_part(3))  
      !Numbers of cells to fill with particles per dimension per processor
      nPcells(1) = 3
      nPcells(2) = 3
      nPcells(3) = 3
      CALL test_CreateParticles(gridDomain,dx_ratio,pdia,dx_part,
     >                          nPcells,part_y,npart_local,nid,nproc)
      part_r       = 0.0D0
      nndistTemp   = 4.0D0*pdia
      filterTemp   = 1.0D-10
      nFilterCells = 2.0
      ! Loop through all cells to find biggest cell dx
      DO j = 1,proc_ncells
        DO i = 4,6
          IF(grid(i,j) > filterTemp(i)) filterTemp(i) = grid(i,j)
        END DO
      END DO
      DO i = 1,3
        filterTemp(i) = nFilterCells*filterTemp(i)
      END DO
      rhop = 7730.0D0 ! steel particles
      DO i = 1,npart_local
        part_y(PPICLF_JVX,i) = 0.0D0
        part_y(PPICLF_JVY,i) = 0.0D0
        part_y(PPICLF_JVZ,i) = 0.0D0
        part_y(PPICLF_JT, i) = 3000.0D0 ! particle temp
        part_y(PPICLF_JOX,i) = 0.0D0
        part_y(PPICLF_JOY,i) = 0.0D0
        part_y(PPICLF_JOZ,i) = 0.0D0
        part_r(PPICLF_R_JRHOP,i) = rhop ! particle density
        part_r(PPICLF_R_JDP,i)   = pdia ! particle diameter
        part_r(PPICLF_R_JVOLP,i) = (4.0D0/3.0D0)*PI
     >                              *(0.5D0*pdia)**3 ! particle volume
        part_r(PPICLF_R_JSPL,i) = 1.0D0 ! Super Particle Loading 
      END DO
      CALL MPI_BARRIER(icomm,ierr)
  
! ppiclF Inputs
!**********************************************************************
      ! Fluid Domain Min/Max
      x_per_min = gridDomain(1,1)
      x_per_max = gridDomain(2,1)
      y_per_min = gridDomain(1,2)
      y_per_max = gridDomain(2,2)
      z_per_min = gridDomain(1,3)
      z_per_max = gridDomain(2,3)
   
      ! Periodicity Setup
      x_per_flag     = 0
      y_per_flag     = 0
      z_per_flag     = 0
      ang_per_flag   = 0
      ang_per_angle  = 0.0D0
      ang_per_xangle = 0.0D0
      ang_per_rin    = 0.0D0
      ang_per_rout   = 0.0D0
      
      ! Setup fluid temperature field
      DO j = 1,proc_ncells
        tpF(j) = 1.0D0 
     >     + 1.0D0*SIN((k/(2*PI))*
     >     (grid(1,j)/(gridDomain(2,1)-gridDomain(1,1))))
     >     + 1.0D0*SIN((k/(2*PI))*
     >     (grid(2,j)/(gridDomain(2,2)-gridDomain(1,2))))
     >     + 1.0D0*SIN((k/(2*PI))*
     >     (grid(3,j)/(gridDomain(2,3)-gridDomain(1,3))))
        END DO
      ! Setup filter(1:3) and nearest neighbor search distance 
      DO i = 1,3
        CALL MPI_Allreduce(filterTemp(i),filter(i),1,MPI_DOUBLE,
     >                                      MPI_MAX,iComm,ierr)
      END DO
      CALL MPI_Allreduce(nndistTemp,nndist,1,MPI_DOUBLE,
     >                                      MPI_MAX,iComm,ierr)
      CALL MPI_BARRIER(icomm,ierr)
      
! Start ppiclF Calls
!**********************************************************************
      IF(npart_local .GT. 0) THEN
        CALL ppiclf_solve_InitParticle(2,3,0,npart_local,
     >                               part_y,part_r,filter,nndist)
      END IF
      CALL ppiclf_solve_Initialize(x_per_flag, x_per_min, x_per_max,
     >                             y_per_flag, y_per_min, y_per_max, 
     >                             z_per_flag, z_per_min, z_per_max, 
     >                             ang_per_flag, ang_per_angle, 
     >                             ang_per_xangle, ang_per_rin,
     >                                                  ang_per_rout)
      CALL ppiclf_comm_InitOverlapMesh(proc_ncells,grid)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JT,tpF)
      CALL ppiclf_solve_InitSolve
      CALL MPI_BARRIER(icomm,ierr)

! Print out results 
!**********************************************************************
      IF(nid .EQ. 0) THEN
      !  CALL test_BinTest()
      END IF
 
      IF(ppiclf_npart .GT. 0) THEN
        DO i = 1,ppiclf_npart
          xp = ppiclf_y(PPICLF_JX,i)
          yp = ppiclf_y(PPICLF_JY,i)
          zp = ppiclf_y(PPICLF_JZ,i) 
          T_truth(i) = 1.0D0 
     >       + 1.0D0*SIN((k/(2*PI))*
     >       (xp/(gridDomain(2,1)-gridDomain(1,1))))
     >       + 1.0D0*SIN((k/(2*PI))*
     >       (yp/(gridDomain(2,2)-gridDomain(1,2))))
     >       + 1.0D0*SIN((k/(2*PI))*
     >       (zp/(gridDomain(2,3)-gridDomain(1,3))))
        END DO
        totErr = 0.0D0
        numErr = 0.0D0
        DO i = 1,ppiclf_npart
          totErr = totErr + ABS(ppiclf_rprop(PPICLF_R_JT,i)-T_truth(i))
     >            /T_truth(i)*100.0D0
          numErr = numErr + 1.0D0
        END DO
      ELSE
        totErr = 0.0D0
        numErr = 0.0D0
      END IF
      CALL MPI_BARRIER(icomm,ierr)
      CALL MPI_Allreduce(totErr,totErr,1,MPI_DOUBLE,
     >                                      MPI_SUM,iComm,ierr)
      CALL MPI_Allreduce(numErr,numErr,1,MPI_DOUBLE,
     >                                      MPI_SUM,iComm,ierr)
      CALL MPI_BARRIER(icomm,ierr)
      ! Print Bin Data
      IF(nid .EQ. 0) THEN
        PRINT*,'Avgeraged Error:',totErr/numErr,'%'
        PRINT*,'Number of x bins:',ppiclf_n_bins(1)
        PRINT*,'Number of y bins:',ppiclf_n_bins(2)
        PRINT*,'Number of z bins:',ppiclf_n_bins(3)
        PRINT*,'Length of x bins:',ppiclf_bins_dx(1)
        PRINT*,'Length of y bins:',ppiclf_bins_dx(2)
        PRINT*,'Length of z bins:',ppiclf_bins_dx(3)
        PRINT*,'Number of total bins:',
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)
      END IF
       CALL MPI_BARRIER(icomm,ierr)
! Close out program
!********************************************************************** 
      CALL MPI_BARRIER(icomm,ierr)
      CALL MPI_FINALIZE(ierr)
      CALL test_PrintBanner(0,nid,nproc)
      END PROGRAM

!----------------------------------------------------------------------

      SUBROUTINE test_PrintBanner(i,nid,nproc)
      ! Input/Output
      ! i  - 1:Start of test, 2:End of test
      ! nid - Processor ID
      ! nproc - Number of Processors
  
      IMPLICIT NONE

      INTEGER*4 i, nid, nproc

      IF(i .EQ. 1) THEN
        IF(nid .EQ. 0) THEN
          PRINT*, ''
          PRINT*, '****************************************************'
          PRINT*, 'ppiclF test run starting'
          PRINT*, 'Number of Processors:',nproc
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
      SUBROUTINE test_CreateGrid(gIn,NCells,ProcID,NProc,gOut,totcells)
      
      ! Input/Output
      ! gIn(1:2,1:3) - fluid domain 1:min,2:max faces in x,y,&z
      ! NCells(3) - x,y,z cells per processor
      ! ProcID - Processor executing this subroutine
      ! NProc - number of processors
      ! gOut 1:3 - centroids, 4:6 - cell lengths, 7 - cell volume
  
      IMPLICIT NONE
      
      INTEGER*4 NCells(3), NProc, ProcID, totcells
      REAL*8    gIn(2,3),gOut(7,PPICLF_LEE),xMin

      ! Local
      INTEGER*4 i, j, k 
      REAL*8    dx(3), ProcDomain(3)

      ! Creates rectangular grid
      DO i = 1,3
        ProcDomain(i) = (gIn(2,i) - gIn(1,i))/REAL(NProc)
        dx(i)         = ProcDomain(i)/NCells(i)
      END DO
      ! Build full y & z domain on each processor
      ! Split x domain by number of processors
      xMin = ProcDomain(1)*ProcID ! This processors Min x
      totcells = 0
      DO i = 1,NCells(1) 
        DO j = 1,NCells(2)*NProc
          DO k = 1,NCells(3)*NProc
            totcells = totcells + 1
            gOut(1,totcells) = xMin     + (i+0.5-1)*dx(1) !x centroid
            gOut(2,totcells) = gIn(1,2) + (j+0.5-1)*dx(2) !y centroid
            gOut(3,totcells) = gIn(1,3) + (k+0.5-1)*dx(3) !z centroid
            gOut(4,totcells) = dx(1)
            gOut(5,totcells) = dx(2)
            gOut(6,totcells) = dx(3)
            gOut(7,totcells) = dx(1)*dx(2)*dx(3)
          END DO !k
        END DO !j
      END DO !i

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE test_CreateParticles(gDom,dxr,dp,dxp,nc,
     >                                part_y,npar,pid,np)
       
      ! Input/Output
      ! gDom(1:2,1:3) - fluid domain min/max in x,y,z
      ! dxr - particle size to grid dx size ratio
      ! dp  - particle diameter
      ! dxp - particle spacing
      ! nc - number of fluid cells to fill with particles
      ! part_y(1:3,1:30,000) - particle centroid x,y,z coordinates
      ! npar - number of particles created on this processor 
      ! pid - Processor ID
      ! np - Number of processors in use
  
      IMPLICIT NONE
      
      REAL*8    dxr, dp, dxp(3), gDom(2,3), 
     >          part_y(PPICLF_LRS,PPICLF_LPART),
     >          procDomain(3)
      INTEGER*4 nc(3), Pcount, PartDomain(3), i, j, k, ii, npar,pid,np

      DO i = 1,3
        ProcDomain(i) = (gDom(2,i)-gDom(1,i))/REAL(np)*REAL(pid)
        PartDomain(i) = INT(dxr)*nc(i)
      END DO

      Pcount = 0
      DO i = 1,PartDomain(1)
        DO j = 1,PartDomain(2)
          DO k = 1,PartDomain(3)
            Pcount = Pcount + 1
            part_y(1,Pcount) = ProcDomain(1) + dxp(1)*(i-0.5)
            part_y(2,Pcount) = ProcDomain(2) + dxp(2)*(j-0.5)
            part_y(3,Pcount) = ProcDomain(3) + dxp(3)*(k-0.5)
          END DO
        END DO
      END DO

      ! Add a single particle in all domain corners besides origin for
      ! periodicity testing
      IF(pid .EQ. 0) THEN
        DO i = 1,2
          DO j = 1,2
            DO k = 1,2
              IF(i .EQ. 1 .AND. j. EQ. 1 .AND. k .EQ. 1) CYCLE
              Pcount = Pcount + 1
              part_y(1,Pcount) = gDom(2,1)*(i-1) - dp/4
              part_y(2,Pcount) = gDom(2,2)*(j-1) - dp/4
              part_y(3,Pcount) = gDom(2,3)*(k-1) - dp/4
              DO ii = 1,3
                IF(part_y(ii,Pcount) .LT. 0.0D0) 
     >                         part_y(ii,Pcount) = dp/4
              END DO
            END DO
          END DO
        END DO
      END IF

      npar = Pcount

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
