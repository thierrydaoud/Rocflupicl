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
     >          filterTemp(3), nFilterCells, tpF(PPICLF_LEE), num_bins,
     >          dx_min(3) 

      ! Particle variables
      REAL*8    part_y(PPICLF_LRS,PPICLF_LPART), pdia, C2Pratio, 
     >          part_r(PPICLF_LRP,PPICLF_LPART),T_truth(PPICLF_LPART),
     >          totErr, numErr, xp, yp, zp,r_npl,r_npt
      INTEGER*4 npart_local,totalParticles

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
      gridDomain(1,1) = 0.0D0 !x domain min
      gridDomain(2,1) = 1.0D0 !x domain max
      nCells(1)       = 10 !Number of x cells in domain

      gridDomain(1,2) = 0.0D0 !y domain min
      gridDomain(2,2) = 1.0D0 !y domain max     
      nCells(2)       = 10 !Number of y cells in domain

      gridDomain(1,3) = 0.0D0 !z domain min
      gridDomain(2,3) = 1.0D0 !z domain max     
      nCells(3)       = 10 !Number of z cells in domain

      CALL test_CreateGrid(gridDomain,nCells,nid,nproc,grid,proc_ncells)
      CALL MPI_BARRIER(icomm, ierr)

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

      ! Find cell filter search distance
      filterTemp   = 1.0D-9 !dummy
      dx_min       = 1.0D9  !dummy
      nFilterCells = 2.0
      DO j = 1,proc_ncells
        DO i = 4,6
          ! Find largest & smallest grid dx, dy, dz
          IF(grid(i,j) > filterTemp(i-3)) filterTemp(i-3) = grid(i,j)
          IF(grid(i,j) < dx_min(i-3)) dx_min(i-3) = grid(i,j)
        END DO
      END DO
      DO i = 1,3
        filterTemp(i) = nFilterCells*filterTemp(i)
      END DO
      
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

      ! Setup filter(1:3) and smallest cell dx across processors 
      DO i = 1,3
        CALL MPI_Allreduce(filterTemp(i),filter(i),1,MPI_DOUBLE,
     >                                      MPI_MAX,iComm,ierr)
        CALL MPI_Allreduce(dx_min(i),dx_min(i),1,MPI_DOUBLE,
     >                                      MPI_MIN,iComm,ierr)
      END DO
      CALL MPI_BARRIER(icomm,ierr)

! Particle Setup   
!********************************************************************** 
      C2Pratio = 3.0 ! C2Pratio = cell dx / particle dx
      pdia     = MIN(dx_min(1)/C2Pratio,
     >               dx_min(2)/C2Pratio,
     >               dx_min(3)/C2Pratio)
      CALL test_CreateParticles(gridDomain,C2Pratio,pdia,dx_min,
     >                        npart_local,part_y,nid,nproc)
      CALL MPI_BARRIER(icomm,ierr)
      r_npl = REAL(npart_local)
      CALL MPI_Allreduce(r_npl,r_npt,1,MPI_DOUBLE,
     >                                      MPI_SUM,iComm,ierr)
      totalParticles = INT(r_npt)
      CALL MPI_BARRIER(icomm,ierr)
      IF(nid .EQ. 0) PRINT*,'Total Particles:',totalParticles 
      part_r = 0.0D0
      rhop   = 7730.0D0 ! steel particles
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
      nndistTemp  = 4.0D0*pdia
      CALL MPI_Allreduce(nndistTemp,nndist,1,MPI_DOUBLE,
     >                                      MPI_MAX,iComm,ierr)
      CALL MPI_BARRIER(icomm,ierr)
  
     
! Start ppiclF Calls
!**********************************************************************
      PPICLF_TEST = .TRUE.
      CALL ppiclf_solve_InitParticle(2,3,0,npart_local,
     >                               part_y,part_r,filter,nndist)
      PPICLF_TEST = .TRUE.
      CALL ppiclf_solve_Initialize(x_per_flag, x_per_min, x_per_max,
     >                             y_per_flag, y_per_min, y_per_max, 
     >                             z_per_flag, z_per_min, z_per_max, 
     >                             ang_per_flag, ang_per_angle, 
     >                             ang_per_xangle, ang_per_rin,
     >                                                  ang_per_rout)
      PPICLF_TEST = .TRUE.
      CALL ppiclf_comm_InitOverlapMesh(proc_ncells,grid)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JT,tpF)
      CALL ppiclf_solve_InitSolve
      CALL MPI_BARRIER(icomm,ierr)

! Print Interpolation results 
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
        PRINT*,'Average Interpolation Error for all particles:'
     >         ,totErr/numErr,'%'
      END IF
      CALL MPI_BARRIER(icomm,ierr)

! Test CreateBin variations
!********************************************************************** 
      IF(nproc .EQ. 1) THEN
        PRINT*,'******************************************************'
        PRINT*,'CreateBin Testing:'
        !DO j = 1,3
          
          DO i = 1,512
            ppiclf_np = i
            num_bins = ppiclf_n_bins(1)*
     >                 ppiclf_n_bins(2)*ppiclf_n_bins(3)
            CALL ppiclf_comm_CreateBin
!            PRINT*,'Num Proc:',i,
!     >           'Num Bins(x,y,z,tot):',ppiclf_n_bins(1),
!     >            ppiclf_n_bins(2),ppiclf_n_bins(3),
!     >            num_bins 
            PRINT*,'Num Bins/Num Proc:', num_bins/ppiclf_np
          END DO
        !END DO
      END IF


! Close out program
!********************************************************************** 
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
      SUBROUTINE test_CreateGrid(gIn,NCells,ProcID,NProc,gOut,procCells)
      
      ! Input/Output
      ! gIn(1:2,1:3) - fluid domain 1:min,2:max faces in x,y,&z
      ! NCells(3) - x,y,z cells in domain
      ! ProcID - Processor executing this subroutine
      ! NProc - number of processors
      ! gOut 1:3 - centroids, 4:6 - cell lengths, 7 - cell volume
  
      IMPLICIT NONE
      
      INTEGER*4 NCells(3), NProc, ProcID, procCells
      REAL*8    gIn(2,3), gOut(7,PPICLF_LEE)

      ! Local
      INTEGER*4 i, j, k, ii, nx_per_proc
      REAL*8    dx(3)

      ! Creates rectangular grid
      DO i = 1,3
        dx(i) = (gIn(2,i) - gIn(1,i))/REAL(NCells(i))
      END DO

      ! Build full y & z domain on each processor
      ! Split x domain by number of processors
      nx_per_proc = CEILING(REAL(NCells(1))/REAL(NProc))+1
      procCells = 0
      DO i = nx_per_proc*ProcID+1, nx_per_proc*(ProcID+1)
        DO j = 1,NCells(2)
          DO k = 1,NCells(3)
            procCells = procCells + 1
            gOut(1,procCells) = gIn(1,1) + (i-0.5)*dx(1) !x centroid
            gOut(2,procCells) = gIn(1,2) + (j-0.5)*dx(2) !y centroid
            gOut(3,procCells) = gIn(1,3) + (k-0.5)*dx(3) !z centroid
            gOut(4,procCells) = dx(1)
            gOut(5,procCells) = dx(2)
            gOut(6,procCells) = dx(3)
            gOut(7,procCells) = dx(1)*dx(2)*dx(3)
            DO ii = 1,3
              IF(gOut(ii,procCells) .GT. gIn(2,ii)) THEN
                procCells = procCells - 1
                EXIT
              END IF
            END DO
          END DO !k
        END DO !j
      END DO !i

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE test_CreateParticles(gDom,dxr,pdia,dx,npar,
     >                                part_y,pid,np)
       
      ! Input/Output
      ! gDom(1:2,1:3)        - fluid domain min/max in x,y,z
      ! dxr                  - ratio of min fluid cell dx and particle dx
      ! pdia                 - particle diameter
      ! dx(1:3)              - minimum fluid cell dx, dy, dz in domain
      ! npar                 - number of particles created on this processor 
      ! part_y(1:3,1:30,000) - particle centroid x,y,z coordinates
      ! pid                  - processor ID
      ! np                   - number of processors in operation
  
      IMPLICIT NONE
      
      REAL*8    dxr, dx(3), gDom(2,3), pdia,
     >          part_y(PPICLF_LRS,PPICLF_LPART), part_dx(3)
      INTEGER*4 Pcount, i, j, k, ii, npar, pid, np, n(3),
     >          nx_perProc 

      
      DO i = 1,3
        part_dx(i) = dx(i)/dxr
        n(i) = INT((gDom(2,i)-gDom(1,i))/part_dx(i)) ! Num particles per dimension
        IF(n(i) .GT. INT(PPICLF_LPART**(1.0/3.0)))
     >    n(i) = INT(PPICLF_LPART**(1.0/3.0))
      END DO

      nx_perProc = INT(REAL(n(1))/REAL(np))+2

      Pcount = 0
      DO i = nx_perProc*pid+1, nx_perProc*(pid+1)
        DO j = 1,n(2)
          DO k = 1,n(3)
            Pcount = Pcount + 1
            part_y(1,Pcount) = gDom(1,1) + pdia + part_dx(1)*(i-1)
            part_y(2,Pcount) = gDom(1,2) + pdia + part_dx(2)*(j-1)
            part_y(3,Pcount) = gDom(1,3) + pdia + part_dx(3)*(k-1)
            ! Ensure particles are within domain
            DO ii = 1,3
              IF(part_y(ii,Pcount) .GT. gDom(2,ii) - pdia) THEN
                Pcount = Pcount - 1
                EXIT
              END IF
            END DO
          END DO
        END DO
      END DO
      IF(pid .EQ. np-1 .AND. nx_perProc*(pid+1) .NE. n(1)) THEN
        DO i = nx_perProc*(pid+1),n(1)
          DO j = 1,n(2)
            DO k = 1,n(3)
              Pcount = Pcount + 1
              part_y(1,Pcount) = gDom(1,1) + pdia + part_dx(1)*(i-1)
              part_y(2,Pcount) = gDom(1,2) + pdia + part_dx(2)*(j-1)
              part_y(3,Pcount) = gDom(1,3) + pdia + part_dx(3)*(k-1)
              ! Ensure particles are within domain
              DO ii = 1,3
                IF(part_y(ii,Pcount) .GT. gDom(2,ii) - pdia) THEN
                  Pcount = Pcount - 1
                  EXIT
                END IF
              END DO
            END DO
          END DO
        END DO  
      END IF
      npar = Pcount

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
