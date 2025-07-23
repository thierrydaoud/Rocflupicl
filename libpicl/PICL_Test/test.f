#include "../../ppiclF/source/PPICLF_USER.h"
#include "../../ppiclF/source/PPICLF_STD.h"

!----------------------------------------------------------------------

      PROGRAM main

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF'

      INTEGER*4 i, j, nproc, nid, icomm, ierr
      REAL*8    nndistTemp, nndist
      ! Grid variables
      INTEGER*4 nCells
      REAL*8    grid(7,PPICLF_LEE), gridDomain(2,3), filter(3), 
     >          filterTemp(3), nFilterCells
      ! Particle variables
      REAL*8    part_y(PPICLF_LRS,PPICLF_LPART), pdia, dx_part(3), 
     >          dx_ratio, part_r(PPICLF_LRP,PPICLF_LPART)
      INTEGER*4 nPcells, npart_local

!     These are already included in PPICLF header
!      INTEGER*4 x_per_flag, y_per_flag, z_per_flag, ang_per_flag 
!      REAL*8    x_per_min, x_per_max, y_per_min, y_per_max,
!     >          z_per_min, z_per_max, ang_per_angle, ang_per_xangle,
!     >          ang_per_rin, ang_per_rout
! 
     
      ! MPI Setup
      CALL MPI_Init(ierr)
      icomm = MPI_COMM_WORLD
      CALL MPI_Comm_size(icomm,nproc,ierr)
      CALL MPI_Comm_rank(icomm,nid,ierr) 
      CALL ppiclf_comm_InitMPI(icomm, nid, nproc)
      CALL test_PrintBanner(1,nid,nproc)

      ! Grid Setup
      ! Create uniform rectangular grid
      DO i = 1,3 !3D uniform grid
        gridDomain(1,i) = 0.0D-2 !Min
        gridDomain(2,i) = 1.0D-2 !Max
      END DO
      ! Fluid Domain Min/Max
      x_per_min = gridDomain(1,1)
      x_per_max = gridDomain(2,1)
      y_per_min = gridDomain(1,2)
      y_per_max = gridDomain(2,2)
      z_per_min = gridDomain(1,3)
      z_per_max = gridDomain(2,3)

      ! Periodicity Setup
      x_per_flag     = 1
      y_per_flag     = 1
      z_per_flag     = 1
      ang_per_flag   = 0
      ang_per_angle  = 0.0D0
      ang_per_xangle = 0.0D0
      ang_per_rin    = 0.0D0
      ang_per_rout   = 0.0D0

      nCells = 8 !Number of cells per processor
      CALL test_CreateGrid(gridDomain,nCells,nid,nproc,grid)

      CALL MPI_BARRIER(icomm,ierr)

      IF(nid .EQ. 0) THEN
        ! Create particle positions
        dx_ratio = 5.0D0 !cell_dx/part_dx
        dx_part(1) = grid(4,1)/dx_ratio 
        dx_part(2) = grid(5,1)/dx_ratio 
        dx_part(3) = grid(6,1)/dx_ratio 
        pdia = MIN(dx_part(1),dx_part(2),dx_part(3)) ! so that particles are not overlapping 
        !Numbers of cells to fill with particles per dimension
        nPcells = 3
        CALL test_CreateParticles(gridDomain,dx_ratio,pdia,dx_part
     >                            ,nPcells,part_y,npart_local)
        part_r       = 0.0D0
        nndistTemp   = 4.0D0*pdia
        filterTemp   = 1.0D-10
        nFilterCells = 1.0
        ! Loop through all cells to find biggest cell dx
        DO j = 1,nCells
          DO i = 4,6
            IF(grid(i,j) > filterTemp(i)) filterTemp(i) = grid(i,j)
          END DO
        END DO
        DO i = 1,3
          filterTemp(i) = nFilterCells*filterTemp(i)
        END DO
      ELSE
        nndistTemp    = 0.0D0
        filterTemp(1) = 0.0D0
        filterTemp(2) = 0.0D0
        filterTemp(3) = 0.0D0
        npart_local   = 0
        part_y        = 0D0
        part_r        = 0.0D0
      END IF
      rhop = 7730.0D0 ! steel particles
      IF(npart_local .GT. 0) THEN
        DO i = 1,npart_local
           part_y(PPICLF_JVX,i) = 0.0D0
           part_y(PPICLF_JVY,i) = 0.0D0
           part_y(PPICLF_JVZ,i) = 0.0D0
           part_y(PPICLF_JT, i) = 300.0D0 ! particle temp
           part_y(PPICLF_JOX,i) = 0.0D0
           part_y(PPICLF_JOY,i) = 0.0D0
           part_y(PPICLF_JOZ,i) = 0.0D0
           part_r(PPICLF_R_JRHOP,i) = rhop ! particle density
           part_r(PPICLF_R_JDP,i)   = pdia ! particle diameter
           part_r(PPICLF_R_JVOLP,i) = (4.0D0/3.0D0)*3.14159265359
     >                                *(0.5D0*pdia)**3 ! particle volume
           ! Super Particle Loading (Real Number of particles = JSPL * number of compuational particles)
           part_r(PPICLF_R_JSPL,i) = 1.0D0
        END DO
      END IF

      CALL MPI_BARRIER(icomm,ierr)
      DO i = 1,3
        CALL MPI_Allreduce(filterTemp(i),filter(i),1,MPI_DOUBLE,
     >                     MPI_MAX,iComm,ierr)
      END DO
      CALL MPI_Allreduce(nndistTemp,nndist,1,MPI_DOUBLE,
     >                   MPI_MAX,iComm,ierr)


      CALL ppiclf_solve_InitParticle(2,3,0,npart_local,
     >                               part_y,part_r,filter,nndist)
 
      CALL ppiclf_solve_Initialize(x_per_flag, x_per_min, x_per_max,
     >                             y_per_flag, y_per_min, y_per_max, 
     >                             z_per_flag, z_per_min, z_per_max, 
     >                             ang_per_flag, ang_per_angle, 
     >                             ang_per_xangle, ang_per_rin,
     >                                                  ang_per_rout)

      CALL ppiclf_comm_InitOverlapMesh(nCells,grid)
      IF(nid .EQ. 0) THEN
        PRINT*,'Number of x bins:',ppiclf_n_bins(1)
        PRINT*,'Number of y bins:',ppiclf_n_bins(2)
        PRINT*,'Number of z bins:',ppiclf_n_bins(3)
        PRINT*,'Length of x bins:',ppiclf_bins_dx(1)
        PRINT*,'Length of y bins:',ppiclf_bins_dx(2)
        PRINT*,'Length of z bins:',ppiclf_bins_dx(3)
        PRINT*,'Number of total bins:',
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)
      END IF


      CALL test_PrintBanner(0,nid,nproc)
      CALL MPI_BARRIER(icomm,ierr)
      CALL MPI_FINALIZE(ierr)

      END PROGRAM

!----------------------------------------------------------------------

      SUBROUTINE test_PrintBanner(i,nid,nproc)

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

      SUBROUTINE test_CreateGrid(gIn,NCells,ProcID,NProc,gOut)
      
      ! Input/Output
      ! gIn(1:2,1:3) - fluid domain 1:min,2:max faces in x,y,&z
      ! NCells - per processor
      ! ProcID - Processor executing this subroutine
      ! NProc - number of processors
      ! gOut 1:3 - centroids, 4:6 - cell lengths, 7 - cell volume
      IMPLICIT NONE
      INTEGER*4 NCells, NProc, ProcID
      REAL*8    gIn(2,3),gOut(7,PPICLF_LEE)

      ! Local
      INTEGER*4 i, j, k, cellCount, NCellsPerDim
      REAL*8    delta(3), xMin(3), dx(3)

      NCellsPerDim = INT(NCells**(1.0/3.0)+0.5)
      ! Create uniform rectangular grid
      DO i = 1,3
        !processors fluid domain length
        delta(i) = (gIn(2,1) - gIn(1,1))/NProc
        xMin(i)  = delta(i)*ProcID ! This processors Min x,y,z
        dx(i)    = delta(i)/NCellsPerDim
      END DO
      cellCount = 0
      DO i = 1,NCellsPerDim 
        DO j = 1,NCellsPerDim
          DO k = 1,NCellsPerDim
            cellCount = cellCount + 1
            IF(cellCount .GT. NCells) EXIT
            gOut(1,cellCount) = xMin(1) + (i+0.5-1)*dx(1) !x centroid
            gOut(2,cellCount) = xMin(2) + (j+0.5-1)*dx(2) !y centroid
            gOut(3,cellCount) = xMin(3) + (k+0.5-1)*dx(3) !z centroid
            gOut(4,cellCount) = dx(1)
            gOut(5,cellCount) = dx(2)
            gOut(6,cellCount) = dx(3)
            gOut(7,cellCount) = dx(1)*dx(2)*dx(3)
          END DO !k
        END DO !j
      END DO !i

      IF(cellCount .NE. NCells) THEN
        PRINT*,'ERROR, wrong NCells generated on proc', NProc
      END IF

      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------

      SUBROUTINE test_CreateParticles(gDom,dxr,dp,dxp,ncell,part_y,np)
       
      ! Input/Output
      ! gDom(1:2,1:3) - fluid domain min/max in x,y,z
      ! dxr - particle size to grid dx size ratio
      ! dp  - particle diameter
      ! dxp - particle spacing
      ! ncell - number of fluid cells to fill with particles
      ! part_y(1:3,1:30,000) - particle centroid x,y,z coordinates
      ! np - number of particles created on this processor 
      IMPLICIT NONE
      
      REAL*8    dxr, dp, dxp(3), gDom(2,3), 
     >          part_y(PPICLF_LRS,PPICLF_LPART)
      INTEGER*4 ncell, Pcount, PartDomain, i, j, k, ii, np

      PartDomain = INT(dxr)*ncell
      Pcount = 0
      DO i = 1,PartDomain
        DO j = 1,PartDomain
          DO k = 1,PartDomain
            Pcount = Pcount + 1
            part_y(1,Pcount) = gDom(1,1) + dxp(1)*(i-1)
            part_y(2,Pcount) = gDom(1,2) + dxp(2)*(j-1)
            part_y(3,Pcount) = gDom(1,3) + dxp(3)*(k-1)
          END DO
        END DO
      END DO

      ! Add a single particle in all domain corners besides origin for
      ! periodicity testing
      DO i = 1,2
        DO j = 1,2
          DO k = 1,2
            IF(i .EQ. 1 .AND. j. EQ. 1 .AND. k .EQ. 1) CYCLE
            Pcount = Pcount + 1
            part_y(1,Pcount) = gDom(2,1)*(i-1) - dp/4
            part_y(2,Pcount) = gDom(2,2)*(j-1) - dp/4
            part_y(3,Pcount) = gDom(2,3)*(k-1) - dp/4
            DO ii = 1,3
              IF(part_y(ii,Pcount) .LT. 0.0D0) part_y(ii,Pcount) = 0.0D0
            END DO
          END DO
        END DO
      END DO

      np = Pcount

      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------

!!*******ADD BELOW TO MAIN TO TEST GRID OUTPUT*******
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


