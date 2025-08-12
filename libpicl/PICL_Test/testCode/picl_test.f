#include "../../ppiclF/source/PPICLF_USER.h"
#include "../../ppiclF/source/PPICLF_STD.h"

!----------------------------------------------------------------------

      PROGRAM main

      IMPLICIT NONE

      INCLUDE 'mpif.h'
      INCLUDE 'PPICLF'

      INTEGER*4 i, j, l, nproc, nid, icomm, ierr, test
      REAL*8    nndistTemp, nndist, PI, k

      ! Grid variables
      INTEGER*4 nCells(3), proc_ncells, ie
      REAL*8    grid(7,PPICLF_LEE), gridDomain(2,3), filter(3), 
     >          filterTemp(3), nFilterCells, tpF(PPICLF_LEE), num_bins,
     >          dx_min(3), feedback(PPICLF_LEE) 

      ! Particle variables
      REAL*8    part_y(PPICLF_LRS,PPICLF_LPART), pdia, C2Pratio, 
     >          part_r(PPICLF_LRP,PPICLF_LPART),T_truth(PPICLF_LPART),
     >          totErr, numErr, xp, yp, zp,r_npl,r_npt, i_err
      INTEGER*4 npart_local,totalParticles, ip

      ! Projection variables
      REAL*8    wsum, dSQl, dSQi, dist, CellVol, GaussianConst,
     >          w(PPICLF_LEE),part_feedbk(PPICLF_LPART),
     >          TrueFeedback(PPICLF_LEE), error

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
      nCells(1)       = 10 !Number of x cells in domain

      gridDomain(1,2) = 0.0D0 !y domain min
      gridDomain(2,2) = 1.0D0 !y domain max     
      nCells(2)       = 10 !Number of y cells in domain

      gridDomain(1,3) = 0.0D0 !z domain min
      gridDomain(2,3) = 1.0D0 !z domain max     
      nCells(3)       = 10 !Number of z cells in domain

      CALL test_CreateGrid(gridDomain,nCells,nid,nproc,grid,proc_ncells)
      CALL MPI_BARRIER(icomm, ierr)

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
     >     + 1.0D0*COS((k/(2*PI))*
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
      C2Pratio = 3.0 ! C2Pratio = cell dx to particle dx ratio
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
! ppiclF Inputs and test case setup
!**********************************************************************

      DO test = 1,9

        ! Periodicity flag Setup
        CALL test_setperiodic(x_per_flag,y_per_flag,z_per_flag,
     >                                            test,testcase)
        ! Will handle angular periodicity separately
        ang_per_flag   = 0
        ang_per_angle  = 0.0D0
        ang_per_xangle = 0.0D0
        ang_per_rin    = 0.0D0
        ang_per_rout   = 0.0D0  

        IF(test .EQ. 1) THEN
          CALL ppiclf_comm_InitMPI(icomm, nid, nproc)
          CALL test_PrintBanner(1,nid,nproc)
          CALL MPI_BARRIER(icomm,ierr)
        END IF
  
! Start ppiclF Calls
!**********************************************************************
        PPICLF_TEST    = .TRUE.
        PPICLF_PERTEST = .TRUE.
        PPICLF_OVERLAP = .FALSE.
        CALL ppiclf_solve_InitParticle(2,3,0,npart_local,
     >                                 part_y,part_r,filter,nndist)
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
        CALL ppiclf_comm_InitOverlapMesh(proc_ncells,grid)
        CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JT,tpF)
        CALL ppiclf_solve_InitSolve
        CALL MPI_BARRIER(icomm,ierr)
        DO ie = 1,proc_ncells
          CALL ppiclf_solve_GetProFld(ie,1,feedback(ie))
        END DO
        CALL MPI_BARRIER(icomm,ierr)

! Print Interpolation results 
!**********************************************************************
 
        IF(ppiclf_npart .GT. 0) THEN
          filename = ''
          filename = TRIM(testcase) // '_' // 'Interpolation_Proc_'
     >                  // TRIM(procString) // '.txt'
          OPEN(UNIT=300,FILE=TRIM(filename), STATUS='REPLACE',
     >                                      ACTION='WRITE')
          WRITE(300,*) 'Particle ID, x, y,',
     >                        'z, Interpolation Error (%).'
          totErr = 0.0D0
          numErr = 0.0D0
          DO i = 1,ppiclf_npart
            xp = ppiclf_y(PPICLF_JX,i)
            yp = ppiclf_y(PPICLF_JY,i)
            zp = ppiclf_y(PPICLF_JZ,i) 
            T_truth(i) = 1.0D0 
     >         + 1.0D0*SIN((k/(2*PI))*
     >         (xp/(gridDomain(2,1)-gridDomain(1,1))))
     >         + 1.0D0*COS((k/(2*PI))*
     >         (yp/(gridDomain(2,2)-gridDomain(1,2))))
     >         + 1.0D0*SIN((k/(2*PI))*
     >         (zp/(gridDomain(2,3)-gridDomain(1,3))))
            i_err = ABS(ppiclf_rprop(PPICLF_R_JT,i)-T_truth(i))
     >                  /T_truth(i)*100.0D0
            totErr = totErr + i_err 
            numErr = numErr + 1.0D0
            WRITE(300,*) i, xp, yp, zp, i_err
          END DO
          CLOSE(UNIT=300)
        ELSE
          totErr = 0.0D0
          numErr = 0.0D0
        END IF
        CALL MPI_BARRIER(icomm,ierr)
        CALL MPI_Allreduce(totErr,totErr,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_Allreduce(numErr,numErr,1,MPI_DOUBLE,
     >                                        MPI_SUM,iComm,ierr)
        CALL MPI_BARRIER(icomm,ierr)
        ! Print Bin Data
        IF(nid .EQ. 0) THEN
          PRINT*,'Average Interpolation Error for all particles:'
     >           ,totErr/numErr,'%'
        END IF
        CALL MPI_BARRIER(icomm,ierr)

! Print Projection Results
!**********************************************************************
        ! Find true projection result
        IF(nproc .EQ. 1) THEN
          TrueFeedback = 0.0D0
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
                  dSQl = MIN( (grid(l,ie) - part_y(l,ip))**2, 
     >             ( (gridDomain(2,l) - gridDomain(1,l)) -
     >                   ABS(grid(l,ie) - part_y(l,ip)) )**2 )
                ELSE
                  dSQl = (grid(l,ie) - part_y(l,ip))**2
                END IF
                dSQi = dSQi + dSQl
              END DO !l

              dist = SQRT(dSQi)
              CellVol = grid(7,ie)
              GaussianConst = 2.305D0
              w(ie) = ABS(CellVol*EXP(-GaussianConst*(dist**2)
     >                  / (CellVol**(2.0D0/3.0D0))))
              wsum = wsum + w(ie)
            END DO !ie

            ! Will range from 0 to 3
            part_feedbk(ip) = part_y(1,ip)
     >              /(gridDomain(2,1) - gridDomain(1,1)) +       
     >              part_y(2,ip)
     >              /(gridDomain(2,2) - gridDomain(1,2)) +       
     >              part_y(3,ip)
     >              /(gridDomain(2,3) - gridDomain(1,3))        
            DO ie = 1,proc_ncells
              TrueFeedback(ie) = TrueFeedback(ie) + 
     >                           w(ie)/wsum*part_feedbk(ip)
            END DO !ie
          END DO !ip

          filename = TRIM(testcase) // '_' //'FeedbackSolution.txt'
          OPEN(UNIT=499,FILE=filename, STATUS='REPLACE',ACTION='WRITE')
          WRITE(499,*) 'Cell ID, x_centroid, y_centroid,',
     >                        'z_centroid, True Feedback Solution'
          DO ie = 1,proc_ncells
            error = ABS(feedback(ie) - TrueFeedback(ie))
     >              / TrueFeedback(ie) * 100.0
            WRITE(499,*) ie, grid(1,ie), grid(2,ie),
     >                           grid(3,ie),TrueFeedback(ie)
          END DO
          CLOSE(UNIT=499)
        END IF


        filename = ''
        CALL MPI_BARRIER(icomm,ierr)
        IF(proc_ncells .GT. 0) THEN
          filename = TRIM(testcase) // '_' // 'Feedback_Proc_' //
     >                              TRIM(procString) // '.txt'
          OPEN(UNIT=400,FILE=filename, STATUS='REPLACE',ACTION='WRITE')
          WRITE(400,*) 'Cell ID, x_centroid, y_centroid,',
     >                        'z_centroid, feedback.'
          DO ie = 1,proc_ncells
            WRITE(400,*) ie, grid(1,ie), grid(2,ie),
     >                           grid(3,ie),feedback(ie)
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
            num_bins = ppiclf_n_bins(1)*
     >                 ppiclf_n_bins(2)*ppiclf_n_bins(3)
            CALL ppiclf_comm_CreateBin
            WRITE(100,*) ppiclf_np, 
     >         ppiclf_n_bins(1), ppiclf_n_bins(2), ppiclf_n_bins(3),
     >         num_bins , num_bins/ppiclf_np*100
          END DO
          CLOSE(UNIT=100)
!        END DO
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
      INTEGER*4    i, j, k, ii, nx_per_proc
      REAL*8       dx(3)
      CHARACTER*50 filename
      LOGICAL      FIRST, GOOD

      ! Creates rectangular grid
      DO i = 1,3
        dx(i) = (gIn(2,i) - gIn(1,i))/REAL(NCells(i))
      END DO

      ! Build full y & z domain on each processor
      ! Split x domain by number of processors
      nx_per_proc = CEILING(REAL(NCells(1))/REAL(NProc))
      IF(NProc .EQ. 1) nx_per_proc = NCells(1)
      procCells = 0
      FIRST = .TRUE.
      DO i = (nx_per_proc*ProcID)+1, nx_per_proc*(ProcID+1)
        DO j = 1,NCells(2)
          DO k = 1,NCells(3)
            GOOD = .TRUE.
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
                GOOD = .FALSE.
                EXIT
              END IF
            END DO
            IF(GOOD) THEN
              IF(FIRST) THEN
                WRITE(filename,'(A,I0,A)') 'Grid_Proc_0',
     >            ProcID, '.txt'
                OPEN(UNIT=200,FILE=filename, STATUS='REPLACE',
     >            ACTION='WRITE')
              END IF
              FIRST = .FALSE.
              WRITE(200,*) procCells, gOut(1,procCells),
     >                        gOut(2,procCells),gOut(3,procCells)
            END IF
          END DO !k
        END DO !j
      END DO !i
      CLOSE(UNIT=200)

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
     >          nx_perProc, istart, iend 

      
      DO i = 1,3
        part_dx(i) = dx(i)/dxr
        n(i) = INT((gDom(2,i)-gDom(1,i))/part_dx(i)) ! Num particles per dimension
        IF(n(i) .GT. INT(PPICLF_LPART**(1.0/3.0)))
     >    n(i) = INT(PPICLF_LPART**(1.0/3.0))
      END DO

      nx_perProc = INT(REAL(n(1))/REAL(np))+2
      istart = nx_perProc*pid + 1
      iend   = nx_perProc*(pid+1)
      IF(pid .EQ. np-1) iend = n(1) !ensures all particles are made
      Pcount = 0
      DO i = istart, iend 
        DO j = 1,n(2)
          DO k = 1,n(3)
            Pcount = Pcount + 1
            part_y(1,Pcount) = gDom(1,1) + pdia + part_dx(1)*(i-1)
            part_y(2,Pcount) = gDom(1,2) + pdia + part_dx(2)*(j-1)
            part_y(3,Pcount) = gDom(1,3) + pdia + part_dx(3)*(k-1)
            DO ii = 1,3
              ! deletes particles inadvertently created outside of domain
              IF(part_y(ii,Pcount) .GT. gDom(2,ii) - pdia) THEN
                Pcount = Pcount - 1
                EXIT
              END IF
            END DO
          END DO
        END DO
      END DO

      ! Add a particle in all 8 corners
      IF(pid .EQ. np-1) THEN
        Pcount = Pcount + 1
        part_y(1,Pcount) = gDom(1,1) + pdia/4
        part_y(2,Pcount) = gDom(1,2) + pdia/4
        part_y(3,Pcount) = gDom(1,3) + pdia/4
        DO i = 1,2
          DO j = 1,2
            DO k = 1,2
              IF(i .EQ. 1 .AND. j .EQ. 1 .AND. k .EQ. 1) CYCLE
              Pcount = Pcount + 1
              part_y(1,Pcount) = gDom(1,1) + gDom(2,1)*(i-1) - pdia/4
              part_y(2,Pcount) = gDom(1,2) + gDom(2,2)*(j-1) - pdia/4
              part_y(3,Pcount) = gDom(1,3) + gDom(2,3)*(k-1) - pdia/4
              DO ii = 1,3
                IF(part_y(ii,Pcount) .LT. 0.0) part_y(ii,Pcount) = 0.0D0
              END DO
            END DO
          END DO
        END DO
      END IF

      npar = Pcount

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
