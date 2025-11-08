!-----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_comm_InitMPI(comm,id,np)
     > bind(C, name="ppiclc_comm_InitMPI")
#else
      SUBROUTINE ppiclf_comm_InitMPI(comm,id,np)
#endif
!
!     This subroutine is called from rocflu/RFLU_InitFlowSolver.F90
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input: 
!
      INTEGER*4 comm
      INTEGER*4 id
      INTEGER*4 np
!
! Code:
!
      ! Ensures a later subroutine init wasn't called out of order
      IF (PPICLF_LINIT .OR. PPICLF_LFILT .OR. PPICLF_OVERLAP)
     >   CALL ppiclf_exittr('InitMPI must be called first$',0.0d0,0)

      ! set ppiclf_processor information
      ppiclf_comm = comm
      ppiclf_nid  = id
      ppiclf_np   = np

      ! GSlib call
      CALL ppiclf_prints('   *Begin InitCrystal$')
         CALL ppiclf_comm_InitCrystal
      CALL ppiclf_prints('    End InitCrystal$')

      ! check to make sure subroutine is called in correct order later
      ! on in the code sequence
      PPICLF_LCOMM = .true.

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_InitCrystal
!
!     This subroutine is called form ppiclf_comm_InitMPI
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input: 
!

!
! Code:
!
      ! GSlib call
      CALL pfgslib_crystal_setup(ppiclf_cr_hndl,ppiclf_comm,ppiclf_np)

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_CreateBin
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal:
!
      INTEGER*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >                            nfacegp, nedgegp, ncornergp
      INTEGER*4 ix, iy, iz, iperiodicx, iperiodicy, iperiodicz, 
     >          npt_total, j, i, idum, jdum, kdum, total_bin, 
     >          sum_value, count, targetTotBin, idealBin(3), iBin(3),
     >          iBinTot, temp,nBinMax,nBinMed,nBinMin, m, l, k,
     >          LBMax,LBMin
      REAL*8 xmin, ymin, zmin, xmax, ymax, zmax, rduml, rdumr, rthresh,
     >       rmiddle, rdiff, binb_length(3),temp1,temp2
      INTEGER*4 ppiclf_iglsum
      EXTERNAL ppiclf_iglsum
      REAL*8 ppiclf_glmin,ppiclf_glmax,ppiclf_glsum
      EXTERNAL ppiclf_glmin,ppiclf_glmax,ppiclf_glsum
!

! face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/
     >                 -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                 -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      IF(ppiclf_ndim .GT. 2) THEN
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      END IF

      ix = 1
      iy = 2
      iz = 1
      IF(ppiclf_ndim .EQ. 3)
     >iz = 3

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)
        
      ! iglsum is integer addition across MPI ranks.
      npt_total = ppiclf_iglsum(ppiclf_npart,1)
      ! compute bin boundaries
      xmin = 1E10
      ymin = 1E10
      zmin = 1E10
      xmax = -1E10
      ymax = -1E10
      zmax = -1E10
      ! Looping through particles on this processor
      DO i=1,ppiclf_npart
         ! Finding min/max particle extremes.
         ! Need to consider filter/neighborwidths
         ! to ensure ppiclf_bins_dx > ppiclf_d2chk(1)
         rduml = ppiclf_y(ix,i) - ppiclf_d2chk(1)
         rdumr = ppiclf_y(ix,i) + ppiclf_d2chk(1)
         IF(rduml .LT. xmin) xmin = rduml
         IF(rdumr .GT. xmax) xmax = rdumr

         rduml = ppiclf_y(iy,i) - ppiclf_d2chk(1)
         rdumr = ppiclf_y(iy,i) + ppiclf_d2chk(1)
         IF(rduml .LT. ymin) ymin = rduml
         IF(rdumr .GT. ymax) ymax = rdumr

         IF(ppiclf_ndim .EQ. 3) THEN
            rduml = ppiclf_y(iz,i) - ppiclf_d2chk(1)
            rdumr = ppiclf_y(iz,i) + ppiclf_d2chk(1)
            IF(rduml .LT. zmin) zmin = rduml
            IF(rdumr .GT. zmax) zmax = rdumr
         END IF
      END DO
      ! Finds global max/mins across MPI ranks
      ppiclf_binb(1) = ppiclf_glmin(xmin,1)
      ppiclf_binb(2) = ppiclf_glmax(xmax,1)
      ppiclf_binb(3) = ppiclf_glmin(ymin,1)
      ppiclf_binb(4) = ppiclf_glmax(ymax,1)
      ppiclf_binb(5) = 0.0d0
      ppiclf_binb(6) = 0.0d0
      IF(ppiclf_ndim .GT. 2) ppiclf_binb(5) = ppiclf_glmin(zmin,1)
      IF(ppiclf_ndim .GT. 2) ppiclf_binb(6) = ppiclf_glmax(zmax,1)

!      if (npt_total .gt. 0) then
!      do i=1,ppiclf_ndim
!         if (ppiclf_bins_balance(i) .eq. 1) then
!            rmiddle = 0.0
!            do j=1,ppiclf_npart
!               rmiddle = rmiddle + ppiclf_y(i,j)
!            enddo
!            rmiddle = ppiclf_glsum(rmiddle,1)
!            rmiddle = rmiddle/npt_total
!
!            rdiff =  max(abs(rmiddle-ppiclf_binb(2*(i-1)+1)),
!     >                   abs(ppiclf_binb(2*(i-1)+2)-rmiddle))
!            ppiclf_binb(2*(i-1)+1) = rmiddle - rdiff
!            ppiclf_binb(2*(i-1)+2) = rmiddle + rdiff
!         endif
!      enddo
!      endif

       if (ppiclf_xdrange(2,1) .lt. ppiclf_binb(2) .or.
     >    ppiclf_xdrange(1,1) .gt. ppiclf_binb(1) .or. 
     >    iperiodicx .eq. 0) then
         ppiclf_binb(1) = ppiclf_xdrange(1,1)
         ppiclf_binb(2) = ppiclf_xdrange(2,1)
       endif

       if (ppiclf_xdrange(2,2) .lt. ppiclf_binb(4) .or.
     >    ppiclf_xdrange(1,2) .gt. ppiclf_binb(3) .or.
     >    iperiodicy .eq. 0) then
         ppiclf_binb(3) = ppiclf_xdrange(1,2)
         ppiclf_binb(4) = ppiclf_xdrange(2,2)
       endif
      
       if (ppiclf_ndim .gt. 2) then
       if (ppiclf_xdrange(2,3) .lt. ppiclf_binb(6) .or.
     >    ppiclf_xdrange(1,3) .gt. ppiclf_binb(5) .or. 
     >    iperiodicz .eq. 0) then
         ppiclf_binb(5) = ppiclf_xdrange(1,3)
         ppiclf_binb(6) = ppiclf_xdrange(2,3)
      endif ! ndim
      endif ! xdrange

      ! End subroutine if no particles present      
      IF(npt_total .LT. 1) RETURN
      LBMax = 0
      LBMin = 0
      temp1 = 1.0D-10
      temp2 = 1.0D10
      ! Find ppiclf bin domain lengths
      ! and Max,Med,Min dimensions
       DO l = 1,3
        binb_length(l) = ppiclf_binb(2*l) -
     >                         ppiclf_binb(2*l-1)
        IF(binb_length(l).GT.temp1) THEN
          temp1 = binb_length(l)
          LBMax = l
        END IF
        IF(binb_length(l).LT.temp2) THEN
          temp2 = binb_length(l)
          LBMin = l
        END IF
       END DO

       IF(ppiclf_ndim .LT. 3)
     >   CALL ppiclf_exittr('CreateBins only supports 3D Grids',0.0D0,0)
      
      ! Update with targetTotBin = ActiveBinNum
      targetTotBin = ppiclf_np

      ! Number of bins calculated based on bin surface
      ! area minimization and bin aspect ratio close to 1
       ppiclf_n_bins(1) = INT((targetTotBin**(1.0D0/3.0D0))*
     >                   (binb_length(1)**(2.0D0/3.0D0))/ 
     >                   ((binb_length(2)**(1.0D0/3.0D0))*
     >                   (binb_length(3))**(1.0D0/3.0D0)))
     
      ppiclf_n_bins(2) = INT((targetTotBin**(1.0D0/3.0D0))*
     >                   (binb_length(2)**(2.0D0/3.0D0))/ 
     >                   ((binb_length(1)**(1.0D0/3.0D0))*
     >                   (binb_length(3))**(1.0D0/3.0D0)))
     
      ppiclf_n_bins(3) = INT((targetTotBin**(1.0D0/3.0D0))*
     >                   (binb_length(3)**(2.0D0/3.0D0))/ 
     >                   ((binb_length(2)**(1.0D0/3.0D0))*
     >                   (binb_length(1))**(1.0D0/3.0D0)))
      ! Since INT trucates, make sure n_bins at least 1 
      DO l = 1,3
        IF(ppiclf_n_bins(l) .LT. 1) ppiclf_n_bins(l) = 1
      END DO

      iBinTot = 0

      ! Filterwidth criteria check.  ppiclf_d2chk(2) automatically
      ! set to be at least 2 fluid cell widths in
      ! PICL_TEMP_InitSolver.F90

       DO l = 1,3
          ! Ensure ppiclf_bin_dx(l) > ppiclf_d2chk(1) 
          if((binb_length(l)/ppiclf_n_bins(l)) .LT. ppiclf_d2chk(1)) 
     >      ppiclf_n_bins(l) = INT(binb_length(l)/ppiclf_d2chk(1))
          if(ppiclf_n_bins(l) .LT. 1) then
            print*, "*** ERROR in ppiclf_CreateBin"
            print*, "ppiclf_n_bins(l) =", ppiclf_n_bins(l), "l =", l
            print*, "ppiclf_d2chk(1) =", ppiclf_d2chk(1)
            CALL ppiclf_exittr('ppiclf_d2chk(1) criteria violated.',
     >      0.0d0,0)
          endif
        idealBin(l) = ppiclf_n_bins(l)
       END DO

      ! Since bin must be an integer, check -1, +0, +1 number of bins for each bin dimension
      ! ideal number of bins will be max value while less than number of total target of bins.
      ! Will not check total bin value (cycle do loop) if
      ! ppiclf_d2chk(1) criteria is violated or ppiclf_n_bins < 1

      total_bin = 0 
      DO ix = 1,3
        iBin(1) = ppiclf_n_bins(1) + (ix-2)
        ppiclf_bins_dx(1) = binb_length(1)/iBin(1)
        IF(ppiclf_bins_dx(1) .LT. ppiclf_d2chk(1) .OR.
     >                           iBin(1) .LT. 1) CYCLE
        DO iy = 1,3
          iBin(2) = ppiclf_n_bins(2) + (iy-2)
          ppiclf_bins_dx(2) = binb_length(2)/iBin(2)
          IF(ppiclf_bins_dx(2) .LT. ppiclf_d2chk(1) .OR.
     >                             iBin(2) .LT. 1) CYCLE
          DO iz = 1,3
            iBin(3) = ppiclf_n_bins(3) + (iz-2)
            ppiclf_bins_dx(3) = binb_length(3)/iBin(3)
            IF(ppiclf_bins_dx(3) .LT. ppiclf_d2chk(1) .OR.
     >                               iBin(3) .LT. 1) CYCLE
            iBinTot = iBin(1)*iBin(2)*iBin(3)
            IF(iBinTot .GT. total_bin .AND.
     >                     iBinTot .LE. targetTotBin) THEN
              total_bin = 1
              DO l = 1,3
                idealBin(l) = iBin(l)
                total_bin = total_bin*idealBin(l)
              END DO
              ! These loops are to make sure the dimension with the longest
              ! ppiclf_binb length gets more bins in the case where two or
              ! more dimensions are within 1 bin division of each other.
              temp = 0
              nBinMax = MAX(idealBin(1),idealBin(2),idealBin(3))
              nBinMin = MIN(idealBin(1),idealBin(2),idealBin(3))
              nBinMed = -99
              DO l = 1,3
                IF(idealBin(l).LT.nBinMax .AND. idealBin(l).GT.nBinMin)
     >             nBinMed = idealBin(l)
              END DO
              IF(nBinMed.EQ. -99) THEN !two number of bins are equal
                DO l = 1,3
                  IF(idealBin(l).EQ.nBinMax) temp = temp + 1
                  IF(idealBin(l).EQ.nBinMin) temp = temp + 10
                END DO
                IF(temp .EQ. 2) THEN
                  nBinMed = nBinMax
                ELSE ! Either two nBinMin or all 3 equal
                  nBinMed = nBinMin
                END IF
              END IF
              DO l = 1,3
                IF(l.EQ.LBMax) THEN
                  idealBin(l)=nBinMax
                ELSE IF(l.EQ.LBMin) THEN
                  idealBin(l)=nBinMin
                ELSE
                  idealBin(l)=nBinMed 
                END IF
              END DO 
            END IF
          END DO !iz
        END DO !iy
      END DO !ix

      ! Set common ppiclf arrays based on above calculation
      DO l = 1,3
        ppiclf_n_bins(l) = idealBin(l)
        ppiclf_bins_dx(l) = binb_length(l)/ppiclf_n_bins(l)
      END DO


      ! Loop to see if we can add one to dimension with largest number of bins
      ! Choose this dimension because it is smallest incremental increase to total bins 
      DO
        IF((total_bin/ppiclf_n_bins(LBMax))*
     >      (ppiclf_n_bins(LBMax)+1) .LT. targetTotBin) THEN
          ! Add a bin and set new bin dx length
          ppiclf_n_bins(LBMax) = ppiclf_n_bins(LBMax)+1
          ppiclf_bins_dx(LBMax) = binb_length(LBMax)/
     >                              ppiclf_n_bins(LBMax)
          IF(ppiclf_bins_dx(LBMax) .LT. ppiclf_d2chk(1)) THEN
            ! If ppiclf_d2chk criteria violated, return to previous bin configuration
            ppiclf_n_bins(LBMax) = ppiclf_n_bins(LBMax)-1
            ppiclf_bins_dx(LBMax) = binb_length(LBMax)/
     >                                ppiclf_n_bins(LBMax)
            EXIT
          END IF
          total_bin = 1
          DO l = 1,3
            total_bin = total_bin*ppiclf_n_bins(l)
          END DO
        ELSE
          EXIT
        END IF
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! David's old binning method left below for now
!      finished(1) = 0
!      finished(2) = 0
!      finished(3) = 0
!      total_bin = 1 
!
!      do i=1,ppiclf_ndim
!         finished(i) = 0
!         exit_1_array(i) = ppiclf_bins_set(i)
!         exit_2_array(i) = 0
!         if (ppiclf_bins_set(i) .ne. 1) ppiclf_n_bins(i) = 1
!         ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                        ppiclf_binb(2*(i-1)+1)  ) / 
!     >                       ppiclf_n_bins(i)
!         ! Make sure exit_2 is not violated by user input
!         if (ppiclf_bins_dx(i) .lt. ppiclf_d2chk(1)) then
!            do while (ppiclf_bins_dx(i) .lt. ppiclf_d2chk(1))
!               ppiclf_n_bins(i) = max(1, ppiclf_n_bins(i)-1)
!               ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                              ppiclf_binb(2*(i-1)+1)  ) / 
!     >                             ppiclf_n_bins(i)
!         WRITE(*,*) "Inf. loop in CreateBin", i, 
!     >              ppiclf_bins_dx(i), ppiclf_d2chk(1)
!         call ppiclf_exittr('Inf. loop in CreateBin$',0.0,0)
!            enddo
!         endif
!         total_bin = total_bin*ppiclf_n_bins(i)
!      enddo
!
!      ! Make sure exit_1 is not violated by user input
!      count = 0
!      do while (total_bin > ppiclf_np)
!          count = count + 1;
!          i = modulo((ppiclf_ndim-1)+count,ppiclf_ndim)+1
!          ppiclf_n_bins(i) = max(ppiclf_n_bins(i)-1,1)
!          ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                         ppiclf_binb(2*(i-1)+1)  ) / 
!     >                        ppiclf_n_bins(i)
!          total_bin = 1
!          do j=1,ppiclf_ndim
!             total_bin = total_bin*ppiclf_n_bins(j)
!          enddo
!          if (total_bin .le. ppiclf_np) exit
!       enddo
!
!       exit_1 = .false.
!       exit_2 = .false.
!
!       do while (.not. exit_1 .and. .not. exit_2)
!          do i=1,ppiclf_ndim
!             if (exit_1_array(i) .eq. 0) then
!                ppiclf_n_bins(i) = ppiclf_n_bins(i) + 1
!                ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                               ppiclf_binb(2*(i-1)+1)  ) / 
!     >                              ppiclf_n_bins(i)
!
!                ! Check conditions
!                ! exit_1
!                total_bin = 1
!                do j=1,ppiclf_ndim
!                   total_bin = total_bin*ppiclf_n_bins(j)
!                enddo
!                if (total_bin .gt. ppiclf_np) then
!                   ! two exit arrays aren't necessary for now, but
!                   ! to make sure exit_2 doesn't slip through, we
!                   ! set both for now
!                   exit_1_array(i) = 1
!                   exit_2_array(i) = 1
!                   ppiclf_n_bins(i) = ppiclf_n_bins(i) - 1
!                   ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                                  ppiclf_binb(2*(i-1)+1)  ) / 
!     >                                  ppiclf_n_bins(i)
!                   exit
!                endif
!                
!                ! exit_2
!                if (ppiclf_bins_dx(i) .lt. ppiclf_d2chk(1)) then
!                   ! two exit arrays aren't necessary for now, but
!                   ! to make sure exit_2 doesn't slip through, we
!                   ! set both for now
!                   exit_1_array(i) = 1
!                   exit_2_array(i) = 1
!                   ppiclf_n_bins(i) = ppiclf_n_bins(i) - 1
!                   ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                                  ppiclf_binb(2*(i-1)+1)  ) / 
!     >                                  ppiclf_n_bins(i)
!                   exit
!                endif
!             endif
!          enddo
!
!          ! full exit_1
!          sum_value = 0
!          do i=1,ppiclf_ndim
!             sum_value = sum_value + exit_1_array(i)
!          enddo
!          if (sum_value .eq. ppiclf_ndim) then
!             exit_1 = .true.
!          endif
!
!          ! full exit_2
!          sum_value = 0
!          do i=1,ppiclf_ndim
!             sum_value = sum_value + exit_2_array(i)
!          enddo
!          if (sum_value .eq. ppiclf_ndim) then
!             exit_2 = .true.
!          endif
!       enddo
!      ! Check for too small bins 
!      rthresh = 1E-12
!      total_bin = 1
!      do i=1,ppiclf_ndim
!         total_bin = total_bin*ppiclf_n_bins(i)
!         if (ppiclf_bins_dx(i) .lt. rthresh) ppiclf_bins_dx(i) = 1.0
!      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! -------------------------------------------------------
! SETUP 3D BACKGROUND GRID PARAMETERS FOR GHOST PARTICLES
! -------------------------------------------------------

!     current box coordinates
      IF(ppiclf_nid .LE. total_bin-1) THEN
         idum = modulo(ppiclf_nid,ppiclf_n_bins(1))
         jdum = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
         kdum = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))
         IF(ppiclf_ndim .LT. 3) kdum = 0
         ppiclf_binx(1,1) = ppiclf_binb(1) + idum    *ppiclf_bins_dx(1)
         ppiclf_binx(2,1) = ppiclf_binb(1) + (idum+1)*ppiclf_bins_dx(1)
         ppiclf_biny(1,1) = ppiclf_binb(3) + jdum    *ppiclf_bins_dx(2)
         ppiclf_biny(2,1) = ppiclf_binb(3) + (jdum+1)*ppiclf_bins_dx(2)
         ppiclf_binz(1,1) = 0.0d0
         ppiclf_binz(2,1) = 0.0d0
         IF(ppiclf_ndim .GT. 2) THEN
            ppiclf_binz(1,1) = ppiclf_binb(5)+kdum    *ppiclf_bins_dx(3)
            ppiclf_binz(2,1) = ppiclf_binb(5)+(kdum+1)*ppiclf_bins_dx(3)
         END IF
      END IF

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateSubBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 nbin, idum, jdum, kdum, ndumx, ndumy, itmp, jtmp, ktmp,
     >          i, j, k
!

      nbin = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)

c     current box coordinates
      if (ppiclf_nid .le. nbin-1) then
         idum = modulo(ppiclf_nid,ppiclf_n_bins(1))
         jdum = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
         kdum = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))
         if (ppiclf_ndim .lt. 3) kdum = 0
         ! interior grid of each bin
         ! +1 for making mesh smaller and +1 since these are vertice counts
         ppiclf_bx = floor(ppiclf_bins_dx(1)/ppiclf_filter) + 1 + 1
         ppiclf_by = floor(ppiclf_bins_dx(2)/ppiclf_filter) + 1 + 1
         ppiclf_bz = 1
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_bz = floor(ppiclf_bins_dx(3)/ppiclf_filter) + 1 + 1

         ppiclf_bx = ppiclf_bx*ppiclf_ngrids
         ppiclf_by = ppiclf_by*ppiclf_ngrids
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_bz = ppiclf_bz*ppiclf_ngrids

         if (ppiclf_bx .gt. PPICLF_BX1)
     >      call ppiclf_exittr('Increase PPICLF_BX1$',0.,ppiclf_bx)
         if (ppiclf_by .gt. PPICLF_BY1)
     >      call ppiclf_exittr('Increase PPICLF_BY1$',0.,ppiclf_by)
         if (ppiclf_bz .gt. PPICLF_BZ1)
     >      call ppiclf_exittr('Increase PPICLF_BZ1$',0.,ppiclf_bz)

         ppiclf_rdx = ppiclf_bins_dx(1)/(ppiclf_bx-1)
         ppiclf_rdy = ppiclf_bins_dx(2)/(ppiclf_by-1)
         ppiclf_rdz = 0
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_rdz = ppiclf_bins_dx(3)/(ppiclf_bz-1)

         ndumx = ppiclf_n_bins(1)*(ppiclf_bx-1) + 1
         ndumy = ppiclf_n_bins(2)*(ppiclf_by-1) + 1
    
         do k=1,ppiclf_bz
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            ppiclf_grid_x(i,j,k) = sngl(ppiclf_binx(1,1) +
     >                                  (i-1)*ppiclf_rdx)
            ppiclf_grid_y(i,j,k) = sngl(ppiclf_biny(1,1) +
     >                                  (j-1)*ppiclf_rdy)
            ppiclf_grid_z(i,j,k) = sngl(ppiclf_binz(1,1) +
     >                                  (k-1)*ppiclf_rdz)

            itmp = idum*(ppiclf_bx-1) + (i-1)
            jtmp = jdum*(ppiclf_by-1) + (j-1)
            ktmp = kdum*(ppiclf_bz-1) + (k-1)
    
            ppiclf_grid_i(i,j,k)  = itmp + ndumx*jtmp + ndumx*ndumy*ktmp

         enddo
         enddo
         enddo

      endif

      return
      end
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_MapOverlapMesh
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE 'mpif.h'
!
! Internal:
!
      INTEGER*4 icalld
      SAVE      icalld
      DATA      icalld /0/
      INTEGER*4 nkey(2), i, j, k,l, ie, iee, ii, jj, kk, ndum, nrank,
     >          nl, nii, njj, nrr, ilow, jlow, klow, nxyz, il,
     >          ihigh, jhigh, khigh, ierr
      INTEGER*4 ix, iy, iz, ixLow, ixHigh, iyLow,
     >          iyHigh, izLow, izHigh 
      REAL*8    rxval, ryval, rzval, EleSizei(3), MaxPoint(3),
     >          MinPoint(3), ppiclf_vlmin, ppiclf_vlmax,
     >          centeri(3)
      LOGICAL   partl, ErrorFound
      EXTERNAL  ppiclf_vlmin, ppiclf_vlmax
!
! Code Start:
!

      nxyz = PPICLF_LEZ*PPICLF_LEY*PPICLF_LEX !Num of vertices per cell
      ppiclf_neltb = 0 !counts number of Rocflu elements on this processor
                       !that are within the ppiclf bounds domain
      DO ie=1,ppiclf_nee ! Number of fluid cells in Rocflu Grid domain
        ! Find fluid cell max x,y,z lengths and centroid
        DO l=1,3
          centeri(l)  = 0.0D0
          MaxPoint(l) = -1000000.0D0
          MinPoint(l) =  1000000.0D0 
          EleSizei(l) =  0.0D0
          DO k=1,PPICLF_LEZ
            DO j=1,PPICLF_LEY
              DO i=1,PPICLF_LEX
                centeri(l) = centeri(l) + ppiclf_xm1bs(i,j,k,l,ie)
                IF (ppiclf_xm1bs(i,j,k,l,ie) .GT. MaxPoint(l)) 
     >              MaxPoint(l) = ppiclf_xm1bs(i,j,k,l,ie)
                IF (ppiclf_xm1bs(i,j,k,l,ie) .LT. MinPoint(l)) 
     >              MinPoint(l) = ppiclf_xm1bs(i,j,k,l,ie)
              END DO !i
            END DO !j
          END DO !k
          centeri(l) = centeri(l)/nxyz

          ! 1.2 times the cell length to ensure that one layers of cells
          ! outside of the ppiclf bin are mapped for interpolation.
          ! This could be changed based on the frequency of ppiclf bin
          ! creation and mapping.
          EleSizei(l) = 1.2*(MaxPoint(l) - MinPoint(l))
          IF(EleSizei(l) .GT. ppiclf_bins_dx(l)) THEN
            PRINT*,'EleSizei > ppiclf_bins_dx in MapOverlapMesh'
            PRINT*, 'Dimension:',l
            PRINT*, "Element size(l) =", EleSizei(l), 
     >              "ppiclf_bins_dx(l) =", ppiclf_bins_dx(l)
            CALL ppiclf_exittr('Error: EleSizei(OverlapMesh)',0.0D0,l)
          END IF
        END DO !l 

        ! Fluid Cell vertex position without additional length
        rxval = centeri(1)
        ryval = centeri(2)
        rzval = 0.0D0
        IF(ppiclf_ndim .GT. 2) rzval = centeri(3)
      
        ! Exits if fluid cell vertex is outside of all bin 
        ! boundaries + Exchange Ghost Fluid Cell Buffer (EleSizei)
        IF (rxval .GT. (ppiclf_binb(2))) CYCLE
        IF (rxval .LT. (ppiclf_binb(1))) CYCLE
        IF (ryval .GT. (ppiclf_binb(4))) CYCLE
        IF (ryval .LT. (ppiclf_binb(3))) CYCLE
        IF (ppiclf_ndim .GT. 2 .AND. rzval .GT. 
     >      (ppiclf_binb(6))) CYCLE
        IF (ppiclf_ndim.GT.2 .AND. rzval .LT.
     >      (ppiclf_binb(5))) CYCLE
 
        ! Determines what bin the fluid cell is nominally mapped to
        ii    = FLOOR((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
        jj    = FLOOR((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
        kk    = FLOOR((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3))

        ! Default is Do loop with ix=iy=iz=2 for fluid cells not near
        ! bin boundary

        ixLow =2
        ixHigh=2
        iyLow =2
        iyHigh=2
        izLow =2
        izHigh=2

        ! These series of if statements check if bin mapping changes
        ! when adding/subtracting multiple of fluid cell length defined
        ! by EleSizei(l). 
        ! This is used to map fluid cells slightly outside of the ppiclf
        ! bin boundary.  If any .NE. 2, then fluid cell is mapped to
        ! multiple ppiclf bins. 
        
        IF (FLOOR((rxval + EleSizei(1) - ppiclf_binb(1))
     >       /ppiclf_bins_dx(1)) .NE. ii)  ixHigh = 3

        IF (FLOOR((rxval - EleSizei(1) - ppiclf_binb(1))
     >       /ppiclf_bins_dx(1)) .NE. ii)  ixLow = 1
        IF (FLOOR((ryval + EleSizei(2) - ppiclf_binb(3))
     >       /ppiclf_bins_dx(2)) .NE. jj)  iyHigh = 3

        IF (FLOOR((ryval - EleSizei(2) - ppiclf_binb(3))
     >       /ppiclf_bins_dx(2)) .NE. jj)  iyLow = 1

        IF (ppiclf_ndim .GT. 2 .AND. FLOOR((rzval + EleSizei(3)
     >    - ppiclf_binb(5))/ppiclf_bins_dx(3)) .NE. kk)  izHigh = 3

        IF (ppiclf_ndim .GT. 2 .AND. FLOOR((rzval - EleSizei(3)
     >    - ppiclf_binb(5))/ppiclf_bins_dx(3)) .NE. kk)  izLow = 1

        DO ix=ixLow,ixHigh
          DO iy=iyLow,iyHigh
            DO iz=izLow,izHigh
              ! Change cell position by EleSizei if ix,iy,or iz NE 2
              rxval = centeri(1) + (ix-2)*EleSizei(1)
              ryval = centeri(2) + (iy-2)*EleSizei(2)
              rzval = 0.0D0
              IF(ppiclf_ndim.GT.2) rzval = centeri(3)
     >                + (iz-2)*EleSizei(3)
              ! Find bin for adjusted rval
              ii    = FLOOR((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
              jj    = FLOOR((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
              kk    = FLOOR((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
              IF (ppiclf_ndim.LT.3) kk = 0
              

              ! This covers ghost exchanged cells for linear periodicity
              ! Maps cells greater than ppiclf bin domain to first bin
              ! Maps cells less than ppiclf bin domain to last bin
              IF (x_per_flag .EQ. 1) THEN
                IF (ii .EQ. ppiclf_n_bins(1)) ii = 0
                IF (ii .EQ. -1) ii = ppiclf_n_bins(1) - 1
              END IF
              IF (y_per_flag .EQ. 1) THEN
                IF (jj .EQ. ppiclf_n_bins(2)) jj = 0
                IF (jj .EQ. -1) jj = ppiclf_n_bins(2) - 1
              END IF
              IF (z_per_flag .EQ. 1) THEN
                IF (kk .EQ. ppiclf_n_bins(3)) kk = 0
                IF (kk .EQ. -1) kk = ppiclf_n_bins(3) - 1
              END IF
              
              ! Ensures duplicate cells don't get sent to same processor
              IF (ii .LT. 0 .OR. ii .GT. ppiclf_n_bins(1)-1) CYCLE
              IF (jj .LT. 0 .OR. jj .GT. ppiclf_n_bins(2)-1) CYCLE
              IF (kk .LT. 0 .OR. kk .GT. ppiclf_n_bins(3)-1) CYCLE


              ! Calculates processor rank
              ndum  = ii + ppiclf_n_bins(1)*jj + 
     >                     ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
              nrank = ndum

                            ppiclf_neltb = ppiclf_neltb + 1
              IF(ppiclf_neltb .GT. PPICLF_LEE) THEN
                PRINT*, '***ERROR*** PPICLF_LEE',PPICLF_LEE, 'in', 
     >           'MapOverlapMesh must be greater than', ppiclf_neltb 
                CALL ppiclf_exittr('Increase PPICLF_LEE$ (MapOverlap)',0.0D0
     >               ,ppiclf_neltb)
              END IF
              ! Stores element to rank mapping.
              ppiclf_er_map(1,ppiclf_neltb) = ie
              ppiclf_er_map(2,ppiclf_neltb) = ppiclf_nid
              ppiclf_er_map(3,ppiclf_neltb) = ndum
              ppiclf_er_map(4,ppiclf_neltb) = nrank
              ppiclf_er_map(5,ppiclf_neltb) = nrank
              ppiclf_er_map(6,ppiclf_neltb) = nrank

!              The loop makes this subroutine 10x slower
!              Replaced with tempCheck below, since cell would 
!              be duplicated in sequential order. It shouldn't happen,
!              so implemented as error vs standard fix in loop.

!              IF (ppiclf_neltb .GT. 1) THEN
!              DO il=1,ppiclf_neltb-1
!                 IF (ppiclf_er_map(1,il) .EQ. ie) THEN
!                 IF (ppiclf_er_map(4,il) .EQ. nrank) THEN
!                    PRINT*, 'AVERY - NELTB Loop remover still used!'
!                    ppiclf_neltb = ppiclf_neltb - 1
!                     CYCLE
!                 END IF
!                 END IF
!              END DO
!              END IF

            END DO !iz
          END DO !iy
        END DO !ix
      END DO !ie
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      DO ie=1,ppiclf_neltb !Number of fluid cells in ppiclf bin domain
       ! These copy all nxyz vertecies since Fortran is column-major
       iee = ppiclf_er_map(1,ie)
       CALL ppiclf_copy(ppiclf_xm1b(1,1,1,1,ie)
     >                 ,ppiclf_xm1bs(1,1,1,1,iee),nxyz)
       CALL ppiclf_copy(ppiclf_xm1b(1,1,1,2,ie)
     >                 ,ppiclf_xm1bs(1,1,1,2,iee),nxyz)
       CALL ppiclf_copy(ppiclf_xm1b(1,1,1,3,ie)
     >                 ,ppiclf_xm1bs(1,1,1,3,iee),nxyz)
      END DO

      ppiclf_neltbb = ppiclf_neltb
      DO ie=1,ppiclf_neltbb
         ! Copies element to rank mapping (integer copy)
         CALL ppiclf_icopy(ppiclf_er_maps(1,ie),ppiclf_er_map(1,ie)
     >             ,PPICLF_LRMAX)
      END DO

      ! GSLIB required info
      ! neltb - number of columns to transfer
      ! PPICLF_LEE - number of columns declared
      ! nl - partl row size (dummy logical variable)
      nl   = 0
      ! nii - ppiclf_er_maps row size declared
      nii  = PPICLF_LRMAX
      ! njj - Row index of ppiclf_er_maps with processor/rank number
      njj  = 6
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      ! nrr - ppiclf_xm1b row size declared
      nrr  = nxyz*3
      ! Defines sorting order
      nkey(1) = 2
      nkey(2) = 1

      CALL pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl,ppiclf_neltb
     >     ,PPICLF_LEE,ppiclf_er_map,nii,partl,nl,ppiclf_xm1b,nrr,njj)
      CALL pfgslib_crystal_tuple_sort    (ppiclf_cr_hndl,ppiclf_neltb
     >     ,ppiclf_er_map,nii,partl,nl,ppiclf_xm1b,nrr,nkey,2)

!*************
! This is only needed for multi-element projection.
! We are currently doing single element projection (hardcoded in rocpicl)
! 
!      DO ie=1,ppiclf_neltb
!      DO k=1,PPICLF_LEZ
!      DO j=1,PPICLF_LEY
!      DO i=1,PPICLF_LEX
!         rxval = ppiclf_xm1b(i,j,k,1,ie)
!         ryval = ppiclf_xm1b(i,j,k,2,ie)
!         rzval = 0.0D0
!         IF(ppiclf_ndim.GT.2) rzval = ppiclf_xm1b(i,j,k,3,ie)
!         
!         ii    = FLOOR((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
!         jj    = FLOOR((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
!         kk    = FLOOR((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
!         IF (ppiclf_ndim.EQ.2) kk = 0
!          IF (ii .EQ. ppiclf_n_bins(1)) ii = ppiclf_n_bins(1) - 1
!          IF (jj .EQ. ppiclf_n_bins(2)) jj = ppiclf_n_bins(2) - 1
!          IF (kk .EQ. ppiclf_n_bins(3)) kk = ppiclf_n_bins(3) - 1
!          IF (ii .EQ. -1) ii = 0
!          IF (jj .EQ. -1) jj = 0
!          IF (kk .EQ. -1) kk = 0
!          ndum  = ii + ppiclf_n_bins(1)*jj + 
!     >                 ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
!
!         ppiclf_modgp(i,j,k,ie,1) = ii
!         ppiclf_modgp(i,j,k,ie,2) = jj
!         ppiclf_modgp(i,j,k,ie,3) = kk
!         ppiclf_modgp(i,j,k,ie,4) = ndum
!   
!      END DO
!      END DO
!      END DO
!      END DO
!**************

      DO ie=1,ppiclf_neltb
         ! Finds minimum and maximum vertex in x,y,z of cell
         ppiclf_xerange(1,1,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,1,ie),nxyz)
         ppiclf_xerange(2,1,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,1,ie),nxyz)
         ppiclf_xerange(1,2,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,2,ie),nxyz)
         ppiclf_xerange(2,2,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,2,ie),nxyz)
         ppiclf_xerange(1,3,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,3,ie),nxyz)
         ppiclf_xerange(2,3,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,3,ie),nxyz)
         
         ! Finds the ppiclf bin that the max/min cell vertex resides in
         ilow  = 
     >     FLOOR((ppiclf_xerange(1,1,ie) - ppiclf_binb(1))/
     >                                             ppiclf_bins_dx(1))
         ihigh = 
     >     FLOOR((ppiclf_xerange(2,1,ie) - ppiclf_binb(1))/
     >                                             ppiclf_bins_dx(1))
         jlow  = 
     >     FLOOR((ppiclf_xerange(1,2,ie) - ppiclf_binb(3))/
     >                                             ppiclf_bins_dx(2))
         jhigh = 
     >     FLOOR((ppiclf_xerange(2,2,ie) - ppiclf_binb(3))/
     >                                             ppiclf_bins_dx(2))
         klow  = 
     >     FLOOR((ppiclf_xerange(1,3,ie) - ppiclf_binb(5))/
     >                                             ppiclf_bins_dx(3))
         khigh = 
     >     FLOOR((ppiclf_xerange(2,3,ie) - ppiclf_binb(5))/
     >                                             ppiclf_bins_dx(3))
         IF (ppiclf_ndim.LT.3) THEN
            klow = 0
            khigh = 0
         END IF

         ! Maps the cell to bin rank range (1,2) and min/max bins in
         ! x,y,z (3-8).  If ppiclf_el_map(1:8,ie) are the same, then 
         ! fluid cell is only in 1 bin.
         ppiclf_el_map(1,ie) = ilow  + ppiclf_n_bins(1)*jlow  
     >                         + ppiclf_n_bins(1)*ppiclf_n_bins(2)*klow
         ppiclf_el_map(2,ie) = ihigh + ppiclf_n_bins(1)*jhigh 
     >                         + ppiclf_n_bins(1)*ppiclf_n_bins(2)*khigh
         ppiclf_el_map(3,ie) = ilow
         ppiclf_el_map(4,ie) = ihigh
         ppiclf_el_map(5,ie) = jlow
         ppiclf_el_map(6,ie) = jhigh
         ppiclf_el_map(7,ie) = klow
         ppiclf_el_map(8,ie) = khigh
      END DO

      IF (icalld .EQ. 0) THEN 
         icalld = icalld + 1
         CALL ppiclf_prints('   *Begin mpi_comm_split$')
            CALL mpi_comm_split(ppiclf_comm
     >                         ,ppiclf_nid
     >                         ,0
     >                         ,ppiclf_comm_nid
     >                         ,ierr)
         CALL ppiclf_prints('    End mpi_comm_split$')
         CALL ppiclf_io_OutputDiagGrid
      END IF

      RETURN
      END 
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_comm_InitOverlapMesh(ncell,lx1,ly1,lz1,
     >                                       xgrid,ygrid,zgrid)
     > bind(C, name="ppiclc_comm_InitOverlapMesh")
#else
      subroutine ppiclf_comm_InitOverlapMesh(ncell,lx1,ly1,lz1,
     >                                       xgrid,ygrid,zgrid)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 ncell
      integer*4 lx1
      integer*4 ly1
      integer*4 lz1
      real*8    xgrid(*)
      real*8    ygrid(*)
      real*8    zgrid(*)
!
! External:
!
      integer*4 nxyz, i, j, ie
      integer*4 k, jj, icont
!
      ppiclf_overlap = .true.

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitOverlap$',0.0d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitOverlap$'
     >                  ,0.0d0,0)

      if (ncell .gt. PPICLF_LEE .or. ncell .lt. 0) then
        PRINT*, '***ERROR*** PPICLF_LEE', PPICLF_LEE, 'in', 
     > 'InitMapOverlapMesh must be greater than', ncell 
        call ppiclf_exittr('Increase LEE in InitOverlap$',0.0d0,ncell)
      endif
      if (lx1 .ne. PPICLF_LEX) 
     >   call ppiclf_exittr('LX1 != LEX in InitOverlap$',0.0d0,ncell)
      if (ly1 .ne. PPICLF_LEY)
     >   call ppiclf_exittr('LY1 != LEY in InitOverlap$',0.0d0,ncell)
      if (lz1 .ne. PPICLF_LEZ)
     >   call ppiclf_exittr('LZ1 != LEZ in InitOverlap$',0.0d0,ncell)

      ppiclf_nee = ncell
      nxyz       = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ

      do ie=1,ppiclf_nee
         ! TLJ changing loop structure
         !do i=1,nxyz
         !   j = (ie-1)*nxyz + i
         !   ppiclf_xm1bs(i,1,1,1,ie) = xgrid(j)
         !   ppiclf_xm1bs(i,1,1,2,ie) = ygrid(j)
         !   ppiclf_xm1bs(i,1,1,3,ie) = zgrid(j)
         !enddo
         icont = 0
         do k=1,PPICLF_LEZ
         do j=1,PPICLF_LEY
         do i=1,PPICLF_LEX
            icont = icont + 1
            jj = (ie-1)*nxyz + icont
            ppiclf_xm1bs(i,j,k,1,ie) = xgrid(jj)
            ppiclf_xm1bs(i,j,k,2,ie) = ygrid(jj)
            ppiclf_xm1bs(i,j,k,3,ie) = zgrid(jj)
         enddo
         enddo
         enddo
      enddo
      
      call ppiclf_solve_InitSolve

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_FindParticle
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 ix, iy, iz, i, ii, jj, kk, ndum, nrank
!
      ix = 1
      iy = 2
      iz = 1
      if (ppiclf_ndim.eq.3)
     >iz = 3

      do i=1,ppiclf_npart
         ! check if particles are greater or less than binb bounds....
         ii  = floor((ppiclf_y(ix,i)-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
         jj  = floor((ppiclf_y(iy,i)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
         kk  = floor((ppiclf_y(iz,i)-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
         if (ppiclf_ndim .lt. 3) kk = 0
         ndum  = ii + ppiclf_n_bins(1)*jj + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
         nrank = ndum

         ppiclf_iprop(8,i)  = ii
         ppiclf_iprop(9,i)  = jj
         ppiclf_iprop(10,i) = kk
         ppiclf_iprop(11,i) = ndum

         ppiclf_iprop(3,i)  = nrank ! where particle is actually moved
         ppiclf_iprop(4,i)  = nrank ! where particle is actually moved
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_MoveParticle
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      logical partl    
      integer*4 lrf
      parameter(lrf = PPICLF_LRS*4 + PPICLF_LRP + PPICLF_LRP2
     >       + PPICLF_LRP3 + PPICLF_LRP4 + PPICLF_LRP5)
      real*8 rwork(lrf,PPICLF_LPART)
      integer*4 i, ic, j0
!

      do i=1,ppiclf_npart
         ic = 1
         call ppiclf_copy(rwork(ic,i),ppiclf_y(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_y1((i-1)*PPICLF_LRS+1)
     >                   ,PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_ydot(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_ydotc(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop(1,i),PPICLF_LRP)
         ic = ic + PPICLF_LRP
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop2(1,i),PPICLF_LRP2)
         ic = ic + PPICLF_LRP2
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop3(1,i),PPICLF_LRP3)
         ic = ic + PPICLF_LRP3
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop4(1,i),PPICLF_LRP4)
         ic = ic + PPICLF_LRP4
         call ppiclf_copy(rwork(ic,i),ppiclf_force(1,i),PPICLF_LRP5)
      enddo

      j0 = 4
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >                                  ,ppiclf_npart,PPICLF_LPART
     >                                  ,ppiclf_iprop,PPICLF_LIP
     >                                  ,partl,0
     >                                  ,rwork,lrf
     >                                  ,j0)

      if (ppiclf_npart .gt. PPICLF_LPART .or. ppiclf_npart .lt. 0)
     >   call ppiclf_exittr('Increase LPART$',0.0d0,ppiclf_npart)

      do i=1,ppiclf_npart
         ic = 1
         call ppiclf_copy(ppiclf_y(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_y1((i-1)*PPICLF_LRS+1),rwork(ic,i)
     >                   ,PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_ydot(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_ydotc(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_rprop(1,i),rwork(ic,i),PPICLF_LRP)
         ic = ic + PPICLF_LRP
         call ppiclf_copy(ppiclf_rprop2(1,i),rwork(ic,i),PPICLF_LRP2)
         ic = ic + PPICLF_LRP2
         call ppiclf_copy(ppiclf_rprop3(1,i),rwork(ic,i),PPICLF_LRP3)
         ic = ic + PPICLF_LRP3
         call ppiclf_copy(ppiclf_rprop4(1,i),rwork(ic,i),PPICLF_LRP4)
         ic = ic + PPICLF_LRP4
         call ppiclf_copy(ppiclf_force(1,i),rwork(ic,i),PPICLF_LRP5)
      enddo
        
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 xdlen,ydlen,zdlen,rxdrng(3),rxnew(3), rfac, rxval, ryval,
     >       rzval, rxl, ryl, rzl, rxr, ryr, rzr, distchk, dist
      integer*4 iadd(3),gpsave(27)
      real*8 map(PPICLF_LRP_PRO)
      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >           nfacegp, nedgegp, ncornergp, iperiodicx, iperiodicy,
     >           iperiodicz, jx, jy, jz, ip, idum, iip, jjp, kkp, ii1,
     >           jj1, kk1, iig, jjg, kkg, iflgx, iflgy, iflgz,
     >           isave, iflgsum, ndumn, nrank, ibctype, i, ifc, ist, j,
     >           k
!

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/ -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                   -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (ppiclf_ndim .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)

! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      jx    = 1
      jy    = 2
      jz    = 3

      ! Thierry - we do not assign the bins to be as big as
      !           the periodic domain in x/y directions anymore. only in z. 
      
      !xdlen = ppiclf_binb(2) - ppiclf_binb(1) ! when bins = periodic domain
      !ydlen = ppiclf_binb(4) - ppiclf_binb(3) ! when bins = periodic domain
      
      ! Thierry - this works whether the bins are as big as periodic domain, or not.
      xdlen = ppiclf_xdrange(2,1) - ppiclf_xdrange(1,1)
      ydlen = ppiclf_xdrange(2,2) - ppiclf_xdrange(1,2)
      
      zdlen = -1.
      if (ppiclf_ndim .gt. 2) 
!     >   zdlen = ppiclf_binb(6) - ppiclf_binb(5)
      ! Thierry - this works whether the bins are as big as periodic domain, or not.
     >   zdlen = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
      if (iperiodicx .ne. 0) xdlen = -1
      if (iperiodicy .ne. 0) ydlen = -1
      if (iperiodicz .ne. 0) zdlen = -1

      rxdrng(1) = xdlen
      rxdrng(2) = ydlen
      rxdrng(3) = zdlen

      ppiclf_npart_gp = 0

      rfac = 1.0d0

      do ip=1,ppiclf_npart

         call ppiclf_user_MapProjPart(map,ppiclf_y(1,ip)
     >         ,ppiclf_ydot(1,ip),ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip),
     >                                              ppiclf_rprop4(1,ip))

c        idum = 1
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
c        idum = 2
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
c        idum = 3
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)

         idum = 0
         do j=1,PPICLF_LRS
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip)
         enddo
         idum = PPICLF_LRS
         do j=1,PPICLF_LRP
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip)
         enddo
         idum = PPICLF_LRS+PPICLF_LRP
         do j=1,PPICLF_LRP_PRO
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = map(j)
         enddo

         rxval = ppiclf_cp_map(1,ip)
         ryval = ppiclf_cp_map(2,ip)
         rzval = 0.0d0
         if (ppiclf_ndim .gt. 2) rzval = ppiclf_cp_map(3,ip)

         iip    = ppiclf_iprop(8,ip)
         jjp    = ppiclf_iprop(9,ip)
         kkp    = ppiclf_iprop(10,ip)

         rxl = ppiclf_binb(1) + ppiclf_bins_dx(1)*iip
         rxr = rxl + ppiclf_bins_dx(1)
         ryl = ppiclf_binb(3) + ppiclf_bins_dx(2)*jjp
         ryr = ryl + ppiclf_bins_dx(2)
         rzl = 0.0d0
         rzr = 0.0d0
         if (ppiclf_ndim .gt. 2) then
            rzl = ppiclf_binb(5) + ppiclf_bins_dx(3)*kkp
            rzr = rzl + ppiclf_bins_dx(3)
         endif

         isave = 0

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) 
            jj1 = jjp + el_face_num(ist+2)
            kk1 = kkp + el_face_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
               iflgx = 1
               iig =modulo(iig,ppiclf_n_bins(1))
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
               iflgy = 1
               jjg =modulo(jjg,ppiclf_n_bins(2))
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
               iflgz = 1  
               kkg =modulo(kkg,ppiclf_n_bins(3))
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 111
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1

            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)

            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  111 continue
         enddo

         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
               iflgx = 1
               iig =modulo(iig,ppiclf_n_bins(1))
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
               iflgy = 1
               jjg =modulo(jjg,ppiclf_n_bins(2))
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
               iflgz = 1  
               kkg =modulo(kkg,ppiclf_n_bins(3))
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 222
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1

            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)

            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  222 continue
         enddo

         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
               iflgx = 1
               iig =modulo(iig,ppiclf_n_bins(1))
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
               iflgy = 1
               jjg =modulo(jjg,ppiclf_n_bins(2))
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
               iflgz = 1  
               kkg =modulo(kkg,ppiclf_n_bins(3))
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 333
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1

            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)

            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  333 continue
         enddo

      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_AngularCreateGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 xdlen,ydlen,zdlen,rxdrng(3),rxnew(3), rfac, rxval, ryval,
     >       rzval, rxl, ryl, rzl, rxr, ryr, rzr, distchk, dist
      integer*4 iadd(3),gpsave(27)
      real*8 map(PPICLF_LRP_PRO)
      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >           nfacegp, nedgegp, ncornergp, iperiodicx, iperiodicy,
     >           iperiodicz, jx, jy, jz, ip, idum, iip, jjp, kkp, ii1,
     >           jj1, kk1, iig, jjg, kkg, iflgx, iflgy, iflgz,
     >           isave, iflgsum, ndumn, nrank, ibctype, i, ifc, ist, j,
     >           k
      ! 08/27/24 - Thierry - added for angular periodicty starts here
      real*8 alpha
      integer*4 xrank, yrank, zrank
      ! 08/27/24 - Thierry - added for angular periodicty ends here
      ! 09/26/24 - Thierry - added for angular periodicty starts here
      real*8 dist1, dist2
      ! 09/26/24 - Thierry - added for angular periodicty ends here
!

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/ -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                   -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (ppiclf_ndim .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)

! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      jx    = 1
      jy    = 2
      jz    = 3

      ! Thierry - we dont use xdlen and ydlen in this algorithm. no need to modify them.
      xdlen = ppiclf_binb(2) - ppiclf_binb(1)
      ydlen = ppiclf_binb(4) - ppiclf_binb(3)
      zdlen = -1.
      if (ppiclf_ndim .gt. 2) 
!     >   zdlen = ppiclf_binb(6) - ppiclf_binb(5)
     >   zdlen = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
      if (iperiodicx .ne. 0) xdlen = -1
      if (iperiodicy .ne. 0) ydlen = -1
      if (iperiodicz .ne. 0) zdlen = -1

      rxdrng(1) = xdlen
      rxdrng(2) = ydlen
      rxdrng(3) = zdlen

      ppiclf_npart_gp = 0

      rfac = 1.0d0

      do ip=1,ppiclf_npart

         call ppiclf_user_MapProjPart(map,ppiclf_y(1,ip)
     >         ,ppiclf_ydot(1,ip),ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip))

c        idum = 1
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
c        idum = 2
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
c        idum = 3
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)

         idum = 0
         do j=1,PPICLF_LRS
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip) ! ppiclf_y(PPICLF_JX/ JY/ JZ/ JVX/ JVY/ JVZ/ JT, ip)
         enddo
         idum = PPICLF_LRS
         do j=1,PPICLF_LRP
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip) ! ppiclf_rprop(PPICLF_R_JRHOP/ R_JRHOF/ .../ R_WDOTZ, ip)
         enddo
         idum = PPICLF_LRS+PPICLF_LRP
         do j=1,PPICLF_LRP_PRO
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = map(j) ! map(PPICLF_P_JPHIP/ JFX/ .../ JPHIPW) - these are found in ppiclf_user_MapProjPart
         enddo

         rxval = ppiclf_cp_map(1,ip) ! ppiclf_y(PPICLF_JX,ip)
         ryval = ppiclf_cp_map(2,ip) ! ppiclf_y(PPICLF_JY,ip)
         rzval = 0.0d0
         if (ppiclf_ndim .gt. 2) rzval = ppiclf_cp_map(3,ip) ! ppiclf_y(PPICLF_JZ,ip)

         iip    = ppiclf_iprop(8,ip) ! ith coordinate of bin
         jjp    = ppiclf_iprop(9,ip) ! jth coordinate of bin
         kkp    = ppiclf_iprop(10,ip) ! kth coordinate of bin

         rxl = ppiclf_binb(1) + ppiclf_bins_dx(1)*iip ! min x of bin
         rxr = rxl + ppiclf_bins_dx(1)                ! max x of bin
         ryl = ppiclf_binb(3) + ppiclf_bins_dx(2)*jjp
         ryr = ryl + ppiclf_bins_dx(2)
         rzl = 0.0d0
         rzr = 0.0d0
         if (ppiclf_ndim .gt. 2) then
            rzl = ppiclf_binb(5) + ppiclf_bins_dx(3)*kkp
            rzr = rzl + ppiclf_bins_dx(3)
         endif

         isave = 0

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) 
            jj1 = jjp + el_face_num(ist+2)
            kk1 = kkp + el_face_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)

            if (ang_case==1) then  ! for wedge geometry

               ! Thierry - I dont think it's code efficient to call this subroutine
               !           for every particle, every ghost face, at every time step
               !           I'm wondering if it's better if we make the plane values 
               !           as global values that are initialized in the beginning 
            
               call ppiclf_solve_InitAngularPlane(ip,
     >                                 ang_per_rin  , ang_per_rout  ,
     >                                 ang_per_angle, ang_per_xangle,
     >                                 dist1, dist2)
               if ((dist .gt. distchk).and.(dist1.gt.distchk)
     >           .and.(dist2.gt.distchk)) cycle
            else
               if (dist .gt. distchk) cycle
            endif

            iflgx = 0
            iflgy = 0
            iflgz = 0
!-----------------------------------------------------------------------
            ! 08/27/24 - Thierry - modification for angular periodicty starts here

               ! angle between particle and x-axis
                alpha = atan2(ppiclf_y(PPICLF_JY,ip), 
     >                        ppiclf_y(PPICLF_JX,ip))
                

                call ppiclf_solve_InvokeAngularPeriodic(ip, 
     >                                                  ang_per_flag,
     >                                                  alpha,         
     >                                                  ang_per_angle,  
     >                                                  ang_per_xangle, 
     >                                                  0)

              ! Thierry - this is how FindParticle implements it
              ! need to find a way to make the code deal with negative xrot values

            xrank = iig ; yrank=jjg; zrank = kkg
            ! Thierry - previously placed before the CheckPeriodicBC call, had to move them for the periodic check
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval ! z-coordinate does not change when angular periodicity is invoked
            
            ! Angular periodicity check in x- and y-directions
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
              iflgx = 1
              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
              if (iperiodicx .ne. 0) cycle
              iig = xrank
              jjg = yrank
            end if
            
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
              iflgy = 1
              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
              if (iperiodicy .ne. 0) cycle
              iig = xrank
              jjg = yrank
            end if
            
            ! Linear periodicity check in z-direction
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
              iflgz = 1
              kkg =modulo(kkg,ppiclf_n_bins(3))
              if (iperiodicz .ne. 0) cycle
              ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
              ! rxdrng(3) = -1.0  if not periodic in Z
              if (rxdrng(3) .gt. 0) then 
                if (iadd(3) .ge. ppiclf_n_bins(3)) then ! particle leaving from max z-face
                  rxnew(3) = rxnew(3) - rxdrng(3)
                elseif (iadd(3) .lt. 0) then ! particle leaving from min z-face
                  rxnew(3) = rxnew(3) + rxdrng(3)
                end if ! iadd
              end if ! rxrdrng
            else ! z-linear periodicity not applicable
              kkg = zrank
            end if ! kkg
            
            iflgsum = iflgx + iflgy + iflgz
            ndumn  = iig + ppiclf_n_bins(1)*jjg + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
             nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            ! 08/27/24 - Thierry - modification for angular periodicty ends here
!-----------------------------------------------------------------------

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 111
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
            
            rxnew(1) = xrot(1)
            rxnew(2) = xrot(2)
            ppiclf_cp_map(4,ip) = vrot(1)
            ppiclf_cp_map(5,ip) = vrot(2)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ! Thierry - we don't need ppiclf_comm_CheckPeriodicBC anymore for the angular periodic ghost algorithm
            !           as this is now taken care of when anticipating where the particle might be when calling
            !           ppiclf_comm_InvokeAngularPeriodic
            !           we only need to assign xr and vr to ppiclf_rprop_gp

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1) ! ppiclf_y(PPICLF_JX, ip) for the periodic ghost particle
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2) ! JY
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3) ! JZ
            
            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  111 continue
         enddo

         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)

            if (ang_case==1) then  ! for wedge geometry

               call ppiclf_solve_InitAngularPlane(ip,
     >                                 ang_per_rin  , ang_per_rout  ,
     >                                 ang_per_angle, ang_per_xangle,
     >                                 dist1, dist2)
               if ((dist .gt. distchk).and.(dist1.gt.distchk)
     >           .and.(dist2.gt.distchk)) cycle
            else
               if (dist .gt. distchk) cycle
            endif

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
!-----------------------------------------------------------------------
            ! 08/27/24 - Thierry - modification for angular periodicty starts here

               ! angle between particle and x-axis
                alpha = atan2(ppiclf_y(PPICLF_JY,ip), 
     >                        ppiclf_y(PPICLF_JX,ip))
                

                call ppiclf_solve_InvokeAngularPeriodic(ip, 
     >                                                  ang_per_flag,
     >                                                  alpha,         
     >                                                  ang_per_angle,  
     >                                                  ang_per_xangle, 
     >                                                  0)

              ! Thierry - this is how FindParticle implements it
              ! need to find a way to make the code deal with negative xrot values

            xrank = iig ; yrank=jjg; zrank = kkg
            ! Thierry - previously placed before the CheckPeriodicBC call, had to move them for the periodic check
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval ! z-coordinate does not change when angular periodicity is invoked
            
            ! Angular periodicity check in x- and y-directions
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
              iflgx = 1
              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
              if (iperiodicx .ne. 0) cycle
              iig = xrank
              jjg = yrank
            end if
            
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
              iflgy = 1
              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
              if (iperiodicy .ne. 0) cycle
              iig = xrank
              jjg = yrank
            end if
            
            ! Linear periodicity check in z-direction
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
              iflgz = 1
              kkg =modulo(kkg,ppiclf_n_bins(3))
              if (iperiodicz .ne. 0) cycle
              ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
              ! rxdrng(3) = -1.0  if not periodic in Z
              if (rxdrng(3) .gt. 0) then ! particle leaving from max z-face
                if (iadd(3) .ge. ppiclf_n_bins(3)) then
                  rxnew(3) = rxnew(3) - rxdrng(3)
                elseif (iadd(3) .lt. 0) then
                  rxnew(3) = rxnew(3) + rxdrng(3)
                end if ! iadd
              end if ! rxrdrng
            else ! z-linear periodicity not applicable
              kkg = zrank
            end if ! kkg
            
            iflgsum = iflgx + iflgy + iflgz
            ndumn  = iig + ppiclf_n_bins(1)*jjg + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
             nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            ! 08/27/24 - Thierry - modification for angular periodicty ends here
!-----------------------------------------------------------------------

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 222
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz

            rxnew(1) = xrot(1)
            rxnew(2) = xrot(2)
            ppiclf_cp_map(4,ip) = vrot(1)
            ppiclf_cp_map(5,ip) = vrot(2)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ! Thierry - we don't need ppiclf_comm_CheckPeriodicBC anymore for the angular periodic ghost algorithm
            !           as this is now taken care of when anticipating where the particle might be when calling
            !           ppiclf_comm_InvokeAngularPeriodic
            !           we only need to assign xr and vr to ppiclf_rprop_gp

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1) ! ppiclf_y(PPICLF_JX, ip) for the periodic ghost particle
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2) ! JY
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3) ! JZ
            
            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  222 continue
         enddo

         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)

            if (ang_case==1) then  ! for wedge geometry
            
               call ppiclf_solve_InitAngularPlane(ip,
     >                                 ang_per_rin  , ang_per_rout  ,
     >                                 ang_per_angle, ang_per_xangle,
     >                                 dist1, dist2)
               if ((dist .gt. distchk).and.(dist1.gt.distchk)
     >           .and.(dist2.gt.distchk)) cycle
            else
               if (dist .gt. distchk) cycle
            endif

            iflgx = 0
            iflgy = 0
            iflgz = 0

!-----------------------------------------------------------------------
            ! 08/27/24 - Thierry - modification for angular periodicty starts here

               ! angle between particle and x-axis
                alpha = atan2(ppiclf_y(PPICLF_JY,ip), 
     >                        ppiclf_y(PPICLF_JX,ip))
                

                call ppiclf_solve_InvokeAngularPeriodic(ip, 
     >                                                  ang_per_flag,
     >                                                  alpha,         
     >                                                  ang_per_angle,  
     >                                                  ang_per_xangle, 
     >                                                  0)

              ! Thierry - this is how FindParticle implements it
              ! need to find a way to make the code deal with negative xrot values

            xrank = iig ; yrank=jjg; zrank = kkg
            ! Thierry - previously placed before the CheckPeriodicBC call, had to move them for the periodic check
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval ! z-coordinate does not change when angular periodicity is invoked
            
            ! Angular periodicity check in x- and y-directions
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
              iflgx = 1
              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
              if (iperiodicx .ne. 0) cycle
              iig = xrank
              jjg = yrank
            end if
            
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
              iflgy = 1
              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
              if (iperiodicy .ne. 0) cycle
              iig = xrank
              jjg = yrank
            end if
            
            ! Linear periodicity check in z-direction
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
              iflgz = 1
              kkg =modulo(kkg,ppiclf_n_bins(3))
              if (iperiodicz .ne. 0) cycle
              ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
              ! rxdrng(3) = -1.0  if not periodic in Z
              if (rxdrng(3) .gt. 0) then ! particle leaving from max z-face
                if (iadd(3) .ge. ppiclf_n_bins(3)) then
                  rxnew(3) = rxnew(3) - rxdrng(3)
                elseif (iadd(3) .lt. 0) then
                  rxnew(3) = rxnew(3) + rxdrng(3)
                end if ! iadd
              end if ! rxrdrng
            else ! z-linear periodicity not applicable
              kkg = zrank
            end if ! kkg
            
            iflgsum = iflgx + iflgy + iflgz
            ndumn  = iig + ppiclf_n_bins(1)*jjg + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
             nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            ! 08/27/24 - Thierry - modification for angular periodicty ends here
!-----------------------------------------------------------------------
            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 333
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz

            rxnew(1) = xrot(1)
            rxnew(2) = xrot(2)
            ppiclf_cp_map(4,ip) = vrot(1)
            ppiclf_cp_map(5,ip) = vrot(2)

            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ! Thierry - we don't need ppiclf_comm_CheckPeriodicBC anymore for the angular periodic ghost algorithm
            !           as this is now taken care of when anticipating where the particle might be when calling
            !           ppiclf_comm_InvokeAngularPeriodic
            !           we only need to assign xr and vr to ppiclf_rprop_gp

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1) ! ppiclf_y(PPICLF_JX, ip) for the periodic ghost particle
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2) ! JY
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3) ! JZ

            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  333 continue
         enddo

      enddo ! ip 

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      real*8 rxdrng(3)
      integer*4 iadd(3)
!
! Input/Output:
!
      real*8 rxnew(3)
!
      ! rxdrng(1) = ppiclf_xdrange(2,1) - ppiclf_xdrange(1,1)
      ! rxdrng(1) = -1.0  if not periodic in X
      ! particle leaving from max x periodic face
      if (rxdrng(1) .gt. 0 ) then
      if (iadd(1) .ge. ppiclf_n_bins(1)) then
         rxnew(1) = rxnew(1) - rxdrng(1)
         goto 123
      endif
      endif
      ! particle leaving from min x periodic face
      if (rxdrng(1) .gt. 0 ) then
      if (iadd(1) .lt. 0) then
         rxnew(1) = rxnew(1) + rxdrng(1)
         goto 123
      endif
      endif

  123 continue    
      ! rxdrng(2) = ppiclf_xdrange(2,2) - ppiclf_xdrange(1,2)
      ! rxdrng(2) = -1.0  if not periodic in Y
      ! particle leaving from max y periodic face
      if (rxdrng(2) .gt. 0 ) then
      if (iadd(2) .ge. ppiclf_n_bins(2)) then
         rxnew(2) = rxnew(2) - rxdrng(2)
         goto 124
      endif
      endif
      if (rxdrng(2) .gt. 0 ) then
      ! particle leaving from min y periodic face
      if (iadd(2) .lt. 0) then
         rxnew(2) = rxnew(2) + rxdrng(2)
         goto 124
      endif
      endif
  124 continue

      if (ppiclf_ndim .gt. 2) then
        ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
        ! rxdrng(3) = -1.0  if not periodic in Z
      ! particle leaving from max z periodic face
         if (rxdrng(3) .gt. 0 ) then
         if (iadd(3) .ge. ppiclf_n_bins(3)) then
            rxnew(3) = rxnew(3) - rxdrng(3)
            goto 125
         endif
         endif
      ! particle leaving from min z periodic face
         if (rxdrng(3) .gt. 0 ) then
         if (iadd(3) .lt. 0) then
            rxnew(3) = rxnew(3) + rxdrng(3)
            goto 125
         endif
         endif
      endif
  125 continue

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_CheckAngularBC(xrank, yrank, zrank)
!
      implicit none
!
      include "PPICLF"
!
! Local:
!
      integer*4 xrank, yrank, zrank
!
! Output:
!

      SELECT CASE (ang_case)
        CASE(1) ! general wedge ; 0 <= angle < 90
!          print*, "Wedge CheckAngularBC"
          xrank  = floor((xrot(1)-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
          yrank  = floor((xrot(2)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
          zrank  = floor((xrot(3)-ppiclf_binb(5))/ppiclf_bins_dx(3))

        CASE(2) ! quarter cylinder ; angle = 90
!          print*, "Quarter Cylinder CheckAngularBC"
          xrank  = floor((abs(xrot(1))-ppiclf_binb(1))
     >                    /ppiclf_bins_dx(1)) 
          yrank  = floor((abs(xrot(2))-ppiclf_binb(3))
     >                   /ppiclf_bins_dx(2)) 
          zrank  = floor((xrot(3)-ppiclf_binb(5))/ppiclf_bins_dx(3))

        CASE(3) ! half cylinder ; angle = 180
          print*, "Half Cylinder CheckAngularBC"

        CASE DEFAULT
            call ppiclf_exittr('Invalid Ghost Rotational Case!$',
     >       0.0d0 ,ppiclf_nid)
          END SELECT

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_MoveGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      logical partl         
!
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >                                  ,ppiclf_npart_gp,PPICLF_LPART_GP
     >                                  ,ppiclf_iprop_gp,PPICLF_LIP_GP
     >                                  ,partl,0
     >                                  ,ppiclf_rprop_gp,PPICLF_LRP_GP
     >                                  ,1)

      return
      end
c----------------------------------------------------------------------
