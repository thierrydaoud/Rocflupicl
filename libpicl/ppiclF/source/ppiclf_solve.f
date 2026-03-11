      SUBROUTINE ppiclf_solve_Initialize(xi1,xpmin,xpmax,
     >        yi1,ypmin,ypmax,zi1,zpmin,zpmax,
     >        ai1,apa,apxa,aprin,aprout)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4 xi1, yi1, zi1, ai1
      REAL*8 xpmin,xpmax,ypmin,ypmax,zpmin,zpmax,
     >       apa,apxa,aprin,aprout
      REAL*8 pi, angled

!
! Code:
!
      ! Called by rocpicl/PICL_TEMP_InitSolver.F90
      ! xdrange adjusts the bin boundaries.  If periodic
      ! in a direction, then ppiclf bin bounds are set to 
      ! fluid domain boundary in that dimension.
      ! User must input minimums and maximums to match fluid
      ! domain boundaries.
!*** future work - can we automate finding fluid domain boundaries?

      ! Linear X-Periodicity
      x_per_flag = xi1
      IF(x_per_flag.EQ.1) THEN
        IF(xpmin .ge. xpmax) CALL ppiclf_exittr('PeriodicX 
     >      must have xmin < xmax$',xpmin,0)
        ppiclf_iperiodic(1) = 0
        x_per_min = xpmin
        x_per_max = xpmax
        ppiclf_xdrange(1,1) = xpmin
        ppiclf_xdrange(2,1) = xpmax
      END IF

      ! Linear Y-Periodicity
      y_per_flag = yi1
      IF(y_per_flag.EQ.1) THEN
        IF(ypmin .ge. ypmax) CALL ppiclf_exittr('PeriodicY 
     >     must have ymin < ymax$',ypmin,0)
        ppiclf_iperiodic(2) = 0
        y_per_min = ypmin
        y_per_max = ypmax
        ppiclf_xdrange(1,2) = ypmin
        ppiclf_xdrange(2,2) = ypmax
      END IF

      ! Linear Z-Periodicity
      z_per_flag = zi1
      IF(z_per_flag.EQ.1) THEN
        IF(zpmin .ge. zpmax) CALL ppiclf_exittr('PeriodicZ 
     >     must have zmin < zmax$',zpmin,0)
        ppiclf_iperiodic(3) = 0
        z_per_min = zpmin
        z_per_max = zpmax
        ppiclf_xdrange(1,3) = zpmin
        ppiclf_xdrange(2,3) = zpmax
      END IF


      ! Angular Periodicity
      ang_per_flag = ai1
      IF(ang_per_flag.EQ.1) THEN
        ppiclf_iperiodic(1) = 0 ! X-Periodicity
        ppiclf_iperiodic(2) = 0 ! Y-Periodicity
        ang_per_angle  = apa
        ang_per_xangle = apxa
        ang_per_rin    = aprin
        ang_per_rout   = aprout
      END IF

      ! User cannot initialize X/Y-Periodicity with Angular Periodicity
      if(((x_per_flag.eq.1).or.(y_per_flag.eq.1))
     >                     .and.(ang_per_flag.eq.1))
     >   call ppiclf_exittr('PPICLF: Invalid Periodicity choice$',0,0)

      ! Thierry - compute ang_case

      pi = ACOS(-1.0)
      angled = ang_per_angle * 180.0d0 / pi ! store angle value in degrees

      IF(ang_per_flag.EQ.0) THEN
         ang_case = 0 ! standard geometry
      ELSE
         IF(angled .lt. 90.0)        ang_case = 1 ! general wedge
         IF(NINT(angled) .EQ. 90.0)  ang_case = 2 ! quarter cylinder
         IF(NINT(angled) .EQ. 180.0) ang_case = 3 ! half cylinder
      END IF

      IF(ppiclf_nid.EQ.0 .AND. ang_case.NE.0) THEN
         PRINT*, " "
         PRINT*, " ======================================="
         PRINT*, " "
         PRINT*, "!!! PPICLF Angular Periodicity Initialized !!!!"
         PRINT*, "  Angular periodicity flag =", ang_per_flag
         !IF(ang_per_flag.EQ.0) THEN
         !   PRINT*, "  Init Angular- ang_case =", ang_case
         !ELSE
         PRINT*, "  Init Angular- angle =", ang_per_angle
         PRINT*, "  Init Angular- angled =", angled
         PRINT*, "  Init Angular- nint(angled) =", NINT(angled)
         PRINT*, "  Init Angular- ang_case =", ang_case
         !END IF
         PRINT*, " "
         PRINT*, " ======================================="
         PRINT*, " "
      END IF

      RETURN
      END
!
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_AddParticles(npart,y,rprop)
     > bind(C, name="ppiclc_solve_AddParticles")
#else
      subroutine ppiclf_solve_AddParticles(npart,y,rprop)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input: 
!
      integer*4  npart
      real*8     y(*)
      real*8     rprop(*)
!
! Internal:
!
      integer*4 ppiclf_iglsum,ntotal
      external ppiclf_iglsum
!

      call ppiclf_prints('   *Begin AddParticles$')

      if (ppiclf_npart+npart .gt. PPICLF_LPART .or. npart .lt. 0)
     >   call ppiclf_exittr('Invalid number of particles$',
     >                      0.0D0,ppiclf_npart+npart)

      call ppiclf_printsi('      -Begin copy particles$',npart)

      ! First, append arrays onto existing arrays
      call ppiclf_copy(ppiclf_y(1,ppiclf_npart+1),
     >                 y,
     >                 npart*PPICLF_LRS)
      call ppiclf_copy(ppiclf_rprop(1,ppiclf_npart+1),
     >                 rprop,
     >                 npart*PPICLF_LRP)
      ppiclf_npart = ppiclf_npart + npart

      call ppiclf_printsi('      -Begin copy particles$',ppiclf_npart)

      if (.not. PPICLF_RESTART) then
         call ppiclf_prints('      -Begin ParticleTag$')
            call ppiclf_solve_SetParticleTag(npart)
         call ppiclf_prints('       End ParticleTag$')
      ENDif

      if (ppiclf_iglsum(ppiclf_npart,1).gt.0) then
         call ppiclf_prints('      -Begin CreateBin$')
            call ppiclf_comm_CreateBin
         call ppiclf_prints('       End CreateBin$')

         call ppiclf_prints('      -Begin FindParticle$')
            call ppiclf_comm_FindParticle
         call ppiclf_prints('       End FindParticle$')

         call ppiclf_prints('      -Begin MoveParticle$')
            call ppiclf_comm_MoveParticle
         call ppiclf_prints('       End MoveParticle$')

      ENDif

      call ppiclf_prints('    End AddParticles$')

      END
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitParticle(imethod,ndim,iendian,npart,y,
     >                                     rprop,filt2,filt3)
     > bind(C, name="ppiclc_solve_InitParticle")
#else
      subroutine ppiclf_solve_InitParticle(imethod,ndim,iendian,npart,y,
     >                                     rprop,filt2,filt3)
#endif
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input: 
!
      integer*4  imethod
      integer*4  ndim
      integer*4  iendian
      integer*4  npart
      integer*4  ierr
      real*8     y(*)
      real*8     rprop(*)
      real*8     filt2,filt3
!
      call mpi_barrier(ppiclf_comm,ierr)
!
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitParticle$',0.0d0
     >   ,ppiclf_nid)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter must be before InitParticle$',0.0d0
     >                  ,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.0d0
     >                  ,0)

      call ppiclf_prints('*Begin InitParticle$')

         call ppiclf_prints('   *Begin InitParam$')
            call ppiclf_solve_InitParam(imethod,ndim,iendian)
            ppiclf_d2chk(2) = filt2
            ppiclf_d2chk(3) = filt3
            ppiclf_d2chk(1) = max( ppiclf_d2chk(2),ppiclf_d2chk(3) )

            ! TLJ added 12/21/2024
            if (ppiclf_nid==0) then
               print*,'TLJ checking d2chk(2) = ',ppiclf_d2chk(2)
               print*,'TLJ checking d2chk(3) = ',ppiclf_d2chk(3)
               print*,'TLJ checking d2chk(1) = ',ppiclf_d2chk(1)
            ENDif

         call ppiclf_prints('    End InitParam$')

         if (.not. PPICLF_RESTART) then
            call ppiclf_prints('   *Begin InitZero$')
               call ppiclf_solve_InitZero
            call ppiclf_prints('   *End InitZero$')
            call ppiclf_prints('   *Begin AddParticles$')
               call ppiclf_solve_AddParticles(npart,y,rprop)
            call ppiclf_prints('   *End AddParticles$')

            ! TLJ - 11/23/2024
            ! The write files at t=0 has been moved to a single
            !   location in ppiclf_solve_WriteVTK
            !call ppiclf_prints('   *Begin WriteParticleVTU$')
            !   call ppiclf_io_WriteParticleVTU('')
            !call ppiclf_prints('    End WriteParticleVTU$')
            !call ppiclf_prints('   *Begin WriteBinVTU$')
            !   call ppiclf_io_WriteBinVTU('')
            !call ppiclf_prints('    End WriteBinVTU$')
         ENDif

      call ppiclf_prints(' End InitParticle$')
!
      call mpi_barrier(ppiclf_comm,ierr)
!

      ! This prints out initial bin information
      call ppiclf_io_OutputDiagGen

      PPICLF_LINIT = .true.

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitParam(imethod,ndim,iendian)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4  imethod
      integer*4  ndim
      integer*4  iendian
!
      if (imethod .eq. 0 .or. imethod .ge. 3 .or. imethod .le. -2)
     >   call ppiclf_exittr('Invalid integration method$',0.0d0,imethod)
      if (ndim .le. 1 .or. ndim .ge. 4)
     >   call ppiclf_exittr('Invalid problem dimension$',0.0d0,ndim)
      if (iendian .lt. 0 .or. iendian .gt. 1)
     >   call ppiclf_exittr('Invalid Endian$',0.0d0,iendian)

      ppiclf_imethod      = imethod
      ppiclf_ndim         = ndim
      ppiclf_iendian      = iendian

      ppiclf_filter       = -1   ! filt default for no projection

      ppiclf_iperiodic(1) = 1    ! periodic in x (== 0) 
      ppiclf_iperiodic(2) = 1    ! periodic in y (==0)
      ppiclf_iperiodic(3) = 1    ! periodic in z (==0)

      ppiclf_cycle  = 0
      ppiclf_iostep = 1
      ppiclf_dt     = 0.0d0
      ppiclf_time   = 0.0d0

      ppiclf_overlap    = .false.
      ppiclf_linit      = .false.
      ppiclf_lfilt      = .false.
      ppiclf_lfiltgauss = .false.
      ppiclf_lfiltbox   = .false.
      ppiclf_lintp      = .false.
      ppiclf_lproj      = .false.
      ppiclf_lsubbin    = .false.
      ppiclf_lsubsubbin = .false.
      if (PPICLF_INTERP .eq. 1)  ppiclf_lintp = .true.
      if (PPICLF_PROJECT .eq. 1) ppiclf_lproj = .true.

      ppiclf_xdrange(1,1) = -1E20
      ppiclf_xdrange(2,1) =  1E20
      ppiclf_xdrange(1,2) = -1E20
      ppiclf_xdrange(2,2) =  1E20
      ppiclf_xdrange(1,3) = -1E20
      ppiclf_xdrange(2,3) =  1E20

      ppiclf_d2chk(1) = 0.0d0
      ppiclf_d2chk(2) = 0.0d0
      ppiclf_d2chk(3) = 0.0d0

      ppiclf_n_bins(1) = 1
      ppiclf_n_bins(2) = 1
      ppiclf_n_bins(3) = 1

      ppiclf_nwall    = 0
      ppiclf_iwallm   = 0

      PPICLF_INT_ICNT = 0

      RETURN
      END
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitNeighborBin(rwidth)
     > bind(C, name="ppiclc_solve_InitNeighborBin")
#else
      subroutine ppiclf_solve_InitNeighborBin(rwidth)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      real*8 rwidth
!
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitNeighborBin$',0.0d0
     >                  ,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitNeighborBin$'
     >                  ,0.0d0,0)

      ppiclf_lsubsubbin = .true.

      ppiclf_d2chk(3) = rwidth

      ! TLJ added 12/21/2024
      if (ppiclf_nid==0) then
         print*,'TLJ checking d2chk(3) = ',ppiclf_d2chk(3)
      ENDif

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_SetNeighborBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
!
      do i=1,ppiclf_npart
         ppiclf_nb_r(1,i) = floor((ppiclf_cp_map(1,i)-ppiclf_binb(1))/
     >                             ppiclf_d2chk(3))
         ppiclf_nb_r(2,i) = floor((ppiclf_cp_map(2,i)-ppiclf_binb(3))/
     >                             ppiclf_d2chk(3))
         ppiclf_nb_r(3,i) = 0
         if (ppiclf_ndim .eq. 3)
     >   ppiclf_nb_r(3,i) = floor((ppiclf_cp_map(3,i)-ppiclf_binb(5))/
     >                             ppiclf_d2chk(3))
      ENDdo

      do i=1,ppiclf_npart_gp
         ppiclf_nb_g(1,i) = floor((ppiclf_rprop_gp(1,i)-ppiclf_binb(1))/
     >                             ppiclf_d2chk(3))
         ppiclf_nb_g(2,i) = floor((ppiclf_rprop_gp(2,i)-ppiclf_binb(3))/
     >                             ppiclf_d2chk(3))
         ppiclf_nb_g(3,i) = 0
         if (ppiclf_ndim .eq. 3)
     >   ppiclf_nb_g(3,i) = floor((ppiclf_rprop_gp(3,i)-ppiclf_binb(5))/
     >                             ppiclf_d2chk(3))
      ENDdo

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitZero
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i,j,ic,k,ie
!
      ic = 0
      do i=1,PPICLF_LPART
      do j=1,PPICLF_LRS
         ic = ic + 1
         ppiclf_y    (j,i) = 0.0d0
         ppiclf_ydot (j,i) = 0.0d0
         ppiclf_ydotc(j,i) = 0.0d0
         ppiclf_y1   (ic ) = 0.0d0
      ENDdo
      do j=1,PPICLF_LRP
         ppiclf_rprop(j,i) = 0.0d0
      ENDdo
      do j=1,PPICLF_LRP2
         ppiclf_rprop2(j,i) = 0.0d0
      ENDdo
      do j=1,PPICLF_LRP3
         ppiclf_rprop3(j,i) = 0.0d0
      ENDdo
      do j=1,PPICLF_LRP4
         ppiclf_rprop4(j,i) = 0.0d0
      ENDdo
      do j=1,PPICLF_LRP5
         ppiclf_force(j,i) = 0.0d0
      ENDdo
      do j=1,PPICLF_LIP
         ppiclf_iprop(j,i) = 0
      enddo
      enddo
      do i=1,PPICLF_LPART_GP
      do j=1,PPICLF_LRP_GP
         ppiclf_rprop_gp(j,i) = 0.0d0
      enddo
      enddo

      ppiclf_npart = 0

      do ie=1,PPICLF_LEE
      do ic=1,PPICLF_LRP_INT
      do k=1,PPICLF_LEZ
      do j=1,PPICLF_LEY
      do i=1,PPICLF_LEX
        ppiclf_int_fld(i,j,k,ic,ie) = 0.0d0
      ENDdo
      ENDdo
      ENDdo
      ENDdo
      ENDdo

      !!call ppiclf_user_InitZero

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_NearestNeighbor(i)
!
      implicit none
!
      include "PPICLF"
      include "omp_lib.h"
! 
! Input:
! 
      integer*4 i
! 
! Internal: 
! 
      real*8 ydum(PPICLF_LRS), rpropdum(PPICLF_LRP)
      real*8 A(3),B(3),C(3),AB(3),AC(3), dist2, xdist2, ydist2,
     >       dist_total
      integer*4 i_iim, i_iip, i_jjm, i_jjp, i_kkm, i_kkp, j, j_ii, j_jj,
     >          j_kk, jp
      real*8 rnx, rny, rnz, area, rpx1, rpy1, rpz1, rpx2, rpy2, rpz2,
     >       rflip, a_sum, rd, rdist, theta, tri_area, rthresh,
     >       ab_dot_ac, ab_mag, ac_mag, zdist2
      integer*4 istride, k, kmax, kp, kkp, kk
! 
      i_iim = ppiclf_nb_r(1,i) - 1
      i_iip = ppiclf_nb_r(1,i) + 1
      i_jjm = ppiclf_nb_r(2,i) - 1
      i_jjp = ppiclf_nb_r(2,i) + 1
      i_kkm = ppiclf_nb_r(3,i) - 1
      i_kkp = ppiclf_nb_r(3,i) + 1

      dist2 = ppiclf_d2chk(3)**2

!$omp parallel do private(i, j)
      do j=1,ppiclf_npart
         if (j .eq. i) cycle

         j_ii = ppiclf_nb_r(1,j)
         j_jj = ppiclf_nb_r(2,j)
         j_kk = ppiclf_nb_r(3,j)

         if (j_ii .gt. i_iip .or. j_ii .lt. i_iim) cycle
         if (j_jj .gt. i_jjp .or. j_jj .lt. i_jjm) cycle
         if (ppiclf_ndim .eq. 3) then
         if (j_kk .gt. i_kkp .or. j_kk .lt. i_kkm) cycle
         ENDif

         xdist2 = (ppiclf_cp_map(1,i)-ppiclf_cp_map(1,j))**2
         if (xdist2 .gt. dist2) cycle
         ydist2 = (ppiclf_cp_map(2,i)-ppiclf_cp_map(2,j))**2
         if (ydist2 .gt. dist2) cycle
         dist_total = xdist2 + ydist2
         if (ppiclf_ndim .eq. 3) then
         zdist2 = (ppiclf_cp_map(3,i)-ppiclf_cp_map(3,j))**2
         if (zdist2 .gt. dist2) cycle
         dist_total = dist_total+zdist2
         ENDif
         if (dist_total .gt. dist2) cycle

         call ppiclf_user_EvalNearestNeighbor(i,j,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ppiclf_cp_map(1,j)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,j))

      ENDdo
!$omp end parallel do

!$omp parallel do private(j)
      do j=1,ppiclf_npart_gp
         j_ii = ppiclf_nb_g(1,j)
         j_jj = ppiclf_nb_g(2,j)
         j_kk = ppiclf_nb_g(3,j)

         if (j_ii .gt. i_iip .or. j_ii .lt. i_iim) cycle
         if (j_jj .gt. i_jjp .or. j_jj .lt. i_jjm) cycle
         if (ppiclf_ndim .eq. 3) then
         if (j_kk .gt. i_kkp .or. j_kk .lt. i_kkm) cycle
         ENDif

         xdist2 = (ppiclf_cp_map(1,i)-ppiclf_rprop_gp(1,j))**2
         if (xdist2 .gt. dist2) cycle
         ydist2 = (ppiclf_cp_map(2,i)-ppiclf_rprop_gp(2,j))**2
         if (ydist2 .gt. dist2) cycle
         dist_total = xdist2 + ydist2
         if (ppiclf_ndim .eq. 3) then
         zdist2 = (ppiclf_cp_map(3,i)-ppiclf_rprop_gp(3,j))**2
         if (zdist2 .gt. dist2) cycle
         dist_total = dist_total+zdist2
         ENDif
         if (dist_total .gt. dist2) cycle

         jp = -1*j
         call ppiclf_user_EvalNearestNeighbor(i,jp,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ppiclf_rprop_gp(1,j)
     >                                 ,ppiclf_rprop_gp(1+PPICLF_LRS,j))

      ENDdo
!$omp end parallel do

      istride = ppiclf_ndim

!$omp parallel do private(j)
      do j=1,ppiclf_nwall

         rnx  = ppiclf_wall_n(1,j)
         rny  = ppiclf_wall_n(2,j)
         rnz  = 0.0d0
         area = ppiclf_wall_n(3,j)
         rpx1 = ppiclf_cp_map(1,i)
         rpy1 = ppiclf_cp_map(2,i)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(1,j)
         rpy2 = ppiclf_wall_c(2,j)
         rpz2 = 0.0d0
         rpx2 = rpx2 - rpx1
         rpy2 = rpy2 - rpy1

         if (ppiclf_ndim .eq. 3) then
            rnz  = ppiclf_wall_n(3,j)
            area = ppiclf_wall_n(4,j)
            rpz1 = ppiclf_cp_map(3,i)
            rpz2 = ppiclf_wall_c(3,j)
            rpz2 = rpz2 - rpz1
         ENDif
    
         rflip = rnx*rpx2 + rny*rpy2 + rnz*rpz2
         if (rflip .gt. 0.0d0) then
            rnx = -1.0d0*rnx
            rny = -1.0d0*rny
            rnz = -1.0d0*rnz
         ENDif


         a_sum = 0.0d0
         kmax = 2
         if (ppiclf_ndim .eq. 3) kmax = 3
         do k=1,kmax 
            kp = k+1
            if (kp .gt. kmax) kp = kp-kmax ! cycle
            
            kk   = istride*(k-1)
            kkp  = istride*(kp-1)
            rpx1 = ppiclf_wall_c(kk+1,j)
            rpy1 = ppiclf_wall_c(kk+2,j)
            rpz1 = 0.0d0
            rpx2 = ppiclf_wall_c(kkp+1,j)
            rpy2 = ppiclf_wall_c(kkp+2,j)
            rpz2 = 0.0d0

            if (ppiclf_ndim .eq. 3) then
               rpz1 = ppiclf_wall_c(kk+3,j)
               rpz2 = ppiclf_wall_c(kkp+3,j)
            ENDif

            rd   = -(rnx*rpx1 + rny*rpy1 + rnz*rpz1)

            rdist = abs(rnx*ppiclf_cp_map(1,i)+rny*ppiclf_cp_map(2,i)
     >                 +rnz*ppiclf_cp_map(3,i)+rd)
            rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

            ! give a little extra room for walls (2x)
            if (rdist .gt. 2.0d0*ppiclf_d2chk(3)) goto 1511

            ydum(1) = ppiclf_cp_map(1,i) - rdist*rnx
            ydum(2) = ppiclf_cp_map(2,i) - rdist*rny
            ydum(3) = 0.0d0

            A(1) = ydum(1)
            A(2) = ydum(2)
            A(3) = 0.0d0

            B(1) = rpx1
            B(2) = rpy1
            B(3) = 0.0d0

            C(1) = rpx2
            C(2) = rpy2
            C(3) = 0.0d0

            AB(1) = B(1) - A(1)
            AB(2) = B(2) - A(2)
            AB(3) = 0.0d0

            AC(1) = C(1) - A(1)
            AC(2) = C(2) - A(2)
            AC(3) = 0.0d0

            if (ppiclf_ndim .eq. 3) then
               ydum(3) = ppiclf_cp_map(3,i) - rdist*rnz
               A(3) = ydum(3)
               B(3) = rpz1
               C(3) = rpz2
               AB(3) = B(3) - A(3)
               AC(3) = C(3) - A(3)

               AB_DOT_AC = AB(1)*AC(1) + AB(2)*AC(2) + AB(3)*AC(3)
               AB_MAG = sqrt(AB(1)**2 + AB(2)**2 + AB(3)**2)
               AC_MAG = sqrt(AC(1)**2 + AC(2)**2 + AC(3)**2)
               theta  = acos(AB_DOT_AC/(AB_MAG*AC_MAG))
               tri_area = 0.5d0*AB_MAG*AC_MAG*sin(theta)
            elseif (ppiclf_ndim .eq. 2) then
               AB_MAG = sqrt(AB(1)**2 + AB(2)**2)
               tri_area = AB_MAG
            ENDif
            a_sum = a_sum + tri_area
         ENDdo

         rthresh = 1.10d0 ! keep it from slipping through crack on edges
         if (a_sum .gt. rthresh*area) cycle

         jp = 0
         call ppiclf_user_EvalNearestNeighbor(i,jp,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ydum
     >                                 ,rpropdum)

 1511 continue
      ENDdo
!$omp end parallel do

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_NearestNeighborSB(i,SBt,SBc,SBm,SBn,iB)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Input:
!
      INTEGER*4 i, SBt, SBn(3), iB(3)  
      INTEGER*4 SBc(0:(SBt-1)),
     >  SBm(0:(SBt-1),(ppiclf_npart+ppiclf_npart_gp))
! 
! Internal: 
! 
      REAL*8 ydum(PPICLF_LRS), rpropdum(PPICLF_LRP), xp(3), bin_xMin(3)
      REAL*8 A(3),B(3),C(3),AB(3),AC(3), dist2, xdist2, ydist2,
     >       dist_total
      INTEGER*4 j, jp, l, iSB, jSB, kSB, loopSB, tempSB, iSBin(3)
      REAL*8 rnx, rny, rnz, area, rpx1, rpy1, rpz1, rpx2, rpy2, rpz2,
     >       rflip, a_sum, rd, rdist, theta, tri_area, rthresh,
     >       ab_dot_ac, ab_mag, ac_mag, zdist2
      INTEGER*4 istride, k, kmax, kp, kkp, kk
! 
      dist2 = ppiclf_d2chk(3)**2
      
      ! find ith particle subbin (tempSB)
      DO l = 1,3
        IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
          xp(l) = ppiclf_y(l,i)
          bin_xMin(l) = ppiclf_binb(2*l-1)+iB(l)*ppiclf_bins_dx(l)
        ELSE
          xp(l) = 0.0
        END IF
      END DO
      DO l = 1,3
        IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
          iSBin(l) = FLOOR((xp(l) - (bin_xMin(l)
     >         - ppiclf_d2chk(3)))/ppiclf_d2chk(3))
        ELSE
          iSBin(l) = 0
        END IF
      END DO
      tempSB = iSBin(1) + iSBin(2)*SBn(1) + iSBin(3)*SBn(1)*SBn(2)
      
      ! Loop through real particles
      DO iSB = 1,3     !to look at -1,current,+1 x-dir subbins
        DO jSB = 1,3   !to look at -1,current,+1 x-dir subbins
          DO kSB = 1,3 !to look at -1,current,+1 x-dir subbins
          ! Loops through 27 adjacent subbins
          loopSB = tempSB + (-2+iSB) + (-2+jSB)*SBn(1) 
     >             + (-2+kSB)*SBn(1)*SBn(2)
          IF (loopSB .GT. -1 .AND. loopSB .LT. SBt) THEN
            DO k = 1,SBc(loopSB) 
              j = SBm(loopSB,k)
              IF (j .GT. 0) THEN ! Real particle
                IF (j .EQ. i) CYCLE
                xdist2 = (ppiclf_cp_map(1,i)-ppiclf_cp_map(1,j))**2
                IF (xdist2 .GT. dist2) CYCLE
                ydist2 = (ppiclf_cp_map(2,i)-ppiclf_cp_map(2,j))**2
                IF (ydist2 .GT. dist2) CYCLE
                dist_total = xdist2 + ydist2
                IF (ppiclf_ndim .EQ. 3) THEN
                  zdist2 = (ppiclf_cp_map(3,i)-ppiclf_cp_map(3,j))**2
                  IF (zdist2 .GT. dist2) CYCLE
                  dist_total = dist_total+zdist2
                END IF
                IF (dist_total .GT. dist2) CYCLE
                CALL ppiclf_user_EvalNearestNeighbor(i,j
     >                                   ,ppiclf_cp_map(1,i)
     >                                   ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                   ,ppiclf_cp_map(1,j)
     >                                   ,ppiclf_cp_map(1+PPICLF_LRS,j))
              ELSE IF (j .LT. 0) THEN ! Ghost Particle
                ! Negative was just use for ghost particle indicator
                ! in subbin mapping array. Need to flip sign
                j = - j                 
                xdist2 = (ppiclf_cp_map(1,i)-ppiclf_rprop_gp(1,j))**2
                IF (xdist2 .GT. dist2) CYCLE
                ydist2 = (ppiclf_cp_map(2,i)-ppiclf_rprop_gp(2,j))**2
                IF (ydist2 .GT. dist2) CYCLE
                dist_total = xdist2 + ydist2
                IF (ppiclf_ndim .EQ. 3) THEN
                zdist2 = (ppiclf_cp_map(3,i)-ppiclf_rprop_gp(3,j))**2
                IF (zdist2 .GT. dist2) CYCLE
                dist_total = dist_total+zdist2
                END IF
                IF (dist_total .GT. dist2) CYCLE
                jp = -1*j
                CALL ppiclf_user_EvalNearestNeighbor(i,jp
     >                                 ,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ppiclf_rprop_gp(1,j)
     >                                 ,ppiclf_rprop_gp(1+PPICLF_LRS,j))
              END IF
            END DO !k
          END IF ! if loopSB is valid
        END DO !kSB
      END DO !jSB
      END DO !iSB


      istride = ppiclf_ndim
      do j=1,ppiclf_nwall

         rnx  = ppiclf_wall_n(1,j)
         rny  = ppiclf_wall_n(2,j)
         rnz  = 0.0d0
         area = ppiclf_wall_n(3,j)
         rpx1 = ppiclf_cp_map(1,i)
         rpy1 = ppiclf_cp_map(2,i)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(1,j)
         rpy2 = ppiclf_wall_c(2,j)
         rpz2 = 0.0d0
         rpx2 = rpx2 - rpx1
         rpy2 = rpy2 - rpy1

         if (ppiclf_ndim .eq. 3) then
            rnz  = ppiclf_wall_n(3,j)
            area = ppiclf_wall_n(4,j)
            rpz1 = ppiclf_cp_map(3,i)
            rpz2 = ppiclf_wall_c(3,j)
            rpz2 = rpz2 - rpz1
         ENDif
    
         rflip = rnx*rpx2 + rny*rpy2 + rnz*rpz2
         if (rflip .gt. 0.0d0) then
            rnx = -1.0d0*rnx
            rny = -1.0d0*rny
            rnz = -1.0d0*rnz
         ENDif


         a_sum = 0.0d0
         kmax = 2
         if (ppiclf_ndim .eq. 3) kmax = 3
         do k=1,kmax 
            kp = k+1
            if (kp .gt. kmax) kp = kp-kmax ! cycle
            
            kk   = istride*(k-1)
            kkp  = istride*(kp-1)
            rpx1 = ppiclf_wall_c(kk+1,j)
            rpy1 = ppiclf_wall_c(kk+2,j)
            rpz1 = 0.0d0
            rpx2 = ppiclf_wall_c(kkp+1,j)
            rpy2 = ppiclf_wall_c(kkp+2,j)
            rpz2 = 0.0d0

            if (ppiclf_ndim .eq. 3) then
               rpz1 = ppiclf_wall_c(kk+3,j)
               rpz2 = ppiclf_wall_c(kkp+3,j)
            ENDif

            rd   = -(rnx*rpx1 + rny*rpy1 + rnz*rpz1)

            rdist = abs(rnx*ppiclf_cp_map(1,i)+rny*ppiclf_cp_map(2,i)
     >                 +rnz*ppiclf_cp_map(3,i)+rd)
            rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

            ! give a little extra room for walls (2x)
            if (rdist .gt. 2.0d0*ppiclf_d2chk(3)) goto 1519

            ydum(1) = ppiclf_cp_map(1,i) - rdist*rnx
            ydum(2) = ppiclf_cp_map(2,i) - rdist*rny
            ydum(3) = 0.0d0

            A(1) = ydum(1)
            A(2) = ydum(2)
            A(3) = 0.0d0

            B(1) = rpx1
            B(2) = rpy1
            B(3) = 0.0d0

            C(1) = rpx2
            C(2) = rpy2
            C(3) = 0.0d0

            AB(1) = B(1) - A(1)
            AB(2) = B(2) - A(2)
            AB(3) = 0.0d0

            AC(1) = C(1) - A(1)
            AC(2) = C(2) - A(2)
            AC(3) = 0.0d0

            if (ppiclf_ndim .eq. 3) then
               ydum(3) = ppiclf_cp_map(3,i) - rdist*rnz
               A(3) = ydum(3)
               B(3) = rpz1
               C(3) = rpz2
               AB(3) = B(3) - A(3)
               AC(3) = C(3) - A(3)

               AB_DOT_AC = AB(1)*AC(1) + AB(2)*AC(2) + AB(3)*AC(3)
               AB_MAG = sqrt(AB(1)**2 + AB(2)**2 + AB(3)**2)
               AC_MAG = sqrt(AC(1)**2 + AC(2)**2 + AC(3)**2)
               theta  = acos(AB_DOT_AC/(AB_MAG*AC_MAG))
               tri_area = 0.5d0*AB_MAG*AC_MAG*sin(theta)
            elseif (ppiclf_ndim .eq. 2) then
               AB_MAG = sqrt(AB(1)**2 + AB(2)**2)
               tri_area = AB_MAG
            ENDif
            a_sum = a_sum + tri_area
         ENDdo

         rthresh = 1.10d0 ! keep it from slipping through crack on edges
         if (a_sum .gt. rthresh*area) cycle

         jp = 0
         call ppiclf_user_EvalNearestNeighbor(i,jp,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ydum
     >                                 ,rpropdum)

 1519 continue
      ENDdo

      RETURN
      END
!-----------------------------------------------------------------------
       subroutine ppiclf_solve_InitWall(xp1,xp2,xp3)
!
      implicit none
!
      include "PPICLF"
! 
! Input:
! 
      real*8 xp1(*)
      real*8 xp2(*)
      real*8 xp3(*)
!
! Internal:
!
      real*8 rpx1, rpy1, rpz1, rpx2, rpy2, rpz2,
     >       a_sum, theta, tri_area, 
     >       ab_dot_ac, ab_mag, ac_mag, rise, run, rmag, 
     >       rpx3, rpy3, rpz3
      integer*4 istride, k, kmax, kp, kkp, kk
      real*8 A(3),B(3),C(3),AB(3),AC(3)
!
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitWall$',0.d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitWall$'
     >                  ,0.d0,0)

      ppiclf_nwall = ppiclf_nwall + 1 

      if (ppiclf_nwall .gt. PPICLF_LWALL)
     >call ppiclf_exittr('Increase LWALL in user file$'
     >                  ,0.d0,ppiclf_nwall)

      istride = ppiclf_ndim
      a_sum = 0.0d0
      kmax = 2
      if (ppiclf_ndim .eq. 3) kmax = 3

      if (ppiclf_ndim .eq. 3) then
         ppiclf_wall_c(1,ppiclf_nwall) = xp1(1)
         ppiclf_wall_c(2,ppiclf_nwall) = xp1(2)
         ppiclf_wall_c(3,ppiclf_nwall) = xp1(3)
         ppiclf_wall_c(4,ppiclf_nwall) = xp2(1)
         ppiclf_wall_c(5,ppiclf_nwall) = xp2(2)
         ppiclf_wall_c(6,ppiclf_nwall) = xp2(3)
         ppiclf_wall_c(7,ppiclf_nwall) = xp3(1)
         ppiclf_wall_c(8,ppiclf_nwall) = xp3(2)
         ppiclf_wall_c(9,ppiclf_nwall) = xp3(3)

         A(1) = (xp1(1) + xp2(1) + xp3(1))/3.0d0
         A(2) = (xp1(2) + xp2(2) + xp3(2))/3.0d0
         A(3) = (xp1(3) + xp2(3) + xp3(3))/3.0d0
      elseif (ppiclf_ndim .eq. 2) then
         ppiclf_wall_c(1,ppiclf_nwall) = xp1(1)
         ppiclf_wall_c(2,ppiclf_nwall) = xp1(2)
         ppiclf_wall_c(3,ppiclf_nwall) = xp2(1)
         ppiclf_wall_c(4,ppiclf_nwall) = xp2(2)

         A(1) = (xp1(1) + xp2(1))/2.0d0
         A(2) = (xp1(2) + xp2(2))/2.0d0
         A(3) = 0.0d0
      ENDif

      ! compute area:
      do k=1,kmax 
         kp = k+1
         if (kp .gt. kmax) kp = kp-kmax ! cycle
         
         kk   = istride*(k-1)
         kkp  = istride*(kp-1)
         rpx1 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy1 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(kkp+1,ppiclf_nwall)
         rpy2 = ppiclf_wall_c(kkp+2,ppiclf_nwall)
         rpz2 = 0.0d0

         B(1) = rpx1
         B(2) = rpy1
         B(3) = 0.0d0
        
         C(1) = rpx2
         C(2) = rpy2
         C(3) = 0.0d0
        
         AB(1) = B(1) - A(1)
         AB(2) = B(2) - A(2)
         AB(3) = 0.0d0
        
         AC(1) = C(1) - A(1)
         AC(2) = C(2) - A(2)
         AC(3) = 0.0d0

         if (ppiclf_ndim .eq. 3) then
             rpz1 = ppiclf_wall_c(kk+3,ppiclf_nwall)
             rpz2 = ppiclf_wall_c(kkp+3,ppiclf_nwall)
             B(3) = rpz1
             C(3) = rpz2
             AB(3) = B(3) - A(3)
             AC(3) = C(3) - A(3)
        
             AB_DOT_AC = AB(1)*AC(1) + AB(2)*AC(2) + AB(3)*AC(3)
             AB_MAG = sqrt(AB(1)**2 + AB(2)**2 + AB(3)**2)
             AC_MAG = sqrt(AC(1)**2 + AC(2)**2 + AC(3)**2)
             theta  = acos(AB_DOT_AC/(AB_MAG*AC_MAG))
             tri_area = 0.5d0*AB_MAG*AC_MAG*sin(theta)
         elseif (ppiclf_ndim .eq. 2) then
             AB_MAG = sqrt(AB(1)**2 + AB(2)**2)
             tri_area = AB_MAG
         ENDif
         a_sum = a_sum + tri_area
      ENDdo
      
      ppiclf_wall_n(ppiclf_ndim+1,ppiclf_nwall) = a_sum

      ! wall normal:
      if (ppiclf_ndim .eq. 2) then

         rise = xp2(2) - xp1(2)
         run  = xp2(1) - xp1(1)

         rmag = sqrt(rise**2 + run**2)
         rise = rise/rmag
         run  = run/rmag
         
         ppiclf_wall_n(1,ppiclf_nwall) = rise
         ppiclf_wall_n(2,ppiclf_nwall) = -run

      elseif (ppiclf_ndim .eq. 3) then

         k  = 1
         kk = istride*(k-1)
         rpx1 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy1 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz1 = ppiclf_wall_c(kk+3,ppiclf_nwall)
         
         k  = 2
         kk = istride*(k-1)
         rpx2 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy2 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz2 = ppiclf_wall_c(kk+3,ppiclf_nwall)
         
         k  = 3
         kk = istride*(k-1)
         rpx3 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy3 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz3 = ppiclf_wall_c(kk+3,ppiclf_nwall)
    
         A(1) = rpx2 - rpx1
         A(2) = rpy2 - rpy1
         A(3) = rpz2 - rpz1

         B(1) = rpx3 - rpx2
         B(2) = rpy3 - rpy2
         B(3) = rpz3 - rpz2

         ppiclf_wall_n(1,ppiclf_nwall) = A(2)*B(3) - A(3)*B(2)
         ppiclf_wall_n(2,ppiclf_nwall) = A(3)*B(1) - A(1)*B(3)
         ppiclf_wall_n(3,ppiclf_nwall) = A(1)*B(2) - A(2)*B(1)

         rmag = sqrt(ppiclf_wall_n(1,ppiclf_nwall)**2 +
     >               ppiclf_wall_n(2,ppiclf_nwall)**2 +
     >               ppiclf_wall_n(3,ppiclf_nwall)**2)

         ppiclf_wall_n(1,ppiclf_nwall) = ppiclf_wall_n(1,ppiclf_nwall)
     >                                  /rmag
         ppiclf_wall_n(2,ppiclf_nwall) = ppiclf_wall_n(2,ppiclf_nwall)
     >                                  /rmag
         ppiclf_wall_n(3,ppiclf_nwall) = ppiclf_wall_n(3,ppiclf_nwall)
     >                                  /rmag

      ENDif

      RETURN
      END
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitPeriodicX(xl,xr)
     > bind(C, name="ppiclc_solve_InitPeriodicX")
#else
      subroutine ppiclf_solve_InitPeriodicX(xl,xr)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8 xl
      real*8 xr
! 
      if (xl .ge. xr)
     >call ppiclf_exittr('PeriodicX must have xl < xr$',xl,0)

      ppiclf_iperiodic(1) = 0

      ppiclf_xdrange(1,1) = xl
      ppiclf_xdrange(2,1) = xr

      call ppiclf_solve_InitSolve

      RETURN
      END
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitPeriodicY(yl,yr)
     > bind(C, name="ppiclc_solve_InitPeriodicY")
#else
      subroutine ppiclf_solve_InitPeriodicY(yl,yr)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8 yl
      real*8 yr
! 
      if (yl .ge. yr)
     >call ppiclf_exittr('PeriodicY must have yl < yr$',yl,0)

      ppiclf_iperiodic(2) = 0

      ppiclf_xdrange(1,2) = yl
      ppiclf_xdrange(2,2) = yr

      call ppiclf_solve_InitSolve

      RETURN
      END
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitPeriodicZ(zl,zr)
     > bind(C, name="ppiclc_solve_InitPeriodicZ")
#else
      subroutine ppiclf_solve_InitPeriodicZ(zl,zr)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8 zl
      real*8 zr
! 
      if (zl .ge. zr)
     >call ppiclf_exittr('PeriodicZ must have zl < zr$',zl,0)
      if (ppiclf_ndim .lt. 3)
     >call ppiclf_exittr('Cannot do PeriodicZ if not 3D$',zl,0)

      ppiclf_iperiodic(3) = 0

      ppiclf_xdrange(1,3) = zl
      ppiclf_xdrange(2,3) = zr

      call ppiclf_solve_InitSolve

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitAngularPeriodic(flag,
     >              rin, rout, angle, xangle)
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      ! Thierry -  07/24/24 - modified code begings here
      ! global variables - user input file
      integer*4 flag
      real*8 rin, rout, angle, xangle
      ! local variables
      real*8 pi, angled

        ! Thierry - 07/24/24 - modified code begins here
        ! Implementation of Angular Periodicity
        ! Just like how Rocflu does it in modflu/RFLU_ModRelatedPatches.F90
        ! this is invoked when particle is leaving x-axis or y-axis
       
        ! sign convention for theta is +ve when measured CCW
        ! switch angle sign when particle is leaving from upper face
        if (rin .ge. rout)
     >   call ppiclf_exittr('Angular Per must have rin < rout$',rout,0)

            ppiclf_iperiodic(1) = 0 ! X-periodic
            ppiclf_iperiodic(2) = 0 ! Y-periodic

            SELECT CASE (ang_case)
              CASE (1) ! general wedge ; 0 <= angle < 90
                if (ppiclf_nid.eq.0) print*,"General Wedge Initialized!"
                ppiclf_xdrange(1,1) = rin  ! Min X-periodic face
                ppiclf_xdrange(2,1) = rout ! Max X-periodic face
                ppiclf_xdrange(1,2) = tan(xangle)*rout ! Min Y-periodic face
                ppiclf_xdrange(2,2) = tan(angle - abs(xangle))*rout ! Max Y-periodic face

              CASE (2) ! quarter cylinder ; angle = 90
                if (ppiclf_nid.eq.0)
     >             print*,"Quarter Cylinder Initialized!"
                ppiclf_xdrange(1,1) = rin  
                ppiclf_xdrange(2,1) = rout 
                ppiclf_xdrange(1,2) = tan(xangle)*rout
                ppiclf_xdrange(2,2) = rout 
              
              CASE (3) ! half cylinder ; angle = 180
                if (ppiclf_nid.eq.0)
     >             print*,"Half Cylinder Initialized!"
                ppiclf_xdrange(1,1) = -1.0*rout
                ppiclf_xdrange(2,1) = rout 
                ppiclf_xdrange(1,2) = tan(xangle)*rout
                ppiclf_xdrange(2,2) = rout 
              
              CASE DEFAULT
                call ppiclf_exittr('Invalid Rotational Case!$',0.0d0
     >             ,ppiclf_nid)
              END SELECT

      call ppiclf_solve_InitSolve
      
      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitAngularPlane(i,rin, rout,
     >                                         angle, xangle,
     >                                         dist1, dist2)
!
      implicit none
!
      include "PPICLF"
! Inputs 
!
      integer*4 i
      real*8 rin, rout, angle, xangle
! Local Variables:
!     
      real*8 p1(3), p2(3), p3(3), p4(3), p5(3), p6(3),
     >       v1(3), v2(3), v3(3), v4(3), n1(3), n2(3),
     >       A, B, C, D, E, F, G, H, zt, xp, yp, zp
!
! Outputs
      real*8 dist1, dist2
!
      xp = ppiclf_y(PPICLF_JX,i)
      yp = ppiclf_y(PPICLF_JY,i)
      zp = ppiclf_y(PPICLF_JZ,i)
      
      zt = ppiclf_binb(6) - ppiclf_binb(5) ! bin thickness in z-direction

      !!! upper plane calculation !!! 
      
      ! plane equation
      ! Ax + By + Cz + D = 0

      ! p1, p2, p3 are 3 points in the upper plane
      p1 = (/rin, tan(angle - abs(xangle))*rin, 0.0d0/) 
      p2 = (/rout, tan(angle - abs(xangle))*rout, 0.0d0/) 
      p3 = (/rin, tan(angle - abs(xangle))*rin, zt/) 

      v1 = p2 - p1 ! vector P1P2
      v2 = p3 - p1 ! vector P1P3
      
      ! upper plane normal vector - n1(A,B,C) = v1 x v2
      ! cross product calculation
      A =  v1(2)*v2(3) - v1(3)*v2(2)
      B = -v1(1)*v2(3) + v1(3)*v2(1)
      C =  v1(1)*v2(2) - v1(2)*v2(1)
      n1(1)=A ; n1(2)=B; n1(3)=C
      
      ! values of either p1, p2, or p3 can be used to calculate D
      D = -A*p1(1) - B*p1(2) - C*p1(3)
      
      ! P(xp, yp, zp) arbitrary point
      ! dist = distance between P and upper plane 
      dist1 = abs(A*xp + B*yp + C*zp + D)
      dist1 = dist1/sqrt(A**2 + B**2 + C**2)
      
      !!! lower plane calculation !!! 
      ! plane equation
      ! Ex + Fy + Gz + H = 0

      ! p4, p5, p6 are 3 points in the lower plane
      p4 = (/rin, -tan(angle - abs(xangle))*rin, 0.0d0/)
      p5 = (/rout, -tan(angle - abs(xangle))*rout, 0.0d0/)
      p6 = (/rin, -tan(angle - abs(xangle))*rin, zt/)
      
      v3 = p5 - p4 ! vector P4P5
      v4 = p6 - p4 ! vector P4P6
      
      ! lower plane normal vector - n2(E,F,G) = v3 x v4
      ! cross product calculation
      E =  v3(2)*v4(3) - v3(3)*v4(2)
      F = -v3(1)*v4(3) + v3(3)*v4(1)
      G =  v3(1)*v4(2) - v3(2)*v4(1)
      n2(1)=E ; n2(2)=F; n2(3)=G
      
      ! values of either p4, p5, or p6 can be used to calculate H
      H = -E*p4(1) - F*p4(2) - G*p4(3)

      ! P(xp, yp, zp) arbitrary point
      ! dist = distance between P and lower plane 
      dist2 = abs(E*xp + F*yp + G*zp + H)
      dist2 = dist2/sqrt(E**2 + F**2 + G**2)

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InvokeLinearPeriodic(i)
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, j, jj(3), jchk, in_part(PPICLF_LPART)
!
      jj(1) = 1 ; jj(2) = 2 ; jj(3) = 3

! Case 1 - Linear Periodicity in any of 3 directions ; NO Anuglar Periodicity
      if(((x_per_flag.eq.1).or.(y_per_flag.eq.1).or.(z_per_flag.eq.1))
     >   .and.(ang_per_flag.eq.0)) then

      do j=0,ppiclf_ndim-1
        jchk = jj(j+1)
        ! particle leaving min. periodic face -> move it relative to 
        !                                         max periodic face
        if(ppiclf_y(jchk,i).lt.ppiclf_xdrange(1,j+1)) then
          ppiclf_y(jchk,i) = ppiclf_xdrange(2,j+1) - 
     >                  abs(ppiclf_xdrange(1,j+1) - ppiclf_y(jchk,i))
          goto 1512
        ENDif

        ! particle leaving max. periodic face -> move it relative to 
        !                                         min periodic face
        if(ppiclf_y(jchk,i).gt.ppiclf_xdrange(2,j+1)) then
          ppiclf_y(jchk,i) = ppiclf_xdrange(1,j+1) + 
     >                  abs(ppiclf_y(jchk,i) - ppiclf_xdrange(2,j+1))
          goto 1512
        ENDif
        
        ! Thierry - I'm not sure what this does but this is how it was implemented
        if (ppiclf_iprop(1,i) .eq. 2) then
             in_part(i) = -1 ! only if periodic check fails it will get here
        ENDif
 1512 continue
        END do ! j=0, ndim-1
      ENDif ! Case 1 

! Case 2 - Linear Periodicity in Z-direction only; WITH Anuglar Periodicity
      if((z_per_flag.eq.1).and.(ang_per_flag.eq.1)) then

        ! particle leaving min. or max. x-face -> delete it
        if((ppiclf_y(1,i).lt.ppiclf_xdrange(1,1)).or.
     >  (ppiclf_y(1,i).gt.ppiclf_xdrange(2,1))) then
          call ppiclf_solve_MarkForRemoval(i)
          goto 1515
        ENDif

        ! particle leaving min. z-periodic face -> move it relative to 
        !                                         max z-periodic face
        if(ppiclf_y(3,i).lt.ppiclf_xdrange(1,3)) then
          ppiclf_y(3,i) = ppiclf_xdrange(2,3) - 
     >                  abs(ppiclf_xdrange(1,3) - ppiclf_y(3,i))
          goto 1515
        ENDif

        ! particle leaving max. z-periodic face -> move it relative to 
        !                                         min z-periodic face
        if(ppiclf_y(3,i).gt.ppiclf_xdrange(2,3)) then
          ppiclf_y(3,i) = ppiclf_xdrange(1,3) + 
     >                  abs(ppiclf_y(3,i) - ppiclf_xdrange(2,3))
          goto 1515
        ENDif
        
        if (ppiclf_iprop(1,i) .eq. 2) then
             in_part(i) = -1 ! only if periodic check fails it will get here
        ENDif
 1515 continue
      ENDif ! Case 2

      RETURN
      END
!-----------------------------------------------------------------------
       subroutine ppiclf_solve_InvokeAngularPeriodic(i,flag,
     >                                              per_alpha,
     >                                              angle, xangle,
     >                                              register)
!
      implicit none
!
      include "PPICLF"
! :
! Input: 
! 
      ! Thierry -  07/24/24 - modified code begins here
      ! global variables
      integer*4 i, flag
      real*8 rin, rout, per_alpha, angle, xangle
      ! local variables
      real*8 ct, st, ex, ey, ez, local_angle
      real*8 rotmat(3,3) , v(3), x(3)
      integer*4 register

        ! Thierry - 07/24/24 - modified code begins here
        ! Implementation of Rotational Periodicity
        ! Just like how Rocflu does it in modflu/RFLU_ModRelatedPatches.F90
        ! this is invoked when particle is leaving x-axis or y-axis
       
!        print*, "!!! Rotational Periodicity Invoked !!!!" 
          
        ! use local angle so the value of angle does not get affected globally
        local_angle = angle
        ! Thierry 
        !    (1) sign convention for theta is +ve when measured CCW
        !           switch angle sign when particle is leaving from 
        !           upper periodic face
        !    (2) 0.5 instead of 1.0 to switch angle for ghost algorithm
        !           since the ghost is being created before the 
        !           particle is leaving domain
        if(per_alpha.gt. 0.5*(xangle+angle)) 
     >    local_angle=-1.0*local_angle
        
        ! Half-cylinder case - particle leaving +ve x-axis 
        !                    - adjust rotation matrix angle accordingly
        if(ang_case .eq. 3) then
          if(per_alpha .lt. xangle) local_angle = 0.0 
        END if

        ! convert from degrees to radians
        ct = cos(local_angle)
        st = sin(local_angle)
        
        SELECT CASE(flag)
          !CASE(1)
          !  ex = 1.0d0
          !  ey = 0.0d0
          !  ez = 0.0d0
          !  print*, "X-Rotational Axis"

          !CASE(2)
          !  ex = 0.0d0
          !  ey = 1.0d0
          !  ez = 0.0d0
          !  print*, "Y-Rotational Axis"

          CASE(1)
            ex = 0.0d0
            ey = 0.0d0
            ez = 1.0d0
!            print*, "Z-Rotational Axis"
          CASE DEFAULT
            call ppiclf_exittr('Invalid Axis of Rotation!$',0.0d0
     >         ,ppiclf_nid)

          END SELECT 
          
          ! Rotation Matrix calculation
          rotmat(1,1) = ct + (1.0d0-ct)*ex*ex
          rotmat(1,2) =      (1.0d0-ct)*ex*ey - st*ez
          rotmat(1,3) =      (1.0d0-ct)*ex*ez + st*ey
          
          rotmat(2,1) =      (1.0d0-ct)*ey*ex + st*ez
          rotmat(2,2) = ct + (1.0d0-ct)*ey*ey
          rotmat(2,3) =      (1.0d0-ct)*ey*ez - st*ex
          
          rotmat(3,1) =      (1.0d0-ct)*ez*ex - st*ey
          rotmat(3,2) =      (1.0d0-ct)*ez*ey + st*ex
          rotmat(3,3) = ct + (1.0d0-ct)*ez*ez

          ! Corrdinates modification
          x(1) = ppiclf_y(PPICLF_JX,i)
          x(2) = ppiclf_y(PPICLF_JY,i)
          x(3) = ppiclf_y(PPICLF_JZ,i)
          
          xrot = MATMUL(rotmat, x)
          
          ! Velocity vector modification

          v(1) = ppiclf_y(PPICLF_JVX,i)
          v(2) = ppiclf_y(PPICLF_JVY,i)
          v(3) = ppiclf_y(PPICLF_JVZ,i)

          vrot = MATMUL(rotmat, v)
          
          ! 08/27/24 - Thierry - we add a register variable to 
          !   choose if we want to register
          !   the angularly modified variables 
          ! register = 1 when called from RemoveParticle -> we want to modify the values
          ! register = 0 when called from AngularCreateGhost -> we don't want to modify values
          
          ! register modified values
          if (register==1) then
            !print*, "Registering values!" 
            ppiclf_y(PPICLF_JX,i) = xrot(1)
            ppiclf_y(PPICLF_JY,i) = xrot(2)
            ppiclf_y(PPICLF_JZ,i) = xrot(3)
            
            ppiclf_y(PPICLF_JVX,i) = vrot(1)
            ppiclf_y(PPICLF_JVY,i) = vrot(2)
            ppiclf_y(PPICLF_JVZ,i) = vrot(3)
            
          END if 
       
      RETURN
      END
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitGaussianFilter(filt,alpha,iwallm)
     > bind(C, name="ppiclc_solve_InitGaussianFilter")
#else
      subroutine ppiclf_solve_InitGaussianFilter(filt,alpha,iwallm)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8    filt
      real*8    alpha
      integer*4 iwallm
! 
! Internal: 
! 
      real*8 rsig
! 
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitFilter$',0.0d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitFilter$',0.0d0
     >                  ,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.0d0
     >                  ,0)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter can only be called once$',0.0d0,0)
      if (iwallm .lt. 0 .or. iwallm .gt. 1)
     >call ppiclf_exittr('0 or 1 must be used to specify filter mirror$'
     >                  ,0.0d0,iwallm)

      ppiclf_filter = filt
      ppiclf_alpha  = alpha 
      ppiclf_iwallm = iwallm

      rsig             = ppiclf_filter/(2.0d0*sqrt(2.0d0*log(2.0d0)))
      ppiclf_d2chk(2)  = rsig*sqrt(-2*log(ppiclf_alpha))

      ! TLJ added 12/21/2024
      if (ppiclf_nid==0) then
         print*,'TLJ recompute d2chk(2) based on Gausian filter = ',
     >     ppiclf_d2chk(2)
      ENDif

      PPICLF_LSUBBIN = .true.
      if (ppiclf_ngrids .eq. 0) PPICLF_LSUBBIN = .false.

      PPICLF_LFILT      = .true.
      PPICLF_LFILTGAUSS = .true.

      ppiclf_ngrids = 0 ! for now leave sub bin off

      RETURN
      END
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitBoxFilter(filt,iwallm,sngl_elem)
     > bind(C, name="ppiclc_solve_InitBoxFilter")
#else
      subroutine ppiclf_solve_InitBoxFilter(filt,iwallm,sngl_elem)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8    filt
      integer*4 iwallm
      integer*4 sngl_elem
! 
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitFilter$',0.0d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitFilter$',0.0d0
     >                   ,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.0d0
     >                   ,0)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter can only be called once$',0.0d0,0)

c     filt = sqrt(1.5d0*filt**2/log(2.0d0) + 1.0d0)

      ppiclf_filter = filt
      ppiclf_iwallm = iwallm

      ppiclf_d2chk(2)  = filt/2.0d0

      ! TLJ added 12/21/2024
      if (ppiclf_nid==0) then
         print*,'TLJ checking d2chk(2) = ',ppiclf_d2chk(2)
      ENDif

      PPICLF_LSUBBIN = .true.
      if (ppiclf_ngrids .eq. 0) PPICLF_LSUBBIN = .false.

      PPICLF_LFILT    = .true.
      PPICLF_LFILTBOX = .true.

      ! option to only use the current element (filter width will be 
      ! ignored)
      ! Note that this assumes the element volume is that of
      ! a cuboid... will need to get a better way for general
      ! hexahedral element eventually
      if ( sngl_elem == 1 ) PPICLF_SNGL_ELEM = .true.

      ppiclf_ngrids = 0 ! for now leave sub bin off

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_SetParticleTag(npart)
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      integer*4 npart
! 
! Internal: 
! 
      integer*4 i
!
      do i=ppiclf_npart-npart+1,ppiclf_npart
         ppiclf_iprop(5,i) = ppiclf_nid 
         ppiclf_iprop(6,i) = ppiclf_cycle
         ppiclf_iprop(7,i) = i
      ENDdo

      RETURN
      END
c----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_WriteVTU(time)
     > bind(C, name="ppiclc_solve_WriteVTU")
#else
      subroutine ppiclf_solve_WriteVTU(time)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8    time
! 
! Internal:
!
      ppiclf_time   = time

      ! already wrote initial conditions
      !if (ppiclf_cycle .ne. 0) then
            call ppiclf_io_WriteParticleVTU('')
            call ppiclf_io_WriteBinVTU('')
      !ENDif

      if (ppiclf_lsubbin)
     >      call ppiclf_io_WriteSubBinVTU('')

      ! Output diagnostics
      call ppiclf_io_OutputDiagAll

      RETURN
      END
c----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_IntegrateParticle(istep,iostep,dt,time)
     > bind(C, name="ppiclc_solve_IntegrateParticle")
#else
      subroutine ppiclf_solve_IntegrateParticle(istep,iostep,dt,time)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      integer*4 istep
      integer*4 iostep
      real*8    dt
      real*8    time
! 
! Internal:
!
      logical iout
!
      ppiclf_cycle  = istep
      ppiclf_iostep = iostep
      ppiclf_dt     = dt
      ppiclf_time   = time

      ! integerate in time
      if (ppiclf_imethod .eq. 1) 
     >   call ppiclf_solve_IntegrateRK3(iout)
      if (ppiclf_imethod .eq. -1) 
     >   call ppiclf_solve_IntegrateRK3s(iout)
      if (ppiclf_imethod .eq. 2)
     >   call ppiclf_solve_IntegrateRK3s_Rocflu(iout)

!      ! output files
!      if (ppiclf_iostep .gt.0)then
!      if (mod(ppiclf_cycle,ppiclf_iostep) .eq. 0 .and. iout) then
!
!         ! already wrote initial conditions
!         if (ppiclf_cycle .ne. 0) then
!            call ppiclf_io_WriteParticleVTU('')
!            call ppiclf_io_WriteBinVTU('')
!         endif
!
!         if (ppiclf_lsubbin)
!     >      call ppiclf_io_WriteSubBinVTU('')
!      endif
!
!      ! Output diagnostics
!      if (mod(ppiclf_cycle,ppiclf_iostep) .eq. 0 .and. iout) then
!         call ppiclf_io_OutputDiagAll
!      endif
!      endif

      RETURN
      END
c----------------------------------------------------------------------
      subroutine ppiclf_solve_IntegrateRK3(iout)
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, ndum, nstage, istage
!
! Output:
!
      logical iout
!
      ! save stage 1 solution
      ndum = PPICLF_NPART*PPICLF_LRS
      do i=1,ndum
         ppiclf_y1(i) = ppiclf_y(i,1)
      ENDdo

      ! get rk3 coeffs
      call ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      do istage=1,nstage

         ! evaluate ydot
         call ppiclf_solve_SetYdot

         ! rk3 integrate
         do i=1,ndum
            ndum = PPICLF_NPART*PPICLF_LRS
            ppiclf_y(i,1) =  ppiclf_rk3coef(1,istage)*ppiclf_y1   (i)
     >                     + ppiclf_rk3coef(2,istage)*ppiclf_y    (i,1)
     >                     + ppiclf_rk3coef(3,istage)*ppiclf_ydot (i,1)
         ENDdo
      ENDdo

      iout = .true.

      RETURN
      END
c----------------------------------------------------------------------
      subroutine ppiclf_solve_IntegrateRK3s(iout)
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, ndum, nstage, istage
      integer*4 icalld
      save      icalld
      data      icalld /0/
!
! Output:
!
      logical iout
!
      icalld = icalld + 1


      ! get rk3 coeffs
      call ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      istage = mod(icalld,nstage)
      if (istage .eq. 0) istage = 3
      iout = .false.
      if (istage .eq. nstage) iout = .true.

      ! save stage 1 solution
      if (istage .eq. 1) then
      ndum = PPICLF_NPART*PPICLF_LRS
      do i=1,ndum
         ppiclf_y1(i) = ppiclf_y(i,1)
      ENDdo
      ENDif

      ! evaluate ydot
      call ppiclf_solve_SetYdot

      ! rk3 integrate
      ndum = PPICLF_NPART*PPICLF_LRS
      do i=1,ndum
         ppiclf_y(i,1) =  ppiclf_rk3coef(1,istage)*ppiclf_y1   (i)
     >                  + ppiclf_rk3coef(2,istage)*ppiclf_y    (i,1)
     >                  + ppiclf_rk3coef(3,istage)*ppiclf_ydot (i,1)
      ENDdo

      RETURN
      END
c----------------------------------------------------------------------
      subroutine ppiclf_solve_IntegrateRK3s_Rocflu(iout)
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, ndum, nstage, istage
      integer*4 icalld
      integer*4 j
      save      icalld
      data      icalld /0/

!
! Output:
!
      logical iout
!
      icalld = icalld + 1

      ! get rk3 coeffs
      call ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      istage = mod(icalld,nstage)
      if (istage .eq. 0) istage = 3
      iout = .false.
      if (istage .eq. nstage) iout = .true.

      ! evaluate ydot
      call ppiclf_solve_SetYdot

      !Zero out for first stage
      ! TLJ this loop is fine
      if (istage .eq. 1) then
        ppiclf_y1 = 0.0d0
        !ndum = PPICLF_NPART*PPICLF_LRS
        !do i=1,ndum
        !  ppiclf_y1(i) = 0.0d0 
        !enddo
      ENDif

      ! TLJ comment Dec 7, 2023
      ! The Rocflu RK3 can be found in equation (7) of:
      ! S. Yu. "Runge-Kutta Methods Combined with Compact
      !   Difference Schemes for the Unsteady Euler Equations".
      !   Center for Modeling of Turbulence and Transition.
      !   Research Briefs, 1991.

      ! TLJ modified loop to prevent -fcheck=all error
      !     must preserve order
      !     ppiclf_y   (PPICLF_LRS, PPICLF_LPART)
      !     ppiclf_ydot(PPICLF_LRS, PPICLF_LPART)
      ndum = 0
      do j=1,PPICLF_NPART
      do i=1,PPICLF_LRS
         ndum = ndum + 1
         ppiclf_y(i,j) =  -ppiclf_rk3coef(1,istage)*ppiclf_y1   (ndum)
     >                   + ppiclf_rk3coef(2,istage)*ppiclf_y    (i,j)
     >                   + ppiclf_rk3coef(3,istage)*ppiclf_ydot (i,j)
      ENDdo
      ENDdo
      !ndum = PPICLF_NPART*PPICLF_LRS
      !do i=1,ndum
      !   ppiclf_y(i,1) =  -ppiclf_rk3coef(1,istage)*ppiclf_y1   (i)
      !>                   + ppiclf_rk3coef(2,istage)*ppiclf_y    (i,1)
      !>                   + ppiclf_rk3coef(3,istage)*ppiclf_ydot (i,1)
      !enddo

!Store Current stage RHS for next stage's use
      ! TLJ modified loop to prevent -fcheck=all error
      !     must preserve order
      ndum = 0
      do j=1,PPICLF_NPART
      do i=1,PPICLF_LRS
         ndum = ndum + 1
         ppiclf_y1(ndum) =  ppiclf_ydot(i,j)
      ENDdo
      ENDdo
      !ndum = PPICLF_NPART*PPICLF_LRS
      !do i=1,ndum
      !   ppiclf_y1(i) = ppiclf_ydot(i,1)
      !enddo

!WARNING: Experimental fix to keep particles unsure where to place this
!          command. Either before or after the storing of the current 
!          storage
        call ppiclf_solve_RemoveParticle      
!End Experimental fix
      RETURN
      END
c----------------------------------------------------------------------
      subroutine ppiclf_solve_SetYdot
!
      implicit none
!
      include "PPICLF"
! 
      call ppiclf_solve_InitSolve
      call ppiclf_user_SetYdot
      call ppiclf_solve_RemoveParticle

      RETURN
      END
c----------------------------------------------------------------------
      subroutine ppiclf_solve_InitSolve
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, j
!
      call ppiclf_comm_CreateBin
      call ppiclf_comm_FindParticle
      call ppiclf_comm_MoveParticle
      if (ppiclf_overlap)
     >   call ppiclf_comm_MapOverlapMesh
      if ((ppiclf_lintp .and. ppiclf_int_icnt .ne. 0) .or.
     >    (ppiclf_lproj .and. ppiclf_sngl_elem))
     >   call ppiclf_solve_InterpParticleGrid
      call ppiclf_solve_RemoveParticle
      if (ppiclf_lsubsubbin .or. ppiclf_lproj) then
         ! Thierry - standard ghost algorithm when no angular periodicity
         if(ang_per_flag.eq.0) then
           call ppiclf_comm_CreateGhost
         elseif(ang_per_flag.eq.1) then
           call ppiclf_comm_AngularCreateGhost
         ENDif
         call ppiclf_comm_MoveGhost
      ENDif

      if (ppiclf_lproj .and. ppiclf_overlap) 
     >   call ppiclf_solve_ProjectParticleGrid
      if (ppiclf_lsubsubbin) 
     >   call ppiclf_solve_SetNeighborBin
      ! Zero 
      do i=1,PPICLF_LPART
      do j=1,PPICLF_LRS
         ppiclf_ydotc(j,i) = 0.0d0
      ENDdo
      ENDdo

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InterpParticleGrid
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 j
!
      call ppiclf_solve_InitInterp
      do j=1,PPICLF_INT_ICNT
         call ppiclf_solve_InterpField(j)
      ENDdo
      call ppiclf_solve_FinalizeInterp

      call ppiclf_solve_LocalInterp

      call ppiclf_solve_PostInterp

      PPICLF_INT_ICNT = 0


      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InterpFieldUser(jp,infld)
!
      implicit none
!
      include "PPICLF"
!
! Input: 
!
      integer*4 jp
      real*8 infld(*)
!
! Internal:
!
      integer*4 n
!
      if (PPICLF_INTERP .eq. 0)
     >call ppiclf_exittr(
     >     'No specified interpolated fields, set PPICLF_LRP_INT$',0.0d0
     >                   ,0)

      PPICLF_INT_ICNT = PPICLF_INT_ICNT + 1

      if (PPICLF_INT_ICNT .gt. PPICLF_LRP_INT)
     >   call ppiclf_exittr('Interpolating too many fields$'
     >                     ,0.0d0,PPICLF_INT_ICNT)
      if (jp .le. 0 .or. jp .gt. PPICLF_LRP)
     >   call ppiclf_exittr('Invalid particle array interp. location$'
     >                     ,0.0d0,jp)

      ! set up interpolation map
      PPICLF_INT_MAP(PPICLF_INT_ICNT) = jp

      ! copy to infld internal storage
      n = PPICLF_NEE*PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      call ppiclf_copy(ppiclf_int_fldu(1,1,1,1,PPICLF_INT_ICNT)
     >                ,infld(1),n)

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitInterp
!
      implicit none
!
      include "PPICLF"
      include 'omp_lib.h'
! 
! Internal: 
! 
      integer*4 n, ie, npt_max, np, ndum
      real*8 tol, bb_t
!
      if (.not.ppiclf_overlap)
     >call ppiclf_exittr('Cannot interpolate unless overlap grid$',0.0d0
     >                   ,0)
      if (.not.ppiclf_lintp) 
     >call ppiclf_exittr('To interpolate, set PPICLF_LRP_PRO to ~= 0$'
     >                   ,0.0d0,0)

      n = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ

!$omp parallel do private(ie)
      do ie=1,ppiclf_neltb
         call ppiclf_copy(ppiclf_xm1bi(1,1,1,ie,1)
     >                   ,ppiclf_xm1b(1,1,1,1,ie),n)
         call ppiclf_copy(ppiclf_xm1bi(1,1,1,ie,2)
     >                   ,ppiclf_xm1b(1,1,1,2,ie),n)
         call ppiclf_copy(ppiclf_xm1bi(1,1,1,ie,3)
     >                   ,ppiclf_xm1b(1,1,1,3,ie),n)
      ENDdo
!$omp end parallel do      

      tol     = 5e-13
      bb_t    = 0.01
      npt_max = 128
      np      = 1
c     ndum    = ppiclf_neltb*n
      ndum    = ppiclf_neltb+2

      ! initiate findpts since mapping can change on next call
      call pfgslib_findpts_setup(ppiclf_fp_hndl
     >                         ,ppiclf_comm_nid
     >                         ,np ! only 1 rank on this comm
     >                         ,ppiclf_ndim
     >                         ,ppiclf_xm1bi(1,1,1,1,1)
     >                         ,ppiclf_xm1bi(1,1,1,1,2)
     >                         ,ppiclf_xm1bi(1,1,1,1,3)
     >                         ,PPICLF_LEX
     >                         ,PPICLF_LEY
     >                         ,PPICLF_LEZ
     >                         ,ppiclf_neltb
     >                         ,2*PPICLF_LEX
     >                         ,2*PPICLF_LEY
     >                         ,2*PPICLF_LEZ
     >                         ,bb_t
     >                         ,ndum
     >                         ,ndum
     >                         ,npt_max
     >                         ,tol)


      ! copy MapOverlapMesh mapping from prior to communicating map
      ppiclf_neltbbb = ppiclf_neltbb
!$omp parallel do private(ie)
      do ie=1,ppiclf_neltbbb
         call ppiclf_icopy(ppiclf_er_mapc(1,ie),ppiclf_er_maps(1,ie)
     >             ,PPICLF_LRMAX)
      ENDdo
!$omp end parallel do      
      !PRINT*, 'Processor ID, ppiclf_neltbbb', ppiclf_nid, ppiclf_neltbbb
      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InterpField(j)
!
      implicit none
!
      include "PPICLF"
      include 'omp_lib.h'
!
! Input: 
!
      integer*4 jp
!
! Internal:
!
      integer*4 n, ie, iee, j
!
      ! use the map to take original grid and map to fld which will be
      ! sent to mapped processors
      n = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
!$omp parallel do private(ie)
      do ie=1,ppiclf_neltbbb
         iee = ppiclf_er_mapc(1,ie)
         call ppiclf_copy(ppiclf_int_fld (1,1,1,j  ,ie)
     >                   ,ppiclf_int_fldu(1,1,1,iee,j ),n)
      ENDdo
!$omp end parallel do      

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_FinalizeInterp
!
      implicit none
!
      include "PPICLF"
      include "omp_lib.h"
!
! Internal: 
!
      real*8 FLD(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE)
      integer*4 nkey(2), nl, nii, njj, nxyz, nrr, ix, iy, iz, i, jp, ie
      logical partl
!
      ! send it all
      nl   = 0
      nii  = PPICLF_LRMAX
      njj  = 6
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      nrr  = nxyz*PPICLF_LRP_INT
      nkey(1) = 2
      nkey(2) = 1
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl,ppiclf_neltbbb
     >      ,PPICLF_LEE,ppiclf_er_mapc,nii,partl,nl,ppiclf_int_fld
     >      ,nrr,njj)
      call pfgslib_crystal_tuple_sort    (ppiclf_cr_hndl,ppiclf_neltbbb
     >       ,ppiclf_er_mapc,nii,partl,nl,ppiclf_int_fld,nrr,nkey,2)

      ! find which cell particle is in locally
      ix = 1
      iy = 2
      iz = 1
      if (ppiclf_ndim .eq. 3)
     >iz = 3

      call pfgslib_findpts(PPICLF_FP_HNDL           !   call pfgslib_findpts( ihndl,
     >        , ppiclf_iprop (1 ,1),PPICLF_LIP        !   $             rcode,1,
     >        , ppiclf_iprop (3 ,1),PPICLF_LIP        !   &             proc,1,
     >        , ppiclf_iprop (2 ,1),PPICLF_LIP        !   &             elid,1,
     >        , ppiclf_rprop2(1 ,1),PPICLF_LRP2       !   &             rst,ndim,
     >        , ppiclf_rprop2(4 ,1),PPICLF_LRP2       !   &             dist,1,
     >        , ppiclf_y     (ix,1),PPICLF_LRS        !   &             pts(    1),1,
     >        , ppiclf_y     (iy,1),PPICLF_LRS        !   &             pts(  n+1),1,
     >        , ppiclf_y     (iz,1),PPICLF_LRS ,PPICLF_NPART) !   &             pts(2*n+1),1,n)

      do i=1,PPICLF_INT_ICNT
         jp = PPICLF_INT_MAP(i)

c!$omp parallel do private(ie)
         do ie=1,ppiclf_neltbbb
            call ppiclf_copy(fld(1,1,1,ie)
     >                      ,ppiclf_int_fld(1,1,1,i,ie),nxyz)
         ENDdo
c!$omp end parallel do         

         ! sam commenting out eval nearest neighbor to use Local Interp instead
         ! leaving findpts call to help with projection, where the element id is
         ! needed. 
         ! interpolate field locally
!         call pfgslib_findpts_eval_local( PPICLF_FP_HNDL
!     >                                  ,ppiclf_rprop (jp,1)
!     >                                  ,PPICLF_LRP
!     >                                  ,ppiclf_iprop (2,1)
!     >                                  ,PPICLF_LIP
!     >                                  ,ppiclf_rprop2(1,1)
!     >                                  ,PPICLF_LRP2
!     >                                  ,PPICLF_NPART
!     >                                  ,fld)

      ENDdo

      ! free since mapping can change on next call
      call pfgslib_findpts_free(PPICLF_FP_HNDL)

      ! Set interpolated fields to zero again
      ! Sam - commenting out for local routine
      !PPICLF_INT_ICNT = 0

      RETURN
      END
!
!-----------------------------------------------------------------------
!
!     Avery's latest version March 27, 2025
!
      SUBROUTINE ppiclf_solve_LocalInterp
      IMPLICIT NONE

      include "PPICLF"
      include "omp_lib.h"

      ! Local Variables
      INTEGER*4 i, j, k, l, ix, iy, iz, ip, ie, iee, nxyz, nnearest, 
     >          inearest(28)
      REAL*8    d2l, d2i, wsum, eps, A(27,4), d2(28), xp(3),  
     >          center(3,28), b(27,1), w(27),binlength(3),  
     >          d2i_EleLen(3), MaxPoint(3), MinPoint(3),d2Max_EleLen(3),
     >          centeri(3,ppiclf_neltbbb), CellLengthMultiplier
      LOGICAL   added, farAway, LinearPerShift(3)
      REAL*8, ALLOCATABLE :: test(:)
      REAL*8    nearest_cell
      !***************************************************************
      IF(ppiclf_neltbbb .EQ. 0 . AND. ppiclf_npart .GT. 0) PRINT*,
     >  'No elements. Num Particles/Proc ID: ', ppiclf_npart, ppiclf_nid
      eps = 1.0e-20 !machine epsilon
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ !number of points in mesh
                                               !element 
      ! Calculate centroid, max cell lengths
c!$omp parallel do private(ie)
      DO ie = 1,ppiclf_neltbbb !Loop fluid cells on this processor
        ! Initialize as zero for each element
        DO l = 1,3
          centeri(l,ie)    =  0.0D0 
          MaxPoint(l)      = -1.0D10 
          MinPoint(l)      =  1.0D10 
          d2i_EleLen(l)    =  0.0D0  
          d2Max_EleLen(l) = 0.0D0  
        ENDDO !l
        ! Add all x,y,z mesh points for centroid and find extremes
        DO l = 1,3
          DO k = 1,PPICLF_LEZ
            DO j = 1,PPICLF_LEY
              DO i = 1,PPICLF_LEX
                centeri(l,ie) = centeri(l,ie) + ppiclf_xm1b(i,j,k,l,ie)
                IF (ppiclf_xm1b(i,j,k,l,ie) .GT. MaxPoint(l)) 
     >            MaxPoint(l) = ppiclf_xm1b(i,j,k,l,ie)  
                IF (ppiclf_xm1b(i,j,k,l,ie) .LT. MinPoint(l))
     >            MinPoint(l) = ppiclf_xm1b(i,j,k,l,ie)
              ENDDO !i
            ENDDO !j
          ENDDO !k
          ! Find individual (mesh element length)^2 in all directions
          d2i_EleLen(l) = (MaxPoint(l)-MinPoint(l))**2
          ! Find max (mesh element length)^2 for all mesh elements in
          ! all directions
          IF (d2i_EleLen(l) .GT. d2Max_EleLen(l)) 
     >      d2Max_EleLen(l) = d2i_EleLen(l)
          ! Divide by number of points in mesh element to find centroid
          centeri(l,ie) = centeri(l,ie) / nxyz
        ENDDO !l
      ENDDO !ie
c!$omp end parallel do      

      ! Find bin lengths for linear periodicity calculations
      DO l = 1,3
        binlength(l) = ppiclf_binb(2*l) - ppiclf_binb((2*l)-1)
        !Zero shift out for non-periodic cases
        LinearPerShift(l) = .FALSE.
      END DO
      ! LinearPerShift created to enable do loop index
      IF(x_per_flag.EQ.1) LinearPerShift(1) = .TRUE.
      IF(y_per_flag.EQ.1) LinearPerShift(2) = .TRUE.
      IF(z_per_flag.EQ.1) LinearPerShift(3) = .TRUE.

      ! Set the multiple of maximum element length in each
      ! dimension that will be used to search for neighboring
      ! fluid elements for interpolation at particle location.
      !
      ! Standard is 1.5*Len, which will give you at least 27 elements
      ! if particle is not near fluid domain boundary.
      CellLengthMultiplier = 1.5D0*1.5D0
      DO l = 1,3 
        d2Max_EleLen(l) = d2Max_EleLen(l)*CellLengthMultiplier
      END DO

c!$omp parallel do private(ip)
      DO ip=1,ppiclf_npart !Loop all particles in this bin
        ! particle centers in all directions
        xp(1) = ppiclf_y(PPICLF_JX, ip)
        xp(2) = ppiclf_y(PPICLF_JY, ip)
        xp(3) = ppiclf_y(PPICLF_JZ, ip)
        nnearest = 0 ! number of nearest elements
        DO ie = 1,28
          inearest(ie) = -1 ! index of nearest elements
          d2(ie) = 1E20 ! distance to center of nearest element
        ENDDO !ie
        DO ie = 1,ppiclf_neltbbb
          ! get distance from particle to center
          d2l     = 0.0
          d2i     = 0.0
          farAway = .FALSE.
          DO l=1,3
            IF(LinearPerShift(l)) THEN
              d2l = MIN((centeri(l,ie) - xp(l))**2, 
     >                (binlength(l)-ABS(centeri(l,ie) - xp(l)))**2)
            ELSE
              d2l = (centeri(l,ie) - xp(l))**2
            END IF
            d2i = d2i + d2l
            IF (d2l > d2Max_EleLen(l)) farAway = .TRUE.
          END DO !l

          ! skip to next fluid cell if greater than 1.5*max cell
          ! distance in respective x,y,z direction.
          IF (farAWAY) CYCLE !ie
          ! Sort closest fluid cell centers
          added = .FALSE.
          DO i=1,27
            j = 27 - i + 1
            IF (d2i .LT. d2(j)) THEN
              d2(j+1) = d2(j)
              inearest(j+1) = inearest(j)
              DO l=1,3
                center(l, j+1) = center(l, j)
              ENDDO
              d2(j) = d2i
              inearest(j) = ie
              DO l=1,3
                center(l, j) = centeri(l,ie)
              END DO
              added = .TRUE.
            ELSE ! If not within closest cell list
              EXIT !i
            END IF
          END DO !i
          IF (added) nnearest = nnearest + 1
        END DO ! ie
        nnearest = min(nnearest, 27)
        IF (nnearest .lt. 1) THEN
          PRINT *, 'Added, T/F:', added
          PRINT *, 'nnearest, neltbbb, num proc, num part, xp(1:3)',
     >              nnearest, ppiclf_neltbbb,
     >              ppiclf_nid, ppiclf_npart, xp
          PRINT *, 'ppiclf_rprop(1:PPICLF_LRP'
          PRINT *, ppiclf_rprop(1:PPICLF_LRP, ip)
          PRINT *, 'ppiclf_y(1:PPICLF_LRS'
          PRINT *, ppiclf_y(1:PPICLF_LRS, ip)
          CALL ppiclf_exittr('Failed to interpolate',0.0d0,nnearest)
        ELSE
            DO i=1,nnearest
              DO j=1,3
                A(i, j) = xp(j) - center(j, i)
              ENDDO !j
              A(i, 4) = 1
            ENDDO !i
            DO i=1,PPICLF_INT_ICNT
              DO k=1,nnearest
                b(k, 1) = 0.0 ! cell averaged properties
                DO iz=1,PPICLF_LEZ
                  DO iy=1,PPICLF_LEY
                    DO ix=1,PPICLF_LEX
                      b(k, 1) = b(k, 1) + ppiclf_int_fld(ix,iy,iz,i,
     >                                                  inearest(k))
                    ENDDO !ix
                  ENDDO !iy
                ENDDO !iz
                b(k, 1) = b(k, 1) / nxyz
              ENDDO ! nnearest
              j = PPICLF_INT_MAP(i)
              ! Inverse Distance Interpolation
              ppiclf_rprop(j, ip) = 0
              wsum = 0
              DO k=1,nnearest
                ! 19/11/2025 - Avery, Thierry, Bala
                ! The inverse interpolation leads to some jump when
                ! particle is moving from cell to another. temporary fix
                ! is to have it inverse distance to power 4. 
                w(k) = 1.0d0 / (SQRT(d2(k)) + eps)**4
                ppiclf_rprop(j, ip) = ppiclf_rprop(j, ip) + w(k)*b(k, 1)
                wsum = wsum + w(k)
              ENDDO ! k
              ppiclf_rprop(j, ip) = ppiclf_rprop(j, ip) / wsum
              IF (isnan(ppiclf_rprop(j,ip))) THEN
                PRINT *, ip,ppiclf_nid,xp,nnearest
              ENDIF
            ENDDO ! i
        ENDIF ! nnearest
      ENDDO ! ip
c!$omp end parallel do      
      RETURN
      END

!
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_LocalInterp_old
      implicit none

      include "PPICLF"

      ! internal
      integer*4 i,j,k,l,ix,iy,iz
      integer*4 ip, ie, iee, inearest(28), nnearest, nxyz, neltgg
      real*8 A(27, 4), d2i, d2(28), xp(3), center(3, 28)
      real*8 U(27, 27), SIG(4), Vt(4, 4), b(27, 1)
      real*8 interp(4, 1) ! for SVD
      integer*4 m, n, lda, ldu, ldvt, lwork, info, ierr
      real*8 w(27), wsum, eps
      real*8 work(5*27+4)
      character jobu, jobv
      logical added 
      integer*4 nl, nii, njj, nkey(2), nrr
      logical partl
      ! Avery added
      real*8 centeri(3,ppiclf_neltbbb), MaxEleSize, EleSizei,
     >       MaxPoint(3), MinPoint(3)
      !integer subbin_part, subbin_ele(piclf_neltbbb,1) 

      eps = 1.0e-12
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      
      ! Avery Add.  Calculate centroid, max cell size, cell subbin
      MaxEleSize = 0.0
      do ie=1,ppiclf_neltbbb
        ! calculate centroid and largest cell diagonal length
        EleSizei = 0.0  
        do l=1,3
          centeri(l,ie) =  0.0
          MaxPoint(l)   = -1.0E10 
          MinPoint(l)   =  1.0E10
        ENDdo !l

        do l=1,3
        do k=1,PPICLF_LEZ
        do j=1,PPICLF_LEY
        do i=1,PPICLF_LEX
          centeri(l,ie) = centeri(l,ie) + ppiclf_xm1b(i, j, k, l, ie)
          if (ppiclf_xm1b(i,j,k,l,ie) .gt. MaxPoint(l)) 
     >        MaxPoint(l) = ppiclf_xm1b(i,j,k,l,ie)  
          if (ppiclf_xm1b(i,j,k,l,ie) .lt. MinPoint(l))
     >        MinPoint(l) = ppiclf_xm1b(i,j,k,l,ie)
        ENDdo !i
        ENDdo !j
        ENDdo !k
        ENDdo !l
        
        do l=1,3
          EleSizei = EleSizei + (MaxPoint(l)-MinPoint(l))**2
        ENDdo

        if (EleSizei .gt. MaxEleSize) MaxEleSize = EleSizei

        do l=1,3
          centeri(l,ie) = centeri(l,ie) / nxyz
        ENDdo !l

      ENDdo !ie
      ! End Avery Added
      
      ! neighbor search O(Nparticles * Nelements)
      ! this should be switched to a KD tree for Log(Nelements) scaling
  
      ! Avery - I'm not sure that KD tree is best for dynamic point
      ! clouds.  Maybe we have multiple "sub-bins" on each processor.
      ! Seems that this is called an Octree in CS speak.
      ! If we do implement some type of tree, we should probably do the
      ! same thing for particle-particle nearest neighbor search as
      ! well.  The current nearest neighbor does the same cycle method
      ! as below.
  
      do ip=1,ppiclf_npart
        ! particle center
        xp(1) = ppiclf_y(PPICLF_JX, ip)
        xp(2) = ppiclf_y(PPICLF_JY, ip)
        xp(3) = ppiclf_y(PPICLF_JZ, ip)
        
        !start particle in subbin
        !subbin_part = !subbin
        !end particle subbin

        nnearest = 0 ! number of nearby elements

        do ie=1,28
          inearest(ie) = -1 ! index of nearest elements
          d2(ie) = 1E20 ! distance to center of nearest element
        ENDdo

        do ie=1,ppiclf_neltbbb
                   ! get distance from particle to center
          d2i = 0
          do l=1,3
            d2i = d2i + (centeri(l,ie) - xp(l))**2
          ENDdo
          ! Avery Added if / cycle
          ! Go to next cell if particle is 1.5*largest grid cell diagonal
          ! direction away from neighboring cell centroid
          if (d2i .gt. 2.25d0*MaxEleSize) cycle !1.5**2 = 2.25

          ! sort
          added = .false.
          do i=1,27
            j = 27 - i + 1

            if (d2i .lt. d2(j)) then
              d2(j+1) = d2(j)
              inearest(j+1) = inearest(j)
              do l=1,3
                center(l, j+1) = center(l, j)
              ENDdo

              d2(j) = d2i
              inearest(j) = ie
              do l=1,3
                center(l, j) = centeri(l,ie)
              ENDdo

              added = .true.
            else ! Avery added else/exit
              exit
            ENDif
          ENDdo !i

          if (added) nnearest = nnearest + 1
          
        ENDdo ! ie
        ! Avery added if
        ! Check for at least 16 due to cases when one layer of 9 cells
        ! isn't available because the particle is near the bin boundary.
        !if (nnearest .lt. 16) then
        !    print *, '***WARNING***: Less than 16 interpolated cells'
        !endif
       
        nnearest = min(nnearest, 27)

        if (nnearest .lt. 1) then
          print *, 'nnearest', nnearest, ip, ppiclf_npart, xp
          print *, ppiclf_rprop(1:PPICLF_LRP, ip)
          print *, ppiclf_y(1:PPICLF_LRS, ip)
          print *, ppiclf_y(1:PPICLF_LRS, ip)
          call ppiclf_exittr('Failed to interpolate',0.0d0,nnearest)
        else
            do i=1,nnearest
              do j=1,3
                A(i, j) = xp(j) - center(j, i)
              END do
              A(i, 4) = 1
            ENDdo
! Avery print
!           write(ppiclf_nid,*) ip, xp, nnearest, ppiclf_nid
            do i=1,PPICLF_INT_ICNT

              do k=1,nnearest
                b(k, 1) = 0.0 ! cell averaged properties
                ! Avery Add
!                if (i ==1) then
!                  write(ppiclf_nid,*) ip,
!     >                     center(1:3,k),indg(k)
!                ENDif
                do iz=1,PPICLF_LEZ
                do iy=1,PPICLF_LEY
                do ix=1,PPICLF_LEX
                  b(k, 1) = b(k, 1) + ppiclf_int_fld(ix,iy,iz,i,
     >                                inearest(k))
                ENDdo
                ENDdo
                ENDdo

                b(k, 1) = b(k, 1) / nxyz
              ENDdo ! nnearest

              j = PPICLF_INT_MAP(i)

              ! Nearest neighbor interpolation
              !ppiclf_rprop(j, ip) = b(1, 1) ! nearest neighbor interpolation

              ! "harmonic" interpolation
              ppiclf_rprop(j, ip) = 0
              wsum = 0
              do k=1,nnearest
                w(k) = 1.0d0 / (sqrt(d2(k)) + eps)
                ppiclf_rprop(j, ip) = ppiclf_rprop(j, ip) + w(k)*b(k, 1)
                wsum = wsum + w(k)
              ENDdo

              ppiclf_rprop(j, ip) = ppiclf_rprop(j, ip) / wsum
            if (isnan(ppiclf_rprop(j,ip))) then
              PRINT *, ip,ppiclf_nid,xp,nnearest
            ENDif
            ENDdo

        ENDif ! nnearest
      ENDdo ! ip

      ! reset here instead of in finalize
      !PPICLF_INT_ICNT = 0

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_PostInterp
!
      implicit none
!
      include "PPICLF"
!
! Internal: 
!
      integer*4 rstride, istride
      parameter(rstride = 7 + PPICLF_LRP_INT)
      parameter(istride = 3)
      real*8 coord(rstride,PPICLF_LPART)
      integer*4 flag(istride,PPICLF_LPART)
      integer*4 fp_handle, i, j, k, npart
      external ppiclf_iglsum
      integer*4 ppiclf_iglsum
      integer*4 npt_max, np, ndum
      real*8 tol, bb_t
      integer*4 copy_back, jp, nxyz, ie
      real*8 xgrid(PPICLF_LEX, PPICLF_LEY, PPICLF_LEZ, PPICLF_LEE)
     >      ,ygrid(PPICLF_LEX, PPICLF_LEY, PPICLF_LEZ, PPICLF_LEE)
     >      ,zgrid(PPICLF_LEX, PPICLF_LEY, PPICLF_LEZ, PPICLF_LEE)
!
      ! Copy not found particles
      npart = 0
      do i=1,ppiclf_npart
         if (ppiclf_iprop(1,i) .eq. 2) then
            npart = npart + 1
            do j=1,ppiclf_ndim
               coord(j,npart) = ppiclf_y(j,i)
            ENDdo
         ENDif
      ENDdo

      if (ppiclf_iglsum(npart,1) .eq. 0) then
         RETURN
      ENDif

      ! Copy grid indexing 
      ! TLJ changing loop structure to prevent -fcheck=all error
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      do ie=1,ppiclf_nee
      !do i=1,nxyz
      !   xgrid(i,1,1,ie) = ppiclf_xm1bs(i,1,1,1,ie)
      !   ygrid(i,1,1,ie) = ppiclf_xm1bs(i,1,1,2,ie)
      !   zgrid(i,1,1,ie) = ppiclf_xm1bs(i,1,1,3,ie)
      !ENDdo
      do k=1,PPICLF_LEZ
      do j=1,PPICLF_LEY
      do i=1,PPICLF_LEX
         xgrid(i,j,k,ie) = ppiclf_xm1bs(i,j,k,1,ie)
         ygrid(i,j,k,ie) = ppiclf_xm1bs(i,j,k,2,ie)
         zgrid(i,j,k,ie) = ppiclf_xm1bs(i,j,k,3,ie)
      ENDdo
      ENDdo
      ENDdo

      ENDdo

      tol     = 5e-13
      bb_t    = 0.01
      npt_max = 128
      np      = ppiclf_np
c     ndum    = ppiclf_nee*n
      ndum    = ppiclf_nee+2

      ! initiate findpts since mapping can change on next call
      call pfgslib_findpts_setup(fp_handle
     >                         ,ppiclf_comm
     >                         ,np 
     >                         ,ppiclf_ndim
     >                         ,xgrid
     >                         ,ygrid
     >                         ,zgrid
     >                         ,PPICLF_LEX
     >                         ,PPICLF_LEY
     >                         ,PPICLF_LEZ
     >                         ,ppiclf_nee
     >                         ,2*PPICLF_LEX
     >                         ,2*PPICLF_LEY
     >                         ,2*PPICLF_LEZ
     >                         ,bb_t
     >                         ,ndum
     >                         ,ndum
     >                         ,npt_max
     >                         ,tol)

      call pfgslib_findpts(fp_handle           !   call pfgslib_findpts( ihndl,
     >        , flag (1 ,1),istride        !   $             rcode,1,
     >        , flag (3 ,1),istride        !   &             proc,1,
     >        , flag (2 ,1),istride        !   &             elid,1,
     >        , coord(4 ,1),rstride       !   &             rst,ndim,
     >        , coord(7 ,1),rstride       !   &             dist,1,
     >        , coord(1,1) ,rstride        !   &             pts(    1),1,
     >        , coord(2,1) ,rstride        !   &             pts(  n+1),1,
     >        , coord(3,1) ,rstride ,npart) !   &             pts(2*n+1),1,n)

      do i=1,PPICLF_LRP_INT

      ! sam - see finalize interp for note
         ! interpolate field (non-local)
!         call pfgslib_findpts_eval( fp_handle
!     >                                  ,coord (7+i,1)
!     >                                  ,rstride
!     >                                  ,flag (1,1)
!     >                                  ,istride
!     >                                  ,flag (3,1)
!     >                                  ,istride
!     >                                  ,flag (2,1)
!     >                                  ,istride
!     >                                  ,coord(4,1)
!     >                                  ,rstride
!     >                                  ,npart
!     >                                  ,ppiclf_int_fldu(1,1,1,1,i))

      ENDdo

      ! free since mapping can change on next call
      call pfgslib_findpts_free(fp_handle)

      ! Now copy particles back (assumes same ordering)
      k = 0
      do i=1,ppiclf_npart
         copy_back = 0
         if (ppiclf_iprop(1,i) .eq. 2) then
            k = k + 1
            if (flag(1,k) .lt. 2) then
               copy_back = 1
            ENDif
         ENDif

         if (copy_back .eq. 1) then
            ppiclf_iprop(1,i) = flag(1,k)
            do j=1,PPICLF_LRP_INT
               jp = PPICLF_INT_MAP(j)
               if(jp .gt. 0) then
                 ppiclf_rprop(jp,i) = coord(7+j,k)
               endif
            ENDdo
         ENDif
      ENDdo

      RETURN
      END
c----------------------------------------------------------------------
      subroutine ppiclf_solve_SetRK3Coeff(dt)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      real*8 dt
!


      if (ppiclf_imethod .eq. 2) then
        !BD:Rocflu's rk3 scheme

        !Folowing form
        !rk3(1,:) = Temp storage i.e. Previous Stage RHS
        !rk3(2,:) = Temp storage i.e. Current Stage iteration
        !rk3(3,:) = Temp storage i.e. Current Stage RHS
        ppiclf_rk3ark(1) = 8.0d0/15.0d0
        ppiclf_rk3ark(2) = 5.0d0/12.0d0
        ppiclf_rk3ark(3) = 0.75d0

        ! 11/20/2025 - These are the time fraction steps between
        ! the different RK stages. Multiply by ppiclf_dt to get the
        ! dt at the current stage
        ppiclf_rk3dtrk(1) = 8.0d0/15.0d0
        ppiclf_rk3dtrk(2) = 2.0d0/15.0d0
        ppiclf_rk3dtrk(3) = 5.0d0/15.0d0

        ppiclf_rk3coef(1,1) = 0.d00
        ppiclf_rk3coef(2,1) = 1.0d0
        ppiclf_rk3coef(3,1) = dt*8.0d0/15.0d0
        ppiclf_rk3coef(1,2) = dt*17.0d0/60.0d0
        ppiclf_rk3coef(2,2) = 1.0d0
        ppiclf_rk3coef(3,2) = dt*5.0d0/12.0d0
        ppiclf_rk3coef(1,3) = dt*5.0d0/12.0d0
        ppiclf_rk3coef(2,3) = 1.0d0
        ppiclf_rk3coef(3,3) = dt*3.0d0/4.0d0
      else
        !BD:Original Code This follows CMT-nek's rk 3 scheme
        ppiclf_rk3coef(1,1) = 0.d00
        ppiclf_rk3coef(2,1) = 1.0d0 
        ppiclf_rk3coef(3,1) = dt
        ppiclf_rk3coef(1,2) = 3.0d0/4.0d0
        ppiclf_rk3coef(2,2) = 1.0d0/4.0d0 
        ppiclf_rk3coef(3,2) = dt/4.0d0
        ppiclf_rk3coef(1,3) = 1.0d0/3.0d0
        ppiclf_rk3coef(2,3) = 2.0d0/3.0d0
        ppiclf_rk3coef(3,3) = dt*2.0d0/3.0d0
        !BD: Original Code END
      END if

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_MarkForRemoval(i)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 i
!
      ppiclf_iprop(1,i) = 3

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_RemoveParticle
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 in_part(PPICLF_LPART), jj(3), iperiodicx, iperiodicy,
     >          iperiodicz,ndim, i, isl, isr, j, jchk, ic
      ! 08/18/24 - Thierry - added for Angular Periodicity begins here
      real*8 per_alpha
      ! 08/18/24 - Thierry - added for Angular Periodicity ends here
!
      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)
      ndim       = ppiclf_ndim

      jj(1) = 1
      jj(2) = 2
      jj(3) = 3

      do i=1,ppiclf_npart
        
         isl = (i -1) * PPICLF_LRS + 1
         in_part(i) = 0
         if (ppiclf_iprop(1,i) .eq. 3) then
            in_part(i) = -1 ! User removed particle
            goto 1513
         ENDif
!-----------------------------------------------------------------------------------
!!!!!!!!!!!!!!!        Rotational Periodicity Starts Here     !!!!!!!!!!!!!!!!!!!!
            ! currently only supports angular rotation around z-axis
            ! and only check theta component and not radial
            ! applying the radial periodicity is very straightforward, but not needed for now

            if(ang_per_flag .eq. 1) then  ! Angular periodicity
           
            ! particle angle w/ x-axis
            ! per_alpha here is obtained in radians
            ! ang_per_angle & ang_per_xangle are transformed 
            !   to radians in PICL_TEMP_InitFlowSolver.F90
            per_alpha = 
     >         atan2(ppiclf_y(PPICLF_JY,i), ppiclf_y(PPICLF_JX,i))

            ! check if particle leaving through lower face or upper face of wedge
            if ((per_alpha .lt. ang_per_xangle) .or. 
     >          (per_alpha .gt. (ang_per_xangle + ang_per_angle))) then
              call ppiclf_solve_InvokeAngularPeriodic(i, ang_per_flag,
     >                                                per_alpha, 
     >                                                ang_per_angle,
     >                                                ang_per_xangle,
     >                                                1)
            ENDif ! per_alpha
           ENDif ! ang_per_flag

        ! Linear Periodicity Invoked
          if((x_per_flag.eq.1) .or. (y_per_flag.eq.1) 
     >                         .or. (z_per_flag.eq.1)) then
!------------------------------------------------------------------------------------------------------
! 10/16/2024 - Thierry - this is now implemented in
! ppiclf_solve_InvokeLinearPeriodic(i). Either delete below or comment
! out
!         do j=0,ndim-1
!            jchk = jj(j+1)
!            ! Thierry - checks if particle is about to leave min. periodic face
!            ! moves it linearly relative to the max. periodic face
!            if (ppiclf_y(jchk,i).lt.ppiclf_xdrange(1,j+1))then
!               if (((iperiodicx.eq.0) .and. (j.eq.0)) .or.   ! periodic
!     >             ((iperiodicy.eq.0) .and. (j.eq.1)) .or.     
!     >             ((iperiodicz.eq.0) .and. (j.eq.2)) ) then
!                   ppiclf_y(jchk,i) = ppiclf_xdrange(2,j+1) - 
!     &                     abs(ppiclf_xdrange(1,j+1) - ppiclf_y(jchk,i))
!!                   ppiclf_y1(isl+j)   = ppiclf_xdrange(2,j+1) +
!!     &                     abs(ppiclf_xdrange(1,j+1) - ppiclf_y1(isl+j))
!                  goto 1512
!                ENDif
!            ENDif
!            ! Thierry - checks if particle is about to leave max. periodic face
!            ! moves it relative to the min. periodic face
!            if (ppiclf_y(jchk,i).gt.ppiclf_xdrange(2,j+1))then
!               if (((iperiodicx.eq.0) .and. (j.eq.0)) .or.   ! periodic
!     >             ((iperiodicy.eq.0) .and. (j.eq.1)) .or.     
!     >             ((iperiodicz.eq.0) .and. (j.eq.2)) ) then
!                   ppiclf_y(jchk,i) = ppiclf_xdrange(1,j+1) +
!     &                     abs(ppiclf_y(jchk,i) - ppiclf_xdrange(2,j+1))
!!                   ppiclf_y1(isl+j)   = ppiclf_xdrange(1,j+1) +
!!     &                     abs(ppiclf_y1(isl+j) - ppiclf_xdrange(2,j+1))
!                  goto 1512
!                ENDif
!            ENDif
!            if (ppiclf_iprop(1,i) .eq. 2) then
!               in_part(i) = -1 ! only if periodic check fails it will get here
!            ENDif
! 1512 continue
!         ENDdo ! j=0,ndim-1
!------------------------------------------------------------------------------------------------------
         call ppiclf_solve_InvokeLinearPeriodic(i)
         ENDif ! x_per_flag
 1513 continue
      ENDdo ! i=1,ppiclf_part

      ic = 0
      do i=1,ppiclf_npart
         if (in_part(i).eq.0) then
            ic = ic + 1 
            if (i .ne. ic) then
               isl = (i -1) * PPICLF_LRS + 1
               isr = (ic-1) * PPICLF_LRS + 1
               call ppiclf_copy
     >              (ppiclf_y     (1,ic),ppiclf_y(1,i)     ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_y1    (isr) ,ppiclf_y1(isl)    ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_ydot  (1,ic),ppiclf_ydot(1,i)  ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_ydotc (1,ic),ppiclf_ydotc(1,i) ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_rprop (1,ic),ppiclf_rprop(1,i) ,PPICLF_LRP)
               call ppiclf_copy
     >              (ppiclf_rprop2(1,ic),ppiclf_rprop2(1,i),PPICLF_LRP2)
               call ppiclf_copy
     >              (ppiclf_rprop3(1,ic),ppiclf_rprop3(1,i),PPICLF_LRP3)
               call ppiclf_copy
     >              (ppiclf_rprop4(1,ic),ppiclf_rprop4(1,i),PPICLF_LRP4)
               call ppiclf_copy
     >              (ppiclf_force(1,ic),ppiclf_force(1,i),PPICLF_LRP5)
               call ppiclf_icopy
     >              (ppiclf_iprop(1,ic) ,ppiclf_iprop(1,i) ,PPICLF_LIP)
            ENDif
         ENDif
      ENDdo
      ppiclf_npart = ic

      RETURN
      END
!----------------------------------------------------------------------
      subroutine ppiclf_solve_FindWallProject(rx2)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
       real*8 rx2(3)
!
! Internal:
!
      real*8 rnx,rny,rnz,rpx1,rpy1,rpz1,rpx2,rpy2,rpz2,rflip,rd,rdist
      integer*4 j, istride
!
      istride = ppiclf_ndim
      ppiclf_nwall_m = 0
      do j = 1,ppiclf_nwall

         rnx  = ppiclf_wall_n(1,j)
         rny  = ppiclf_wall_n(2,j)
         rnz  = 0.0d0
         rpx1 = rx2(1)
         rpy1 = rx2(2)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(1,j)
         rpy2 = ppiclf_wall_c(2,j)
         rpz2 = 0.0d0
         rpx2 = rpx2 - rpx1
         rpy2 = rpy2 - rpy1

         if (ppiclf_ndim .eq. 3) then
            rnz  = ppiclf_wall_n(3,j)
            rpz1 = rx2(3)
            rpz2 = ppiclf_wall_c(3,j)
            rpz2 = rpz2 - rpz1
         ENDif
    
         rflip = rnx*rpx2 + rny*rpy2 + rnz*rpz2
         if (rflip .gt. 0.0d0) then
            rnx = -1.0d0*rnx
            rny = -1.0d0*rny
            rnz = -1.0d0*rnz
         ENDif

         rpx1 = ppiclf_wall_c(1,j)
         rpy1 = ppiclf_wall_c(2,j)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(istride+1,j)
         rpy2 = ppiclf_wall_c(istride+2,j)
         rpz2 = 0.0d0

         if (ppiclf_ndim .eq. 3) then
            rpz1 = ppiclf_wall_c(3,j)
            rpz2 = ppiclf_wall_c(istride+3,j)
         ENDif

         rd   = -(rnx*rpx1 + rny*rpy1 + rnz*rpz1)

         rdist = abs(rnx*rx2(1)+rny*rx2(2)
     >              +rnz*rx2(3)+rd)
         rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)
         rdist = rdist*2.0d0 ! Mirror

         if (rdist .gt. ppiclf_d2chk(2)) goto 1511

         ppiclf_nwall_m = ppiclf_nwall_m + 1

         ppiclf_xyz_mirror(1,ppiclf_nwall_m) = rx2(1) - rdist*rnx
         ppiclf_xyz_mirror(2,ppiclf_nwall_m) = rx2(2) - rdist*rny
         ppiclf_xyz_mirror(3,ppiclf_nwall_m) = 0.0d0
         if (ppiclf_ndim .eq. 3) then
            ppiclf_xyz_mirror(3,ppiclf_nwall_m) = rx2(3) - rdist*rnz
         ENDif

 1511 continue
      ENDdo

      RETURN
      END
c----------------------------------------------------------------------
      subroutine ppiclf_solve_ProjectParticleGrid
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8    rproj(1+PPICLF_LRP_GP,PPICLF_LPART+PPICLF_LPART_GP)
      integer*4 iproj(4,PPICLF_LPART+PPICLF_LPART_GP)
      integer*4 ppiclf_jxgp,ppiclf_jygp,ppiclf_jzgp
      logical partl, if3d
      integer*4 nkey(2), nxyz, nxyzdum, i, j, k, idum, ic, ip, iip, jjp,
     >          kkp, ilow, ihigh, jlow, jhigh, klow, khigh, ie, jj, j1,
     >          neltbc, ndum, nl, nii, njj, nrr, nlxyzep, iee, ndumdum
      integer*4 jp
      real*8 pi, d2chk2_sq, rdum, multfci, rsig, rdist2, rexp, rx2(3),
     >       rx22, ry22, rz22, rtmp2, evol

      ! Sam - variables for general hex volume calculation
      real*8 v1(3), v2(3), v3(3), cross(3), centroid(3), voltet
      real*8 face(2,2,3) ! ifacex, ifacey, xyz
      integer*4 face_map(3,2,2,2,3) ! xyz,front/back,ifacex,ifacey,ixyz
      integer*4 ix,ix2,iface,inode,ia,ib
!
      if3d = .false.
      if (ppiclf_ndim .eq. 3) if3d = .true.

      PI=4.D0*DATAN(1.D0)

      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ

      nxyzdum = nxyz*PPICLF_LRP_PRO*PPICLF_LEE
      ! TLJ changed loop structure to prevent -fcheck=all error
      ppiclf_pro_fldb = 0.0d0
      !do i=1,nxyzdum
      !   ppiclf_pro_fldb(i,1,1,1,1) = 0.0d0
      !enddo
      !do jp=1,PPICLF_LEE
      !do ip=1,PPICLF_LRP_PRO
      !do k=1,PPICLF_LEZ
      !do j=1,PPICLF_LEY
      !do i=1,PPICLF_LEX
      !   ppiclf_pro_fldb(i,j,k,ip,jp) = 0.0d0
      !enddo
      !enddo
      !enddo
      !enddo
      !enddo

      d2chk2_sq = ppiclf_d2chk(2)**2

      ! real particles
      ppiclf_jxgp  = 1
      ppiclf_jygp  = 2
      ppiclf_jzgp  = 1
      if (if3d) ppiclf_jzgp  = 3

      rdum = 0.0d0
      if (ppiclf_lfiltgauss) then
         rsig    = ppiclf_filter/(2.0d0*sqrt(2.0d0*log(2.0d0)))
         multfci = 1.0d0/(sqrt(2.0d0*pi)**2 * rsig**2) ! in 2D
         if (if3d) multfci = multfci**(1.5d0) ! in 3D
         rdum   = 1.0d0/(-2.0d0*rsig**2)
      ENDif

      if (ppiclf_lfiltbox) then
         if (ppiclf_sngl_elem) then
           multfci = 1.0d0
           rdum = multfci
         else
           multfci = 1.0d0/(PI/4.0d0*ppiclf_filter**2)
           if (if3d) multfci = multfci/(1.0d0/1.5d0*ppiclf_filter)
           rdum = multfci
         ENDif
      ENDif

      ! real particles
      do ip=1,ppiclf_npart

         rproj(1 ,ip) = rdum
         rproj(2 ,ip) = ppiclf_cp_map(ppiclf_jxgp,ip)
         rproj(3 ,ip) = ppiclf_cp_map(ppiclf_jygp,ip)
         if (if3d) 
     >   rproj(4 ,ip) = ppiclf_cp_map(ppiclf_jzgp,ip)

         idum = PPICLF_LRS+PPICLF_LRP
         ic = 4
         ! TLJ modifed loop to remove out of bounds in first index
         !do j=idum+1,idum+PPICLF_LRP_GP
         do j=idum+1,PPICLF_LRP_GP
            ic = ic + 1
            rproj(ic,ip) = ppiclf_cp_map(j,ip)*multfci
         ENDdo

         iproj(1,ip)  = ppiclf_iprop(8,ip)
         iproj(2,ip)  = ppiclf_iprop(9,ip)
         if (if3d)
     >   iproj(3,ip)  = ppiclf_iprop(10,ip)
         iproj(4,ip)  = ppiclf_iprop(11,ip)
      ENDdo

      if (.not. ppiclf_sngl_elem) then
  
        ! ghost particles
        do ip=1,ppiclf_npart_gp
  
           rproj(1 ,ip+ppiclf_npart) = rdum
           rproj(2 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jxgp,ip)
           rproj(3 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jygp,ip)
           if (if3d) 
     >     rproj(4 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jzgp,ip)
  
           idum = PPICLF_LRS+PPICLF_LRP
           ic = 4
           ! TLJ modified loop to remove out of bounds in first index
           !do j=idum+1,idum+PPICLF_LRP_GP
           do j=idum+1,PPICLF_LRP_GP
              ic = ic + 1
              rproj(ic,ip+ppiclf_npart) = ppiclf_rprop_gp(j,ip)*multfci
           ENDdo
                      
           iproj(1,ip+ppiclf_npart)  = ppiclf_iprop_gp(2,ip)
           iproj(2,ip+ppiclf_npart)  = ppiclf_iprop_gp(3,ip)
           if (if3d)
     >     iproj(3,ip+ppiclf_npart)  = ppiclf_iprop_gp(4,ip)
           iproj(4,ip+ppiclf_npart)  = ppiclf_iprop_gp(5,ip)
        ENDdo
  
        ndum = ppiclf_npart+ppiclf_npart_gp
  
        do ip=1,ndum
           iip      = iproj(1,ip)
           jjp      = iproj(2,ip)
           if (if3d)
     >     kkp      = iproj(3,ip)
           ndumdum  = iproj(4,ip)
  
           ilow  = iip-1
           ihigh = iip+1
           jlow  = jjp-1
           jhigh = jjp+1
           if (if3d) then
              klow  = kkp-1
              khigh = kkp+1
           ENDif
  
           ! Find if particle near wall and should mirror itself
           if (ppiclf_iwallm .eq. 1) then
              rx2(1) = rproj(2,ip)
              rx2(2) = rproj(3,ip)
              rx2(3) = rproj(4,ip)
              call ppiclf_solve_FindWallProject(rx2)
           ENDif
  
           do ie=1,ppiclf_neltb
  
                 if (ppiclf_el_map(1,ie) .gt. ndumdum) exit
                 if (ppiclf_el_map(2,ie) .lt. ndumdum) cycle 
           
                 if (ppiclf_el_map(3,ie) .gt. ihigh) cycle
                 if (ppiclf_el_map(4,ie) .lt. ilow)  cycle
                 if (ppiclf_el_map(5,ie) .gt. jhigh) cycle
                 if (ppiclf_el_map(6,ie) .lt. jlow)  cycle
                 if (if3d) then
                    if (ppiclf_el_map(7,ie) .gt. khigh) cycle
                    if (ppiclf_el_map(8,ie) .lt. klow)  cycle
                 ENDif
  
           do k=1,PPICLF_LEZ
           do j=1,PPICLF_LEY
           do i=1,PPICLF_LEX
              if (ppiclf_modgp(i,j,k,ie,4).ne.ndumdum) cycle
  
              rdist2  = (ppiclf_xm1b(i,j,k,1,ie) - rproj(2,ip))**2 +
     >                  (ppiclf_xm1b(i,j,k,2,ie) - rproj(3,ip))**2
              if(if3d) rdist2 = rdist2 +
     >                  (ppiclf_xm1b(i,j,k,3,ie) - rproj(4,ip))**2
  
              if (rdist2 .gt. d2chk2_sq) cycle
  
              rexp = 1.0d0  ! for box filter
              if (ppiclf_lfiltgauss)
     >           rexp = exp(rdist2*rproj(1,ip))
  
              ! add wall effects
              if (ppiclf_iwallm .eq. 1) then
                 do jj=1,ppiclf_nwall_m
                    rx22 = (ppiclf_xm1b(i,j,k,1,ie) 
     >                     -ppiclf_xyz_mirror(1,jj))**2
                    ry22 = (ppiclf_xm1b(i,j,k,2,ie)
     >                     -ppiclf_xyz_mirror(2,jj))**2
                    rtmp2 = rx22 + ry22
                    if (if3d) then
                       rz22 = (ppiclf_xm1b(i,j,k,3,ie)
     >                        -ppiclf_xyz_mirror(3,jj))**2
                       rtmp2 = rtmp2 + rz22
                    ENDif
                    if (ppiclf_lfiltgauss) then
                       rexp = rexp + exp(rtmp2*rproj(1,ip))
                    else
                       rexp = rexp + 1.0d0
                    ENDif
                 ENDdo
              ENDif
  
              
              do jj=1,PPICLF_LRP_PRO
                 j1 = jj+4
                 ppiclf_pro_fldb(i,j,k,jj,ie) = 
     >                           ppiclf_pro_fldb(i,j,k,jj,ie) 
     >                         + rproj(j1,ip)*rexp
              ENDdo
           ENDdo
           ENDdo
           ENDdo
           ENDdo
        ENDdo
      ! sngl elem
      else

!BD: This is where rproj gets stored in fld, this is one of the steps
!that should be tracked

        ! Sam - build map to face coordinates for 3D general hex volume
        ! calculation
        do ix=1,3
          ia = max(3-ix,1)
          ib = min(5-ix,3)
          do iface=1,2
            do j=1,2
            do i=1,2
              face_map(ix,iface,i,j,ix) = iface ! constant
              face_map(ix,iface,i,j,ia) = i ! ix
              face_map(ix,iface,i,j,ib) = j ! iy
            END do
            END do
          END do
        END do

        do ip=1,ppiclf_npart
           !do ie=1,ppiclf_neltb
  
             !if (ie .ne. ppiclf_iprop(2,ip)+1) cycle
             ie = ppiclf_iprop(2, ip) + 1
             if ((ie .lt. 1) .or. (ie .gt. ppiclf_neltb)) cycle

             ! Sam - general hexahedron volume calculation
             if (if3d) then
               ! get centroid of hexahedron
               do ix=1,3
                 centroid(ix) = 0.0
               END do

               do ix=1,3
               do k=1,PPICLF_LEZ
               do j=1,PPICLF_LEY
               do i=1,PPICLF_LEX
                 centroid(ix) = centroid(ix) + ppiclf_xm1b(i,j,k,ix,ie)
               END do
               END do
               END do
               END do

               do ix=1,3
                 centroid(ix) = centroid(ix) / 8.0
               END do


               ! calculate volume based on two contributions from each
               ! face as tetrahedrons
               evol = 0.0
               do ix=1,3
                 do iface=1,2

                   ! get face coordinates
                   do j=1,2
                   do i=1,2
                   do ix2=1,3
                     face(i,j,ix2) = ppiclf_xm1b(
     >                                 face_map(ix,iface,i,j,1),
     >                                 face_map(ix,iface,i,j,2),
     >                                 face_map(ix,iface,i,j,3),
     >                                 ix2,ie)
                   END do
                   END do
                   END do

                   do ix2=1,3
                     v1(ix2) = face(1,2,ix2) - face(2,1,ix2)
                     v2(ix2) = centroid(ix2) - face(2,1,ix2)
                   END do ! ix2

                   ! take cross product
                   cross(1) = v1(2)*v2(3) - v1(3)*v2(2)
                   cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
                   cross(3) = v1(1)*v2(2) - v1(2)*v2(1)

                   ! get contriubtions to volume from each tetrahedron
                   do inode=1,2
                   do ix2=1,3
                     v3(ix2) = face(inode,inode,ix2) - face(2,1,ix2)
                   END do ! ix2

                   ! really 6 times the volume of the tet, but we can
                   ! save an operation by dividing at the end
                   voltet = 0.0
                   do ix2=1,3
                     voltet = voltet + v3(ix2)*cross(ix2)
                   END do ! ix2
                   evol = evol + abs(voltet)
                   END do ! inode
                   
                 END do ! iface
              END do ! ix
               evol = evol / 6.0
             else
               ! Sam - default to naive solution for 2D. ASSUMES
               ! rectangular elements. This will
               ! probably never get used, but if it does throw an error
               ! so the user is absolutely sure of what they're doing.
!               call ppiclf_exittr('Single element projection only
!     >          supported in 3D for general hex elements. Comment and
!     >          ignore this error if your elements are perfect
!     >          rectangles. $',0.0d0,0)

               evol = (ppiclf_xm1b(PPICLF_LEX,1,1,1,ie) 
     >               - ppiclf_xm1b(1,1,1,1,ie))
               evol = evol
     >              * (ppiclf_xm1b(1,PPICLF_LEY,1,2,ie) 
     >               - ppiclf_xm1b(1,1,1,2,ie))
!               if (if3d) evol = evol
!     >              * (ppiclf_xm1b(1,1,PPICLF_LEZ,3,ie) 
!     >               - ppiclf_xm1b(1,1,1,3,ie))
            END if ! if3d

             rexp = 1.0 / evol
           do k=1,PPICLF_LEZ
           do j=1,PPICLF_LEY
           do i=1,PPICLF_LEX
              do jj=1,PPICLF_LRP_PRO
                 j1 = jj+4
                 ppiclf_pro_fldb(i,j,k,jj,ie) = 
     >                           ppiclf_pro_fldb(i,j,k,jj,ie) 
     >                         + rproj(j1,ip)*rexp
              ENDdo
           ENDdo
           ENDdo
           ENDdo
           ENDdo
        !ENDdo ! ppiclf_neltb
      ENDif ! ppiclf_npart

      ! now send xm1b to the processors in nek that hold xm1

      neltbc = ppiclf_neltb
      ndum = PPICLF_LRMAX*neltbc
      call ppiclf_icopy(ppiclf_er_mapc,ppiclf_er_map,ndum)
      do ie=1,neltbc
         ppiclf_er_mapc(5,ie) = ppiclf_er_mapc(2,ie)
         ppiclf_er_mapc(6,ie) = ppiclf_er_mapc(2,ie)
      ENDdo
      nl = 0
      nii = PPICLF_LRMAX
      njj = 6
      nrr = nxyz*PPICLF_LRP_PRO
      nkey(1) = 2
      nkey(2) = 1
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl,neltbc,
     >   PPICLF_LEE,ppiclf_er_mapc,nii,partl,nl,ppiclf_pro_fldb,nrr,njj)
      call pfgslib_crystal_tuple_sort    (ppiclf_cr_hndl,neltbc
     >       ,ppiclf_er_mapc,nii,partl,nl,ppiclf_pro_fldb,nrr,nkey,2)

      ! add the fields from the bins to ptw array
      nlxyzep = nxyz*PPICLF_LEE*PPICLF_LRP_PRO
      ! TLJ changed looping to prevent -fcheck=all error
      ppiclf_pro_fld = 0.0d0
      !do i=1,nlxyzep
      !   ppiclf_pro_fld(i,1,1,1,1) = 0.0d0
      !enddo
      !do jp=1,PPICLF_LRP_PRO
      !do ip=1,PPICLF_LEE
      !do k=1,PPICLF_LEZ
      !do j=1,PPICLF_LEY
      !do i=1,PPICLF_LEX
      !   ppiclf_pro_fld(i,j,k,ip,jp) = 0.0d0
      !enddo
      !enddo
      !enddo
      !enddo
      !enddo


      do ie=1,neltbc
         iee = ppiclf_er_mapc(1,ie)
         do ip=1,PPICLF_LRP_PRO
         do k=1,PPICLF_LEZ
         do j=1,PPICLF_LEY
         do i=1,PPICLF_LEX
           ppiclf_pro_fld(i,j,k,iee,ip) = ppiclf_pro_fld(i,j,k,iee,ip) +
     >                                    ppiclf_pro_fldb(i,j,k,ip,ie)

         ENDdo
         ENDdo
         ENDdo
         ENDdo
      ENDdo

      RETURN
      END
c----------------------------------------------------------------------
      subroutine ppiclf_solve_ProjectParticleSubBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8    rproj(1+PPICLF_LRP_GP,PPICLF_LPART+PPICLF_LPART_GP)
      integer*4 iproj(4,PPICLF_LPART+PPICLF_LPART_GP)
      integer*4 ppiclf_jxgp,ppiclf_jygp,ppiclf_jzgp, nxyz, nxyzdum,
     >          idum, jdum, kdum, ic, i, j, k, ip, ndum, il, ir, jl, jr,
     >          kl, kr, jj, j1, iip, jjp, kkp
      logical if3d
      real*8 pi, d2chk2_sq, rdum, rsig, multfci, rexp, rdist2
!

      if3d = .false.
      if (ppiclf_ndim .eq. 3) if3d = .true.

      PI=4.D0*DATAN(1.D0)

      nxyz = PPICLF_BX1*PPICLF_BY1*PPICLF_BZ1

      nxyzdum = nxyz*PPICLF_LRP_PRO
      do i=1,nxyzdum
         ppiclf_grid_fld(i,1,1,1) = 0.0d0
      ENDdo

      d2chk2_sq = ppiclf_d2chk(2)**2

      ! real particles
      ppiclf_jxgp  = 1
      ppiclf_jygp  = 2
      ppiclf_jzgp  = 1
      if (if3d)
     >ppiclf_jzgp  = 3

      rdum = 0.0d0
      if (ppiclf_lfiltgauss) then
         rsig    = ppiclf_filter/(2.0d0*sqrt(2.0d0*log(2.0d0)))
         multfci = 1.0d0/(sqrt(2.0d0*pi)**2 * rsig**2) 
         if (if3d) multfci = multfci**(1.5d0)
         rdum   = 1.0d0/(-2.0d0*rsig**2)
      ENDif

      if (ppiclf_lfiltbox) then
         multfci = 1.0d0/(PI/4.0d0*ppiclf_filter**2)
         if (if3d) multfci = multfci/(1.0d0/1.5d0*ppiclf_filter)
      ENDif

      ! real particles
      do ip=1,ppiclf_npart

         rproj(1 ,ip) = rdum
         rproj(2 ,ip) = ppiclf_cp_map(ppiclf_jxgp,ip)
         rproj(3 ,ip) = ppiclf_cp_map(ppiclf_jygp,ip)
         if (if3d)
     >   rproj(4 ,ip) = ppiclf_cp_map(ppiclf_jzgp,ip)

         idum = PPICLF_LRS+PPICLF_LRP
         ic = 4
         ! TLJ modifed loop to remove out of bounds in first index
         !do j=idum+1,idum+PPICLF_LRP_GP
         do j=idum+1,PPICLF_LRP_GP
            ic = ic + 1
            rproj(ic,ip) = ppiclf_cp_map(j,ip)*multfci
         ENDdo

         iproj(1,ip) = 
     >       floor( (rproj(2,ip) - ppiclf_binx(1,1))/ppiclf_rdx)
         iproj(2,ip) = 
     >       floor( (rproj(3,ip) - ppiclf_biny(1,1))/ppiclf_rdy)
         if (if3d)
     >   iproj(3,ip) = 
     >       floor( (rproj(4,ip) - ppiclf_binz(1,1))/ppiclf_rdz)
      ENDdo

      ! ghost particles
      do ip=1,ppiclf_npart_gp

         rproj(1 ,ip+ppiclf_npart) = rdum
         rproj(2 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jxgp,ip)
         rproj(3 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jygp,ip)
         if (if3d)
     >   rproj(4 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jzgp,ip)

         idum = PPICLF_LRS+PPICLF_LRP
         ic = 4
         ! TLJ modifed loop to remove out of bounds in first index
         !do j=idum+1,idum+PPICLF_LRP_GP
         do j=idum+1,PPICLF_LRP_GP
            ic = ic + 1
            rproj(ic,ip+ppiclf_npart) = ppiclf_rprop_gp(j,ip)*multfci
         ENDdo
                    
         iproj(1,ip+ppiclf_npart) = 
     >     floor((rproj(2,ip+ppiclf_npart)-ppiclf_binx(1,1))/ppiclf_rdx)
         iproj(2,ip+ppiclf_npart) = 
     >     floor((rproj(3,ip+ppiclf_npart)-ppiclf_biny(1,1))/ppiclf_rdy)
         if (if3d)
     >   iproj(3,ip+ppiclf_npart) = 
     >     floor((rproj(4,ip+ppiclf_npart)-ppiclf_binz(1,1))/ppiclf_rdz)
      ENDdo

      ndum = ppiclf_npart+ppiclf_npart_gp

      if (ppiclf_lfiltgauss) then
         idum = floor(ppiclf_filter/2.0d0/ppiclf_rdx
     >    *sqrt(-log(ppiclf_alpha)/log(2.0d0)))+1
         jdum = floor(ppiclf_filter/2.0d0/ppiclf_rdy
     >    *sqrt(-log(ppiclf_alpha)/log(2.0d0)))+1
         kdum = 999999999
         if (if3d)
     >   kdum = floor(ppiclf_filter/2.0d0/ppiclf_rdz
     >    *sqrt(-log(ppiclf_alpha)/log(2.0d0)))+1
      ENDif

      if (ppiclf_lfiltbox) then
         idum = ppiclf_ngrids/2+1
         jdum = ppiclf_ngrids/2+1
         kdum = 999999999
         if (if3d)
     >   kdum = ppiclf_ngrids/2+1
      ENDif

      do ip=1,ndum
         iip = iproj(1,ip)
         jjp = iproj(2,ip)
         if (if3d)
     >   kkp = iproj(3,ip)

         il  = max(1     ,iip-idum)
         ir  = min(ppiclf_bx,iip+idum)
         jl  = max(1     ,jjp-jdum)
         jr  = min(ppiclf_by,jjp+jdum)
         kl  = 1
         kr  = 1
         if (if3d) then
         kl  = max(1     ,kkp-kdum)
         kr  = min(ppiclf_bz,kkp+kdum)
         ENDif

c        do k=kl,kr
c        do j=jl,jr
c        do i=il,ir
         do k=1,ppiclf_bz
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            rdist2  = (ppiclf_grid_x(i,j,k) - rproj(2,ip))**2 +
     >                (ppiclf_grid_y(i,j,k) - rproj(3,ip))**2
            if(if3d) rdist2 = rdist2 +
     >                (ppiclf_grid_z(i,j,k) - rproj(4,ip))**2

            if (rdist2 .gt. d2chk2_sq) cycle

            rexp = 1.0d0
            if (ppiclf_lfiltgauss)
     >         rexp = exp(rdist2*rproj(1,ip))

            do jj=1,PPICLF_LRP_PRO
               j1 = jj+4
               ppiclf_grid_fld(i,j,k,jj) = 
     >                         ppiclf_grid_fld(i,j,k,jj) 
     >                       + sngl(rproj(j1,ip)*rexp)
            ENDdo
         ENDdo
         ENDdo
         ENDdo
      ENDdo

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_GetProFldIJKEF(i,j,k,e,m,fld)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 i,j,k,e,m
!
! Output:
!
      real*8 fld
!
      fld = ppiclf_pro_fld(i,j,k,e,m)

      RETURN
      END
!-----------------------------------------------------------------------
