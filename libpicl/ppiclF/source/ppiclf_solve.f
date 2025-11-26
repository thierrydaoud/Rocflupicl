      SUBROUTINE ppiclf_solve_Initialize(xi1,xpmin,xpmax,
     >           yi1,ypmin,ypmax,zi1,zpmin,zpmax,
     >           ai1,apa,apxa,aprin,aprout)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4 xi1, yi1, zi1, ai1
      REAL*8    xpmin,xpmax,ypmin,ypmax,zpmin,zpmax,
     >          apa,apxa,aprin,aprout, pi, angled

      ! Called by rocpicl/PICL_TEMP_InitSolver.F90
      ! xdrange adjusts the bin boundaries to ensure they aren't 
      ! larger than the cartesian fluid domain extremes.

      ! Linear X-Periodicity
      ppiclf_xdrange(1,1) = xpmin
      ppiclf_xdrange(2,1) = xpmax
      x_per_flag = xi1
      IF(x_per_flag.EQ.1) THEN
        IF(xpmin .ge. xpmax) THEN
          PRINT*,'Failed in Initialize.'
          PRINT*,'xpMin > xpMax'
          CALL ppiclf_exittr('PeriodicX must have xmin < xmax$',xpmin,0)
        END IF
        ppiclf_linperiodic(1) = .TRUE.
      END IF

      ! Linear Y-Periodicity
      y_per_flag = yi1
      ppiclf_xdrange(1,2) = ypmin
      ppiclf_xdrange(2,2) = ypmax
      IF(y_per_flag.EQ.1) THEN
        IF(ypmin .ge. ypmax) THEN
          PRINT*,'Failed in Initialize.'
          PRINT*,'ypMin > ypMax'
          CALL ppiclf_exittr('PeriodicX must have ymin < ymax$',ypmin,0)
        END IF
        ppiclf_linperiodic(2) = .TRUE.
      END IF

      ! Linear Z-Periodicity
      z_per_flag = zi1
      ppiclf_xdrange(1,3) = zpmin
      ppiclf_xdrange(2,3) = zpmax
      IF(z_per_flag.EQ.1) THEN
        IF(zpmin .ge. zpmax) THEN
          PRINT*,'Failed in Initialize.'
          PRINT*,'ypMin > ypMax'
          CALL ppiclf_exittr('PeriodicZ must have zmin < zmax$',zpmin,0)
        END IF
        ppiclf_linperiodic(3) = .TRUE.
      END IF

!*** THIS CODE WILL CHANGE
!      ! Angular Periodicity
!      ang_per_flag = ai1
!      IF(ang_per_flag.EQ.1) THEN
!        ppiclf_linperiodic(1) = .TRUE. ! X-Periodicity
!        ppiclf_linperiodic(2) = .TRUE. ! Y-Periodicity
!        ang_per_angle  = apa
!        ang_per_xangle = apxa
!        ang_per_rin    = aprin
!        ang_per_rout   = aprout
!      END IF
!
!      ! User cannot initialize X/Y-Periodicity with Angular Periodicity
!      IF(((x_per_flag.EQ.1).OR.(y_per_flag.EQ.1))
!     >                     .AND.(ang_per_flag.EQ.1))
!     >   CALL ppiclf_exittr('PPICLF: Invalid Periodicity choice$',0,0)
!
!      ! Thierry - compute ang_case
!
!      pi = ACOS(-1.0)
!      angled = ang_per_angle * 180.0d0 / pi ! store angle value in degrees
!
!      IF(ang_per_flag.EQ.0) THEN
!         ang_case = 0 ! standard geometry
!      ELSE
!         IF(angled .lt. 90.0)        ang_case = 1 ! general wedge
!         IF(NINT(angled) .EQ. 90.0)  ang_case = 2 ! quarter cylinder
!         IF(NINT(angled) .EQ. 180.0) ang_case = 3 ! half cylinder
!      END IF
!
!      IF(ppiclf_nid.EQ.0 .AND. ang_case.NE.0) THEN
!         PRINT*, " "
!         PRINT*, " ======================================="
!         PRINT*, " "
!         PRINT*, "!!! PPICLF Angular Periodicity Initialized !!!!"
!         PRINT*, "  Angular periodicity flag =", ang_per_flag
!         PRINT*, "  Init Angular- angle =", ang_per_angle
!         PRINT*, "  Init Angular- angled =", angled
!         PRINT*, "  Init Angular- nint(angled) =", NINT(angled)
!         PRINT*, "  Init Angular- ang_case =", ang_case
!         PRINT*, " "
!         PRINT*, " ======================================="
!         PRINT*, " "
!      END IF
!! *** END CHANGE ***
!
      RETURN
      END
!
!-----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_solve_InitParticle(imethod,ndim,iendian,npart,y,
     >                                     rprop,filt2,filt3)
     > bind(C, name="ppiclc_solve_InitParticle")
#else
      SUBROUTINE ppiclf_solve_InitParticle(imethod,ndim,iendian,npart,y,
     >                                     rprop,filt2,filt3)
#endif
!
! Called from rocpicl/PICL_TEMP_InitSolver.F90

      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE 'mpif.h'
!
! Input: 
!
      INTEGER*4  imethod ! From rocpicl: 2 (Same RK3 as Rocflu)
      INTEGER*4  ndim    ! From rocpicl: 3
      INTEGER*4  iendian ! From rocpicl: 0
      INTEGER*4  npart
      INTEGER*4  l
      REAL*8     y(*)
      REAL*8     rprop(*)
      REAL*8     filt2(3)
      REAL*8     filt3

      IF(.NOT. PPICLF_LCOMM)
     >CALL ppiclf_exittr('InitMPI must be before InitParticle$',0.0d0
     >   ,ppiclf_nid)
      IF(PPICLF_OVERLAP)
     >CALL ppiclf_exittr('InitFilter must be before InitOverlap$',0.0d0
     >                  ,0)

      CALL ppiclf_prints('*Begin InitParticle$')
      CALL ppiclf_prints('   *Begin InitParam$')

      CALL ppiclf_solve_InitParam(imethod,ndim,iendian)
      DO l = 1,3
        ppiclf_filter(l) = filt2(l)
      END DO
      ppiclf_nndist = filt3
      
      CALL ppiclf_prints('    End InitParam$')

      IF(.NOT. PPICLF_RESTART) THEN
        CALL ppiclf_prints('   *Begin InitZero$')
        CALL ppiclf_solve_InitZero
        CALL ppiclf_prints('   *End InitZero$')
        CALL ppiclf_prints('   *Begin AddParticles$')
        CALL ppiclf_solve_AddParticles(npart,y,rprop)
        CALL ppiclf_prints('   *End AddParticles$')
      END IF

      CALL ppiclf_prints(' End InitParticle$')
!
      PPICLF_LINIT = .TRUE.

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InitParam(imethod,ndim,iendian)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4  imethod
      INTEGER*4  ndim
      INTEGER*4  iendian
!
      IF(imethod .EQ. 0 .OR. imethod .GE. 3 .OR. imethod .LE. -2)
     >   CALL ppiclf_exittr('Invalid integration method$',0.0d0,imethod)
      IF(ndim .NE. 3) THEN
        PRINT*, 'ERROR: Rocpicl is written for 3D grids.'
        CALL ppiclf_exittr('Invalid problem dimension$',0.0d0,ndim)
      END IF
      IF(iendian .LT. 0 .OR. iendian .GT. 1)
     >   CALL ppiclf_exittr('Invalid Endian$',0.0d0,iendian)

      ppiclf_imethod      = imethod
      ppiclf_ndim         = ndim
      ppiclf_iendian      = iendian

      ppiclf_linperiodic(1) = .FALSE.    
      ppiclf_linperiodic(2) = .FALSE.   
      ppiclf_linperiodic(3) = .FALSE. 

      ppiclf_cycle  = 0
      ppiclf_iostep = 1
      ppiclf_dt     = 0.0d0
      ppiclf_time   = 0.0d0

!      ppiclf_readytosolve = .FALSE.
      ppiclf_overlap      = .FALSE.
      ppiclf_linit        = .FALSE.
      ppiclf_lintp        = .FALSE.
      ppiclf_lproj        = .FALSE.
      ppiclf_binchanged   = .TRUE.
      ppiclf_printbinvtu  = .TRUE.
      IF(PPICLF_INTERP .EQ. 1)  ppiclf_lintp = .TRUE.
      IF(PPICLF_PROJECT .EQ. 1) ppiclf_lproj = .TRUE.

      ppiclf_xdrange(1,1) = -1E20
      ppiclf_xdrange(2,1) =  1E20
      ppiclf_xdrange(1,2) = -1E20
      ppiclf_xdrange(2,2) =  1E20
      ppiclf_xdrange(1,3) = -1E20
      ppiclf_xdrange(2,3) =  1E20

      ppiclf_filter(1)        = 0.0D0
      ppiclf_filter(2)        = 0.0D0
      ppiclf_filter(3)        = 0.0D0
      ppiclf_nndist           = 0.0D0
      ppiclf_interp_dchk(1) = 0.0D0
      ppiclf_interp_dchk(2) = 0.0D0
      ppiclf_interp_dchk(3) = 0.0D0
 
      ppiclf_n_bins(1) = 1
      ppiclf_n_bins(2) = 1
      ppiclf_n_bins(3) = 1

      ppiclf_binb(1) = 0.0
      ppiclf_binb(2) = 0.0
      ppiclf_binb(3) = 0.0
      ppiclf_binb(4) = 0.0
      ppiclf_binb(5) = 0.0
      ppiclf_binb(6) = 0.0

      ppiclf_previousbinb(1) =  1.0E9
      ppiclf_previousbinb(2) = -1.0E9
      ppiclf_previousbinb(3) =  1.0E9
      ppiclf_previousbinb(4) = -1.0E9
      ppiclf_previousbinb(5) =  1.0E9
      ppiclf_previousbinb(6) = -1.0E9

      ppiclf_bins_dx(1) = 1.0
      ppiclf_bins_dx(2) = 1.0
      ppiclf_bins_dx(3) = 1.0     

      ppiclf_nwall    = 0
      ppiclf_iwallm   = 0

      PPICLF_INT_ICNT = 0

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InitZero
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal:
!
      INTEGER*4 i,j,ie
!
      ! zero'ing real particle properties
      DO i=1,PPICLF_LPART 
        DO j=1,PPICLF_LRS
          ppiclf_y    (j,i) = 0.0D0
          ppiclf_ydot (j,i) = 0.0D0
          ppiclf_ydotc(j,i) = 0.0D0
          ppiclf_y1   (j,i) = 0.0D0
        END DO
        DO j=1,PPICLF_LRP
           ppiclf_rprop(j,i) = 0.0D0
        END DO
        DO j=1,PPICLF_LRP2
           ppiclf_rprop2(j,i) = 0.0D0
        END DO
        DO j=1,PPICLF_LRP3
           ppiclf_rprop3(j,i) = 0.0D0
        END DO
        DO j=1,PPICLF_LRP4
           ppiclf_rprop4(j,i) = 0.0D0
        END DO
        DO j=1,PPICLF_LRP5
           ppiclf_rprop5(j,i) = 0.0D0
        END DO
        DO j=1,PPICLF_LIP
           ppiclf_iprop(j,i) = 0
        END DO
        DO j = 1,PPICLF_LRP_PRO
          ppiclf_feedbk(j,i) = 0.0D0
        END DO
      END DO
      ! zero'ing ghost particle properties
      DO i=1,PPICLF_LPART_GP
        DO j =1,PPICLF_LIP_GP
          ppiclf_iprop_gp(j,i) = 0
        END DO
        DO j=1,PPICLF_LRP_GP
           ppiclf_rprop_gp(j,i) = 0.0D0
        END DO
      END DO

      ppiclf_npart = 0
      ! zero'ing grid properties for interpolation
      DO ie=1,PPICLF_LEE
        DO j=1,PPICLF_LRP_INT
          ppiclf_int_fld(j,ie) = 0.0D0
        END DO
      END DO

      CALL ppiclf_user_InitZero

      RETURN
      END
!-----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_solve_AddParticles(npart,y,rprop)
     > bind(C, name="ppiclc_solve_AddParticles")
#else
      SUBROUTINE ppiclf_solve_AddParticles(npart,y,rprop)
#endif
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input: 
!
      INTEGER*4  npart
      REAL*8     y(*)
      REAL*8     rprop(*)
!
! Internal:
!
      INTEGER*4 ppiclf_iglsum,ntotal,i,j
      external ppiclf_iglsum
!

      CALL ppiclf_prints('   *Begin AddParticles$')

      IF(ppiclf_npart+npart .gt. PPICLF_LPART .or. npart .lt. 0)
     >   CALL ppiclf_exittr('Invalid number of particles$',
     >                      0.0D0,ppiclf_npart+npart)

      CALL ppiclf_printsi('      -Begin copy particles$',npart)

      CALL ppiclf_copy(ppiclf_y(1,ppiclf_npart+1),
     >                 y,
     >                 npart*PPICLF_LRS)
      CALL ppiclf_copy(ppiclf_rprop(1,ppiclf_npart+1),
     >                 rprop, npart*PPICLF_LRP)

      ppiclf_npart = ppiclf_npart + npart

      CALL ppiclf_printsi('      -Begin copy particles$',ppiclf_npart)

      IF(.NOT. PPICLF_RESTART) THEN
         CALL ppiclf_prints('      -Begin ParticleTag$')
            CALL ppiclf_solve_SetParticleTag(npart)
         CALL ppiclf_prints('       End ParticleTag$')
      END IF

      CALL ppiclf_prints('    End AddParticles$')

      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_SetParticleTag(npart)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Input: 
! 
      INTEGER*4 npart
! 
! Internal: 
! 
      INTEGER*4 i
!
      DO i=ppiclf_npart-npart+1,ppiclf_npart
         ppiclf_iprop(1,i) = i
         ppiclf_iprop(2,i) = ppiclf_nid
         ppiclf_iprop(3,i) = ppiclf_cycle
      END DO

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
      REAL*8 ydum(PPICLF_LRS), rpropdum(PPICLF_LRP), xp(3), bin_xMin(3),
     >       A(3),B(3),C(3),AB(3),AC(3), distSQ, xdistSQ, ydistSQ,
     >       dist_total, rnx, rny, rnz, area, rpx1, rpy1, rpz1, rpx2,
     >       rpy2, rpz2, rflip, a_sum, rd, rdist, theta, tri_area,
     >       ab_dot_ac, ab_mag, ac_mag, zdistSQ, rthresh
      INTEGER*4 istride, k, kmax, kp, kkp, kk, j, jp, l, iSB, jSB, kSB,
     >          loopSB, tempSB, iSBin(3)
! 
      distSQ = ppiclf_nndist**2
      
      ! find ith particle subbin (tempSB)
      DO l = 1,3
        xp(l) = ppiclf_y(l,i)
        bin_xMin(l) = ppiclf_bin_pos(1,l)
      END DO
      DO l = 1,3
        iSBin(l) = FLOOR((xp(l) - (bin_xMin(l)
     >         - ppiclf_nndist))/ppiclf_nndist)
      END DO
      tempSB = iSBin(1) + iSBin(2)*SBn(1) + iSBin(3)*SBn(1)*SBn(2)
#ifdef TEST
      PARTICLE_NN(i) = 0 
      PPICLF_TOTNNDIST(i) = 0.0D0
#endif
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
                IF (j .EQ. i) CYCLE ! Same particle
                xdistSQ = (ppiclf_cp_map(1,i)-ppiclf_cp_map(1,j))**2
                IF (xdistSQ .GT. distSQ) CYCLE
                ydistSQ = (ppiclf_cp_map(2,i)-ppiclf_cp_map(2,j))**2
                IF (ydistSQ .GT. distSQ) CYCLE
                dist_total = xdistSQ + ydistSQ
                zdistSQ = (ppiclf_cp_map(3,i)-ppiclf_cp_map(3,j))**2
                IF (zdistSQ .GT. distSQ) CYCLE
                dist_total = dist_total+zdistSQ
                IF (dist_total .GT. distSQ) CYCLE
#ifdef TEST
                PARTICLE_NN(i) = PARTICLE_NN(i) + 1
                PPICLF_TOTNNDIST(i) = PPICLF_TOTNNDIST(i) + dist_total
                CYCLE !Don't want to call EvalNN. Just testing
                      ! nneighbor search
#endif
                CALL ppiclf_user_EvalNearestNeighbor(i,j
     >                                   ,ppiclf_cp_map(1,i)
     >                                   ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                   ,ppiclf_cp_map(1,j)
     >                                   ,ppiclf_cp_map(1+PPICLF_LRS,j))
              ELSE IF (j .LT. 0) THEN ! Ghost Particle
                ! Negative was just use for ghost particle indicator
                ! in subbin mapping array. Need to flip sign
                j = - j                 
                xdistSQ = (ppiclf_cp_map(1,i)-ppiclf_rprop_gp(1,j))**2
                IF (xdistSQ .GT. distSQ) CYCLE
                ydistSQ = (ppiclf_cp_map(2,i)-ppiclf_rprop_gp(2,j))**2
                IF (ydistSQ .GT. distSQ) CYCLE
                dist_total = xdistSQ + ydistSQ
                IF (ppiclf_ndim .EQ. 3) THEN
                zdistSQ = (ppiclf_cp_map(3,i)-ppiclf_rprop_gp(3,j))**2
                IF (zdistSQ .GT. distSQ) CYCLE
                dist_total = dist_total+zdistSQ
                END IF
                IF (dist_total .GT. distSQ) CYCLE
#ifdef TEST
                PARTICLE_NN(i) = PARTICLE_NN(i) + 1
                PPICLF_TOTNNDIST(i) = PPICLF_TOTNNDIST(i) + dist_total
                CYCLE !Don't want to call EvalNN. Just testing
                      ! nneighbor search
#endif
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
         rnz  = ppiclf_wall_n(3,j)
         area = ppiclf_wall_n(4,j)
         rpz1 = ppiclf_cp_map(3,i)
         rpz2 = ppiclf_wall_c(3,j)
         rpz2 = rpz2 - rpz1
    
         rflip = rnx*rpx2 + rny*rpy2 + rnz*rpz2
         IF(rflip .GT. 0.0d0) THEN
            rnx = -1.0d0*rnx
            rny = -1.0d0*rny
            rnz = -1.0d0*rnz
         END IF


         a_sum = 0.0d0
         kmax = 3
         DO k=1,kmax 
            kp = k+1
            IF(kp .GT. kmax) kp = kp-kmax ! cycle
            
            kk   = istride*(k-1)
            kkp  = istride*(kp-1)
            rpx1 = ppiclf_wall_c(kk+1,j)
            rpy1 = ppiclf_wall_c(kk+2,j)
            rpz1 = 0.0d0
            rpx2 = ppiclf_wall_c(kkp+1,j)
            rpy2 = ppiclf_wall_c(kkp+2,j)
            rpz2 = 0.0d0

            rpz1 = ppiclf_wall_c(kk+3,j)
            rpz2 = ppiclf_wall_c(kkp+3,j)

            rd   = -(rnx*rpx1 + rny*rpy1 + rnz*rpz1)

            rdist = abs(rnx*ppiclf_cp_map(1,i)+rny*ppiclf_cp_map(2,i)
     >                 +rnz*ppiclf_cp_map(3,i)+rd)
            rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

            ! give a little extra room for walls (2x)
            IF(rdist .GT. 2.0d0*ppiclf_nndist) GOTO 1519

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
            a_sum = a_sum + tri_area
         END DO

         rthresh = 1.10d0 ! keep it from slipping through crack on edges
         IF(a_sum .GT. rthresh*area) CYCLE

         jp = 0
         CALL ppiclf_user_EvalNearestNeighbor(i,jp,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ydum
     >                                 ,rpropdum)

 1519 continue
      ENDdo

      RETURN
      END
!-----------------------------------------------------------------------
       SUBROUTINE ppiclf_solve_InitWall(xp1,xp2,xp3)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Input:
! 
      REAL*8 xp1(*)
      REAL*8 xp2(*)
      REAL*8 xp3(*)
!
! Internal:
!
      REAL*8 rpx1, rpy1, rpz1, rpx2, rpy2, rpz2,
     >       a_sum, theta, tri_area, 
     >       ab_dot_ac, ab_mag, ac_mag, rise, run, rmag, 
     >       rpx3, rpy3, rpz3
      INTEGER*4 istride, k, kmax, kp, kkp, kk
      REAL*8 A(3),B(3),C(3),AB(3),AC(3)
!
      if (.not.PPICLF_LCOMM)
     >CALL ppiclf_exittr('InitMPI must be before InitWall$',0.d0,0)
      if (.not.PPICLF_LINIT)
     >CALL ppiclf_exittr('InitParticle must be before InitWall$'
     >                  ,0.d0,0)

      ppiclf_nwall = ppiclf_nwall + 1 

      if (ppiclf_nwall .gt. PPICLF_LWALL)
     >CALL ppiclf_exittr('Increase LWALL in user file$'
     >                  ,0.d0,ppiclf_nwall)

      istride = ppiclf_ndim
      a_sum = 0.0d0
      kmax = 2
      if (ppiclf_ndim .EQ. 3) kmax = 3

      if (ppiclf_ndim .EQ. 3) then
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
      elseif (ppiclf_ndim .EQ. 2) then
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

         if (ppiclf_ndim .EQ. 3) then
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
         elseif (ppiclf_ndim .EQ. 2) then
             AB_MAG = sqrt(AB(1)**2 + AB(2)**2)
             tri_area = AB_MAG
         ENDif
         a_sum = a_sum + tri_area
      ENDdo
      
      ppiclf_wall_n(ppiclf_ndim+1,ppiclf_nwall) = a_sum

      ! wall normal:
      if (ppiclf_ndim .EQ. 2) then

         rise = xp2(2) - xp1(2)
         run  = xp2(1) - xp1(1)

         rmag = sqrt(rise**2 + run**2)
         rise = rise/rmag
         run  = run/rmag
         
         ppiclf_wall_n(1,ppiclf_nwall) = rise
         ppiclf_wall_n(2,ppiclf_nwall) = -run

      elseif (ppiclf_ndim .EQ. 3) then

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

      END IF

      RETURN
      END
c----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_solve_WriteVTU(time)
     > bind(C, name="ppiclc_solve_WriteVTU")
#else
      SUBROUTINE ppiclf_solve_WriteVTU(time)
#endif
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Input: 
! 
      REAL*8    time
! 
! Internal:
!
      ppiclf_time   = time

      CALL ppiclf_io_WriteParticleVTU('')
      CALL ppiclf_io_WriteBinVTU('')
      ! Output diagnostics
      CALL ppiclf_io_OutputDiagAll

      RETURN
      END
c----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_solve_IntegrateParticle(istep,iostep,dt,time)
     > bind(C, name="ppiclc_solve_IntegrateParticle")
#else
      SUBROUTINE ppiclf_solve_IntegrateParticle(istep,iostep,dt,time)
#endif
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Input: 
! 
      INTEGER*4 istep
      INTEGER*4 iostep
      REAL*8    dt
      REAL*8    time
! 
! Internal:
!
      LOGICAL iout
!
      ppiclf_cycle  = istep
      ppiclf_iostep = iostep
      ppiclf_dt     = dt
      ppiclf_time   = time

      ! integerate in time
!************************************************
! NOT USING THESE in rocpicl!
!      if (ppiclf_imethod .EQ. 1) 
!     >   CALL ppiclf_solve_IntegrateRK3(iout)
!      if (ppiclf_imethod .EQ. -1) 
!     >   CALL ppiclf_solve_IntegrateRK3s(iout)
!************************************************
      IF(ppiclf_imethod .EQ. 2) THEN
        CALL ppiclf_solve_IntegrateRK3s_Rocflu(iout)
      ELSE
        PRINT*, 'ERROR: Wrong RK selected for rocpicl. Use RK3!'
        CALL ppiclf_exittr('Wrong RK for rocpicl',0.D0,0)
      END IF

      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_IntegrateRK3s_Rocflu(iout)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Internal: 
! 
      INTEGER*4 i, ndum, nstage, istage
      INTEGER*4 icalld
      INTEGER*4 j
      save      icalld
      data      icalld /0/

      EXTERNAL  ppiclf_glsum
      REAL*8    ppiclf_glsum
      REAL*8    fxsum, fysum, fzsum, fxabssum, fyabssum, fzabssum

!
! Output:
!
      LOGICAL iout
!
      icalld = icalld + 1

      ! get rk3 coeffs
      CALL ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      istage = MOD(icalld,nstage)
      if (istage .EQ. 0) istage = 3
      iout = .FALSE.
      if (istage .EQ. nstage) iout = .TRUE.

      ! evaluate ydot
      CALL ppiclf_solve_SetYdot

      !Zero out for first stage
      if (istage .EQ. 1) then
        ppiclf_y1 = 0.0D0
      END IF

      ! The Rocflu RK3 can be found in equation (7) of:
      ! S. Yu. "Runge-Kutta Methods Combined with Compact
      !   Difference Schemes for the Unsteady Euler Equations".
      !   Center for Modeling of Turbulence and Transition.
      !   Research Briefs, 1991.

      DO i = 1,PPICLF_NPART
        DO j = 1,PPICLF_LRS
          ppiclf_y(j,i) =  - ppiclf_rk3coef(1,istage)*ppiclf_y1   (j,i)
     >                     + ppiclf_rk3coef(2,istage)*ppiclf_y    (j,i)
     >                     + ppiclf_rk3coef(3,istage)*ppiclf_ydot (j,i)
        END DO
      END DO
      
      IF(ppiclf_linperiodic(1) .OR. ppiclf_linperiodic(2) .OR.
     >                             ppiclf_linperiodic(3)) THEN
        CALL ppiclf_solve_PeriodicParticleShift
      END IF

      !Store Current stage RHS for next stage's use
      DO i= 1,PPICLF_NPART
        DO j= 1,PPICLF_LRS
          ppiclf_y1(j,i) =  ppiclf_ydot(j,i)
        END DO
      END DO
      
!      PPICLF_READYTOSOLVE = .FALSE.
!      CALL ppiclf_solve_PostSolve

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_PeriodicParticleShift
!
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal:
!
      INTEGER*4 i  
      REAL*8 per_alpha
!
      DO i=1,ppiclf_npart
!!!!!!!!!!!!!!!!        Rotational Periodicity Starts Here     !!!!!!!!!!!!!!!!!!!!
!            ! currently only supports angular rotation around z-axis
!            ! and only check theta component and not radial
!            ! applying the radial periodicity is very straightforward, but not needed for now
!
!         IF(ang_per_flag .EQ. 1) THEN  ! Angular periodicity
!           
!           ! particle angle w/ x-axis
!           ! per_alpha here is obtained in radians
!           ! ang_per_angle & ang_per_xangle are transformed 
!           !   to radians in PICL_TEMP_InitFlowSolver.F90
!           per_alpha = 
!     >         atan2(ppiclf_y(PPICLF_JY,i), ppiclf_y(PPICLF_JX,i))
!
!           ! check if particle leaving through lower face or upper face of wedge
!           IF((per_alpha .LT. ang_per_xangle) .OR. 
!     >         (per_alpha .GT. (ang_per_xangle + ang_per_angle))) THEN
!             CALL ppiclf_solve_InvokeAngularPeriodic(i, ang_per_flag,
!     >                                                per_alpha, 
!     >                                                ang_per_angle,
!     >                                                ang_per_xangle,
!     >                                                1)
!           END IF ! per_alpha
!         END IF ! ang_per_flag
!
         ! Linear Periodicity Invoked
         IF(ppiclf_linperiodic(1) .OR. ppiclf_linperiodic(2) 
     >                            .OR. ppiclf_linperiodic(3)) THEN
           CALL ppiclf_solve_InvokeLinearPeriodic(i)
         END IF 
      END DO ! i=1,ppiclf_part

      RETURN
      END

!----------------------------------------------------------------------- 

      SUBROUTINE ppiclf_solve_InvokeLinearPeriodic(i)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Internal: 
! 
      INTEGER*4 i, j
!

! Case 1 - Linear Periodicity in any of 3 directions ; NO Anuglar Periodicity
      IF(((x_per_flag.EQ.1).OR.(y_per_flag.EQ.1).OR.(z_per_flag.EQ.1))
     >   .AND.(ang_per_flag.EQ.0)) THEN

        DO j= 1,3
          ! particle leaving min. periodic face -> move it relative to 
          !                                         max periodic face
          IF(ppiclf_y(j,i) .LT. ppiclf_xdrange(1,j)) THEN
             ppiclf_y(j,i) = ppiclf_xdrange(2,j) - 
     >                    ABS(ppiclf_xdrange(1,j) - ppiclf_y(j,i))
            CYCLE
          END IF

          ! particle leaving max. periodic face -> move it relative to 
          !                                         min periodic face
          IF(ppiclf_y(j,i).GT.ppiclf_xdrange(2,j)) THEN
             ppiclf_y(j,i) = ppiclf_xdrange(1,j) + 
     >                    ABS(ppiclf_y(j,i) - ppiclf_xdrange(2,j))
          END IF
        END DO ! j
      END IF 
 
      RETURN
      END
!-----------------------------------------------------------------------

      SUBROUTINE ppiclf_solve_SetRK3Coeff(dt)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      REAL*8 dt
!
      IF(ppiclf_imethod .EQ. 2) THEN
        !BD:Rocflu's rk3 scheme
        !Folowing form
        !rk3(1,:) = Temp storage i.e. Previous Stage RHS
        !rk3(2,:) = Temp storage i.e. Current Stage iteration
        !rk3(3,:) = Temp storage i.e. Current Stage RHS
        ppiclf_rk3ark(1) = 8.0d0/15.0d0
        ppiclf_rk3ark(2) = 5.0d0/12.0d0
        ppiclf_rk3ark(3) = 0.75d0

        ppiclf_rk3coef(1,1) = 0.d00
        ppiclf_rk3coef(2,1) = 1.0d0
        ppiclf_rk3coef(3,1) = dt*8.0d0/15.0d0
        ppiclf_rk3coef(1,2) = dt*17.0d0/60.0d0
        ppiclf_rk3coef(2,2) = 1.0d0
        ppiclf_rk3coef(3,2) = dt*5.0d0/12.0d0
        ppiclf_rk3coef(1,3) = dt*5.0d0/12.0d0
        ppiclf_rk3coef(2,3) = 1.0d0
        ppiclf_rk3coef(3,3) = dt*3.0d0/4.0d0
!      ELSE
!        !BD:Original Code This follows CMT-nek's rk 3 scheme
!        ppiclf_rk3coef(1,1) = 0.d00
!        ppiclf_rk3coef(2,1) = 1.0d0 
!        ppiclf_rk3coef(3,1) = dt
!        ppiclf_rk3coef(1,2) = 3.0d0/4.0d0
!        ppiclf_rk3coef(2,2) = 1.0d0/4.0d0 
!        ppiclf_rk3coef(3,2) = dt/4.0d0
!        ppiclf_rk3coef(1,3) = 1.0d0/3.0d0
!        ppiclf_rk3coef(2,3) = 2.0d0/3.0d0
!        ppiclf_rk3coef(3,3) = dt*2.0d0/3.0d0
!        !BD: Original Code END
      END IF

      RETURN
      END

!----------------------------------------------------------------------

      SUBROUTINE ppiclf_solve_SetYdot
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
!      IF(.NOT. PPICLF_READYTOSOLVE)
!     >  CALL ppiclf_solve_InitSolve
      CALL ppiclf_solve_InitSolve
      CALL ppiclf_user_SetYdot

      RETURN
      END
!----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InitSolve
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
! 
! Internal: 
! 
      INTEGER*4 :: i, j
      ! ppiclf_binchanged set in CreateBin
      ! ppiclf_binchanged .TRUE. means
      ! bin coordinates changed
      CALL ppiclf_comm_CreateBin

      ! ppiclf_particleMoved set in FindParticle
      ! ppiclf_particleMoved .EQ. 0 means all particles
      ! stayed in same bin as previous RK Stage.
      CALL ppiclf_comm_FindParticle
      IF(ppiclf_particleMoved .NE. 0 .OR.
     >              ppiclf_binchanged) THEN
        CALL ppiclf_comm_MoveParticle
      END IF

      IF(ppiclf_overlap .AND. ppiclf_binchanged) THEN
        CALL ppiclf_comm_MapOverlapGrid
      END IF

      IF(ppiclf_overlap) THEN
        ! Interpolate fluid solver grid to particle
        CALL ppiclf_solve_InterpParticleGrid
        ! Project particle feedback to fluid solver grid
        CALL ppiclf_solve_ProjectParticleGrid
      END IF

!      IF(ppiclf_gprequired) THEN
      ! Ghost particles are needed 
        CALL ppiclf_comm_CreateGhost
        CALL ppiclf_comm_MoveGhost
!      END IF

      ! Zero collisions 
      ppiclf_ydotc = 0.0D0

      RETURN
      END

!----------------------------------------------------------------------

!      SUBROUTINE ppiclf_solve_PostSolve
!!
!      IMPLICIT NONE
!!
!      INCLUDE "PPICLF"
!! 
!! Internal: 
!! 
!      INTEGER*4 :: i, j, ierr
!
!      ! ppiclf_binchanged set in CreateBin
!      ! ppiclf_binchanged .TRUE. means
!      ! bin coordinates changed
!      CALL ppiclf_comm_CreateBin
!   
!      CALL MPI_BARRIER(ppiclf_comm,ierr)
!
!      ! ppiclf_particleMoved set in FindParticle
!      ! ppiclf_particleMoved .EQ. 0 means all particles
!      ! stayed in same bin as previous RK Stage.
!      CALL ppiclf_comm_FindParticle
!
!      IF(ppiclf_particleMoved .NE. 0 .OR.
!     >              ppiclf_binchanged) THEN
!        CALL ppiclf_comm_MoveParticle
!      END IF
!
!      IF(ppiclf_overlap .AND. ppiclf_binchanged) THEN
!        CALL ppiclf_comm_MapOverlapGrid
!      END IF
!
!      IF(ppiclf_overlap) THEN
!        ! Interpolate fluid solver grid to particle
!        CALL ppiclf_solve_InterpParticleGrid
!        ! Project particle feedback to fluid solver grid
!        CALL ppiclf_solve_ProjectParticleGrid
!      END IF
!
!!      IF(ppiclf_gprequired) THEN
!      ! Ghost particles are needed 
!        CALL ppiclf_comm_CreateGhost
!        CALL ppiclf_comm_MoveGhost
!!      END IF
!
!      ! Zero collisions 
!      ppiclf_ydotc = 0.0D0
!
!      RETURN
!      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InterpParticleGrid
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal:
!
      INTEGER*4 j
!
      ! Copies Grid Cell ID for all Rocflu elements that map
      ! to ppiclf domain for GSLIB Transfer.  This copy is from
      ! MapOverlapGrid.
      CALL ppiclf_solve_InitInterp

      ! Makes array (ppiclf_int_fld_input) of all rprop data
      ! for grid cellss that map to ppiclf domain.
      DO j=1,PPICLF_INT_ICNT
         CALL ppiclf_solve_InterpField(j)
      END DO
      
      ! Transfers ppiclf_er_mapc & ppiclf_int_fld for all Rocflu Grid
      ! cells that map to ppiclf domain.
      CALL ppiclf_solve_InterpTupleTransfer

      ! Maps up to 27 closest cell centers to particle
      ! Includes: CellID, total dist, x dist, y dist, z dist
      CALL ppiclf_solve_SBParticleToCellMap

      ! Interpolates rprop data for ppiclf domain cells in this bin
      CALL ppiclf_solve_Interpolate

      PPICLF_INT_ICNT = 0


      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InterpFieldUser(jp,infld)
!
! This is called by rocpicl/PICL_TEMP_Runge.F90
! There is a call for each quantity that should be interpolated
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input: 
!

      INTEGER*4 jp,i !rprop index
      REAL*8 infld(*) !value to set rprop to
!
! Internal:
!
      INTEGER*4 n
!
      IF(PPICLF_INTERP .EQ. 0)
     >CALL ppiclf_exittr(
     >     'No specified interpolated fields, set PPICLF_LRP_INT$',0.0d0
     >                   ,0)

      PPICLF_INT_ICNT = PPICLF_INT_ICNT + 1

      IF(PPICLF_INT_ICNT .GT. PPICLF_LRP_INT)
     >   CALL ppiclf_exittr('Interpolating too many fields$'
     >                     ,0.0d0,PPICLF_INT_ICNT)
      IF(jp .LE. 0 .OR. jp .GT. PPICLF_LRP)
     >   CALL ppiclf_exittr('Invalid particle array interp. location$'
     >                     ,0.0d0,jp)

      ! set up interpolation map
      PPICLF_INT_MAP(PPICLF_INT_ICNT) = jp

      ! copy to infld internal storage
      n = ppiclf_nFVCells
      CALL ppiclf_copy(ppiclf_int_fld_input(1,PPICLF_INT_ICNT)
     >                ,infld(1),n)

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InitInterp
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Internal: 
! 
      INTEGER*4 ie
!
      IF(.NOT.ppiclf_overlap)
     >CALL ppiclf_exittr('Cannot interpolate unless overlap grid$',0.0d0
     >                   ,0)
      IF(.NOT.ppiclf_lintp) 
     >CALL ppiclf_exittr('To interpolate, set PPICLF_LRP_PRO to ~= 0$'
     >                   ,0.0d0,0)
      ppiclf_nCells_Interp = ppiclf_nCells_FV2PICL_Orig
      DO ie=1,ppiclf_nCells_Interp
        CALL ppiclf_icopy(ppiclf_cell_map_interp(1,ie),
     >        ppiclf_cell_map_Orig(1,ie), PPICLF_LRMAX)
      END DO
      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InterpField(j)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input: 
!
      INTEGER*4 jp
!
! Internal:
!
      INTEGER*4 n, ie, iee, j
!
      ! use the map to take original grid and map to fld which will be
      ! sent to mapped processors
      DO ie=1,ppiclf_nCells_Interp
         ! iee is the Rocflu element ID from previous MapOverlapGrid
         ! subroutine
         ! j is the rprop index
         iee = ppiclf_cell_map_interp(1,ie) 
         CALL ppiclf_copy(ppiclf_int_fld (j,ie)
     >                   ,ppiclf_int_fld_input(iee,j),1)
      END DO

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InterpTupleTransfer
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal: 
!
      REAL*8 FLD(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE),
     >       Max_CellLen(3)
      INTEGER*4 nkey(2), nl, nii, njj, nrr, ie, l 
      LOGICAL partl
!
      ! send it all
      nl   = 0
      nii  = PPICLF_LRMAX
      njj  = 3
      nrr  = PPICLF_LRP_INT
      nkey(1) = 2
      nkey(2) = 1

      CALL pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl ! Setup
     >      ,ppiclf_nCells_Interp, PPICLF_LEE ! Amount of columns to transfer
     >      ,ppiclf_cell_map_interp, nii      ! Integer communication
     >      ,partl, nl                        ! Logical communication
     >      ,ppiclf_int_fld, nrr              ! Real communication
     >      ,njj)                             ! Proc index to send to
      CALL pfgslib_crystal_tuple_sort(ppiclf_cr_hndl ! Setup
     >      ,ppiclf_nCells_Interp             ! Amount of columns to sort
     >      ,ppiclf_cell_map_interp,nii       ! Integer data
     >      ,partl,nl                         ! Logical data
     >      ,ppiclf_int_fld,nrr               ! Real data
     >      ,nkey,2)                          ! Sorting order

      ! Find distance check for interpolation.
      ! This is 1.5*MaxCellLength to ensure that at least
      ! 27 neighboring cells are mapped.
      Max_CellLen(1) = 0.0D0
      Max_CellLen(2) = 0.0D0
      Max_CellLen(3) = 0.0D0
      DO ie = 1,ppiclf_nCells_Interp ! Loop through cells mapped to bin
        DO l = 1,3
          ! Find max cell lengths in all dimensions
          IF(ppiclf_picl_grid(3+l,ie) .GT. Max_CellLen(l))
     >      Max_CellLen(l) = ppiclf_picl_grid(3+l,ie)
        END DO !l
      END DO !ie

      DO l = 1,3
        ! Multiply by 1.5 so particle near face will
        ! find center one cell over in farthest direction
        ppiclf_interp_dchk(l) = Max_CellLen(l)*1.5D0
      END DO

      RETURN
      END

!
!______________________________________________________________________
!
      SUBROUTINE ppiclf_solve_SBParticleToCellMap

      IMPLICIT NONE

      INCLUDE "PPICLF"

      ! Local Variables
      INTEGER*4 i, j, k, l, ix, iy, iz, ip, ie, iee, nxyz, nnearest, 
     >          CellID_nearest(28), partCount
      REAL*8    dSQl, dSQi, dSQ(28), xp(3),  
     >          CellCenter(3,28), w(27), binblength(3),  
     >          Max_CellLen(3), Max_CellLenSQ(3), dSQchk(3)
      LOGICAL   added, farAway, alreadyMapped
 
      INTEGER*4  SBin_map( 0 : (
     > (FLOOR((ppiclf_bins_dx(1) + 2*ppiclf_interp_dchk(1))
     >                       /ppiclf_interp_dchk(1)) + 1) *
     > (FLOOR((ppiclf_bins_dx(2) + 2*ppiclf_interp_dchk(2))
     >                       /ppiclf_interp_dchk(2)) + 1) *
     > (FLOOR((ppiclf_bins_dx(3) + 2*ppiclf_interp_dchk(3))
     >                       /ppiclf_interp_dchk(3)) + 1) - 1)
     > ,ppiclf_nCells_Interp)

      INTEGER*4  SBin_counter( 0 : (
     > (FLOOR((ppiclf_bins_dx(1) + 2*ppiclf_interp_dchk(1))
     >                       /ppiclf_interp_dchk(1)) + 1) *
     > (FLOOR((ppiclf_bins_dx(2) + 2*ppiclf_interp_dchk(2))
     >                       /ppiclf_interp_dchk(2)) + 1) *
     > (FLOOR((ppiclf_bins_dx(3) + 2*ppiclf_interp_dchk(3))
     >                       /ppiclf_interp_dchk(3)) + 1) - 1))
       INTEGER*4 n_SBin(3), tot_SBin, i_SBin(3), iTemp_SBin(3),
     >           temp_SBin, iSB, jSB, kSB, loopSB, i_count,
     >           firstSB(3), lastSB(3)  
      REAL*8    bin_Min(3), x_range(3), size_SBin(3)
      !***************************************************************

      IF(ppiclf_npart .LT. 1) RETURN
      IF(ppiclf_nCells_Interp .EQ. 0 . AND. ppiclf_npart .GT. 0) THEN
        PRINT*,'ERROR: ',ppiclf_npart, 'Particles mapped to bin:'
     >         ,ppiclf_nid
        PRINT*,'No cells mapped to bin for Interpolation/Projection.'
        CALL ppiclf_exittr('Failure in particle to cell mapping',0.D0,0)
      END IF
 
      DO l = 1,3
        binblength(l) = ppiclf_binb(2*l) - ppiclf_binb((2*l)-1)
        bin_Min(l) = ppiclf_bin_pos(1,l) - ppiclf_interp_dchk(l)
        n_SBin(l) = FLOOR((ppiclf_bins_dx(l) + 2*ppiclf_interp_dchk(l))
     >                                     / ppiclf_interp_dchk(l)) + 1 
        dSQchk(l) = (ppiclf_interp_dchk(l))**2
        ! SB at bin min boundary
        firstSB(l) = FLOOR((ppiclf_bin_pos(1,l) - bin_Min(l))
     >               / ppiclf_interp_dchk(l))
        ! SB at bin miax boundary
        lastSB(l)  =  FLOOR((ppiclf_bin_pos(2,l) - bin_Min(l))
     >               / ppiclf_interp_dchk(l))
      END DO
 
      SBin_Counter = 0
      tot_SBin = n_SBin(1)*n_SBin(2)*n_SBin(3)
     
      ! Loop through all elements to map to subbins.
      ! Particles don't need to be mapped, since the particle
      ! subbin is determined in following loop.
      DO ie = 1,ppiclf_nCells_Interp  
        DO l = 1,3
          i_SBin(l) = FLOOR((ppiclf_picl_grid(l,ie) - 
     >                bin_Min(l)) / ppiclf_interp_dchk(l))
        END DO
        ! In the i,j,k loops below, 0 takes care of non-periodic mapping
        ! and 1 takes care of periodic mapping.  If a cell is in corner,
        ! it'll be mapped to 2*2*2=8 subbins.
!*** need to adjust for angular periodicity ***
        DO i = 0,1
          IF(i .EQ. 0) THEN
            iTemp_SBin(1) = i_SBin(1)
            IF(iTemp_SBin(1) .LT. 0) 
     >        iTemp_SBin(1) = 0
            IF(iTemp_SBin(1) .GT. n_SBin(1) - 1) 
     >        iTemp_SBin(1) = n_SBin(1) - 1
          ELSE ! i .EQ. 1
            IF(ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)) THEN
              IF(i_SBin(1) .LE. firstSB(1) + 1) THEN
                iTemp_SBin(1) = lastSB(1) - 1
                IF(iTemp_SBin(1) .EQ. i_SBin(1)) CYCLE 
              ELSE IF(i_SBin(1) .GE. lastSB(1) - 1) THEN
                iTemp_SBin(1) = firstSB(1) + 1
                IF(iTemp_SBin(1) .EQ. i_SBin(1)) CYCLE
              ELSE
                CYCLE
              END IF
            ELSE 
              CYCLE
            END IF
          END IF
          DO j = 0,1
            IF(j .EQ. 0) THEN
              iTemp_SBin(2) = i_SBin(2)
              IF(iTemp_SBin(2) .LT. 0) 
     >          iTemp_SBin(2) = 0
               IF(iTemp_SBin(2) .GT. n_SBin(2) - 1) 
     >          iTemp_SBin(2) = n_SBin(2) - 1
            ELSE ! j .EQ. 1
              ! This takes care of periodicity for single processor
              IF(ppiclf_linperiodic(2) .AND. ppiclf_EqualDomain(2)) THEN
                IF(i_SBin(2) .LE. firstSB(2) + 1) THEN
                  iTemp_SBin(2) = lastSB(2) - 1
                  IF(iTemp_SBin(2) .EQ. i_SBin(2)) CYCLE
                ELSE IF(i_SBin(2) .GE. lastSB(2) - 1) THEN
                  iTemp_SBin(2) = firstSB(2) + 1
                  IF(iTemp_SBin(2) .EQ. i_SBin(2)) CYCLE
                ELSE
                  CYCLE
                END IF
              ELSE 
                CYCLE
              END IF
            END IF
            DO k = 0,1
              IF(k .EQ. 0) THEN
                iTemp_SBin(3) = i_SBin(3)
                IF(iTemp_SBin(3) .LT. 0) 
     >            iTemp_SBin(3) = 0
                IF(iTemp_SBin(3) .GT. n_SBin(3) - 1) 
     >            iTemp_SBin(3) = n_SBin(3) - 1
              ELSE ! k .EQ. 1
                ! This takes care of periodicity for single processor
                IF(ppiclf_linperiodic(3) .AND. 
     >                                      ppiclf_EqualDomain(3)) THEN 
                  IF(i_SBin(3) .LE. firstSB(3) + 1) THEN
                    iTemp_SBin(3) = lastSB(3) - 1
                    IF(iTemp_SBin(3) .EQ. i_SBin(3)) CYCLE
                  ELSE IF(i_SBin(3) .GE. lastSB(3) - 1) THEN
                    iTemp_SBin(3) = firstSB(3) + 1
                    IF(iTemp_SBin(3) .EQ. i_SBin(3)) CYCLE
                  ELSE
                    CYCLE
                  END IF
                ELSE 
                  CYCLE
                END IF
              END IF
              ! Finally, add the cell to a subbin 
              temp_SBin = iTemp_SBin(1) + iTemp_SBin(2)*n_SBin(1) +
     >                    iTemp_SBin(3)*n_SBin(1)*n_SBin(2)
              SBin_Counter(temp_SBin) = SBin_Counter(temp_SBin) + 1
              IF(SBin_Counter(temp_SBin) .GT. ppiclf_nCells_Interp)
     >          PRINT*, 'counter more than interp cells. SB:',
     >          temp_SBin, SBin_Counter(temp_SBin), ppiclf_nCells_Interp
              SBin_Map(temp_SBin,SBin_Counter(temp_SBin)) = ie
            END DO !k
          END DO !j 
        END DO !i
      END DO !ie
      partCount = 0
      DO ip=1,ppiclf_npart !Loop all particles in this bin
        nnearest = 0 ! number of nearest elements
        DO ie = 1,28
          CellID_nearest(ie) = -1 ! index of nearest elements
          dSQ(ie) = 1E20 ! distance to center of nearest element
        END DO !ie
        ! particle centers in all directions
        xp(1) = ppiclf_y(PPICLF_JX, ip)
        xp(2) = ppiclf_y(PPICLF_JY, ip)
        xp(3) = ppiclf_y(PPICLF_JZ, ip)
        DO l = 1,3
          i_SBin(l) = FLOOR((xp(l) - bin_Min(l)) 
     >                /ppiclf_interp_dchk(l))
        END DO
        temp_SBin = i_SBin(1) + i_SBin(2)*n_SBin(1) +
     >              i_SBin(3)*n_SBin(1)*n_SBin(2)
        DO iSB = 1,3      ! -1,+0,+1 subbin in x-dir
          DO jSB = 1,3    ! -1,+0,+1 subbin in y-dir
            DO kSB = 1,3  ! -1,+0,+1 subbin in z-dir
              loopSB = temp_SBin + (-2+iSB) + (-2+jSB)*n_SBin(1)
     >                 + (-2 + kSB)*n_SBin(1)*n_SBin(2)  
              IF(loopSB .GT. (-1) .AND. loopSB .LT. tot_SBin) THEN
                DO i_count = 1,SBin_Counter(loopSB)
                  ie = SBin_Map(loopSB,i_count) 
                  ! get distance from particle to center
                  dSQi    = 0.0
                  dSQl    = 0.0
                  farAway = .FALSE.
                  DO l=1,3
                    IF(ppiclf_linperiodic(l) .AND.
     >                                      ppiclf_EqualDomain(l)) THEN
                      dSQl = MIN((ppiclf_picl_grid(l,ie) - xp(l))**2, 
     >                       (binblength(l)-ABS(ppiclf_picl_grid(l,ie)
     >                        - xp(l)))**2)
                    ELSE
                      dSQl = (ppiclf_picl_grid(l,ie) - xp(l))**2
                    END IF
                    dSQi = dSQi + dSQl
                    IF (dSQl .GT. dSQchk(l)) farAway = .TRUE.
                  END DO !l
                  ! skip to next fluid cell if greater than 1.5*max cell
                  ! distance in respective x,y,z direction.
                  IF (farAWAY) CYCLE !i_count
                  ! Sort closest fluid cell centers
                  added = .FALSE.
                  alreadyMapped = .FALSE.
                  DO i = 1,27
                    IF(ie .EQ. CellID_nearest(i)) alreadyMapped = .TRUE.
                  END DO
                  DO i = 1,27
                    IF(alreadyMapped) EXIT ! go to next cell in SB
                    j = 27 - i + 1
                    IF (dSQi .LT. dSQ(j)) THEN
                      dSQ(j+1) = dSQ(j)
                      CellID_nearest(j+1) = CellID_nearest(j)
                      DO l=1,3
                        CellCenter(l, j+1) = CellCenter(l, j)
                      END DO
                      dSQ(j) = dSQi
                      CellID_nearest(j) = ie
                      DO l=1,3
                        CellCenter(l,j) = ppiclf_picl_grid(l,ie)
                      END DO
                      added = .TRUE.
                    ELSE ! If not within closest cell list
                      EXIT !i
                    END IF
                  END DO !i
                  IF (added) nnearest = nnearest + 1  
                END DO ! i_count
              END IF !SB out of domain
            END DO !kSB
          END DO !jSB
        END DO !iSB
        nnearest = MIN(nnearest,27)
        IF (nnearest .lt. 1) THEN
          ! Particle is outside of fluid domain.
          ! iprop(8,ip) set to -1 means it will be removed
          ! from ppiclf_y & ppiclf_rprop, rprop2, rprop3, rprop4, rprop5
          ppiclf_iprop(9,ip) = -1
          ppiclf_remove_particle = .TRUE.
          PRINT*, 'part # on proc # removed.', ppiclf_iprop(1,ip),
     >            ppiclf_nid
        ELSE
          partCount = partCount + 1
          ! use partCount since ip includes possible removed particles
          ppiclf_nPart2Cell(partCount) = nnearest
          DO i = 1,nnearest
            ppiclf_Part2Cell_map(partCount,i) = CellID_nearest(i) ! Cell ID
            ! Particle center to cell center distance
            ppiclf_Part2Cell_dist(partCount,i) = SQRT(dSQ(i)) 
          END DO !i
        END IF !nnearest
      END DO !ip

      IF(ppiclf_remove_particle) THEN
        ! Delete particles that are outside of fluid grid
        CALL ppiclf_solve_RemoveParticle
        ppiclf_remove_particle = .FALSE.
      END IF

      RETURN
      END
!
!______________________________________________________________________
!
      SUBROUTINE ppiclf_solve_ParticleToCellMap

      IMPLICIT NONE

      INCLUDE "PPICLF"

      ! Local Variables
      INTEGER*4 i, j, k, l, ix, iy, iz, ip, ie, iee, nxyz, nnearest, 
     >          CellID_nearest(28), partCount
      REAL*8    dSQl, dSQi, dSQ(28), xp(3), dSQchk(3), 
     >          CellCenter(3,28), w(27),binblength(3),  
     >          Max_CellLen(3),Max_CellLenSQ(3)
      LOGICAL   added, farAway

      !***************************************************************

      IF(ppiclf_nCells_Interp .EQ. 0 . AND. ppiclf_npart .GT. 0) THEN
        PRINT*,'No cells mapped to ppiclf bin. Num Particles/Proc ID:',
     >  ppiclf_npart, ppiclf_nid
        CALL ppiclf_exittr('Failure in particle to cell mapping',0.D0,0)
      END IF
 

      ! Find bin lengths for linear periodicity calculations
      DO l = 1,3
        binblength(l) = ppiclf_binb(2*l) - ppiclf_binb((2*l)-1)
        dSQchk(l) = (ppiclf_interp_dchk(l))**2
      END DO

      partCount = 0
      DO ip=1,ppiclf_npart !Loop all particles in this bin
        ! particle centers in all directions
        xp(1) = ppiclf_y(PPICLF_JX, ip)
        xp(2) = ppiclf_y(PPICLF_JY, ip)
        xp(3) = ppiclf_y(PPICLF_JZ, ip)
        nnearest = 0 ! number of nearest elements
        DO ie = 1,28
          CellID_nearest(ie) = -1 ! index of nearest elements
          dSQ(ie) = 1D20 ! distance to center of nearest element
        ENDDO !ie
        DO ie = 1,ppiclf_nCells_Interp
          ! get distance from particle to center
          dSQl     = 0.0D0
          dSQi     = 0.0D0
          farAway = .FALSE.
          DO l=1,3
            IF(ppiclf_linperiodic(l) .AND. ppiclf_EqualDomain(l)) THEN
              dSQl = MIN((ppiclf_picl_grid(l,ie) - xp(l))**2, 
     >           (binblength(l)-ABS(ppiclf_picl_grid(l,ie) - xp(l)))**2)
            ELSE
              dSQl = (ppiclf_picl_grid(l,ie) - xp(l))**2
            END IF
            dSQi = dSQi + dSQl
            IF (dSQl .GT. dSQchk(l)) farAway = .TRUE.
          END DO !l

          ! skip to next fluid cell if greater than 1.5*max cell
          ! distance in respective x,y,z direction.
          IF (farAWAY) CYCLE !ie
          ! Sort closest fluid cell centers
          ! *** Slow and should be updated.  
          ! No need to have closest 27 cells sorted.
          ! just need to exclude cells farther than 27.***
          added = .FALSE.
          DO i=1,27
            j = 27 - i + 1
            IF (dSQi .LT. dSQ(j)) THEN
              dSQ(j+1) = dSQ(j)
              CellID_nearest(j+1) = CellID_nearest(j)
              DO l=1,3
                CellCenter(l, j+1) = CellCenter(l, j)
              END DO
              dSQ(j) = dSQi
              CellID_nearest(j) = ie
              DO l=1,3
                CellCenter(l,j) = ppiclf_picl_grid(l,ie)
              END DO
              added = .TRUE.
            ELSE ! If not within closest cell list
              EXIT !i
            END IF
          END DO !i
          IF (added) nnearest = nnearest + 1
        END DO ! ie
        nnearest = MIN(nnearest, 27)
        IF (nnearest .LT. 1) THEN
          ! Particle is outside of fluid domain.
          ! iprop(8,ip) set to -1 means it will be removed
          ! from ppiclf_y & ppiclf_rprop, rprop2, rprop3, rprop4, rprop5
          ppiclf_iprop(9,ip) = -1
          ppiclf_remove_particle = .TRUE.
          PRINT*, 'part # on proc # removed with postion of:',
     >            ppiclf_iprop(1,ip),ppiclf_nid,xp(1),xp(2),xp(3) 
        ELSE
          partCount = partCount + 1
          ! use partCount since ip includes possible removed particles
          ppiclf_nPart2Cell(partCount) = nnearest
          IF(nnearest .GT. 27) PRINT*,'nn > 27!',partCount
          DO i = 1,nnearest
            ppiclf_Part2Cell_map(partCount,i) = CellID_nearest(i) ! Cell ID
            ! Particle center to cell center distance
            ppiclf_Part2Cell_dist(partCount,i) = SQRT(dSQ(i)) 
           END DO
        END IF 
      END DO !ip

      IF(ppiclf_remove_particle) THEN
        ! Delete particles that are outside of fluid grid
        CALL ppiclf_solve_RemoveParticle
        ppiclf_remove_particle = .FALSE.
      END IF

      RETURN
      END
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ppiclf_solve_Interpolate

      IMPLICIT NONE

      INCLUDE "PPICLF"

      ! Local Variables
      INTEGER*4 i, j, k, ip, nnearest,cellID 
      REAL*8    wsum, eps, dist, a(27), w(27)  

      IF(ppiclf_npart .LT. 1) RETURN

      eps = 1.0D-60 ! Machine epsilon to avoid dividing by zero
      DO ip = 1,ppiclf_npart
        nnearest = ppiclf_nPart2Cell(ip)
        w = 0.0D0
        wsum = 0.0D0
        DO k = 1,nnearest
          ! Interpolation Weighting: 1/(distance^3)
          dist = ppiclf_Part2Cell_dist(ip,k)**3 + eps
          w(k) = 1.0d0 / dist
          wsum = w(k) + wsum
        END DO ! k
        DO i = 1,PPICLF_INT_ICNT
          j = PPICLF_INT_MAP(i)
          ppiclf_rprop(j, ip) = 0.0D0
          ! Inverse Distance Interpolation
          DO k = 1,nnearest
            cellID = ppiclf_Part2Cell_map(ip,k) 
            a(k) = ppiclf_int_fld(i,cellID)
            ppiclf_rprop(j, ip) = ppiclf_rprop(j, ip) + w(k)*a(k)/wsum
          END DO ! k
          IF (isnan(ppiclf_rprop(j,ip))) THEN
            PRINT *, 'INTERP NAN: Particle, processor id, nnearest', ip,
     >                                   ppiclf_nid,nnearest
            PRINT*, 'Index:',j, 'Value:',ppiclf_rprop(j,ip)
            CALL ppiclf_exittr('rprop NaN in Interpolate',0.D0,0)
          END IF
        END DO ! i
      END DO ! ip

      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_RemoveParticle
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"

! Internal:
!
      INTEGER*4 i, icount
      icount = 0
      DO i=1,ppiclf_npart
         IF(ppiclf_iprop(9,i) .NE. -1) THEN
            ! Keep particle - copy the column with particle information
            icount = icount + 1 
            IF(i .NE. icount) THEN
               CALL ppiclf_copy
     >          (ppiclf_y     (1,icount),     ppiclf_y(1,i), PPICLF_LRS)
               CALL ppiclf_copy
     >          (ppiclf_y1    (1,icount),    ppiclf_y1(1,i), PPICLF_LRS)
               CALL ppiclf_copy
     >          (ppiclf_ydot  (1,icount),  ppiclf_ydot(1,i), PPICLF_LRS)
               CALL ppiclf_copy
     >          (ppiclf_ydotc (1,icount), ppiclf_ydotc(1,i), PPICLF_LRS)
               CALL ppiclf_copy
     >          (ppiclf_rprop (1,icount), ppiclf_rprop(1,i), PPICLF_LRP)
               IF(PPICLF_LRP2 .GT. 1) THEN
                 CALL ppiclf_copy
     >          (ppiclf_rprop2(1,icount),ppiclf_rprop2(1,i),PPICLF_LRP2)
               END IF
               IF(PPICLF_LRP3 .GT. 1) THEN
                 CALL ppiclf_copy
     >          (ppiclf_rprop3(1,icount),ppiclf_rprop3(1,i),PPICLF_LRP3)
               END IF
               IF(PPICLF_LRP4 .GT. 1) THEN
                 CALL ppiclf_copy
     >          (ppiclf_rprop4(1,icount),ppiclf_rprop4(1,i),PPICLF_LRP4)
               END IF
               IF(PPICLF_LRP5 .GT. 1) THEN
                 CALL ppiclf_copy
     >          (ppiclf_rprop5(1,icount),ppiclf_rprop5(1,i),PPICLF_LRP5)
               END IF
               CALL ppiclf_copy(ppiclf_feedbk(1,icount), 
     >                          ppiclf_feedbk(1,i), PPICLF_LRP_PRO)
               CALL ppiclf_icopy
     >          (ppiclf_iprop(1,icount) , ppiclf_iprop(1,i), PPICLF_LIP)
            END IF
         ! Else - don't copy particle column if marked for removal
         ! Particles marked for removal if outside fluid domain, which
         ! is found in the particle to cell mapping during interpolation
         END IF
      END DO

      ppiclf_npart = icount

      RETURN
      END

!----------------------------------------------------------------------

      SUBROUTINE ppiclf_solve_ProjectParticleGrid
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"

      ! Internal:
      INTEGER*4 i, j, ip, ie, nCellProj, CellID, nl, nii, njj,
     >          nrr, nkey(2), iee
      REAL*8    CellVol, GaussianConst, dist, w(27), wsum,
     >          x_norm, y_norm, z_norm, PI, eps
      LOGICAL   partl 
 
      PI = 4*ATAN(1.0D0)
      GaussianConst = 2.305D0 ! Distribution over 2 cell widths
      ppiclf_pro_fld_picl = 0.0d0
      eps = 1.0D-60
      DO ip=1,ppiclf_npart
        ! Update volume fraction for feedback - important for 1st RK
        ! step at time = 0.0
        ppiclf_feedbk(PPICLF_P_JPHIP,ip) =
     >  ppiclf_rprop(PPICLF_R_JVOLP,ip) * ppiclf_rprop(PPICLF_R_JSPL,ip)
        nCellProj = ppiclf_nPart2Cell(ip)
        wsum = 0.0D0
        ! Loop to find individual cell weightings
        DO i = 1,nCellProj
          CellID = ppiclf_Part2Cell_map(ip,i) 
          dist = ppiclf_Part2Cell_dist(ip,i) + eps
          CellVol = ppiclf_picl_grid(7,CellID)
          w(i) = ABS(CellVol*EXP(-GaussianConst*(dist**2)
     >              / (CellVol**(2.0D0/3.0D0))))
          wsum = wsum + w(i)
        END DO !i
#ifdef TEST
        ! These are same feedback equations used in unit testing
        x_norm = (ppiclf_y(PPICLF_JX,ip) - ppiclf_binb(1))
     >           / (ppiclf_binb(2) - ppiclf_binb(1))
        y_norm = (ppiclf_y(PPICLF_JY, ip) - ppiclf_binb(3))
     >           / (ppiclf_binb(4) - ppiclf_binb(3))
        z_norm = (ppiclf_y(PPICLF_JZ, ip) - ppiclf_binb(5))
     >           / (ppiclf_binb(6) - ppiclf_binb(5))

        ppiclf_feedbk(1,ip) = 1
        ppiclf_feedbk(2,ip) = SIN(2*PI*x_norm) + 
     >                        SIN(2*PI*y_norm) + SIN(2*PI*z_norm)
#endif     
        DO j=1,PPICLF_LRP_PRO
          ! Loop through cells to apply feedback     
          DO i = 1,nCellProj
            CellID = ppiclf_Part2Cell_map(ip,i)
            ppiclf_pro_fld_picl(j,CellID) = 
     >         ppiclf_pro_fld_picl(j,CellID) 
     >         + ppiclf_feedbk(j,ip)*w(i)/wsum
          END DO !i
        END DO !j
      END DO !ip

      ! Now send feedback information to processor that contains 
      ! the cell for the fluid solver

      ppiclf_nCells_Proj = ppiclf_nCells_Interp
      DO i = 1,ppiclf_nCells_Proj
        CALL ppiclf_icopy(ppiclf_cell_map_proj(1,i),
     >         ppiclf_cell_map_interp(1,i),PPICLF_LRMAX)
      END DO

      nl = 0
      nii = PPICLF_LRMAX
      njj = 2 ! original processor with cell for fluid grid
      nrr = PPICLF_LRP_PRO
      nkey(1) = 2
      nkey(2) = 1
      CALL pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl ! Setup
     >      ,ppiclf_nCells_Proj, PPICLF_LEE ! Amount of columns to transfer
     >      ,ppiclf_cell_map_proj, nii      ! Integer communication
     >      ,partl, nl                      ! Logical communication
     >      ,ppiclf_pro_fld_picl, nrr       ! Real communication
     >      ,njj)                           ! Proc index to send to
      CALL pfgslib_crystal_tuple_sort(ppiclf_cr_hndl ! Setup
     >      ,ppiclf_nCells_Proj             ! Amount of columns to sort
     >      ,ppiclf_cell_map_proj,nii       ! Integer data
     >      ,partl,nl                       ! Logical data
     >      ,ppiclf_pro_fld_picl,nrr        ! Real data
     >      ,nkey,2)                        ! Sorting order

      ppiclf_pro_fld = 0.0d0
      DO ie=1,ppiclf_nCells_Proj
         iee = ppiclf_cell_map_Proj(1,ie)
         DO j=1,PPICLF_LRP_PRO
           ! Mapped to the fluid solver domain
           ppiclf_pro_fld(iee,j) = ppiclf_pro_fld(iee,j) +
     >                                    ppiclf_pro_fld_picl(j,ie)
         END DO
      END DO

      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_GetProFld(e,m,fld)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4 e,m
! e - fluid_grid element number
! m - projection property index

!
! Output:
!
      REAL*8 fld
!
      fld = ppiclf_pro_fld(e,m)

      RETURN
      END

!-----------------------------------------------------------------------
!       SUBROUTINE ppiclf_solve_InvokeAngularPeriodic(i,flag,
!     >                                              per_alpha,
!     >                                              angle, xangle,
!     >                                              register)
!!
!      IMPLICIT NONE
!!
!      INCLUDE "PPICLF"
!! :
!! Input: 
!! 
!      ! Thierry -  07/24/24 - modified code begins here
!      ! global variables
!      INTEGER*4 i, flag
!      REAL*8 rin, rout, per_alpha, angle, xangle
!      ! local variables
!      REAL*8 ct, st, ex, ey, ez, local_angle
!      REAL*8 rotmat(3,3) , v(3), x(3)
!      INTEGER*4 register
!
!        ! Thierry - 07/24/24 - modified code begins here
!        ! Implementation of Rotational Periodicity
!        ! Just like how Rocflu does it in modflu/RFLU_ModRelatedPatches.F90
!        ! this is invoked when particle is leaving x-axis or y-axis
!       
!!        print*, "!!! Rotational Periodicity Invoked !!!!" 
!          
!        ! use local angle so the value of angle does not get affected globally
!        local_angle = angle
!        ! Thierry 
!        !    (1) sign convention for theta is +ve when measured CCW
!        !           switch angle sign when particle is leaving from 
!        !           upper periodic face
!        !    (2) 0.5 instead of 1.0 to switch angle for ghost algorithm
!        !           since the ghost is being created before the 
!        !           particle is leaving domain
!        if(per_alpha.gt. 0.5*(xangle+angle)) 
!     >    local_angle=-1.0*local_angle
!        
!        ! Half-cylinder case - particle leaving +ve x-axis 
!        !                    - adjust rotation matrix angle accordingly
!        if(ang_case .EQ. 3) then
!          if(per_alpha .lt. xangle) local_angle = 0.0 
!        END if
!
!        ! convert from degrees to radians
!        ct = cos(local_angle)
!        st = sin(local_angle)
!        
!        SELECT CASE(flag)
!          !CASE(1)
!          !  ex = 1.0d0
!          !  ey = 0.0d0
!          !  ez = 0.0d0
!          !  print*, "X-Rotational Axis"
!
!          !CASE(2)
!          !  ex = 0.0d0
!          !  ey = 1.0d0
!          !  ez = 0.0d0
!          !  print*, "Y-Rotational Axis"
!
!          CASE(1)
!            ex = 0.0d0
!            ey = 0.0d0
!            ez = 1.0d0
!!            print*, "Z-Rotational Axis"
!          CASE DEFAULT
!            CALL ppiclf_exittr('Invalid Axis of Rotation!$',0.0d0
!     >         ,ppiclf_nid)
!
!          END SELECT 
!          
!          ! Rotation Matrix calculation
!          rotmat(1,1) = ct + (1.0d0-ct)*ex*ex
!          rotmat(1,2) =      (1.0d0-ct)*ex*ey - st*ez
!          rotmat(1,3) =      (1.0d0-ct)*ex*ez + st*ey
!          
!          rotmat(2,1) =      (1.0d0-ct)*ey*ex + st*ez
!          rotmat(2,2) = ct + (1.0d0-ct)*ey*ey
!          rotmat(2,3) =      (1.0d0-ct)*ey*ez - st*ex
!          
!          rotmat(3,1) =      (1.0d0-ct)*ez*ex - st*ey
!          rotmat(3,2) =      (1.0d0-ct)*ez*ey + st*ex
!          rotmat(3,3) = ct + (1.0d0-ct)*ez*ez
!
!          ! Corrdinates modification
!          x(1) = ppiclf_y(PPICLF_JX,i)
!          x(2) = ppiclf_y(PPICLF_JY,i)
!          x(3) = ppiclf_y(PPICLF_JZ,i)
!          
!          xrot = MATMUL(rotmat, x)
!          
!          ! Velocity vector modification
!
!          v(1) = ppiclf_y(PPICLF_JVX,i)
!          v(2) = ppiclf_y(PPICLF_JVY,i)
!          v(3) = ppiclf_y(PPICLF_JVZ,i)
!
!          vrot = MATMUL(rotmat, v)
!          
!          ! 08/27/24 - Thierry - we add a register variable to 
!          !   choose if we want to register
!          !   the angularly modified variables 
!          ! register = 1 when called from RemoveParticle -> we want to modify the values
!          ! register = 0 when called from AngularCreateGhost -> we don't want to modify values
!          
!          ! register modified values
!          if (register==1) then
!            !print*, "Registering values!" 
!            ppiclf_y(PPICLF_JX,i) = xrot(1)
!            ppiclf_y(PPICLF_JY,i) = xrot(2)
!            ppiclf_y(PPICLF_JZ,i) = xrot(3)
!            
!            ppiclf_y(PPICLF_JVX,i) = vrot(1)
!            ppiclf_y(PPICLF_JVY,i) = vrot(2)
!            ppiclf_y(PPICLF_JVZ,i) = vrot(3)
!            
!          END if 
!       
!      RETURN
!      END
!-----------------------------------------------------------------------
!      SUBROUTINE ppiclf_solve_InitAngularPlane(i,rin, rout,
!     >                                         angle, xangle,
!     >                                         dist1, dist2)
!!
!      IMPLICIT NONE
!!
!      INCLUDE "PPICLF"
!! Inputs 
!!
!      INTEGER*4 i
!      REAL*8 rin, rout, angle, xangle
!! Local Variables:
!!     
!      REAL*8 p1(3), p2(3), p3(3), p4(3), p5(3), p6(3),
!     >       v1(3), v2(3), v3(3), v4(3), n1(3), n2(3),
!     >       A, B, C, D, E, F, G, H, zt, xp, yp, zp
!!
!! Outputs
!      REAL*8 dist1, dist2
!!
!      xp = ppiclf_y(PPICLF_JX,i)
!      yp = ppiclf_y(PPICLF_JY,i)
!      zp = ppiclf_y(PPICLF_JZ,i)
!      
!      zt = ppiclf_binb(6) - ppiclf_binb(5) ! bin thickness in z-direction
!
!      !!! upper plane calculation !!! 
!      
!      ! plane equation
!      ! Ax + By + Cz + D = 0
!
!      ! p1, p2, p3 are 3 points in the upper plane
!      p1 = (/rin, tan(angle - abs(xangle))*rin, 0.0d0/) 
!      p2 = (/rout, tan(angle - abs(xangle))*rout, 0.0d0/) 
!      p3 = (/rin, tan(angle - abs(xangle))*rin, zt/) 
!
!      v1 = p2 - p1 ! vector P1P2
!      v2 = p3 - p1 ! vector P1P3
!      
!      ! upper plane normal vector - n1(A,B,C) = v1 x v2
!      ! cross product calculation
!      A =  v1(2)*v2(3) - v1(3)*v2(2)
!      B = -v1(1)*v2(3) + v1(3)*v2(1)
!      C =  v1(1)*v2(2) - v1(2)*v2(1)
!      n1(1)=A ; n1(2)=B; n1(3)=C
!      
!      ! values of either p1, p2, or p3 can be used to calculate D
!      D = -A*p1(1) - B*p1(2) - C*p1(3)
!      
!      ! P(xp, yp, zp) arbitrary point
!      ! dist = distance between P and upper plane 
!      dist1 = abs(A*xp + B*yp + C*zp + D)
!      dist1 = dist1/sqrt(A**2 + B**2 + C**2)
!      
!      !!! lower plane calculation !!! 
!      ! plane equation
!      ! Ex + Fy + Gz + H = 0
!
!      ! p4, p5, p6 are 3 points in the lower plane
!      p4 = (/rin, -tan(angle - abs(xangle))*rin, 0.0d0/)
!      p5 = (/rout, -tan(angle - abs(xangle))*rout, 0.0d0/)
!      p6 = (/rin, -tan(angle - abs(xangle))*rin, zt/)
!      
!      v3 = p5 - p4 ! vector P4P5
!      v4 = p6 - p4 ! vector P4P6
!      
!      ! lower plane normal vector - n2(E,F,G) = v3 x v4
!      ! cross product calculation
!      E =  v3(2)*v4(3) - v3(3)*v4(2)
!      F = -v3(1)*v4(3) + v3(3)*v4(1)
!      G =  v3(1)*v4(2) - v3(2)*v4(1)
!      n2(1)=E ; n2(2)=F; n2(3)=G
!      
!      ! values of either p4, p5, or p6 can be used to calculate H
!      H = -E*p4(1) - F*p4(2) - G*p4(3)
!
!      ! P(xp, yp, zp) arbitrary point
!      ! dist = distance between P and lower plane 
!      dist2 = abs(E*xp + F*yp + G*zp + H)
!      dist2 = dist2/sqrt(E**2 + F**2 + G**2)
!
!      RETURN
!      END

!-----------------------------------------------------------------------
! David's old nearest neighbor method
!
!      SUBROUTINE ppiclf_solve_NearestNeighbor(i)
!!
!      IMPLICIT NONE
!!
!      INCLUDE "PPICLF"
!! 
!! Input:
!! 
!      INTEGER*4 i
!! 
!! Internal: 
!! 
!      REAL*8 ydum(PPICLF_LRS), rpropdum(PPICLF_LRP)
!      REAL*8 A(3),B(3),C(3),AB(3),AC(3), dist2, xdist2, ydist2,
!     >       dist_total
!      INTEGER*4 i_iim, i_iip, i_jjm, i_jjp, i_kkm, i_kkp, j, j_ii, j_jj,
!     >          j_kk, jp
!      REAL*8 rnx, rny, rnz, area, rpx1, rpy1, rpz1, rpx2, rpy2, rpz2,
!     >       rflip, a_sum, rd, rdist, theta, tri_area, rthresh,
!     >       ab_dot_ac, ab_mag, ac_mag, zdist2
!      INTEGER*4 istride, k, kmax, kp, kkp, kk
!! 
!      i_iim = ppiclf_nb_r(1,i) - 1
!      i_iip = ppiclf_nb_r(1,i) + 1
!      i_jjm = ppiclf_nb_r(2,i) - 1
!      i_jjp = ppiclf_nb_r(2,i) + 1
!      i_kkm = ppiclf_nb_r(3,i) - 1
!      i_kkp = ppiclf_nb_r(3,i) + 1
!
!      dist2 = ppiclf_nndist**2
!
!      do j=1,ppiclf_npart
!         if (j .EQ. i) cycle
!
!         j_ii = ppiclf_nb_r(1,j)
!         j_jj = ppiclf_nb_r(2,j)
!         j_kk = ppiclf_nb_r(3,j)
!
!         if (j_ii .gt. i_iip .or. j_ii .lt. i_iim) cycle
!         if (j_jj .gt. i_jjp .or. j_jj .lt. i_jjm) cycle
!         if (ppiclf_ndim .EQ. 3) then
!         if (j_kk .gt. i_kkp .or. j_kk .lt. i_kkm) cycle
!         ENDif
!
!         xdist2 = (ppiclf_cp_map(1,i)-ppiclf_cp_map(1,j))**2
!         if (xdist2 .gt. dist2) cycle
!         ydist2 = (ppiclf_cp_map(2,i)-ppiclf_cp_map(2,j))**2
!         if (ydist2 .gt. dist2) cycle
!         dist_total = xdist2 + ydist2
!         if (ppiclf_ndim .EQ. 3) then
!         zdist2 = (ppiclf_cp_map(3,i)-ppiclf_cp_map(3,j))**2
!         if (zdist2 .gt. dist2) cycle
!         dist_total = dist_total+zdist2
!         ENDif
!         if (dist_total .gt. dist2) cycle
!
!         CALL ppiclf_user_EvalNearestNeighbor(i,j,ppiclf_cp_map(1,i)
!     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
!     >                                 ,ppiclf_cp_map(1,j)
!     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,j))
!
!      ENDdo
!
!      do j=1,ppiclf_npart_gp
!         j_ii = ppiclf_nb_g(1,j)
!         j_jj = ppiclf_nb_g(2,j)
!         j_kk = ppiclf_nb_g(3,j)
!
!         if (j_ii .gt. i_iip .or. j_ii .lt. i_iim) cycle
!         if (j_jj .gt. i_jjp .or. j_jj .lt. i_jjm) cycle
!         if (ppiclf_ndim .EQ. 3) then
!         if (j_kk .gt. i_kkp .or. j_kk .lt. i_kkm) cycle
!         ENDif
!
!         xdist2 = (ppiclf_cp_map(1,i)-ppiclf_rprop_gp(1,j))**2
!         if (xdist2 .gt. dist2) cycle
!         ydist2 = (ppiclf_cp_map(2,i)-ppiclf_rprop_gp(2,j))**2
!         if (ydist2 .gt. dist2) cycle
!         dist_total = xdist2 + ydist2
!         if (ppiclf_ndim .EQ. 3) then
!         zdist2 = (ppiclf_cp_map(3,i)-ppiclf_rprop_gp(3,j))**2
!         if (zdist2 .gt. dist2) cycle
!         dist_total = dist_total+zdist2
!         ENDif
!         if (dist_total .gt. dist2) cycle
!
!         jp = -1*j
!         CALL ppiclf_user_EvalNearestNeighbor(i,jp,ppiclf_cp_map(1,i)
!     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
!     >                                 ,ppiclf_rprop_gp(1,j)
!     >                                 ,ppiclf_rprop_gp(1+PPICLF_LRS,j))
!
!      ENDdo
!
!      istride = ppiclf_ndim
!      do j=1,ppiclf_nwall
!
!         rnx  = ppiclf_wall_n(1,j)
!         rny  = ppiclf_wall_n(2,j)
!         rnz  = 0.0d0
!         area = ppiclf_wall_n(3,j)
!         rpx1 = ppiclf_cp_map(1,i)
!         rpy1 = ppiclf_cp_map(2,i)
!         rpz1 = 0.0d0
!         rpx2 = ppiclf_wall_c(1,j)
!         rpy2 = ppiclf_wall_c(2,j)
!         rpz2 = 0.0d0
!         rpx2 = rpx2 - rpx1
!         rpy2 = rpy2 - rpy1
!
!         if (ppiclf_ndim .EQ. 3) then
!            rnz  = ppiclf_wall_n(3,j)
!            area = ppiclf_wall_n(4,j)
!            rpz1 = ppiclf_cp_map(3,i)
!            rpz2 = ppiclf_wall_c(3,j)
!            rpz2 = rpz2 - rpz1
!         ENDif
!    
!         rflip = rnx*rpx2 + rny*rpy2 + rnz*rpz2
!         if (rflip .gt. 0.0d0) then
!            rnx = -1.0d0*rnx
!            rny = -1.0d0*rny
!            rnz = -1.0d0*rnz
!         ENDif
!
!
!         a_sum = 0.0d0
!         kmax = 2
!         if (ppiclf_ndim .EQ. 3) kmax = 3
!         do k=1,kmax 
!            kp = k+1
!            if (kp .gt. kmax) kp = kp-kmax ! cycle
!            
!            kk   = istride*(k-1)
!            kkp  = istride*(kp-1)
!            rpx1 = ppiclf_wall_c(kk+1,j)
!            rpy1 = ppiclf_wall_c(kk+2,j)
!            rpz1 = 0.0d0
!            rpx2 = ppiclf_wall_c(kkp+1,j)
!            rpy2 = ppiclf_wall_c(kkp+2,j)
!            rpz2 = 0.0d0
!
!            if (ppiclf_ndim .EQ. 3) then
!               rpz1 = ppiclf_wall_c(kk+3,j)
!               rpz2 = ppiclf_wall_c(kkp+3,j)
!            ENDif
!
!            rd   = -(rnx*rpx1 + rny*rpy1 + rnz*rpz1)
!
!            rdist = abs(rnx*ppiclf_cp_map(1,i)+rny*ppiclf_cp_map(2,i)
!     >                 +rnz*ppiclf_cp_map(3,i)+rd)
!            rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)
!
!            ! give a little extra room for walls (2x)
!            if (rdist .gt. 2.0d0*ppiclf_nndist) goto 1511
!
!            ydum(1) = ppiclf_cp_map(1,i) - rdist*rnx
!            ydum(2) = ppiclf_cp_map(2,i) - rdist*rny
!            ydum(3) = 0.0d0
!
!            A(1) = ydum(1)
!            A(2) = ydum(2)
!            A(3) = 0.0d0
!
!            B(1) = rpx1
!            B(2) = rpy1
!            B(3) = 0.0d0
!
!            C(1) = rpx2
!            C(2) = rpy2
!            C(3) = 0.0d0
!
!            AB(1) = B(1) - A(1)
!            AB(2) = B(2) - A(2)
!            AB(3) = 0.0d0
!
!            AC(1) = C(1) - A(1)
!            AC(2) = C(2) - A(2)
!            AC(3) = 0.0d0
!
!            if (ppiclf_ndim .EQ. 3) then
!               ydum(3) = ppiclf_cp_map(3,i) - rdist*rnz
!               A(3) = ydum(3)
!               B(3) = rpz1
!               C(3) = rpz2
!               AB(3) = B(3) - A(3)
!               AC(3) = C(3) - A(3)
!
!               AB_DOT_AC = AB(1)*AC(1) + AB(2)*AC(2) + AB(3)*AC(3)
!               AB_MAG = sqrt(AB(1)**2 + AB(2)**2 + AB(3)**2)
!               AC_MAG = sqrt(AC(1)**2 + AC(2)**2 + AC(3)**2)
!               theta  = acos(AB_DOT_AC/(AB_MAG*AC_MAG))
!               tri_area = 0.5d0*AB_MAG*AC_MAG*sin(theta)
!            elseif (ppiclf_ndim .EQ. 2) then
!               AB_MAG = sqrt(AB(1)**2 + AB(2)**2)
!               tri_area = AB_MAG
!            ENDif
!            a_sum = a_sum + tri_area
!         ENDdo
!
!         rthresh = 1.10d0 ! keep it from slipping through crack on edges
!         if (a_sum .gt. rthresh*area) cycle
!
!         jp = 0
!         CALL ppiclf_user_EvalNearestNeighbor(i,jp,ppiclf_cp_map(1,i)
!     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
!     >                                 ,ydum
!     >                                 ,rpropdum)
!
! 1511 continue
!      ENDdo
!
!      RETURN
!      END
!!-----------------------------------------------------------------------
!
!
!!
!!
!! Maybe this one can be deleted??? - AVERY ***
!!
!!
!      SUBROUTINE ppiclf_solve_InitAngularPeriodic(flag,
!     >              rin, rout, angle, xangle)
!!
!      IMPLICIT NONE
!!
!      INCLUDE "PPICLF"
!! 
!! Input: 
!! 
!      ! Thierry -  07/24/24 - modified code begings here
!      ! global variables - user input file
!      INTEGER*4 flag
!      REAL*8 rin, rout, angle, xangle
!      ! local variables
!      REAL*8 pi, angled
!
!        ! Thierry - 07/24/24 - modified code begins here
!        ! Implementation of Angular Periodicity
!        ! Just like how Rocflu does it in modflu/RFLU_ModRelatedPatches.F90
!        ! this is invoked when particle is leaving x-axis or y-axis
!       
!        ! sign convention for theta is +ve when measured CCW
!        ! switch angle sign when particle is leaving from upper face
!        if (rin .ge. rout)
!     >   CALL ppiclf_exittr('Angular Per must have rin < rout$',rout,0)
!
!            ppiclf_iperiodic(1) = 0 ! X-periodic
!            ppiclf_iperiodic(2) = 0 ! Y-periodic
!
!            SELECT CASE (ang_case)
!              CASE (1) ! general wedge ; 0 <= angle < 90
!                if (ppiclf_nid.EQ.0) print*,"General Wedge Initialized!"
!                ppiclf_xdrange(1,1) = rin  ! Min X-periodic face
!                ppiclf_xdrange(2,1) = rout ! Max X-periodic face
!                ppiclf_xdrange(1,2) = tan(xangle)*rout ! Min Y-periodic face
!                ppiclf_xdrange(2,2) = tan(angle - abs(xangle))*rout ! Max Y-periodic face
!
!              CASE (2) ! quarter cylinder ; angle = 90
!                if (ppiclf_nid.EQ.0)
!     >             print*,"Quarter Cylinder Initialized!"
!                ppiclf_xdrange(1,1) = rin  
!                ppiclf_xdrange(2,1) = rout 
!                ppiclf_xdrange(1,2) = tan(xangle)*rout
!                ppiclf_xdrange(2,2) = rout 
!              
!              CASE (3) ! half cylinder ; angle = 180
!                if (ppiclf_nid.EQ.0)
!     >             print*,"Half Cylinder Initialized!"
!                ppiclf_xdrange(1,1) = -1.0*rout
!                ppiclf_xdrange(2,1) = rout 
!                ppiclf_xdrange(1,2) = tan(xangle)*rout
!                ppiclf_xdrange(2,2) = rout 
!              
!              CASE DEFAULT
!                CALL ppiclf_exittr('Invalid Rotational Case!$',0.0d0
!     >             ,ppiclf_nid)
!              END SELECT
!
!      
!      RETURN
!      END
!!-----------------------------------------------------------------------
!      SUBROUTINE ppiclf_solve_IntegrateRK3(iout)
!!
!      IMPLICIT NONE
!!
!      INCLUDE "PPICLF"
!! 
!! Internal: 
!! 
!      INTEGER*4 i, nstage, istage
!!
!! Output:
!!
!      LOGICAL iout
!!
!      ! save stage 1 solution
!      DO i = 1, ppiclf_npart
!        DO j = 1,PPICLF_LRS
!          ppiclf_y1(j,i) = ppiclf_y(j,i)
!        END DO
!      END DO
!
!      ! get rk3 coeffs
!      CALL ppiclf_solve_SetRK3Coeff(ppiclf_dt)
!
!      nstage = 3
!      DO istage=1,nstage
!         ! evaluate ydot
!         CALL ppiclf_solve_SetYdot
!
!         ! rk3 integrate
!         DO i= 1,npart
!           DO j = 1,PPICLF_LRS
!            ppiclf_y(j,i) =  ppiclf_rk3coef(1,istage)*ppiclf_y1   (j,i)
!     >                     + ppiclf_rk3coef(2,istage)*ppiclf_y    (j,i)
!     >                     + ppiclf_rk3coef(3,istage)*ppiclf_ydot (j,i)
!         END DO
!      END DO
!
!      iout = .TRUE.
!
!      RETURN
!      END
!c----------------------------------------------------------------------
!      SUBROUTINE ppiclf_solve_IntegrateRK3s(iout)
!!
!      IMPLICIT NONE
!!
!      INCLUDE "PPICLF"
!! 
!! Internal: 
!! 
!      INTEGER*4 i, ndum, nstage, istage
!      INTEGER*4 icalld
!      save      icalld
!      data      icalld /0/
!!
!! Output:
!!
!      LOGICAL iout
!!
!      icalld = icalld + 1
!
!
!      ! get rk3 coeffs
!      CALL ppiclf_solve_SetRK3Coeff(ppiclf_dt)
!
!      nstage = 3
!      istage = mod(icalld,nstage)
!      IF(istage .EQ. 0) istage = 3
!      iout = .FALSE.
!      IF(istage .EQ. nstage) iout = .TRUE.
!
!      ! save stage 1 solution
!      IF(istage .EQ. 1) THEN
!        DO i = 1,ppiclf_npart
!          DO j = 1,PPICLF_LRS
!           ppiclf_y1(j,i) = ppiclf_y(j,i)
!          END DO
!        END DO
!      END IF
!
!      ! evaluate ydot
!      CALL ppiclf_solve_SetYdot
!
!      ! rk3 integrate
!      DO i = 1,npart
!        DO j = 1,PPICLF_LRS
!         ppiclf_y(j,i) =  ppiclf_rk3coef(1,istage)*ppiclf_y1   (j,i)
!     >                  + ppiclf_rk3coef(2,istage)*ppiclf_y    (j,i)
!     >                  + ppiclf_rk3coef(3,istage)*ppiclf_ydot (j,i)
!        END DO
!      END DO
!
!      RETURN
!      END
