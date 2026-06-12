      SUBROUTINE ppiclf_solve_Initialize(PP, xi1,xpmin,xpmax,
     >           yi1,ypmin,ypmax,zi1,zpmin,zpmax,
     >           ai1,apa,apxa,aprin,aprout)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      LOGICAL   PP
      INTEGER*4 xi1, yi1, zi1, ai1
      REAL*8    xpmin,xpmax,ypmin,ypmax,zpmin,zpmax,
     >          apa,apxa,aprin,aprout, pi, angled

      ! Called by rocpicl/PICL_TEMP_InitSolver.F90
      ! xdrange adjusts the bin boundaries to ensure they aren't 
      ! larger than the cartesian fluid domain extremes.
      
      ! Establishes if particle-particle interactions are included.
      ! Includes collisions, Added mass, pseudo-turbulence, and qs
      ! fluctuations
      PPICLF_PPInteractions = PP

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

!*** The start for angular periodicity below.
!    The subbin NN search will make angular periodicity a little more 
!    complicated to implement...
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
      
      PPICLF_PPInteractions = .FALSE.

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
      ppiclf_printbinvtu  = .FALSE.
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
         ! 9 is set to -1 when particle should be removed
         ppiclf_iprop(9,i) = 0
      END DO

      RETURN
      END

!
!-----------------------------------------------------------------------
!
      SUBROUTINE ppiclf_solve_NearestNeighborSB(ip)
!
      USE ppiclf_DynamicAllocation

      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Internal: 
! 
      REAL*8    ydum(PPICLF_LRS), rpropdum(PPICLF_LRP), xp(3), 
     >          A(3),B(3),C(3),AB(3),AC(3), distSQ, xdistSQ, ydistSQ,
     >          dist_total, rnx, rny, rnz, area, rpx1, rpy1, rpz1, rpx2,
     >          rpy2, rpz2, rflip, a_sum, rd, rdist, theta, tri_area,
     >          ab_dot_ac, ab_mag, ac_mag, zdistSQ,rthresh,bin_xMin(3),
     >          cross_x, cross_y, cross_z, unx, uny, unz
      INTEGER*4 i,k, kmax, kp, kkp, kk, j, jp, l, iSB, jSB, kSB,
     >          loopSB, tempSB, iSBin(3), istride, ip, ii, jj
      INTEGER*4 n_SBin1x2, i_SBin(3)

!
!
      distSQ = ppiclf_nndist**2

      n_SBin1x2 = ppiclf_nSBin(1)*ppiclf_nSBin(2)
      ! particle centers in all directions
      xp(1) = ppiclf_y(PPICLF_JX, ip)
      xp(2) = ppiclf_y(PPICLF_JY, ip)
      xp(3) = ppiclf_y(PPICLF_JZ, ip)
      i_SBin(1) = ppiclf_iprop(5,ip) - ppiclf_binOffset(1)
      i_SBin(2) = ppiclf_iprop(6,ip) - ppiclf_binOffset(2)
      i_SBin(3) = ppiclf_iprop(7,ip) - ppiclf_binOffset(3)  
      tempSB = i_SBin(1) + i_SBin(2)*ppiclf_nSBin(1) 
     >         + i_SBin(3)*n_SBin1x2
#ifdef TEST
      PARTICLE_NN(ip) = 0 
      PPICLF_TOTNNDIST(ip) = 0.0D0
#endif
      ! Loop through real particles
      DO iSB = -1,1     !to look at -1,current,+1 x-dir subbins
        ii = i_SBin(1) + iSB
        IF(ii .LT. 0 .OR. ii .GT. ppiclf_nSBin(1)-1) CYCLE
        DO jSB = -1,1   !to look at -1,current,+1 x-dir subbins
          jj = i_SBin(2) + jSB
          IF(jj .LT. 0 .OR. jj .GT. ppiclf_nSBin(2)-1) CYCLE
          DO kSB = -1,1 !to look at -1,current,+1 x-dir subbins
            kk = i_SBin(3) + kSB
            IF(kk .LT. 0 .OR. kk .GT. ppiclf_nSBin(3)-1) CYCLE
            ! Loops through 27 adjacent subbins
            loopSB = ii + jj*ppiclf_nSBin(1) + kk*n_SBin1x2 
            DO k = 1,ppiclf_binPartCount(loopSB) 
              j = ppiclf_binPartList(loopSB,k)
              IF (j .GT. 0) THEN ! Real particle
                ! Cycle when same particle         
                IF(ppiclf_iprop(1,j) .EQ. ppiclf_iprop(1,ip) .AND.
     >             ppiclf_iprop(2,j) .EQ. ppiclf_iprop(2,ip) .AND.
     >             ppiclf_iprop(3,j) .EQ. ppiclf_iprop(3,ip)) THEN
                  CYCLE
                END IF
 
                xdistSQ = (ppiclf_y(1,ip)-ppiclf_y(1,j))**2
                IF (xdistSQ .GE. distSQ) CYCLE
                ydistSQ = (ppiclf_y(2,ip)-ppiclf_y(2,j))**2
                IF (ydistSQ .GE. distSQ) CYCLE
                dist_total = xdistSQ + ydistSQ
                zdistSQ = (ppiclf_y(3,ip)-ppiclf_y(3,j))**2
                IF (zdistSQ .GE. distSQ) CYCLE
                dist_total = dist_total+zdistSQ
                IF (dist_total .GE. distSQ) CYCLE

#ifdef TEST
                PARTICLE_NN(ip) = PARTICLE_NN(ip) + 1
                PPICLF_TOTNNDIST(ip) = PPICLF_TOTNNDIST(ip)+dist_total
                CYCLE !Don't want to call EvalNN. Just testing
                      ! nneighbor search
#endif
                CALL ppiclf_user_EvalNearestNeighbor(ip,j
     >                                 ,ppiclf_y(1,ip)
     >                                 ,ppiclf_rprop(1,ip)
     >                                 ,ppiclf_y(1,j)
     >                                 ,ppiclf_rprop(1,j))
              ELSE IF (j .LT. 0) THEN ! Ghost Particle
                ! Negative was just use for ghost particle indicator
                ! in subbin mapping array. Need to flip sign
                j = - j        
                IF(ppiclf_iprop_gp(1,j) .EQ. ppiclf_iprop(1,ip) .AND.
     >             ppiclf_iprop_gp(2,j) .EQ. ppiclf_iprop(2,ip) .AND.
     >             ppiclf_iprop_gp(3,j) .EQ. ppiclf_iprop(3,ip)) THEN
                  CYCLE
                END IF

                xdistSQ =(ppiclf_y(1,ip)-ppiclf_rprop_gp(1,j))**2
                IF (xdistSQ .GE. distSQ) CYCLE
                ydistSQ =(ppiclf_y(2,ip)-ppiclf_rprop_gp(2,j))**2
                IF (ydistSQ .GE. distSQ) CYCLE
                dist_total = xdistSQ + ydistSQ
                zdistSQ =(ppiclf_y(3,ip)-ppiclf_rprop_gp(3,j))**2
                IF (zdistSQ .GE. distSQ) CYCLE
                dist_total = dist_total+zdistSQ
                IF (dist_total .GE. distSQ) CYCLE

#ifdef TEST
                PARTICLE_NN(ip) = PARTICLE_NN(ip) + 1
                PPICLF_TOTNNDIST(ip) = PPICLF_TOTNNDIST(ip)+dist_total
                CYCLE !Don't want to call EvalNN. Just testing
                      ! nneighbor search
#endif
                jp = j
                CALL ppiclf_user_EvalNearestNeighbor(ip,jp
     >                             ,ppiclf_y(1,ip)
     >                             ,ppiclf_rprop(1,ip)
     >                             ,ppiclf_rprop_gp(1:PPICLF_LRS,j)
     >      ,ppiclf_rprop_gp(1+PPICLF_LRS:PPICLF_LRP+PPICLF_LRS,j))
              END IF
            END DO !k
          END DO !kSB
        END DO !jSB
      END DO !iSB
      
      ! --- PARTICLE-WALL INTERACTION LOOP ---
      istride = ppiclf_ndim
      
      wall_loop: DO j=1,ppiclf_nwall
         rnx  = ppiclf_wall_n(1,j)
         rny  = ppiclf_wall_n(2,j)
         rnz  = ppiclf_wall_n(3,j)
         area = ppiclf_wall_n(4,j)

         rpx1 = ppiclf_y(1,ip)
         rpy1 = ppiclf_y(2,ip)
         rpz1 = ppiclf_y(3,ip)

         rpx2 = ppiclf_wall_c(1,j) - rpx1
         rpy2 = ppiclf_wall_c(2,j) - rpy1
         rpz2 = ppiclf_wall_c(3,j) - rpz1
    
         rflip = rnx*rpx2 + rny*rpy2 + rnz*rpz2
         IF(rflip .GT. 0.0d0) THEN
            rnx = -1.0d0*rnx
            rny = -1.0d0*rny
            rnz = -1.0d0*rnz
         END IF
 
         unx = rnx
         uny = rny
         unz = rnz
 
         ! Equation of the plane using the first vertex of the wall
         rpx1 = ppiclf_wall_c(1,j)
         rpy1 = ppiclf_wall_c(2,j)
         rpz1 = ppiclf_wall_c(3,j)
         rd   = -(rnx*rpx1 + rny*rpy1 + rnz*rpz1)

         ! Orthogonal distance from particle to the infinite plane
         rdist = abs(rnx*ppiclf_y(1,ip) + rny*ppiclf_y(2,ip)
     >               + rnz*ppiclf_y(3,ip) + rd)
         
         ! If the particle is far from the plane, 
         ! it cannot interact with ANY facet. Skip wall entirely.
         IF(rdist .GT. 2.0d0*ppiclf_nndist) CYCLE wall_loop

         ! Mathematically correct projection using unit normal vectors
         ydum(1) = ppiclf_y(1,ip) - rdist*unx
         ydum(2) = ppiclf_y(2,ip) - rdist*uny
         ydum(3) = ppiclf_y(3,ip) - rdist*unz
         
         A(1:3) = ydum(1:3)
         
         a_sum = 0.0d0
         kmax = 3
         
         facet_loop: DO k=1,kmax 
            kp = k+1
            IF(kp .GT. kmax) kp = kp-kmax 
            
            kk   = istride*(k-1)
            kkp  = istride*(kp-1)
            
            B(1) = ppiclf_wall_c(kk+1,j)
            B(2) = ppiclf_wall_c(kk+2,j)
            B(3) = ppiclf_wall_c(kk+3,j)
            
            C(1) = ppiclf_wall_c(kkp+1,j)
            C(2) = ppiclf_wall_c(kkp+2,j)
            C(3) = ppiclf_wall_c(kkp+3,j)
            
            AB(1:3) = B(1:3) - A(1:3)
            AC(1:3) = C(1:3) - A(1:3)

            ! Cross Product for sub-triangle area
            cross_x = AB(2)*AC(3) - AB(3)*AC(2)
            cross_y = AB(3)*AC(1) - AB(1)*AC(3)
            cross_z = AB(1)*AC(2) - AB(2)*AC(1)

            tri_area = 0.5d0 * sqrt(cross_x**2 
     >                              + cross_y**2 + cross_z**2)
            a_sum = a_sum + tri_area
         END DO facet_loop
         
         rthresh = 1.10d0 
         IF(a_sum .GT. rthresh*area) CYCLE wall_loop

         jp = 0
         CALL ppiclf_user_EvalNearestNeighbor(ip,jp,ppiclf_cp_map(1,ip)
     >                                 ,ppiclf_rprop(1,ip)
     >                                 ,ydum
     >                                 ,rpropdum)
      END DO wall_loop


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
      REAL*8 rmag
      REAL*8 AB(3), AC(3), cross_product(3)
!
      IF(.NOT.PPICLF_LCOMM)
     >CALL ppiclf_exittr('InitMPI must be before InitWall$',0.d0,0)
      IF(.NOT.PPICLF_LINIT)
     >CALL ppiclf_exittr('InitParticle must be before InitWall$'
     >                  ,0.d0,0)

      ppiclf_nwall = ppiclf_nwall + 1 
      IF(ppiclf_nwall .gt. PPICLF_LWALL)
     >CALL ppiclf_exittr('Increase LWALL in user file$'
     >                  ,0.d0,ppiclf_nwall)
      IF(ppiclf_ndim .LT. 3) THEN
        PRINT*, 'ERROR: Expected 3D Input for Boundary Wall Geometry'
        CALL ppiclf_exittr('',0.0D0,0)
      END IF
      
      ! Store vertices
      ppiclf_wall_c(1,ppiclf_nwall) = xp1(1)
      ppiclf_wall_c(2,ppiclf_nwall) = xp1(2)
      ppiclf_wall_c(3,ppiclf_nwall) = xp1(3)
      ppiclf_wall_c(4,ppiclf_nwall) = xp2(1)
      ppiclf_wall_c(5,ppiclf_nwall) = xp2(2)
      ppiclf_wall_c(6,ppiclf_nwall) = xp2(3)
      ppiclf_wall_c(7,ppiclf_nwall) = xp3(1)
      ppiclf_wall_c(8,ppiclf_nwall) = xp3(2)
      ppiclf_wall_c(9,ppiclf_nwall) = xp3(3)

      AB(1) = xp2(1) - xp1(1)
      AB(2) = xp2(2) - xp1(2)
      AB(3) = xp2(3) - xp1(3)
      
      AC(1) = xp3(1) - xp1(1)
      AC(2) = xp3(2) - xp1(2)
      AC(3) = xp3(3) - xp1(3)
      
      ! Cross Product (AB x AC) gives the normal vector
      ! scaled by 2*Area
      cross_product(1) = AB(2)*AC(3) - AB(3)*AC(2)
      cross_product(2) = AB(3)*AC(1) - AB(1)*AC(3)
      cross_product(3) = AB(1)*AC(2) - AB(2)*AC(1)
      
      ! Magnitude of the cross product
      rmag = SQRT(cross_product(1)**2 + cross_product(2)**2
     >            + cross_product(3)**2)
      
      ! Store the Area (1/2 the magnitude of the cross product)
      ! Note: ppiclf_ndim must be 3, so ndim+1 is 4
      ppiclf_wall_n(4, ppiclf_nwall) = 0.5D0 * rmag 

      IF (rmag > 1.0D-12) THEN
         ppiclf_wall_n(1,ppiclf_nwall) = cross_product(1) / rmag
         ppiclf_wall_n(2,ppiclf_nwall) = cross_product(2) / rmag
         ppiclf_wall_n(3,ppiclf_nwall) = cross_product(3) / rmag
      ELSE
         ! Safety catch for degenerate/zero-area triangles
         ppiclf_wall_n(1:3,ppiclf_nwall) = 0.0D0
      END IF

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
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
      INCLUDE "mpif.h"
! 
! Input: 
! 
      REAL*8    time
#ifdef PERF
      REAL *8 tstart,tfinal     
#endif
! 
! Internal:
!
#ifdef PERF
      tstart = MPI_WTIME()
#endif
      ppiclf_time   = time

      CALL ppiclf_io_WriteParticleVTU('')
      !CALL ppiclf_io_WriteBinVTU('')
      ! Output diagnostics
      CALL ppiclf_io_OutputDiagAll

#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TWriteSolution = tfinal - tstart
#endif

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
      INCLUDE "mpif.h"
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

#ifdef PERF
      REAL*8    tstart, tfinal
      REAL*8    tsPeriodic, tfPeriodic

      ! Reset timers at start of each RK stage
      PPICLF_TCreateBin              = 0.0D0
      PPICLF_TSendParticles          = 0.0D0
      PPICLF_TSendGridOverlap        = 0.0D0               
      PPICLF_TSendFluidFields        = 0.0D0             
      PPICLF_TParticleParticleModels = 0.0D0       
      PPICLF_TFluidParticleModels    = 0.0D0     
      PPICLF_TSendGhostParticles     = 0.0D0     
      PPICLF_TMapParticlesCells      = 0.0D0    
      PPICLF_TInterpolation          = 0.0D0    
      PPICLF_TProjection             = 0.0D0     
      PPICLF_TWriteSolution          = 0.0D0                  
      PPICLF_TIntegration            = 0.0D0    
      PPICLF_TPeriodicity            = 0.0D0       
      PPICLF_TDataTransfers          = 0.0D0        
      PPICLF_TTotal                  = 0.0D0

      tstart = MPI_WTIME()
#endif

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
#ifdef PERF
      tsPeriodic = MPI_WTIME()
#endif
      IF(ppiclf_linperiodic(1) .OR. ppiclf_linperiodic(2) .OR.
     >                             ppiclf_linperiodic(3)) THEN
        CALL ppiclf_solve_PeriodicParticleShift
      END IF
#ifdef PERF
      tfPeriodic = MPI_WTIME()
      PPICLF_TPeriodicity = tfPeriodic - tsPeriodic
#endif
      CALL ppiclf_solve_PostTimeStepPartLB
#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TTotal = tfinal - tstart
      CALL ppiclf_io_WritePerformance()
#endif
      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_IntegrateRK3s_Rocflu(iout)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
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
#ifdef PERF
      REAL *8 tstart,tfinal     
#endif
!
#ifdef PERF
      tstart = MPI_WTIME()
#endif
      icalld = icalld + 1

      ! get rk3 coeffs
      CALL ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      istage = MOD(icalld,nstage)
      if (istage .EQ. 0) istage = 3
      iout = .FALSE.
      if (istage .EQ. nstage) iout = .TRUE.

#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TIntegration = tfinal - tstart
#endif
      ! evaluate ydot
      CALL ppiclf_solve_SetYdot

#ifdef PERF
      tstart = MPI_WTIME()
#endif

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
      
      !Store Current stage RHS for next stage's use
      IF(istage .NE. 3) THEN
        DO i= 1,PPICLF_NPART
          DO j= 1,PPICLF_LRS
            ppiclf_y1(j,i) =  ppiclf_ydot(j,i)
          END DO
        END DO
      END IF

#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TIntegration = PPICLF_TIntegration + (tfinal - tstart)
#endif

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_PeriodicParticleShift
!
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
!
! Internal:
!
      INTEGER*4 i  
      REAL*8 per_alpha
#ifdef PERF
      REAL *8 tstart,tfinal     
#endif
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
      IF( ( (x_per_flag.EQ.1) .OR. (y_per_flag.EQ.1)
     >                        .OR. (z_per_flag.EQ.1) ) ) THEN
!     >   .AND.(ang_per_flag.EQ.0)) THEN

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
      END IF

      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_SetYdot
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
!
      INTEGER*4 j, ierr
#ifdef PERF
      REAL *8 tstart,tfinal     
#endif
! 
! Assumes cells have already been mapped to particles
!

#ifdef PERF
      tstart = MPI_WTIME()     
#endif

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

#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TSendFluidFields = tfinal - tstart    
      tstart = MPI_WTIME()
#endif

      !CALL MPI_BARRIER(ppiclf_comm,ierr)
      ! Interpolates rprop data for ppiclf domain cells in this bin
      CALL ppiclf_solve_Interpolate

      ! Reset for next iteration. Input from rocpicl/PICL_TEMP_Runge
      PPICLF_INT_ICNT = 0

#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TInterpolation = tfinal - tstart    
#endif

#ifdef PERF
        tfinal = MPI_WTIME()
        PPICLF_TSendGridOverlap = tfinal - tstart
#endif
#ifdef PERF
      PPICLF_TSendGhostParticles = 0.0D0
#endif
      IF(PPICLF_PPInteractions) THEN
#ifdef PERF
      tstart = MPI_WTIME()
#endif
        ! Ghost particles are needed 
        ! They're made here after interpolated
        ! values are updated.
        CALL ppiclf_comm_CreateGhostPartLB
        CALL ppiclf_comm_MoveGhostPartLB
        CALL ppiclf_comm_subbinGhostParticleMap
        ! Zero collisions 
        ppiclf_ydotc = 0.0D0
#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TSendGhostParticles = tfinal-tstart
#endif
      END IF
      CALL ppiclf_user_SetYdot

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_InitSolvePartLB
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
! This is called during the initialization. It forms the first bin and
! paritcle-cell mappings.
! 
! Internal: 
! 
      INTEGER*4 :: j,ierr
      CALL ppiclf_comm_CreateBinPartLB
      CALL ppiclf_comm_FindParticlePartLB
      CALL ppiclf_comm_PartLoadBalance
      CALL ppiclf_comm_MoveParticlePartLB
      CALL ppiclf_comm_subbinRealParticleMap
      CALL ppiclf_comm_MapOverlapGridPartLB

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

      IF(PPICLF_PPInteractions) THEN
        ! Ghost particles are needed 
        CALL ppiclf_comm_CreateGhostPartLB
        CALL ppiclf_comm_MoveGhostPartLB
        CALL ppiclf_comm_subbinGhostParticleMap
        ! Zero collisions 
        ppiclf_ydotc = 0.0D0
      END IF

      ! Maps up to 27 closest cell centers to particle
      ! Includes: CellID, total dist, x dist, y dist, z dist
      CALL ppiclf_comm_subbinCellMap
      CALL ppiclf_solve_SBParticleToCellMap

      ! Interpolates rprop data for ppiclf domain cells in this bin
      CALL ppiclf_solve_Interpolate
      ! Reset for next iteration. Input from rocpicl/PICL_TEMP_Runge
      PPICLF_INT_ICNT = 0

      ! Project particle feedback to fluid solver grid
      CALL ppiclf_solve_ProjectParticleGrid
 
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_PrintQuantities

      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      INTEGER*4  maxGP, maxRP, maxOC, ierr

      CALL MPI_ALLREDUCE(ppiclf_npart,maxRP,1
     >                   ,MPI_INTEGER4, MPI_MAX
     >                   ,ppiclf_comm, ierr)
      CALL MPI_ALLREDUCE(ppiclf_npart_gp,maxGP,1
     >                   ,MPI_INTEGER4, MPI_MAX
     >                   ,ppiclf_comm, ierr)
      CALL MPI_ALLREDUCE(ppiclf_nCells_Interp,maxOC,1
     >                   ,MPI_INTEGER4, MPI_MAX
     >                   ,ppiclf_comm, ierr)

      IF(ppiclf_nid .EQ. 0) THEN
        PRINT*, '-----------------------------------------------------'
        PRINT*, 'Max number of Particles per rank:      ',maxRP 
        PRINT*, 'Max number of Ghost Particles per rank:',maxGP 
        PRINT*, 'Max number of Overlap Cells per rank:  ',maxOC 
        PRINT*, '-----------------------------------------------------'
        PRINT*, ''
      END IF

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------

      SUBROUTINE ppiclf_solve_InitSolve
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
! This is called during the initialization. It forms the first bin and
! paritcle-cell mappings.
! 
! Internal: 
! 
      INTEGER*4 :: j
      ! ppiclf_binchanged set in CreateBin
      ! ppiclf_binchanged .TRUE. means
      ! bin coordinates changed
      CALL ppiclf_comm_CreateBin

!changed to T/F      ! ppiclf_particleMoved set in FindParticle
      ! ppiclf_particleMoved .EQ. 0 means all particles
      ! stayed in same bin as previous RK Stage.
      CALL ppiclf_comm_FindParticle
!      IF(ppiclf_particleMoved .NE. 0 .OR.
!     >              ppiclf_binchanged) THEN
        CALL ppiclf_comm_MoveParticle
!      END IF

      IF(ppiclf_overlap .AND. ppiclf_binchanged) THEN
        CALL ppiclf_comm_MapOverlapGrid
      END IF

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

      IF(PPICLF_PPInteractions) THEN
        ! Ghost particles are needed 
        CALL ppiclf_comm_CreateGhost
        CALL ppiclf_comm_MoveGhost
        ! Zero collisions 
        ppiclf_ydotc = 0.0D0
      END IF

      ! Maps up to 27 closest cell centers to particle
      ! Includes: CellID, total dist, x dist, y dist, z dist
      CALL ppiclf_solve_ParticleToCellMap

      ! Interpolates rprop data for ppiclf domain cells in this bin
      CALL ppiclf_solve_Interpolate
      ! Reset for next iteration. Input from rocpicl/PICL_TEMP_Runge
      PPICLF_INT_ICNT = 0

      ! Project particle feedback to fluid solver grid
      CALL ppiclf_solve_ProjectParticleGrid


      RETURN
      END

!----------------------------------------------------------------------

      SUBROUTINE ppiclf_solve_PostTimeStepPartLB
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
! 

! Internal: 
! 
      INTEGER*4 :: i, j,ierr
#ifdef PERF
      REAL *8 tstart,tfinal     
      tstart = MPI_WTIME()
#endif

! There is some error in the binchanged, particleMoved, or rebalance
! logic below.  Commenting out for now, but could be an efficiency
! gain later after troubleshooting.

      ! ppiclf_binchanged set in CreateBin
      ! ppiclf_binchanged .TRUE. means
      ! bin coordinates changed
      CALL ppiclf_comm_CreateBinPartLB
#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TCreateBin = tfinal - tstart
      tstart = MPI_WTIME()
#endif
      CALL ppiclf_comm_FindParticlePartLB
!      IF(ppiclf_rebalance) THEN
         CALL ppiclf_comm_PartLoadBalance
!      END IF
!      IF(ppiclf_particleMoved) THEN
!        IF(.NOT. ppiclf_rebalance) THEN
          ! Already called when ppiclf_rebalance=.TRUE.
!          CALL ppiclf_comm_setEmptyIndicator
!          CALL ppiclf_comm_setInterfaceIndicator
!        END IF
        CALL ppiclf_comm_MoveParticlePartLB
!      END IF
      CALL ppiclf_comm_subbinRealParticleMap
#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TSendParticles = tfinal - tstart
      PPICLF_TSendGridOverlap = 0.0D0
#endif
!      IF(ppiclf_binchanged) THEN
#ifdef PERF
        tstart = MPI_WTIME()
#endif
        CALL ppiclf_comm_MapOverlapGridPartLB
!      END IF
      CALL ppiclf_comm_subbinCellMap
#ifdef PERF
      tstart= MPI_WTIME()
#endif
      ! Maps up to 27 closest cell centers to particle
      ! Includes: CellID, total dist, x dist, y dist, z dist
      CALL ppiclf_solve_SBParticleToCellMap
#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TMapParticlesCells = tfinal - tstart
      tstart = MPI_WTIME()
#endif
      ! Project particle feedback to fluid solver grid
      CALL ppiclf_solve_ProjectParticleGrid
#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TProjection = tfinal - tstart
#endif
      RETURN
      END

!----------------------------------------------------------------------

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
      INCLUDE "mpif.h"
!
! Internal: 
!
      REAL*8 FLD(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE)
      INTEGER*4 nkey(2), nl, nii, njj, nrr  
      LOGICAL partl
#ifdef PERF
      REAL *8 tstart,tfinal     
#endif
!
      ! send it all
      nl   = 0
      nii  = PPICLF_LRMAX
      njj  = 3
      nrr  = PPICLF_LRP_INT
      nkey(1) = 2
      nkey(2) = 1

#ifdef PERF
      tstart = MPI_WTIME()
#endif

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

#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TDataTransfers = PPICLF_TDataTransfers + (tfinal - tstart)
#endif

      RETURN
      END

!
!______________________________________________________________________
!
      SUBROUTINE ppiclf_solve_SBParticleToCellMap

      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      ! --- Local Variables ---
      INTEGER*4  i, j, k, l, ip, ie, nnearest, partCount
      INTEGER*4  n_SBin1x2
      INTEGER*4  i_SBin(3), i_count, ii, jj, kk, loopSB
      INTEGER*4  CellID_nearest(27), idx_worst
      REAL*8     dSQl, dSQi, dSQ(27), xp(3), binblength(3), dSQchk(3)
      REAL*8     dSQ_worst
      LOGICAL, ALLOCATABLE, SAVE :: cell_is_in_list(:)
      LOGICAL    farAway
      LOGICAL    wrap_x, wrap_y, wrap_z

#ifdef PERF
      REAL*8    tstart,tfinal
#endif
      
      IF(ppiclf_npart < 1) RETURN

      ! --- Initialization ---
      ! Allocate the hash table only if needed
      IF (.NOT. ALLOCATED(cell_is_in_list)) THEN
          ALLOCATE(cell_is_in_list(1:PPICLF_LEE))
      END IF

      DO i = 1,3
        binblength(i) = ppiclf_binb(2*i) - ppiclf_binb(2*i-1)
        dSQchk(i) = ppiclf_interp_dchk(i)**2
      END DO

      n_SBin1x2 = ppiclf_nSBin(1)*ppiclf_nSBin(2)

      wrap_x = ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)
      wrap_y = ppiclf_linperiodic(2) .AND. ppiclf_EqualDomain(2)
      wrap_z = ppiclf_linperiodic(3) .AND. ppiclf_EqualDomain(3)
      cell_is_in_list = .FALSE. ! Reset hash table
      partCount = 0
      ! --- Main Particle Loop ---
      DO ip=1,ppiclf_npart
        ! Initialize arrays for this particle's search
        nnearest = 0
        dSQ(1:27) = 1.0E20
        CellID_nearest(1:27) = -1


        xp(1) = ppiclf_y(PPICLF_JX, ip)
        xp(2) = ppiclf_y(PPICLF_JY, ip)
        xp(3) = ppiclf_y(PPICLF_JZ, ip)

        i_SBin(1) = ppiclf_iprop(5,ip) - ppiclf_binOffset(1)
        i_SBin(2) = ppiclf_iprop(6,ip) - ppiclf_binOffset(2)
        i_SBin(3) = ppiclf_iprop(7,ip) - ppiclf_binOffset(3)

        ! --- Neighbor Sub-Bin Search Loop ---
        DO k = -1, 1
          kk = i_SBin(3) + k
          IF(kk < 0 .OR. kk >= ppiclf_nSBin(3)) CYCLE
          DO j = -1, 1
            jj = i_SBin(2) + j
            IF(jj < 0 .OR. jj >= ppiclf_nSBin(2)) CYCLE
            DO i = -1, 1
              ii = i_SBin(1) + i
              IF(ii < 0 .OR. ii >= ppiclf_nSBin(1)) CYCLE
              
              loopSB = ii + jj*ppiclf_nSBin(1) + kk*n_SBin1x2
              IF(loopSB < 0 .OR. loopSB >= ppiclf_total_SBin) CYCLE

              DO i_count = 1, ppiclf_binCellCount(loopSB)
                ie = ppiclf_binCellList(loopSB, i_count)
                
                IF (cell_is_in_list(ie)) CYCLE

                dSQi = 0.0D0
                farAway = .FALSE.
                DO l=1,3
                  IF( (l==1.AND.wrap_x) .OR.
     >                (l==2.AND.wrap_y) .OR. (l==3.AND.wrap_z) ) THEN
                    dSQl = MIN((ppiclf_picl_grid(l,ie) - xp(l))**2, 
     >                     (binblength(l)-
     >                      ABS(ppiclf_picl_grid(l,ie) - xp(l)))**2)
                  ELSE
                    dSQl = (ppiclf_picl_grid(l,ie) - xp(l))**2
                  END IF
                  dSQi = dSQi + dSQl
                  IF (dSQl > dSQchk(l)) farAway = .TRUE.
                END DO ! l
                
                IF (farAway) CYCLE

                ! --- "Replace-the-worst" Neighbor Finding ---
                IF (nnearest < 27) THEN
                  nnearest = nnearest + 1
                  dSQ(nnearest) = dSQi
                  CellID_nearest(nnearest) = ie
                  cell_is_in_list(ie) = .TRUE.
                  IF (nnearest .EQ. 27) THEN
                    ! List is now full, find the initial worst neighbor
                    dSQ_worst = dSQ(1)
                    idx_worst = 1
                    DO l = 2, 27
                      IF (dSQ(l) > dSQ_worst) THEN
                        dSQ_worst = dSQ(l)
                        idx_worst = l
                      END IF
                    END DO
                  END IF
                ELSE IF (dSQi < dSQ_worst) THEN
                  ! New cell is better than the worst one, replace it
                  cell_is_in_list(CellID_nearest(idx_worst)) = .FALSE.
                  dSQ(idx_worst) = dSQi
                  CellID_nearest(idx_worst) = ie
                  cell_is_in_list(ie) = .TRUE.
                  ! Find the new worst in the updated list
                  dSQ_worst = dSQ(1)
                  idx_worst = 1
                  DO l = 2, 27
                    IF (dSQ(l) > dSQ_worst) THEN
                      dSQ_worst = dSQ(l)
                      idx_worst = l
                    END IF
                  END DO ! l
                END IF
              END DO ! i_count
            END DO ! i
          END DO ! j
        END DO ! k
        
        ! --- Finalize mapping for this particle ---
        IF (nnearest < 1) THEN
          ppiclf_iprop(9,ip) = -1
          ppiclf_remove_particle = .TRUE.
          PRINT*, 'part # on proc # & bin # removed.'
     >            ,ppiclf_iprop(1,ip),ppiclf_nid, ppiclf_iprop(8,ip)
        ELSE
          partCount = partCount + 1
          ppiclf_nPart2Cell(partCount) = nnearest
          DO i = 1, nnearest
            ppiclf_Part2Cell_map(partCount,i) = CellID_nearest(i)
            ppiclf_Part2Cell_dist(partCount,i) = SQRT(dSQ(i))
            cell_is_in_list(CellID_nearest(i)) = .FALSE.
          END DO
        END IF
      END DO ! ip

      ! --- Cleanup removed particles ---
      IF(ppiclf_remove_particle) THEN
        CALL ppiclf_solve_RemoveParticle
        ppiclf_remove_particle = .FALSE.
      END IF

      RETURN
      END SUBROUTINE

!
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ppiclf_solve_Interpolate

      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

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
      INCLUDE "mpif.h"

! Internal:
!
      INTEGER*4 i, icount
!
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
               IF(PPICLF_LRP2 .GT. 0) THEN
                 CALL ppiclf_copy
     >          (ppiclf_rprop2(1,icount),ppiclf_rprop2(1,i),PPICLF_LRP2)
               END IF
               IF(PPICLF_LRP3 .GT. 0) THEN
                 CALL ppiclf_copy
     >          (ppiclf_rprop3(1,icount),ppiclf_rprop3(1,i),PPICLF_LRP3)
               END IF
               IF(PPICLF_LRP4 .GT. 0) THEN
                 CALL ppiclf_copy
     >          (ppiclf_rprop4(1,icount),ppiclf_rprop4(1,i),PPICLF_LRP4)
               END IF
               IF(PPICLF_LRP5 .GT. 0) THEN
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
!
!----------------------------------------------------------------------
! 
      SUBROUTINE ppiclf_solve_ProjectParticleGrid
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      ! Internal:
      INTEGER*4 i, j, ip, ie, nCellProj, CellID, nl, nii, njj,
     >          nrr, nkey(2), iee
      REAL*8    CellVol, GaussianConst, dist, w(27), wsum,
     >          x_norm, y_norm, z_norm, PI, eps, avg_vol
      LOGICAL   partl 

#ifdef PERF
      REAL*8    tstart, tfinal
#endif
! 
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
        avg_vol = 0.0D0
        DO i = 1,nCellProj
          CellID = ppiclf_Part2Cell_map(ip,i) 
          avg_vol = avg_vol + ppiclf_picl_grid(7,CellID)
        END DO
        avg_vol = avg_vol/DBLE(nCellProj)
        
        ! Loop to find individual cell weightings
        DO i = 1,nCellProj
          CellID = ppiclf_Part2Cell_map(ip,i) 
          dist = ppiclf_Part2Cell_dist(ip,i) + eps
          CellVol = ppiclf_picl_grid(7,CellID)
          avg_vol = CellVol !seeing if this is unit test difference
          ! Multiply CellVol by ppiclf_feedbk(PPICLF_P_JPHIP,ip) below.  
          ! Need to think through way to do it if VF not defined.
          w(i) = ABS(CellVol*EXP(-GaussianConst*(dist**2)
     >              / (avg_vol**(2.0D0/3.0D0))))
          wsum = wsum + w(i)
        END DO !i

        IF (wsum .GT. eps) THEN
          DO i = 1, nCellProj
            w(i) = w(i) / wsum
          END DO
        ELSE
          DO i = 1, nCellProj
            w(i) = 0.0D0
          END DO
        END IF

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

        ! Loop through cells to apply feedback     
        ! DO j is outer to match column-major 
        ! ppiclf_pro_fld_picl(j, CellID)
        DO j=1,PPICLF_LRP_PRO
          DO i = 1,nCellProj
            CellID = ppiclf_Part2Cell_map(ip,i)
            ppiclf_pro_fld_picl(j,CellID) = 
     >         ppiclf_pro_fld_picl(j,CellID) 
     >         + ppiclf_feedbk(j,ip)*w(i)
          END DO !i
        END DO !j
      END DO !ip

      ! Now send feedback information to processor that contains 
      ! the cell for the fluid solver
      ppiclf_nCells_Proj = ppiclf_nCells_FV2PICL
      DO i = 1,ppiclf_nCells_Proj
        CALL ppiclf_icopy(ppiclf_cell_map_proj(1,i),
     >         ppiclf_cell_map(1,i),PPICLF_LRMAX)
      END DO

      nl = 0
      nii = PPICLF_LRMAX
      njj = 2 ! original processor with cell for fluid grid
      nrr = PPICLF_LRP_PRO
      nkey(1) = 2
      nkey(2) = 1

#ifdef PERF
      tstart = MPI_WTIME()
#endif
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
#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TDataTransfers = PPICLF_TDataTransfers + (tfinal - tstart)
#endif

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
      END SUBROUTINE

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
      
!______________________________________________________________________
!

      SUBROUTINE OLDppiclf_solve_PostTimeStep
!     This is the 1 bin per rank version
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
! 
! Internal: 
! 
      INTEGER*4 :: i, j,ierr
#ifdef PERF
      REAL *8 tstart,tfinal     
      tstart = MPI_WTIME()
#endif
      ! ppiclf_binchanged set in CreateBin
      ! ppiclf_binchanged .TRUE. means
      ! bin coordinates changed
      CALL ppiclf_comm_CreateBin
#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TCreateBin = tfinal - tstart
      tstart = MPI_WTIME()
#endif
! Changed to T/F      ! ppiclf_particleMoved set in FindParticle
      ! ppiclf_particleMoved .EQ. 0 means all particles
      ! stayed in same bin as previous RK Stage.
      CALL ppiclf_comm_FindParticle
!      IF(ppiclf_particleMoved .NE. 0 .OR.
!     >              ppiclf_binchanged) THEN
        CALL ppiclf_comm_MoveParticle
!      END IF
#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TSendParticles = tfinal - tstart
      PPICLF_TSendGridOverlap = 0.0D0
#endif
      IF(ppiclf_overlap .AND. ppiclf_binchanged) THEN
#ifdef PERF
        tstart = MPI_WTIME()
#endif
        CALL ppiclf_comm_MapOverlapGrid
#ifdef PERF
        tfinal = MPI_WTIME()
        PPICLF_TSendGridOverlap = tfinal - tstart
#endif
      END IF
#ifdef PERF
      PPICLF_TSendGhostParticles = 0.0D0
#endif
      IF(PPICLF_PPInteractions) THEN
#ifdef PERF
      tstart = MPI_WTIME()
#endif
        ! Ghost particles are needed 
        CALL ppiclf_comm_CreateGhost
        CALL ppiclf_comm_MoveGhost
        ! Zero collisions 
        ppiclf_ydotc = 0.0D0
#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TSendGhostParticles = tfinal-tstart
#endif
      END IF
#ifdef PERF
      tstart= MPI_WTIME()
#endif
      CALL MPI_BARRIER(ppiclf_comm,ierr)
      ! Maps up to 27 closest cell centers to particle
      ! Includes: CellID, total dist, x dist, y dist, z dist
      CALL ppiclf_solve_SBParticleToCellMap
#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TMapParticlesCells = tfinal - tstart
      tstart = MPI_WTIME()
#endif
      ! Project particle feedback to fluid solver grid
      CALL ppiclf_solve_ProjectParticleGrid
#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TProjection = tfinal - tstart
#endif
      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_ParticleToCellMap
!     This one is O(NPART**2) - much slower than subbin NN check
!     OLD SUBROUTINE!!!
      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      ! Local Variables
      INTEGER*4 i, j, k, l, ix, iy, iz, ip, ie, iee, nxyz, nnearest, 
     >          CellID_nearest(28), partCount
      REAL*8    dSQl, dSQi, dSQ(28), xp(3), dSQchk(3), 
     >          CellCenter(3,28), w(27),binblength(3),  
     >          Max_CellLen(3),Max_CellLenSQ(3)
      LOGICAL   added, farAway
      !***************************************************************

      IF(ppiclf_nCells_FV2PICL .EQ. 0 . AND. ppiclf_npart .GT. 0) THEN
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
        DO ie = 1,ppiclf_nCells_FV2PICL
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
