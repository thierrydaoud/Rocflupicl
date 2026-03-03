#:include "PPICLF_PARTMACROS.fypp"
#include "PPICLF_STD.h"
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine to find nearest neighbor for particle-particle and
!     particle-wall interactions. Subroutine also includes
!     the new added-mass binary terms, developed by Sam Briney.
!
! Added: if collisional_flag = 1  F = Fn
!                            = 2  F = Fn + Ft + Tt
!                            = 3  F = Fn + Ft + Tt + Th + Tr
! where Tt = collisional torque
!       Th = hydrodynamic torque
!       Tr = rolling torque
!
! Note that the collisional and rolling torques are due to
!     particle-particle interactions and thus are evaulated
!     in ppiclf_user_EvalNearestNeighbor. Therefore, only the
!     hydrodynamic torque is left to be calculated.
!-----------------------------------------------------------------------
!
submodule (ppiclf_user) ppiclf_user_EvalNearestNeighbor_imp
    ! particle data
    use ppiclf_data, only: ppiclf_npart
    use ppiclf_m_particledata, only: @{USEMODVAR(PPICLF_t_particle, ppiclf_parts)}@
    ! grid data
    use ppiclf_data, only:
    use ppiclf_data, only:
    use ppiclf_data, only:
    ! particle options variables
    use ppiclf_data, only:
    use ppiclf_data, only: ppiclf_ndim
    use ppiclf_data, only: ppiclf_nndist, ppiclf_dt, ppiclf_time, ppiclf_rk3ark, ppiclf_filter
    ! use ppiclf_data, only:
    ! comm variables
    use ppiclf_data, only: ppiclf_nid
    ! binning variables
    use ppiclf_data, only: ppiclf_n_bins, ppiclf_bins_dx
    ! ghost particle variables
    use ppiclf_data, only: ppiclf_npart_gp
    ! wall support variables
    use ppiclf_data, only:
    ! AngularPeriodic variables (?)(SEE NOTE IN ppiclf_data)
    use ppiclf_data, only:


    use ppiclf_m_user_data
    use ppiclf_m_user_RFLUdata

    use ppiclf_user_AM_functions, only: resistance_pair
    use ppiclf_op, only: ppiclf_exittr
    implicit none
    contains

module procedure ppiclf_user_EvalNearestNeighbor
    !
    ! Input:
    !
    !      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
    !     >   collisional_flag, heattransfer_flag, feedback_flag,
    !     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
    !     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
    !     >   qs_fluct_filter_adapt_flag,
    !     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData, ppiclf_nTimeBH,
    !     >   sbNearest_flag, burnrate_flag, flow_model
    !      real*8 :: rmu_ref, tref, suth, ksp, erest
    !      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
    !     >   collisional_flag, heattransfer_flag, feedback_flag,
    !     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
    !     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
    !     >   qs_fluct_filter_adapt_flag, ksp, erest,
    !     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData, ppiclf_nTimeBH,
    !     >   sbNearest_flag, burnrate_flag, flow_model

    ! @{USEMODVAR(PPICLF_t_ghostParticle, neighbor)}@
    !
    ! Internal:
    !
    real*8 pi2, rthresh, rxdiff, rydiff, rzdiff, rdiff, rm1, rm2,   &
        rmult, eta_n, rbot, rn_12x, rn_12y, rn_12z, rdelta12,       &
        rv12_mag, rv12_mage, rksp_max, rnmag, rksp_wall, rextra,    &
        JDP_i,JDP_j
    real*8 eps

    ! Sam - from TLJ for box filter
    integer*4 ifilt
    real*8 adptfilter, dpl, phip, dist2, xdist2, ydist2, zdist2
    real*8 dist, rsig
    real*8 sig2, gkern, pi
    ! TODO: Fix this
    parameter(pi = 3.1415926535)

    ! 06/06/2024 - Thierry - Added Mass code
    integer*4 k, l, kk, ll
    real*8 alpha_local, rad
    real*8 rxdiff1, rydiff1, rzdiff1

    ! 07/16/2024 - TLJ added tangential component
    ! - does not take into account angular velocity
    real*8 unx, uny, unz, utx, uty, utz, ut_mag, rn_mag
    real*8 rt_12x, rt_12y, rt_12z
    real*8 Fn_mag, Ftmin
    real*8 eta_t, mu_c
    real*8 A12x, A12y, A12z 
    real*8 rad1, rad2
    real*8 u12x, u12y, u12z
    real*8 tcx, tcy, tcz
    real*8 trx, try, trz
    real*8 thetar, dp1, dp2, r12
    real*8 omgrx, omgry, omgrz, omgr_mag
    real*8 Ftx, Fty, Ftz
    
    ! 04/03/2025 - TLJ added for spring stiffness coefficient
    real*8 nu1, nu2
    real*8 E1, E2, Estar
    real*8 r1, r2, Rstar 
    real*8 ksp1, ksp2, ksp_min

    ! 01/29/2025 - Thierry - added for particle collision with conical
    !                         wall domain
    real*8 rp, yp, zp, vp, wp, rp_new, yp_new, zp_new, vp_new, wp_new, rbound, urp, thetap

    !
    ! Code:
    !
    pi2  = rpi*rpi

    ! other particles
    if (j .ne. 0) then
        !Added spload and radius factor

        ! Compute mean particle diameter between i and j; delta_{ij}
        rthresh  = 0.5d0*(@{USEPARTICLE(ppiclf_parts(i)%rprop%DP)}@ + @{USEPARTICLE(neighbor%rprop%DP)}@)

        ! Compute vector components and distance between 
        !    centers of particles i and j; D_{ij}
        rxdiff = (@{USEPARTICLE(neighbor%y%pos%X)}@) - (@{USEPARTICLE(ppiclf_parts(i)%y%pos%X)}@)
        rydiff = (@{USEPARTICLE(neighbor%y%pos%Y)}@) - (@{USEPARTICLE(ppiclf_parts(i)%y%pos%Y)}@)
        rzdiff = (@{USEPARTICLE(neighbor%y%pos%Z)}@) - (@{USEPARTICLE(ppiclf_parts(i)%y%pos%Z)}@)
        
        rdiff = rxdiff**2 + rydiff**2 + rzdiff**2
        rdiff = sqrt(rdiff)

        !-----------------------------------------------------------------------
        !
        ! For binary added-mass for Briney model

        ! 06/06/2024 - Thierry - Added Mass code continues here
        ! 07/09/2024 - TLJ - Updated
        ! 07/14/2024 - Thierry - Updated overlapping particles if statement

        ! Filter widths are set to be equal to 2*cell length in x,y,z
        ! directions (1:3)
        ! ppiclf_nndist is neighbor width - user defined
        !   = max(NEIGHBORWIDTH,4*Dp)
    
        if (am_flag == 2 .and. rdiff <= ppiclf_nndist) then
            ! Do not overwrite rxdiff, rydiff, rzdiff
            rxdiff1 = rxdiff
            rydiff1 = rydiff
            rzdiff1 = rzdiff

            ! Check if particles are overlapping, replace value
            !   if yes, since we get crazy resistance values
            if (rdiff < rthresh) then
                rxdiff1 = rthresh
                rydiff1 = rthresh
                rzdiff1 = rthresh 
            endif

            ! Model only valid for local volume fraction
            ! less than 0.4, so we limit it here without
            ! over riding rphip
            ! limit alpha to mitigate misuse
            alpha_local = min(0.4, rphip) 
            
            ! Compute the resistance matrix
            ! Only valid for monodispersed particles
            rad = 0.5d0*(@{USEPARTICLE(ppiclf_parts(i)%rprop%DP)}@)

            call resistance_pair(rxdiff1, rydiff1, rzdiff1, alpha_local, rad, R_pair)
                
            ! accumulate number of neighbors
            nneighbors = nneighbors + 1
                
            do k=1,3
                do l=1,3
                    ! added mass
                    Fam(k) = Fam(k) + R_pair(k,l)   * @{USEPARTICLE(ppiclf_parts(i)%rprop%WDOT, skipIndex)}@(l)
                    ! induced added mass
                    Fam(k) = Fam(k) + R_pair(k,l+3) * @{USEPARTICLE(neighbor%rprop%WDOT, skipIndex)}@(l)
                end do ! l-loop

                ! accumulate neighbor acceleration
                Wdot_neighbor_mean(k) = Wdot_neighbor_mean(k)  + @{USEPARTICLE(neighbor%rprop%WDOT, skipIndex)}@(k)
            end do ! k-loop
        end if ! am_flag==2 .and. rdiff <= ppiclf_nndist

        !-----------------------------------------------------------------------
        !
        ! For particle-particle collision

        ! Cycle if rdiff > rthresh
        eps = 0.0d0
        if (rdiff .lt. rthresh+eps) then

            ! Compute spring stiffness constant dynamically.
            ! The number of collision timesteps (ksp) is set by the user
            ! k1 = k_{n,limit}
            ksp1 = rmass*rpi*rpi/((ksp*ppiclf_dt)**2)
            ! k2 = k_{hertzian}
            E1  = 1.0d9  ! Assumed value for Young's modulus
            E2  = 1.0d9  ! Assumed value for Young's modulus
            nu1 = 0.35d0 ! Assumed value for Poisson's ratio
            nu2 = 0.35d0 ! Assumed value for Poisson's ratio
            Estar = (1.0d0-nu1*nu1)/E1 + (1.0d0-nu2*nu2)/E2
            Estar = 1.0d0/Estar
            r1 = 0.5d0*(@{USEPARTICLE(ppiclf_parts(i)%rprop%DP)}@)
            r2 = 0.5d0*(@{USEPARTICLE(neighbor%rprop%DP)}@)
            Rstar = r1*r2/(r1+r2)
            ksp2 = (4.0d0/3.0d0)*Estar*sqrt(Rstar)
            ksp2 = ksp2*sqrt(abs(rdiff-rthresh))
            ! kn = min(k1,k2)
            ksp_min = min(ksp1,ksp2)

            rm1 = (@{USEPARTICLE(ppiclf_parts(i)%rprop%RHOP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@)
            rm2 = (@{USEPARTICLE(neighbor%rprop%RHOP)}@) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
            
            rmult = (rm1*rm2)/(rm1+rm2)
            eta_n = -2.0d0*sqrt(ksp_min)*log(erest)/sqrt(log(erest)**2+pi2)*sqrt(rmult)

            !            print*,'COLLS: ',i,j,ksp1,ksp2,ksp_min,
            !     >              eta_n,rdiff-rthresh,vmag

            ! Compute unit normal vector along line of contact 
            !   pointing from particle i to particle j
            rbot = 1.0d0/rdiff
            rn_12x = rxdiff*rbot
            rn_12y = rydiff*rbot
            rn_12z = rzdiff*rbot
            rn_mag = rdiff
         
            ! Relative velocity in normal direction
            u12x = (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%X)}@) - (@{USEPARTICLE(neighbor%y%Vel%X)}@)
            u12y = (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%Y)}@) - (@{USEPARTICLE(neighbor%y%Vel%Y)}@)
            u12z = (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%Z)}@) - (@{USEPARTICLE(neighbor%y%Vel%Z)}@)

            if (collisional_flag>=2) then
               ! Add contribution from angular velocity
               rad1 = 0.5d0*(@{USEPARTICLE(ppiclf_parts(i)%rprop%DP)}@)
               rad2 = 0.5d0*(@{USEPARTICLE(neighbor%rprop%DP)}@)
               A12x = rad1 * (@{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%X)}@) + rad2 * (@{USEPARTICLE(neighbor%y%ang_vel%X)}@)
               A12y = rad1 * (@{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%Y)}@) + rad2 * (@{USEPARTICLE(neighbor%y%ang_vel%Y)}@)
               A12z = rad1 * (@{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%Z)}@) + rad2 * (@{USEPARTICLE(neighbor%y%ang_vel%Z)}@)

               u12x = u12x + (A12y*rn_12z - A12z*rn_12y)
               u12y = u12y + (A12z*rn_12x - A12x*rn_12z)
               u12z = u12z + (A12x*rn_12y - A12y*rn_12x)
            endif

            ! Compute (u_ij \cdot n_ij)
            rv12_mag = u12x*rn_12x + u12y*rn_12y + u12z*rn_12z
         
            ! Compute delta_12 and normal parameters
            rdelta12 = rthresh - rdiff
            rksp_max  = ksp_min*rdelta12
            rv12_mage = rv12_mag*eta_n
            rnmag     = -rksp_max - rv12_mage

            ! Normal collision force Fn = -rnmag*n_{ij}
            ! Scalar magnitude |Fn| = abs(rnmag)
            Fn_mag = abs(rnmag)

            ! Compute tangential unit vector
            unx = rv12_mag*rn_12x
            uny = rv12_mag*rn_12y
            unz = rv12_mag*rn_12z
            utx = u12x - unx
            uty = u12y - uny
            utz = u12z - unz
            ut_mag = sqrt(utx*utx + uty*uty + utz*utz)
            ut_mag = max(ut_mag,1.0d-8)
            rt_12x = utx/ut_mag
            rt_12y = uty/ut_mag
            rt_12z = utz/ut_mag

            ! Compute tangential collision force
            Ftmin = 0.0d0
            if (collisional_flag>=2) then ! Tangential component
                if (ut_mag > 0) then
                    mu_c  = 0.4d0  ! Dimensionless; Coulomb
                    eta_t = eta_n  ! Set to normal; damping
                    Ftmin  = -min(mu_c*Fn_mag,eta_t*ut_mag)  
                endif
            endif

            ! Compute contributions to angular velocities
            tcx = 0.0d0; tcy = 0.0d0; tcz = 0.0d0;
            trx = 0.0d0; try = 0.0d0; trz = 0.0d0;

            if (collisional_flag>=2) then

               ! Tangential force and Collision torque contributions
               Ftx = Ftmin*rt_12x
               Fty = Ftmin*rt_12y
               Ftz = Ftmin*rt_12z
               rad1 = 0.5d0* (@{USEPARTICLE(ppiclf_parts(i)%rprop%DP)}@) ! rpropi(PPICLF_R_JDP)
               tcx = rad1*(rn_12y*Ftz - rn_12z*Fty)
               tcy = rad1*(rn_12z*Ftx - rn_12x*Ftz)
               tcz = rad1*(rn_12x*Fty - rn_12y*Ftx)

               if (collisional_flag>=3) then
                  ! Add Rolling torque contribution
                  thetar = 0.06  ! Needs to be calibrated
                  dp1 = @{USEPARTICLE(ppiclf_parts(i)%rprop%DP)}@
                  dp2 = @{USEPARTICLE(neighbor%rprop%DP)}@
                  r12 = 0.5d0*(dp1*dp2)/(dp1+dp2)
                  omgrx = @{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%X)}@ - @{USEPARTICLE(neighbor%y%ang_vel%X)}@
                  omgry = @{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%Y)}@ - @{USEPARTICLE(neighbor%y%ang_vel%Y)}@
                  omgrz = @{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%Z)}@ - @{USEPARTICLE(neighbor%y%ang_vel%Z)}@
                  omgr_mag = sqrt(omgrx*omgrx+omgry*omgry+omgrz*omgrz)
                  omgr_mag = max(omgr_mag,1.d-8)
                  trx = -thetar*Fn_mag*r12*omgrx/omgr_mag
                  try = -thetar*Fn_mag*r12*omgry/omgr_mag
                  trz = -thetar*Fn_mag*r12*omgrz/omgr_mag
               endif
            endif


            ! Now update that part of the RHS of equations 
            !   that involve nearest neighbors

            ! Particle velocities
            @{USEPARTICLE(ppiclf_parts(i)%ydotc%vel%X)}@ = @{USEPARTICLE(ppiclf_parts(i)%ydotc%vel%X)}@ + rnmag*rn_12x + Ftmin*rt_12x
            @{USEPARTICLE(ppiclf_parts(i)%ydotc%vel%Y)}@ = @{USEPARTICLE(ppiclf_parts(i)%ydotc%vel%Y)}@ + rnmag*rn_12y + Ftmin*rt_12y
            @{USEPARTICLE(ppiclf_parts(i)%ydotc%vel%Z)}@ = @{USEPARTICLE(ppiclf_parts(i)%ydotc%vel%Z)}@ + rnmag*rn_12z + Ftmin*rt_12z

            ! Particle angular velocities
            @{USEPARTICLE(ppiclf_parts(i)%ydotc%ang_vel%X)}@ = @{USEPARTICLE(ppiclf_parts(i)%ydotc%ang_vel%X)}@ + tcx + trx
            @{USEPARTICLE(ppiclf_parts(i)%ydotc%ang_vel%Y)}@ = @{USEPARTICLE(ppiclf_parts(i)%ydotc%ang_vel%Y)}@ + tcy + try
            @{USEPARTICLE(ppiclf_parts(i)%ydotc%ang_vel%Z)}@ = @{USEPARTICLE(ppiclf_parts(i)%ydotc%ang_vel%Z)}@ + tcz + trz

        end if ! rdiff lt rthresh

        !-----------------------------------------------------------------------
        !
        ! Feedback fluctuation mean

        dist2 = MAX(ppiclf_filter(1),ppiclf_filter(2),ppiclf_filter(3))

        ! Box filter half-width dist2
        IF(qs_fluct_filter_adapt_flag.NE.0) THEN
            ! Adaptive filter defined wrt particle i
            ! Used for adaptive box or gaussian
            dpl = @{USEPARTICLE(ppiclf_parts(i)%rprop%DP)}@
            phip = @{USEPARTICLE(ppiclf_parts(i)%rprop%PHIP)}@
            adptfilter = ( 10.*(dpl**3)/max(1.e-4,phip) )**(1./3.)
            adptfilter = adptfilter/2.0
            IF(adptfilter .GT. dist2) dist2 = adptfilter
        END IF

        ! Check if particle lies inside box or gaussian filter
        xdist2 = abs((@{USEPARTICLE(ppiclf_parts(i)%y%pos%X)}@) - (@{USEPARTICLE(neighbor%y%pos%X)}@))
        if (xdist2 .gt. dist2) return

        ydist2 = abs((@{USEPARTICLE(ppiclf_parts(i)%y%pos%Y)}@) - (@{USEPARTICLE(neighbor%y%pos%Y)}@))
        if (ydist2 .gt. dist2) return

        if (ppiclf_ndim .eq. 3) then
            zdist2 = abs((@{USEPARTICLE(ppiclf_parts(i)%y%pos%Z)}@) - (@{USEPARTICLE(neighbor%y%pos%Z)}@))
            if (zdist2 .gt. dist2) return
        endif

        !
        ! The mean is calcuated according to Lattanzi etal,
        !   Physical Review Fluids, 2022.
        !
        if (j.ne.0) then
            if (qs_fluct_filter_flag==0) then
                upmean   = upmean + (@{USEPARTICLE(neighbor%y%vel%X)}@)
                vpmean   = vpmean + (@{USEPARTICLE(neighbor%y%vel%Y)}@)
                wpmean   = wpmean + (@{USEPARTICLE(neighbor%y%vel%Z)}@)
                u2pmean  = u2pmean + (@{USEPARTICLE(neighbor%y%vel%X)}@)**2
                v2pmean  = v2pmean + (@{USEPARTICLE(neighbor%y%vel%Y)}@)**2
                w2pmean  = w2pmean + (@{USEPARTICLE(neighbor%y%vel%Z)}@)**2
                icpmean  = icpmean + 1
            else if (qs_fluct_filter_flag==1) then
                ! See https://dpzwick.github.io/ppiclF-doc/algorithms/overlap_mesh.html
                dist = sqrt(xdist2**2 + ydist2**2 + zdist2**2)
                gkern = sqrt(pi*dist2**2/(4.0d0*log(2.0d0)))**(-ppiclf_ndim) * exp(-dist**2/(dist2**2/(4.0d0*log(2.0d0))))

                phipmean = phipmean + gkern*(@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                upmean   = upmean + gkern * (@{USEPARTICLE(neighbor%y%vel%X)}@) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                vpmean   = vpmean + gkern * (@{USEPARTICLE(neighbor%y%vel%Y)}@) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                wpmean   = wpmean + gkern * (@{USEPARTICLE(neighbor%y%vel%Z)}@) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                u2pmean  = u2pmean + gkern * ((@{USEPARTICLE(neighbor%y%vel%X)}@)**2) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                v2pmean  = v2pmean + gkern * ((@{USEPARTICLE(neighbor%y%vel%Y)}@)**2) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                w2pmean  = w2pmean + gkern * ((@{USEPARTICLE(neighbor%y%vel%Z)}@)**2) * (@{USEPARTICLE(neighbor%rprop%VOLP)}@)
                icpmean = icpmean + 1
            end if
        end if


        !-----------------------------------------------------------------------
        !
        ! boundaries
    elseif (j .eq. 0) then

        !rksp_wall = ksp
        !rksp_wall = 1000

        ! give a bit larger collision threshold for walls
        rextra   = 0.05d0 !
        ! add sploading and radius factor 
        rthresh  = (0.5d0+rextra)* (@{USEPARTICLE(ppiclf_parts(i)%rprop%DP)}@)
        
        rxdiff = @{USEPARTICLE(neighbor%y%pos%X)}@ - @{USEPARTICLE(ppiclf_parts(i)%y%pos%X)}@
        rydiff = @{USEPARTICLE(neighbor%y%pos%Y)}@ - @{USEPARTICLE(ppiclf_parts(i)%y%pos%Y)}@
        rzdiff = @{USEPARTICLE(neighbor%y%pos%Z)}@ - @{USEPARTICLE(ppiclf_parts(i)%y%pos%Z)}@
        
        rdiff = rxdiff**2 + rydiff**2 + rzdiff**2
        rdiff = sqrt(rdiff)
        
        if (rdiff .gt. rthresh) return

        rm1 = (@{USEPARTICLE(ppiclf_parts(i)%rprop%RHOP)}@)*(@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@)

        ! Compute spring stiffness constant dynamically, 
        !   which overrides the user defined value
        ! Need to make sure this formula is valid for a wall
        ! k1 = k_{n,limit}
        ksp1 = rm1*rpi*rpi/((ksp*ppiclf_dt)**2)
        ! k2 = k_{hertzian}
        E1  = 1.0d9  ! Assumed value for Young's modulus
        nu1 = 0.35d0 ! Assumed value for Poisson's ratio
        Estar = E1/(1.0d0-nu1*nu1)
        r1 = 0.5d0*@{USEPARTICLE(ppiclf_parts(i)%rprop%DP)}@
        r2 = r1
        Rstar = r1*r2/(r1+r2)
        ksp2 = (2.0d0/3.0d0)*Estar*sqrt(Rstar)
        ksp2 = ksp2*sqrt(abs(rdiff-rthresh))
        ! kn = min(k1,k2)
        rksp_wall = min(ksp1,ksp2)
        
        rmult = sqrt(rm1)
        eta_n = 2.0d0*sqrt(rksp_wall)*log(erest)/sqrt(log(erest)**2+pi2)*rmult
         
        rbot = 1.0d0/rdiff
        rn_12x = rxdiff*rbot
        rn_12y = rydiff*rbot
        rn_12z = rzdiff*rbot
    
        rdelta12 = rthresh - rdiff
        
        rv12_mag = - (@{USEPARTICLE(ppiclf_parts(i)%y%vel%X)}@) * rn_12x- (@{USEPARTICLE(ppiclf_parts(i)%y%vel%Y)}@) * rn_12y- (@{USEPARTICLE(ppiclf_parts(i)%y%vel%Z)}@) * rn_12z

        rv12_mage = rv12_mag*eta_n
        rksp_max  = rksp_wall*rdelta12
        rnmag     = -rksp_max - rv12_mage

         
        ! @{USEPARTICLE(ppiclf_parts(i)%ydotc%Vel%X)}@ = @{USEPARTICLE(ppiclf_parts(i)%ydotc%Vel%X)}@ + rnmag*rn_12x
        ! @{USEPARTICLE(ppiclf_parts(i)%ydotc%Vel%Y)}@ = @{USEPARTICLE(ppiclf_parts(i)%ydotc%Vel%Y)}@ + rnmag*rn_12y
        ! @{USEPARTICLE(ppiclf_parts(i)%ydotc%Vel%Z)}@ = @{USEPARTICLE(ppiclf_parts(i)%ydotc%Vel%Z)}@ + rnmag*rn_12z



        ! Particles leavind the domain with wall collisions
        ! Simple fix for a conical geometry
        yp = @{USEPARTICLE(ppiclf_parts(i)%y%pos%y)}@ !yi(PPICLF_JY)
        zp = @{USEPARTICLE(ppiclf_parts(i)%y%pos%z)}@ !yi(PPICLF_JZ)
        vp = @{USEPARTICLE(ppiclf_parts(i)%y%vel%y)}@ !yi(PPICLF_JVY)
        wp = @{USEPARTICLE(ppiclf_parts(i)%y%vel%z)}@ !yi(PPICLF_JVZ)

        ! rbound = sqrt(yj(PPICLF_JY)**2 + yj(PPICLF_JZ)**2)
        rbound = sqrt((@{USEPARTICLE(neighbor%y%pos%y)}@)**2 + (@{USEPARTICLE(neighbor%y%pos%z)}@)**2)
        rp = sqrt(yp**2 + zp**2)
        urp = sqrt(vp**2 + wp**2)
        thetap = atan2(zp, yp)

        if(rp > rbound) then
        
            rp_new = rp - (rp - rbound)
            yp_new = rp_new * cos(thetap)
            zp_new = rp_new * sin(thetap)

            urp = - urp
            
            vp_new = urp * cos(thetap)
            wp_new = urp * sin(thetap)

            ! ppiclf_y(PPICLF_JY,i) = yp_new
            ! ppiclf_y(PPICLF_JZ,i) = zp_new

            ! ppiclf_y(PPICLF_JVY,i) = vp_new
            ! ppiclf_y(PPICLF_JVZ,i) = wp_new
            @{USEPARTICLE(ppiclf_parts(i)%y%pos%Y)}@ = yp_new
            @{USEPARTICLE(ppiclf_parts(i)%y%pos%Z)}@ = zp_new

            @{USEPARTICLE(ppiclf_parts(i)%y%Vel%Y)}@ = vp_new
            @{USEPARTICLE(ppiclf_parts(i)%y%Vel%Z)}@ = wp_new
        endif
        
        !write(*,*) "Wall NEAR",i,ppiclf_ydotc(PPICLF_JVY,i)  
    endif


    return
end procedure ppiclf_user_EvalNearestNeighbor

end submodule ppiclf_user_EvalNearestNeighbor_imp

