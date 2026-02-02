#include "PPICLF_USER.h"
#include "PPICLF_STD.h"
!-----------------------------------------------------------------------
!
! Sam - Code for solving ydot = F(y,t)
!
! Set internal flags
!   These flags allows the user to turn on or
!      off various force terms, or to select
!      between different versions of each force model
!   To turn off all particle motion, set stationary = 0
!   To turn off a force set the flag to zero
!   For two-way coupling set collisional_flag = 0
!   For four-way coupling set collisional_flag = 1
!   To turn off feedback force set feedback_flag = 0
!   To use fluctuations turn corresponding flag = 1
!
!
!   stationary = 0 if 1, do not move particles but do
!         calculate drag forces; feedback force can also
!         be turned on
!
!   qs_flag = 2  ! none = 0; Parmar = 1; Osnes = 2
!   am_flag = 1  ! Parmar = 1; Briney = 2
!   pg_flag = 1
!   collisional_flag = 1 ! two way coupled = 0 ; four-way = 1
!
!   heattransfer_flag = 1
!
!   ViscousUnsteady_flag = 0 no viscous unsteady drag
!                        = 1 history kernal for visc. unsteady drag
!
!   feedback_flag = 1
!
!   qs_fluct_flag = 1 ! None = 0 ; Lattanzi = 1 ; Osnes = 2 
!
!
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_SetYdot
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 :: ppiclf_rcp_part, ppiclf_p0
      integer :: ppiclf_moveparticle, ierr
      CHARACTER(12) :: ppiclf_matname
      common /RFLU_ppiclf_misc01/ ppiclf_rcp_part
      common /RFLU_ppiclf_misc02/ ppiclf_matname
      common /RFLU_ppiclf_misc03/ ppiclf_p0, ppiclf_moveparticle
      real*8 fqsx, fqsy, fqsz
      real*8 fqsforce
      real*8 fqs_fluct(3)
      real*8 xi_par, xi_perp, xi_T
      real*8 famx, famy, famz 
      real*8 fdpdx, fdpdy, fdpdz
      real*8 fdpvdx, fdpvdy, fdpvdz
      real*8 fcx, fcy, fcz
      real*8 fbx, fby, fbz 
      real*8 fvux, fvuy, fvuz

      real*8 ug, vg, wg

      real*8 beta,cd

      real*8 factor, rcp_fluid, rmass_add

      real*8 gkern

! Needed for Added mass calculation
      integer*4 j, l
      real*8 SDrho, maxFilter

      real*8 vgradrhog
      integer*4 i, n, ic, k
      integer*4 store_forces

! Needed for heat transfer
      real*8 qq, rmass_therm, temp

! Needed for reactive particles
      integer*4 burnrate_model
      real*8    mdot_me, mdot_ox
      real*8    Pres

! Needed for angular velocity
      real*8 taux, tauy, tauz, rmass_omega,
     >       taux_hydro, tauy_hydro, tauz_hydro 
      real*8 tau
      real*8 liftx, lifty, liftz
      real*8 lift

! Finite Diff Material derivative Variables
      integer*4 nstage, istage
      integer*4 icallb
      save      icallb
      data      icallb /0/
      integer*4 idebug
      save      idebug
      data      idebug /0/

! Print Data to file
      LOGICAL I_EXIST 
      Character(LEN=25) str 
      integer*4 f_dump
      save      f_dump  
      data      f_dump /1/

      logical exist_file
!
!-----------------------------------------------------------------------
!   

      ! Avery added 10/10/2024 for subbin nearest neighbor search
      
      INTEGER*4 SBin_map( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_nndist)/ppiclf_nndist) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_nndist)/ppiclf_nndist)
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_nndist)/ppiclf_nndist) 
     >       + 1) - 1), (ppiclf_npart+ppiclf_npart_gp))
      INTEGER*4  SBin_counter( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_nndist)/ppiclf_nndist) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_nndist)/ppiclf_nndist)
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_nndist)/ppiclf_nndist) 
     >       + 1) - 1))
      INTEGER*4 i_Bin(3), n_SBin(3), tot_SBin


! Unit Test only Code:
!-----------------------------------------------------------------------
!
#ifdef TEST
      ! unit_test only tests the SB nearest neighbor search in this
      ! subroutine.  The full subroutine is called to ensure that
      ! the array initialization is correct.
      CALL ppiclf_user_subbinMap(i_Bin, n_SBin, tot_SBin 
     >                           ,SBin_counter ,SBin_map)
      DO i = 1,ppiclf_npart
        CALL ppiclf_solve_NearestNeighborSB(
     >         i,tot_SBin,SBin_counter,SBin_map,n_SBin,i_Bin)
      END DO
      RETURN
#endif

! 
!-----------------------------------------------------------------------
!
!
! Code:
!

      icallb = icallb + 1
      nstage = 3
      istage = mod(icallb,nstage)
      if (istage .eq. 0) istage = 3  

      ! Count every iStage=1 for debug output
      if (iStage .eq. 1) idebug = idebug + 1

      ! Print dt and time every time step
      if (ppiclf_nid==0) then
        if (istage .eq. 1) then
          write(6,'(a,2x,2(1pe14.6),2x,i3)') '*** PPICLF dt, time = ',
     >      ppiclf_dt,ppiclf_time
        endif
      endif

      burnrate_model = 0
      if (burnrate_flag .gt. 0) then
         if ( TRIM(ppiclf_matname)=='AL'.or. 
     >        TRIM(ppiclf_matname)=='Al' ) then
            burnrate_model = 1
         elseif ( TRIM(ppiclf_matname)=='Mg' ) then
            burnrate_model = 2
         elseif ( TRIM(ppiclf_matname)=='C' ) then
            burnrate_model = 3
         else
            print*,'Error: no burn rate model'
            stop
         endif
      endif

      rpi        = acos(-1.0d0)
      rcp_part   = ppiclf_rcp_part
      rpr        = 0.70d0
      rcp_fluid  = 1004.64d0

      fac = ppiclf_rk3ark(iStage)*ppiclf_dt
      if (1==2) then
         if (ppiclf_nid==0) print*,'dt,fac=',
     >      istage,ppiclf_dt,fac,
     >      stationary, qs_flag, am_flag, pg_flag,
     >      collisional_flag, heattransfer_flag, feedback_flag,
     >      qs_fluct_flag, ppiclf_debug, ppiclf_nTimeBH,
     >      ppiclf_nUnsteadyData
      endif

      OneThird = 1.0d0/3.0d0

!
!-----------------------------------------------------------------------
!
! Reapply axi-sym collision correction
! Right now hard coding smallest radius  
      do i=1,ppiclf_npart
        ppiclf_rprop(PPICLF_R_JDPe,i) = 
     > (0.00005/ppiclf_rprop(PPICLF_R_JSPT,i))
     > * ppiclf_rprop(PPICLF_R_JDP,i)  
      end do 
!
!-----------------------------------------------------------------------
!
! Reset arrays for Viscous Unsteady Force
!
      if (ViscousUnsteady_flag>=1) then
         call ppiclf_user_prop2plag
      endif
!
!-----------------------------------------------------------------------
!
! Avery added 10/10/2024 - Map particles to subbins if collisional force, 
! Briney Added Mass force, QS fluctation force, or pseudo turbulence is flagged
!
      if((am_flag==2).or.(collisional_flag>=1)
     >                  .or.(qs_fluct_flag>=1)
     >                .or.(pseudoTurb_flag==1)) then
        call ppiclf_user_subbinMap(i_Bin, n_SBin, tot_SBin 
     >                               ,SBin_counter ,SBin_map)
      endif ! Collisions, QS Fluct, Briney AM, or pseudoTurb flags on
!
!-----------------------------------------------------------------------
!
      ! Set initial max values - must be done npart loop
      if (ppiclf_debug >= 1) then
         phimax    = 0.d0
         fqsx_max  = 0.d0; fqsy_max  = 0.d0; fqsz_max  = 0.d0
         famx_max  = 0.d0; famy_max  = 0.d0; famz_max  = 0.d0
         fdpdx_max = 0.d0; fdpdy_max = 0.d0; fdpdz_max = 0.d0
         fcx_max   = 0.d0; fcy_max   = 0.d0; fcz_max   = 0.d0
         fvux_max  = 0.d0; fvuy_max  = 0.d0; fvuz_max  = 0.d0
         qq_max    = 0.d0;
         fqsx_fluct_max = 0.d0; fqsy_fluct_max = 0.d0
         fqsz_fluct_max = 0.d0
         fqsx_total_max = 0.d0; fqsy_total_max = 0.d0
         fqsz_total_max = 0.d0
         fqs_mag = 0.0d0; fam_mag = 0.0d0; fdp_mag = 0.0d0
         fc_mag  = 0.0d0
         umean_max = 0.d0; vmean_max = 0.d0; wmean_max = 0.d0
      endif


!
!-----------------------------------------------------------------------
!
!
! Evaluate ydot, the rhs of the equations of motion
! for the particles
!

      do i=1,ppiclf_npart

         ! Choose viscosity law
         if (rmu_flag==rmu_fixed_param) then
            ! Constant viscosity law
            rmu = rmu_ref
         elseif (rmu_flag==rmu_suth_param) then
            ! Sutherland law
            temp    = ppiclf_rprop(PPICLF_R_JT,i)
            rmu     = rmu_ref*sqrt(temp/tref)
     >                   *(1.0d0+suth/tref)/(1.0d0+suth/temp)
         else
             call ppiclf_exittr('Unknown viscosity law$', 0.0d0, 0)
         endif
         rkappa = rcp_fluid*rmu/rpr


         ! Useful values
         rmass  = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >              *ppiclf_rprop(PPICLF_R_JRHOP,i)
         vx     = ppiclf_rprop(PPICLF_R_JUX,i) - ppiclf_y(PPICLF_JVX,i)
         vy     = ppiclf_rprop(PPICLF_R_JUY,i) - ppiclf_y(PPICLF_JVY,i)
         vz     = ppiclf_rprop(PPICLF_R_JUZ,i) - ppiclf_y(PPICLF_JVZ,i)
         vmag   = sqrt(vx*vx + vy*vy + vz*vz)
         rhof   = ppiclf_rprop(PPICLF_R_JRHOF,i)
         dp     = ppiclf_rprop(PPICLF_R_JDP,i)
         rep    = vmag*dp*rhof/rmu
         rphip  = ppiclf_rprop(PPICLF_R_JPHIP,i)
         rphif  = 1.0d0-ppiclf_rprop(PPICLF_R_JPHIP,i)
         asndf  = ppiclf_rprop(PPICLF_R_JCS,i)
         rmachp = vmag/asndf
         rhop   = ppiclf_rprop(PPICLF_R_JRHOP,i)

         ! TLJ - 04/03/2025; Do not calculate forces if vmag = 0
         !       Otherwise the particles might move before the 
         !       shock arrives
         if (vmag <= 1.d-8) cycle

         ! 08/08/2025 - Thierry  - 1.d-8 is very small
         if (vmag <= 1.d-3) cycle
         !***
         ! Avery - This is planar-shock curtain specific.  Shouldn't be in main code
         !***
 
         ! Thierry - at initial times when rmachp and vmag are very
         ! small, we would get very large CD (>100)
         ! This leads to NaN projected force when PseudoTurbulence is 
         ! enabled. One fix I'm trying is to cycle when rmachp and vmag
         ! are small

         ! Thierry - seeing if that fixes the issue of very large CD
         ! initially 
         if(vmag .lt. 1.0 .or. rmachp .lt. 1.d-3) cycle
 
         ! TLJ - redefined rprop(PPICLF_R_JSPT,i) to be the particle
         !   velocity magnitude for plotting purposes - 01/03/2025
 
         ! ***This is bad practice and leads to hard to debug code -Avery***
         ppiclf_rprop(PPICLF_R_JSPT,i) = sqrt(
     >       ppiclf_y(PPICLF_JVX,i)**2 +
     >       ppiclf_y(PPICLF_JVY,i)**2 +
     >       ppiclf_y(PPICLF_JVZ,i)**2)

         rep = max(0.1d0,rep)

         ! Redefine volume fractions
         ! Need to make sure phi_p + phi_f = 1
         rphip = ppiclf_rprop(PPICLF_R_JPHIP,i)
         rphip = min(rphip,0.62d0)
         rphif = 1.0d0-rphip
         !*** AVERY - we should rethink this limit of 62% pVF ***

         ! TLJ: Needed for viscous unsteady force
         !      Using same nomenclature as rocinteract subroutines
         reyL = dp*vmag*rhof/rmu
         rnu = rmu/rhof

         ! Zero out for each particle i
         famx = 0.0d0; famy = 0.0d0; famz = 0.0d0; rmass_add = 0.0d0;
         Fam(1) = 0.0d0; Fam(2) = 0.0d0; Fam(3) = 0.0d0
         FamUnary(1)=0.0d0;FamUnary(2)=0.0d0;FamUnary(3)=0.0d0;
         FamBinary(1)=0.0d0;FamBinary(2)=0.0d0;FamBinary(3)=0.0d0;
         Wdot_neighbor_mean(1) = 0.0d0; Wdot_neighbor_mean(2) = 0.0d0;
         Wdot_neighbor_mean(3) = 0.0d0; nneighbors = 0.0d0
         fqsx = 0.0d0; fqsy = 0.0d0; fqsz = 0.0d0; beta = 0.0d0;
         fqs_fluct(1)=0.0d0;fqs_fluct(2)=0.0d0;fqs_fluct(3)=0.0d0;
         fdpdx = 0.0d0; fdpdy = 0.0d0; fdpdz = 0.0d0;
         fcx = 0.0d0; fcy = 0.0d0; fcz = 0.0d0;
         taux = 0.0d0; tauy = 0.0d0; tauz = 0.0d0;
         liftx = 0.0d0; lifty = 0.0d0; liftz = 0.0d0;
         fvux = 0.0d0; fvuy = 0.0d0; fvuz = 0.0d0;
         qq=0.0d0
         mdot_me = 0.0d0; mdot_ox = 0.0d0;
         upmean = 0.0; vpmean = 0.0; wpmean = 0.0;
         u2pmean = 0.0; v2pmean = 0.0; w2pmean = 0.0;
         fdpvdx = 0.0d0; fdpvdy = 0.0d0; fdpvdz = 0.0d0;
         !--- Added for PseudoTurbulence
         Rsg = 0.0d0; T_par = 0.0d0

!
! Step 1a: New Added-Mass model of Briney
!
         ! 07/15/2024 - If am_flag = 2, then we need to call
         !   the Unary term before we make any calls to nearest
         !   neighbor
         ! 06/05/2024 - Thierry - For each particle i, initialize
         ! variables to be used in nearest neighbors to zero
         ! before looping over particle j (j neq i)
         ! Briney Added Mass flag
         if (am_flag == 2) then 
            ! 07/14/24 - Thierry - If Briney Algorithm flag and fluct_flag
            !   are ON -> evaluate added-mass unary term before evaluating
            !   neighbor-induced acceleration in EvalNearestNeighbor
            call ppiclf_user_AM_Briney_Unary(i,iStage,
     >           famx,famy,famz,rmass_add)
         endif ! end am_flag = 2


!
! Step 1b: Call NearestNeighbor if particles i and j interact
!
         if ((am_flag==2).or.(collisional_flag>=1)
     >          .or.(qs_fluct_flag>=1)
     >          .or.(pseudoTurb_flag==1)) then

         !AVERY - we should fix vmag ~0 bug and remove conditional check
         if ((qs_fluct_flag>=1) .and. (vmag .gt. 1.d-8)) then
            ! Compute mean for particle i
            !    add neighbor particle j afterward
            ! Box filter is used if qs_fluct_filter_flag=0
            !   The box filter used here is a simple cube centered
            !     at particle i with half-width dist2 (see
            !     ppiclf_user_EvalNearestNeighbor.f for definition)
            !   We use a simple arithmetic mean
            !   phipmean is not used
            ! Gaussian filter is used if qs_fluct_filter_flag=1
            !   We use the value of the Gaussian times the volume
            !     of the particle to get the filtered particle volume

            if (qs_fluct_filter_flag==0) then
               ! box filter
               !***AVERY - I think phipmean is using wrong index ***
               ! Also, none of the below are means...
               phipmean = ppiclf_rprop(PPICLF_R_JVOLP,i)
               upmean   = ppiclf_y(PPICLF_JVX,i)
               vpmean   = ppiclf_y(PPICLF_JVY,i)
               wpmean   = ppiclf_y(PPICLF_JVZ,i)
               u2pmean  = upmean**2
               v2pmean  = vpmean**2
               w2pmean  = wpmean**2
               icpmean  = 1
            else if (qs_fluct_filter_flag==1) then
               ! gaussian kernel
               ! r = 0
               !***AVERY - verify max filter dimension should be used***
               maxFilter = MAX(ppiclf_filter(1),ppiclf_filter(2),
     >                         ppiclf_filter(3))
               gkern = sqrt(rpi*maxFilter**2/
     >                (4.0d0*log(2.0d0)))**(-ppiclf_ndim)
               !***AVERY - I think phipmean is using wrong index ***
               ! Also, none of the below are means...
               phipmean = gkern*ppiclf_rprop(PPICLF_R_JVOLP,i)
               upmean   = gkern*ppiclf_y(PPICLF_JVX,i)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               vpmean   = gkern*ppiclf_y(PPICLF_JVY,i)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               wpmean   = gkern*ppiclf_y(PPICLF_JVZ,i)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               u2pmean  = gkern*(ppiclf_y(PPICLF_JVX,i)**2)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               v2pmean  = gkern*(ppiclf_y(PPICLF_JVY,i)**2)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               w2pmean  = gkern*(ppiclf_y(PPICLF_JVZ,i)**2)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               icpmean = 1
            end if
         end if

         ! add neighbors
         ! AVERY - always use subbins
         IF ( 1 .EQ. 1 ) THEN !sbNearest_flag .EQ. 1) THEN
            CALL ppiclf_solve_NearestNeighborSB(
     >           i,tot_SBin,SBin_counter,SBin_map,n_SBin,i_Bin)
!         ELSE
!             CALL ppiclf_solve_NearestNeighbor(i)
         END IF

         end if ! end Step 1b; nearestneighbor


!
! Step 2: Force component quasi-steady
!
         if (qs_flag==1) then 
           call ppiclf_user_QS_Parmar(i,beta,cd)
         else if (qs_flag==2) then 
           call ppiclf_user_QS_Osnes (i,beta,cd)
         else if (qs_flag==3) then 
           call ppiclf_user_QS_ModifiedParmar(i,beta,cd)
         else if (qs_flag==4) then 
           call ppiclf_user_QS_Gidaspow(i,beta,cd)
         else
           print*, "***PPICLF: Error in QS Model Selection!"   
           call ppiclf_exittr('Wrong QS Model Choice$', 0.0, 0)
         endif

         fqsx = beta*vx
         fqsy = beta*vy
         fqsz = beta*vz

!
! Step 3: Force fluctuation for quasi-steady force
!
         ! Note: QS fluctuations needs nearest neighbors,
         !   and is called above in Step 1b
         if (qs_fluct_flag==1) then
            call ppiclf_user_QS_fluct_Lattanzi(i,iStage,fqs_fluct)
         elseif (qs_fluct_flag==2 .or. pseudoTurb_flag==1) then
            call ppiclf_user_QS_fluct_Osnes(i,iStage,fqs_fluct,
     >                                     xi_par,xi_perp,xi_T,
     >                                      fqsx, fqsy, fqsz)
         endif

         ! Add fluctuation part to quasi-steady force
         fqsx = fqsx + fqs_fluct(1)
         fqsy = fqsy + fqs_fluct(2)
         fqsz = fqsz + fqs_fluct(3)

         ! Store quasi-steady fluctuating force
         ppiclf_rprop(PPICLF_R_FLUCTFX,i) = fqs_fluct(1)
         ppiclf_rprop(PPICLF_R_FLUCTFY,i) = fqs_fluct(2)
         ppiclf_rprop(PPICLF_R_FLUCTFZ,i) = fqs_fluct(3)
         
         ! Store normally distributed random variables xi for PseudoTurbulence
         ppiclf_rprop(PPICLF_R_XIPAR,i)  = xi_par
         ppiclf_rprop(PPICLF_R_XIPERP,i) = xi_perp
         ppiclf_rprop(PPICLF_R_XIT,i)    = xi_T


!
! Step 4: Force component added mass
!
         if (am_flag == 1) then 
            call ppiclf_user_AM_Parmar(i,iStage,
     >           famx,famy,famz,rmass_add)

!-----------------------------------------------------------------------
!Thierry - Added Mass code continues here
         
         elseif (am_flag == 2) then 

            ! Thierry - binary_model.f90 evaluates the terms
            ! in the folllowing order:
            !   (1) Unary Term
            !   (2) Evaluates Neighbor Acceleration
            !   (3) Binary Term
            ! We replicate that here by calling them in the same order
            ! Unary and Binary calculations are now under 
            !   two separate subroutines
            ! Thierry - need to make sure NearestNeighbor is called
            !    if fluct_flag = 0 (ie, no QS fluctuations)

            ! Binary subroutine only valid when number of neighbors .gt. 0
            if (nneighbors .gt. 0) then
               call ppiclf_user_AM_Briney_Binary(i,iStage,
     >              famx,famy,famz,rmass_add)
               FamBinary(1) = famx - FamUnary(1)
               FamBinary(2) = famy - FamUnary(2)
               FamBinary(3) = famz - FamUnary(3)
            else
            ! if particle has no neighbors, need to multiply added mass forces
            ! by volume, as this is taken care of in Binary subroutine
               famx = famx*ppiclf_rprop(PPICLF_R_JVOLP,i)
               famy = famy*ppiclf_rprop(PPICLF_R_JVOLP,i)
               famz = famz*ppiclf_rprop(PPICLF_R_JVOLP,i)
            endif
         endif


!-----------------------------------------------------------------------

!
! Step 5: Force component pressure gradient
!
         if (pg_flag == 1) then
            fdpdx = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >               ppiclf_rprop(PPICLF_R_JDPDX,i)
            fdpdy = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >               ppiclf_rprop(PPICLF_R_JDPDY,i)
            fdpdz = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >               ppiclf_rprop(PPICLF_R_JDPDZ,i)

            if (flow_model == 1) then ! Navier-Stokes Flow Model
               fdpvdx = ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >                  ppiclf_rprop(PPICLF_R_JDPVDX,i)
               fdpvdy = ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >                  ppiclf_rprop(PPICLF_R_JDPVDY,i)
               fdpvdz = ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >                  ppiclf_rprop(PPICLF_R_JDPVDZ,i)
            endif ! flow_model

            fdpdx = fdpdx + fdpvdx
            fdpdy = fdpdy + fdpvdy
            fdpdz = fdpdz + fdpvdz
         endif ! end pg_flag = 1

!
! Step 6: Force component collisional force, ie, particle-particle
!
         if (collisional_flag >= 1) then
            ! Collision force:
            !  A discrete numerical model for granular assemblies
            !  - Cundall and Strack (1979)
            !  - Geotechnique

            ! Sam - STILL NEED TO VALIDATE COLLISION FORCE
            ! Sam - Step 1b already calls nearest neighbor
            
            fcx  = ppiclf_ydotc(PPICLF_JVX,i)
            fcy  = ppiclf_ydotc(PPICLF_JVY,i)
            fcz  = ppiclf_ydotc(PPICLF_JVZ,i) 

         endif ! collisional_flag >= 1


!
! Step 7: Viscous unsteady force with history kernel
!
         if (ViscousUnsteady_flag==1) then
            call ppiclf_user_VU_Rocflu(i,iStage,fvux,fvuy,fvuz)
         endif


!
! Step 8a: Combustion model for reactive particles
!
         rmass_therm = rmass*rcp_part
         qq = 0.0d0

         if (burnrate_flag >= 1) then
            call ppiclf_user_BR_driver(i,iStage,
     >         burnrate_model,qq,mdot_me,mdot_ox)
         endif
!
! Step 8b: Heat transfer model
!
         if (heattransfer_flag >= 1) then
            call ppiclf_user_HT_driver(i,qq)
         endif ! heattransfer_flag >= 1


!
! Step 9a: Angular velocity model
!
         rmass_omega = rmass*dp*dp/10.0d0

         if (collisional_flag >= 2) then
            taux  = ppiclf_ydotc(PPICLF_JOX,i)
            tauy  = ppiclf_ydotc(PPICLF_JOY,i)
            tauz  = ppiclf_ydotc(PPICLF_JOZ,i) 
            call ppiclf_user_Torque_driver(i,iStage,taux,tauy,tauz,
     >                            taux_hydro,tauy_hydro,tauz_hydro)
         endif ! collisional_flag >= 2

!
! Step 9b: Saffman and Magnus Lift models
!          Lift models requires gas-phase vorticity and
!          particle angular velocity
!
         if (collisional_flag == 4) then
            call ppiclf_user_Lift_driver(i,iStage,liftx,lifty,liftz)
         endif ! collisional_flag == 4


!
! Step 10: Set ydot for all PPICLF_SLN number of equations
!
         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JZ, i) = ppiclf_y(PPICLF_JVZ,i)
         ppiclf_ydot(PPICLF_JVX,i) = (fqsx+famx+fdpdx+fvux+liftx+fcx)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JVY,i) = (fqsy+famy+fdpdy+fvuy+lifty+fcy)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JVZ,i) = (fqsz+famz+fdpdz+fvuz+liftz+fcz)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JT,i)  = qq/rmass_therm
         ppiclf_ydot(PPICLF_JOX,i) = taux/rmass_omega
         ppiclf_ydot(PPICLF_JOY,i) = tauy/rmass_omega
         ppiclf_ydot(PPICLF_JOZ,i) = tauz/rmass_omega
         ppiclf_ydot(PPICLF_JMETAL,i)  = mdot_me
         ppiclf_ydot(PPICLF_JOXIDE,i)  = mdot_ox

!
! Update data for viscous unsteady case
!
         if (ViscousUnsteady_flag>=1) then
            call ppiclf_user_UpdatePlag(i)
         endif

!
! Step 11: Feed Back force to the gas phase
!
         !    Note that Rocflu uses a negative of the RHS, and
         !    so ppiclf must respect this odd convention.
         !
         ! Project work done by hydrodynamic forces:
         !   Inter-phase heat transfer and energy coupling in turbulent 
         !   dispersed multiphase flows
         !   - Ling et al. (2016)
         !   - Physics of Fluids
         ! See also for more details
         !   Explosive dispersal of particles in high speed environments
         !   - Durant et al. (2022)
         !   - Journal of Applied Physics


        IF(feedback_flag==0) THEN
          ppiclf_feedbk(PPICLF_P_JFX,i) = 0.0d0 
          ppiclf_feedbk(PPICLF_P_JFY,i) = 0.0d0 
          ppiclf_feedbk(PPICLF_P_JFZ,i) = 0.0d0 
          ppiclf_feedbk(PPICLF_P_JE,i)  = 0.0d0
        END IF

        IF(feedback_flag==1) THEN
          ! Momentum equations feedback terms
          ppiclf_feedbk(PPICLF_P_JFX,i) = ppiclf_rprop(PPICLF_R_JSPL,i)*
     >      (ppiclf_ydot(PPICLF_JVX,i)*rmass - fcx)
          ppiclf_feedbk(PPICLF_P_JFY,i) = ppiclf_rprop(PPICLF_R_JSPL,i)*
     >      (ppiclf_ydot(PPICLF_JVY,i)*rmass - fcy)
          ppiclf_feedbk(PPICLF_P_JFZ,i) = ppiclf_rprop(PPICLF_R_JSPL,i)*
     >      (ppiclf_ydot(PPICLF_JVZ,i)*rmass - fcz)

          ! Energy equation feedback term
          ! 09/19/2025 - Thierry - Added Lift force
          ! Still need to add Torue \cdot angular velocity
          ppiclf_feedbk(PPICLF_P_JE,i) = ppiclf_rprop(PPICLF_R_JSPL,i)
     >     * ( (fqsx+fvux+liftx)*ppiclf_y(PPICLF_JVX,i) + 
     >         (fqsy+fvuy+lifty)*ppiclf_y(PPICLF_JVY,i) + 
     >         (fqsz+fvuz+liftz)*ppiclf_y(PPICLF_JVZ,i) +
     >                famx*ppiclf_rprop(PPICLF_R_JUX,i) +
     >                famy*ppiclf_rprop(PPICLF_R_JUY,i) +
     >                famz*ppiclf_rprop(PPICLF_R_JUZ,i) +
     >                taux_hydro*ppiclf_y(PPICLF_JOX,i) +
     >                tauy_hydro*ppiclf_y(PPICLF_JOY,i) +
     >                tauz_hydro*ppiclf_y(PPICLF_JOZ,i) +
     >         qq )
          IF(pseudoTurb_flag==1) THEN
            ! 09/02/2025 -  Addition of PTKE to Rocflu's Energy Equation
            ppiclf_feedbk(PPICLF_P_JE,i) = ppiclf_rprop(PPICLF_R_JSPL,i)
     >       * ( (fqsx+fvux+famx+liftx)*ppiclf_y(PPICLF_JVX,i) + 
     >           (fqsy+fvuy+famy+lifty)*ppiclf_y(PPICLF_JVY,i) + 
     >           (fqsz+fvuz+famz+liftz)*ppiclf_y(PPICLF_JVZ,i) +
     >           qq )
          ELSE
            Rsg   = 0.0D0
            T_par = 0.0D0
          END IF ! pseudoTurb_flag
          ! 07/21/2025 - Thierry - Added Reynolds Subgrid Stress Feedback
          ppiclf_feedbk(PPICLF_P_JRSG11,i) = Rsg(1,1) 
     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
          ppiclf_feedbk(PPICLF_P_JRSG12,i) = Rsg(1,2) 
     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
          ppiclf_feedbk(PPICLF_P_JRSG13,i) = Rsg(1,3) 
     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
          ppiclf_feedbk(PPICLF_P_JRSG21,i) = Rsg(2,1) 
     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
          ppiclf_feedbk(PPICLF_P_JRSG22,i) = Rsg(2,2) 
     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
          ppiclf_feedbk(PPICLF_P_JRSG23,i) = Rsg(2,3) 
     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
          ppiclf_feedbk(PPICLF_P_JRSG31,i) = Rsg(3,1) 
     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
          ppiclf_feedbk(PPICLF_P_JRSG32,i) = Rsg(3,2) 
     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
          ppiclf_feedbk(PPICLF_P_JRSG33,i) = Rsg(3,3)  
     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)

          ppiclf_feedbk(PPICLF_P_JTSG1,i) = T_par(1) 
     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
          ppiclf_feedbk(PPICLF_P_JTSG2,i) = T_par(2) 
     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
          ppiclf_feedbk(PPICLF_P_JTSG3,i) = T_par(3) 
     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)

        END IF ! Feedback flag

        ! Update volume fraction feedback quantities with feedback on or off
        ppiclf_feedbk(PPICLF_P_JPHIP,i) = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >   *ppiclf_rprop(PPICLF_R_JSPL,i)
        ppiclf_feedbk(PPICLF_P_JPHIPD,i) =
     >   ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JRHOP,i)
        ppiclf_feedbk(PPICLF_P_JPHIPU,i) = 
     >   ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JVX,i)
        ppiclf_feedbk(PPICLF_P_JPHIPV,i) =
     >   ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JVY,i)
        ppiclf_feedbk(PPICLF_P_JPHIPW,i) = 
     >   ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JVZ,i)
        ppiclf_feedbk(PPICLF_P_JPHIPT,i) =
     >   ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JT,i)


! Step 12: If stationary, don't move particles. Feedback can still be on
! though.
!
         if (stationary .gt. 0) then
            if (stationary==1) then
               ppiclf_ydot(PPICLF_JX ,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JY ,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JZ, i)  = 0.0d0
               ppiclf_ydot(PPICLF_JVX,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JVY,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JVZ,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JT,i)   = 0.0d0
               ppiclf_ydot(PPICLF_JOX,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JOY,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JOZ,i)  = 0.0d0
            else
               call ppiclf_exittr('Unknown stationary flag$', 0.0d0, 0)
            endif
         elseif(stationary .lt. 0) then
            call ppiclf_user_unit_tests(i,iStage,famx,famy,famz)
            ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
            ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
            ppiclf_ydot(PPICLF_JZ, i) = ppiclf_y(PPICLF_JVZ,i)
            ppiclf_ydot(PPICLF_JVX,i) = ppiclf_ydot(PPICLF_JVX,i)+famx
            ppiclf_ydot(PPICLF_JVY,i) = ppiclf_ydot(PPICLF_JVY,i)+famy
            ppiclf_ydot(PPICLF_JVZ,i) = ppiclf_ydot(PPICLF_JVZ,i)+famz
            ppiclf_ydot(PPICLF_JT,i)  = 0.0d0
            ppiclf_ydot(PPICLF_JOX,i) = 0.0d0
            ppiclf_ydot(PPICLF_JOY,i) = 0.0d0
            ppiclf_ydot(PPICLF_JOZ,i) = 0.0d0
         endif



!
! Step 13: Store forces

         IF(PPICLF_LRP5 .EQ. 19) THEN
           ppiclf_rprop5(PPICLF_R_FQSX,i)  = fqsx
           ppiclf_rprop5(PPICLF_R_FQSY,i)  = fqsy
           ppiclf_rprop5(PPICLF_R_FQSZ,i)  = fqsz
           ppiclf_rprop5(PPICLF_R_FAMX,i)  = famx - rmass_add*
     >                                       ppiclf_ydot(PPICLF_JVX,i)
           ppiclf_rprop5(PPICLF_R_FAMY,i)  = famy - rmass_add*
     >                                       ppiclf_ydot(PPICLF_JVY,i)
           ppiclf_rprop5(PPICLF_R_FAMZ,i)  = famz - rmass_add*
     >                                       ppiclf_ydot(PPICLF_JVZ,i)
           ppiclf_rprop5(PPICLF_R_FAMBX,i) = FamBinary(1)
           ppiclf_rprop5(PPICLF_R_FAMBY,i) = FamBinary(2)
           ppiclf_rprop5(PPICLF_R_FAMBZ,i) = FamBinary(3)
           ppiclf_rprop5(PPICLF_R_FCX,i)   = fcx
           ppiclf_rprop5(PPICLF_R_FCY,i)   = fcy
           ppiclf_rprop5(PPICLF_R_FCZ,i)   = fcz
           ppiclf_rprop5(PPICLF_R_FVUX,i)  = fvux
           ppiclf_rprop5(PPICLF_R_FVUY,i)  = fvuy
           ppiclf_rprop5(PPICLF_R_FVUZ,i)  = fvuz
           ppiclf_rprop5(PPICLF_R_QQ,i)    = qq
           ppiclf_rprop5(PPICLF_R_FPGX,i)  = fdpdx
           ppiclf_rprop5(PPICLF_R_FPGY,i)  = fdpdy
           ppiclf_rprop5(PPICLF_R_FPGZ,i)  = fdpdz
         END IF
!
! Step 14: If debug mode is ON, calculate and print the max values.
!          The user should not have this ON for production runs.
!
         if (ppiclf_debug .ge. 1) then
            phimax = max(phimax,abs(rphip))

            fqsx_max = max(fqsx_max,abs(fqsx))
            fqsy_max = max(fqsy_max,abs(fqsy))
            fqsz_max = max(fqsz_max,abs(fqsz))
            fqs_mag  = max(fqs_mag,
     >                 sqrt(fqsx*fqsx+fqsy*fqsy+fqsz*fqsz))

            fqsx_fluct_max = max(fqsx_fluct_max, abs(fqs_fluct(1)))
            fqsy_fluct_max = max(fqsy_fluct_max, abs(fqs_fluct(2)))
            fqsz_fluct_max = max(fqsz_fluct_max, abs(fqs_fluct(3)))

            fqsx_total_max = max(fqsx_total_max, abs(fqsx))
            fqsy_total_max = max(fqsy_total_max, abs(fqsy))
            fqsz_total_max = max(fqsz_total_max, abs(fqsz))

            umean_max = max(umean_max, abs(upmean))
            vmean_max = max(vmean_max, abs(vpmean))
            wmean_max = max(wmean_max, abs(wpmean))

            famx_max = max(famx_max,abs(famx))
            famy_max = max(famy_max,abs(famy))
            famz_max = max(famz_max,abs(famz))
            fam_mag  = max(fam_mag,
     >                 sqrt(famx*famx+famy*famy+famz*famz))

            fdpdx_max = max(fdpdx_max,abs(fdpdx))
            fdpdy_max = max(fdpdy_max,abs(fdpdy))
            fdpdz_max = max(fdpdz_max,abs(fdpdz))
            fdp_mag   = max(fdp_mag,sqrt(fdpdx*fdpdx+fdpdy*fdpdy
     >                  +fdpdz*fdpdz))

            fcx_max = max(fcx_max, abs(fcx))
            fcy_max = max(fcy_max, abs(fcy))
            fcz_max = max(fcz_max, abs(fcz))
            fc_mag  = max(fc_mag,sqrt(fcx*fcx+fcy*fcy+fcz*fcz))

            fvux_max = max(fvux_max, abs(fvux))
            fvuy_max = max(fvuy_max, abs(fvuy))
            fvuz_max = max(fvuz_max, abs(fvuz))

            qq_max = max(qq_max, abs(qq))
 
            tau = sqrt(taux*taux + tauy*tauy + tauz*tauz)
            tau_max = max(tau_max, abs(tau))

            lift = sqrt(liftx**2 + lifty**2 + liftz**2)
            lift_max = max(lift_max,lift)

            if (ppiclf_debug.eq.2 .and. ppiclf_nid.eq.0) then
               if (iStage==3) then
                  if (i==1) then
                     write(7010,*) i,ppiclf_time,rmass,vmag,rhof,dp,
     >                rep,rphip,rphif,rmachp,rhop,rhoMixt,reyL,
     >             rmu,rnu,rkappa
                  endif
                  if (i==ppiclf_npart) then
                     write(7011,*) i,ppiclf_time,rmass,vmag,rhof,dp,
     >                rep,rphip,rphif,rmachp,rhop,rhoMixt,reyL,
     >                rmu,rnu,rkappa
                  endif
               endif
            endif

         endif ! ppiclf_debug .ge. 1
          
         ! write out for debug
         if (ppiclf_debug==3) then
         if (ppiclf_nid==0 .and. iStage==1) then
         if (mod(idebug,1)==0) then
            if (i<=5) then
               write(7020+i,*) i, ppiclf_time, rhof,
     >             ppiclf_rprop(PPICLF_R_JSDRX,i),                   
     >             ppiclf_rprop(PPICLF_R_JSDRY,i), 
     >             ppiclf_rprop(PPICLF_R_JSDRZ,i),
     >             ppiclf_ydot(PPICLF_JVX,i),
     >             ppiclf_ydot(PPICLF_JVY,i),
     >             ppiclf_ydot(PPICLF_JVZ,i),
     >             ppiclf_y(PPICLF_JVX,i),
     >             ppiclf_y(PPICLF_JVY,i),
     >             ppiclf_y(PPICLF_JVZ,i),
     >             ppiclf_y(PPICLF_JOX,i),
     >             ppiclf_y(PPICLF_JOY,i),
     >             ppiclf_y(PPICLF_JOZ,i)

               write(7040+i,*) i, ppiclf_time, 
     >              ppiclf_rprop(PPICLF_R_JSDRX:PPICLF_R_JSDRZ,i), ! Du/Dt
     >              ppiclf_rprop(PPICLF_R_JSDOX:PPICLF_R_JSDOZ,i)  ! DOmega/Dt

               write(7050+i,*) i, ppiclf_time, 
     >              fqs_mag,fam_mag,fdp_mag,fc_mag,tau_max

               write(7060+i,*) i, ppiclf_time, 
     >              fcx,fcy,fcz,
     >              liftx,lifty,liftz,
     >              taux,tauy,tauz
            endif
         endif
         endif
         endif



      enddo ! do i=1,ppiclf_npart

!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!06/05/2024 - Thierry - Store density-weighted acceleration
!
      ! Briney Added Mass flag
      if (am_flag==2) then 
         do i=1,ppiclf_npart
     
            ! Substantial derivative of density - how rocflu does it 
            SDrho = ppiclf_rprop(PPICLF_R_JRHSR,i)
     >         + ppiclf_y(PPICLF_JVX,i) * ppiclf_rprop(PPICLF_R_JPGCX,i)
     >         + ppiclf_y(PPICLF_JVY,i) * ppiclf_rprop(PPICLF_R_JPGCY,i)
     >         + ppiclf_y(PPICLF_JVZ,i) * ppiclf_rprop(PPICLF_R_JPGCZ,i)
            
            ! material derivative is phi weighted in Rocflu
            ! drho/dt
            SDrho = SDrho / (rphif)  
            vgradrhog = vx * ppiclf_rprop(PPICLF_R_JRHOGX,i) +
     >                  vy * ppiclf_rprop(PPICLF_R_JRHOGY,i) +
     >                  vz * ppiclf_rprop(PPICLF_R_JRHOGZ,i)
      
            ! Fluid density
            rhof   = ppiclf_rprop(PPICLF_R_JRHOF,i)

            vx = ppiclf_rprop(PPICLF_R_JUX,i) - ppiclf_y(PPICLF_JVX,i)
            vy = ppiclf_rprop(PPICLF_R_JUY,i) - ppiclf_y(PPICLF_JVY,i)
            vz = ppiclf_rprop(PPICLF_R_JUZ,i) - ppiclf_y(PPICLF_JVZ,i)
            ug = ppiclf_rprop(PPICLF_R_JUX,i)
            vg = ppiclf_rprop(PPICLF_R_JUY,i)
            wg = ppiclf_rprop(PPICLF_R_JUZ,i)
            ! Unary added mass solves rho^g d(u^p)/dt implicitly
            ! Binary added mass solves it explicitly and not implicitly
            ! WDOTX = D(rho^g u^g)/Dt - d(rho^g u^p)/dt)
            ! X-acceleration
            ppiclf_rprop(PPICLF_R_WDOTX,i) =
     >                vx*SDrho + rhof*ppiclf_rprop(PPICLF_R_JSDRX,i)
     >              + ug*vgradrhog
     >              - rhof*ppiclf_ydot(PPICLF_JVX,i)
          
            ! Y-acceleration
            ppiclf_rprop(PPICLF_R_WDOTY,i) =
     >                vy*SDrho + rhof*ppiclf_rprop(PPICLF_R_JSDRY,i)
     >              + vg*vgradrhog 
     >              - rhof*ppiclf_ydot(PPICLF_JVY,i)

            ! Z-acceleration
            ppiclf_rprop(PPICLF_R_WDOTZ,i) =
     >                vz*SDrho + rhof*ppiclf_rprop(PPICLF_R_JSDRZ,i)
     >              + wg*vgradrhog
     >              - rhof*ppiclf_ydot(PPICLF_JVZ,i)
          
            ! write out for debug
            if (ppiclf_debug==2) then
            if (ppiclf_nid==0 .and. iStage==1) then
            if (mod(idebug,10)==0) then
               if (i<=3) then
                  write(7020+i,*) i, ppiclf_time, rhof,
     >                ppiclf_rprop(PPICLF_R_JSDRX,i),                   
     >                ppiclf_rprop(PPICLF_R_JSDRY,i), 
     >                ppiclf_rprop(PPICLF_R_JSDRZ,i),
     >                ppiclf_ydot(PPICLF_JVX,i),
     >                ppiclf_ydot(PPICLF_JVY,i),
     >                ppiclf_ydot(PPICLF_JVZ,i),
     >                ppiclf_y(PPICLF_JVX,i),
     >                ppiclf_y(PPICLF_JVY,i),
     >                ppiclf_y(PPICLF_JVZ,i)

                  write(7030+i,*) i, ppiclf_time, 
     >              ppiclf_rprop(PPICLF_R_WDOTX:PPICLF_R_WDOTZ,i)

               endif
            endif
            endif
            endif

         enddo
      endif

!
! ----------------------------------------------------------------------
!

      ! Use ppiclf ALLREDUCE to compute values across processors
      ! Note that ALLREDUCE uses MPI_BARRIER, which is cpu expensive
      ! Print out every 10th iStage=1 counts
      if (ppiclf_debug   .ge. 1) then
      if (iStage         .eq. 1) then
      if (mod(idebug,10) .eq. 0) then
         call ppiclf_user_debug
      endif
      endif
      endif

!
! ----------------------------------------------------------------------
!
!

      !
      ! Reset arrays for Viscous Unsteady Force
      !
      if (ViscousUnsteady_flag>=1) then
         if (iStage==3) call ppiclf_user_ShiftUnsteadyData
         call ppiclf_user_plag2prop
      endif


! ----------------------------------------------------------------------

      return
      end
!-----------------------------------------------------------------------
!
! Created Oct. 18, 2024
!
! Subroutine to map both real and ghost particles to subbins
! for nearest neighbor search
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ppiclf_user_subbinMap(i_Bin, n_SBin, tot_SBin,
     >                                  SBin_counter, SBin_map)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4  SBin_map( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_nndist)/ppiclf_nndist) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_nndist)/ppiclf_nndist)
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_nndist)/ppiclf_nndist) 
     >       + 1) - 1), (ppiclf_npart+ppiclf_npart_gp))
      INTEGER*4  SBin_counter( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_nndist)/ppiclf_nndist) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_nndist)/ppiclf_nndist)
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_nndist)/ppiclf_nndist) 
     >       + 1) - 1))
      INTEGER*4  i_Bin(3), n_SBin(3), tot_SBin
     >          
!
! Internal:
!
      REAL*8    xp(3), bin_xMin(3)
      INTEGER*4 temp_SBin, i_SBin(3),ibinTemp, i, j, k, l 

!
! Code:
!
        ! Determine ppiclf bin in each dimension for this processor
        ! All real particles are in the same bin.  Look at 1st r particle
        DO l = 1,3
            i_Bin(l) = FLOOR((ppiclf_y(l,1) - ppiclf_binb(2*l-1))
     >                 /ppiclf_bins_dx(l))
            bin_xMin(l) = ppiclf_binb(2*l-1)+i_Bin(l)*ppiclf_bins_dx(l) 
        END DO 
        ! Determine the number of subbins in each dimension
        DO l = 1,3
          IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
            n_SBin(l) = FLOOR((ppiclf_bins_dx(l)+2*ppiclf_nndist)
     >                       /ppiclf_nndist) + 1
          ELSE
            n_SBin(l) = 0
          END IF
        END DO

        ! Determine total number of subbins
        tot_SBin = n_SBin(1)*n_SBin(2)
        IF (ppiclf_ndim .EQ. 3) THEN
          tot_SBin = tot_SBin*n_SBin(3)
        END IF
        ! Assign Subbin counters to 0
        SBin_counter = 0

        ! Map each real particle to a subbin
        DO i = 1 , ppiclf_npart
           DO l = 1,3
              IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
                 xp(l) = ppiclf_y(l,i)
              ELSE
                 xp(l) = 0.0
              END IF
           END DO 
           ! Determine subbin
           DO l = 1,3
              IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
                  i_SBin(l) = FLOOR((xp(l) - (bin_xMin(l) 
     >            - ppiclf_nndist))/ppiclf_nndist) 
              ELSE
                 i_SBin(l) = 0
              END IF
           END DO
           temp_SBin = i_SBin(1) + n_SBin(1)*i_SBin(2) +
     >                n_SBin(1)*n_SBin(2)*i_SBin(3)
           SBin_counter(temp_SBin) = SBin_counter(temp_SBin) + 1
           SBin_map(temp_SBin,SBin_counter(temp_SBin)) = i
        END DO ! real particle loop


        ! Map each ghost particle to a subbin
        DO i = 1 , ppiclf_npart_gp
          DO l = 1,3
            IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
              xp(l) = ppiclf_rprop_gp(l,i)
            ELSE
              xp(l) = 0.0
            END IF
          END DO
          ! Only map ghost particles within one neighborwidth
          ! from bin edge to subbins. All others are outside
          ! of collision search distance.
          IF (xp(1) .GT. (bin_xMin(1)-ppiclf_nndist)
     >  .AND. xp(2) .GT. (bin_xMin(2)-ppiclf_nndist)
     >  .AND. xp(3) .GT. (bin_xMin(3)-ppiclf_nndist)
     >  .AND. xp(1) .LT. (bin_xMin(1)+ppiclf_bins_dx(1)+ppiclf_nndist)
     >  .AND. xp(2) .LT. (bin_xMin(2)+ppiclf_bins_dx(2)+ppiclf_nndist)
     >  .AND. xp(3) .LT. (bin_xMin(3)+ppiclf_bins_dx(3)+ppiclf_nndist)
     >        ) THEN
            ! Determine subbin
            DO l = 1,3
              IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
                i_SBin(l) = FLOOR((xp(l) - (bin_xMin(l) 
     >          - ppiclf_nndist))/ppiclf_nndist) 
              ELSE
                i_SBin(l) = 0
              END IF
            END DO
            temp_SBin = i_SBin(1) + n_SBin(1)*i_SBin(2) +
     >                 n_SBin(1)*n_SBin(2)*i_SBin(3)
            SBin_counter(temp_SBin) = SBin_counter(temp_SBin) + 1
            ! negative in subbin map means it is ghost particle
            SBin_map(temp_SBin,SBin_counter(temp_SBin)) = -i
          END IF 
        END DO ! gp loop

      RETURN
      END 
!-----------------------------------------------------------------------
!
! Created Spet. 14, 2024
!
! Subroutine for unit test problems
!      
! stationary = -1: azimuthal velocity only
!            = -2: radial velocity only
!            = -3: radial + azimuthal velocity
!
! NOTE: Time step ppiclf_dt is hardcoded in ppiclf_solve_IntegrateRK3s_Rocflu
!       subroutine to be ppiclf_dt = 5.0000000000000004E-008      
!-----------------------------------------------------------------------
      subroutine ppiclf_user_unit_tests(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz
!
! Code:
!
      if (stationary==-1) then
         call ppiclf_user_unit01(i,iStage,famx,famy,famz)
      elseif (stationary==-2) then
         call ppiclf_user_unit02(i,iStage,famx,famy,famz)
      elseif (stationary==-3) then
         call ppiclf_user_unit03(i,iStage,famx,famy,famz)
      endif

      return
      end
!      
!-----------------------------------------------------------------------
!
! Azimuthal(theta) velocity only
!
      subroutine ppiclf_user_unit01(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz
      real*8 pi, omg00, omg02, tau

!
! Code:
!
      pi = acos(-1.0d0)
      tau = 0.1d-3 ! time for 2pi revolution. Can be max simulation time
      omg00 = (2.0d0*pi)/tau
      omg02 = omg00*omg00

      if(ppiclf_time .eq. 0.0) then
        ppiclf_y(PPICLF_JVX,i) = -omg00*ppiclf_y(PPICLF_JY,i)
        ppiclf_y(PPICLF_JVY,i) =  omg00*ppiclf_y(PPICLF_JX,i)
        ppiclf_y(PPICLF_JVZ,i) =  0.0d0
      endif
      
      famx = -omg02*ppiclf_y(PPICLF_JX,i)
      famy = -omg02*ppiclf_y(PPICLF_JY,i)
      famz = 0.0d0

      return
      end
!
!-----------------------------------------------------------------------
!
! Radial velocity only
!
      subroutine ppiclf_user_unit02(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz
      real*8 drdt00, drdt02, theta
!
! Code:
!
      ! dx/dt = dr/dt * cos theta
      ! dy/dt = dr/dt * sin theta
      drdt00 = 79.928163327808903
      drdt02 = drdt00*drdt00
      
      ! angle between particle and +ve x-axis
      theta = atan2(ppiclf_y(PPICLF_JY,i), ppiclf_y(PPICLF_JX,i))

      if(ppiclf_time .eq. 0.0) then
        ppiclf_y(PPICLF_JVX,i) =  -drdt00*cos(theta)
        ppiclf_y(PPICLF_JVY,i) =  -drdt00*sin(theta)
        ppiclf_y(PPICLF_JVZ,i) =  0.0d0
      endif
      
      famx = -drdt02*cos(theta)
      famy = -drdt02*sin(theta)
      famz = 0.0d0


      return
      end
!
!-----------------------------------------------------------------------
!
! Radial + Azimuthal velocity
!
      subroutine ppiclf_user_unit03(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz
      real*8 pi, tau, omg00, omg02, drdt00, theta
!
! Code:
!
      ! dx/dt = dr/dt * cos theta
      ! dy/dt = dr/dt * sin theta
      
      pi = acos(-1.0d0)

      tau = 0.1d-3
      omg00 = (2.0d0*pi)/tau
      omg02 = omg00*omg00
      drdt00 = 79.928163327808903
      
      ! angle between particle and +ve x-axis
      theta = atan2(ppiclf_y(PPICLF_JY,i), ppiclf_y(PPICLF_JX,i))

      if(ppiclf_time .eq. 0.0) then
        ppiclf_y(PPICLF_JVX,i) = -omg00*ppiclf_y(PPICLF_JY,i)
     >                           -drdt00*cos(theta)
        ppiclf_y(PPICLF_JVY,i) =  omg00*ppiclf_y(PPICLF_JX,i) 
     >                           -drdt00*sin(theta)
        ppiclf_y(PPICLF_JVZ,i) =  0.0d0
      endif
      
      famx = 2.0d0*drdt00*omg00*sin(theta) 
     >       -omg02*ppiclf_y(PPICLF_JX,i)
      famy = -2.0d0*drdt00*omg00*cos(theta) 
     >       -omg02*ppiclf_y(PPICLF_JY,i)
      famz = 0.0d0

      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for quasi-steady force
!
! Quasi-steady force (Re_p and Ma_p corrections):
!   Improved Drag Correlation for Spheres and Application 
!   to Shock-Tube Experiments 
!   - Parmar et al. (2010)
!   - AIAA Journal
!
! Quasi-steady force (phi corrections):
!   The Added Mass, Basset, and Viscous Drag Coefficients 
!   in Nondilute Bubbly Liquids Undergoing Small-Amplitude 
!   Oscillatory Motion
!   - Sangani et al. (1991)
!   - Phys. Fluids A
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_Parmar(i,beta)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 gamma
      real*8 rcd1,rmacr,rcd_mcr,rcd_std,rmach_rat,rcd_M1,
     >   rcd_M2,C1,C2,C3,f1M,f2M,f3M,lrep,factor,cd,beta,phi_corr

!
! Code:
!
      gamma = 1.4d0
      mp  = dmax1(rmachp,0.01d0)
      phi = dmax1(rphip,0.0001d0)
      re  = dmax1(rep,0.1d0)

      if(re .lt. 1E-14) then
         rcd1 = 1.0
      else 
         rmacr= 0.6 ! Critical rmachp no.
         rcd_mcr = (1.+0.15*re**(0.684)) + 
     >                (re/24.0)*(0.513/(1.+483./re**(0.669)))
       if (mp .le. rmacr) then
          rcd_std = (1.+0.15*re**(0.687)) + 
     >                (re/24.0)*(0.42/(1.+42500./re**(1.16)))
          rmach_rat = mp/rmacr
          rcd1 = rcd_std + (rcd_mcr - rcd_std)*rmach_rat
       else if (mp .le. 1.0) then
         rcd_M1 = (1.0+0.118*re**0.813) +
     >                (re/24.0)*0.69/(1.0+3550.0/re**.793)
         C1 =  6.48
         C2 =  9.28
         C3 = 12.21
         f1M = -1.884 +8.422*mp -13.70*mp**2 +8.162*mp**3
         f2M = -2.228 +10.35*mp -16.96*mp**2 +9.840*mp**3
         f3M =  4.362 -16.91*mp +19.84*mp**2 -6.296*mp**3
         lrep = log(re)
         factor = f1M*(lrep-C2)*(lrep-C3)/((C1-C2)*(C1-C3))
     >              +f2M*(lrep-C1)*(lrep-C3)/((C2-C1)*(C2-C3))
     >              +f3M*(lrep-C1)*(lrep-C2)/((C3-C1)*(C3-C2)) 
         rcd1 = rcd_mcr + (rcd_M1-rcd_mcr)*factor
       else if (mp .lt. 1.75) then
         rcd_M1 = (1.0+0.118*re**0.813) +
     >              (re/24.0)*0.69/(1.0+3550.0/re**.793)
         rcd_M2 = (1.0+0.107*re**0.867) +
     >              (re/24.0)*0.646/(1.0+861.0/re**.634)
         C1 =  6.48
         C2 =  8.93
         C3 = 12.21
         f1M = -2.963 +4.392*mp -1.169*mp**2 -0.027*mp**3
     >             -0.233*exp((1.0-mp)/0.011)
         f2M = -6.617 +12.11*mp -6.501*mp**2 +1.182*mp**3
     >             -0.174*exp((1.0-mp)/0.010)
         f3M = -5.866 +11.57*mp -6.665*mp**2 +1.312*mp**3
     >             -0.350*exp((1.0-mp)/0.012)
         lrep = log(re)
         factor = f1M*(lrep-C2)*(lrep-C3)/((C1-C2)*(C1-C3))
     >              +f2M*(lrep-C1)*(lrep-C3)/((C2-C1)*(C2-C3))
     >              +f3M*(lrep-C1)*(lrep-C2)/((C3-C1)*(C3-C2)) 
         rcd1 = rcd_M1 + (rcd_M2-rcd_M1)*factor
       else
         rcd1 = (1.0+0.107*re**0.867) +
     >                (re/24.0)*0.646/(1.0+861.0/re**.634)
       end if ! mp
      endif    ! re

      ! Sangani's volume fraction correction for dilute random arrays
      ! Capping volume fraction at 0.5
      phi_corr = (1.0+5.94*min(rphip,0.5))
      
      cd = (24.0/re)*rcd1*phi_corr

      beta = rcd1*3.0*rpi*rmu*dp

      beta = beta*phi_corr

      return
      end
!-----------------------------------------------------------------------
! Created Feb. 1, 2024
! Modified July 1, 2025
!
! Subroutine for quasi-steady force
!
! Quasi-steady force (Re_p and Ma_p corrections):
!   Improved Drag Correlation for Spheres and Application 
!   to Shock-Tube Experiments 
!   - Parmar et al. (2010)
!   - AIAA Journal
!      
! Quasi-steady force (phi corrections):
!   Sangani et al. (1991) volume fraction correction overshoots 
!   the drag coefficient. 
!      
!   We adopt instead Osnes et al. (2023) volume fraction correction
!   based on Tenneti et al. with one extra term. 
!      
!   At Mach=0, the drag coefficient from this subroutine matches very
!   well with the one calculated using the Osnes subroutine, for various
!   Reynolds numbers and volume fractions. 
!
!-----------------------------------------------------------------------
      subroutine ppiclf_user_QS_ModifiedParmar(i,beta)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 gamma
      real*8 rcd1,rmacr,rcd_mcr,rcd_std,rmach_rat,rcd_M1,
     >   rcd_M2,C1,C2,C3,f1M,f2M,f3M,lrep,factor,cd,beta,phi_corr,
     >   b1,b2,b3

!
! Code:
!
      gamma = 1.4d0
      mp  = dmax1(rmachp,0.01d0)
      phi = dmax1(rphip,0.0001d0)
      re  = dmax1(rep,0.1d0)

      if(re .lt. 1E-14) then
         rcd1 = 1.0
      else 
         rmacr= 0.6 ! Critical rmachp no.
         rcd_mcr = (1.+0.15*re**(0.684)) + 
     >                (re/24.0)*(0.513/(1.+483./re**(0.669)))
       if (mp .le. rmacr) then
          rcd_std = (1.+0.15*re**(0.687)) + 
     >                (re/24.0)*(0.42/(1.+42500./re**(1.16)))
          rmach_rat = mp/rmacr
          rcd1 = rcd_std + (rcd_mcr - rcd_std)*rmach_rat
       else if (mp .le. 1.0) then
         rcd_M1 = (1.0+0.118*re**0.813) +
     >                (re/24.0)*0.69/(1.0+3550.0/re**.793)
         C1 =  6.48
         C2 =  9.28
         C3 = 12.21
         f1M = -1.884 +8.422*mp -13.70*mp**2 +8.162*mp**3
         f2M = -2.228 +10.35*mp -16.96*mp**2 +9.840*mp**3
         f3M =  4.362 -16.91*mp +19.84*mp**2 -6.296*mp**3
         lrep = log(re)
         factor = f1M*(lrep-C2)*(lrep-C3)/((C1-C2)*(C1-C3))
     >              +f2M*(lrep-C1)*(lrep-C3)/((C2-C1)*(C2-C3))
     >              +f3M*(lrep-C1)*(lrep-C2)/((C3-C1)*(C3-C2)) 
         rcd1 = rcd_mcr + (rcd_M1-rcd_mcr)*factor
       else if (mp .lt. 1.75) then
         rcd_M1 = (1.0+0.118*re**0.813) +
     >              (re/24.0)*0.69/(1.0+3550.0/re**.793)
         rcd_M2 = (1.0+0.107*re**0.867) +
     >              (re/24.0)*0.646/(1.0+861.0/re**.634)
         C1 =  6.48
         C2 =  8.93
         C3 = 12.21
         f1M = -2.963 +4.392*mp -1.169*mp**2 -0.027*mp**3
     >             -0.233*exp((1.0-mp)/0.011)
         f2M = -6.617 +12.11*mp -6.501*mp**2 +1.182*mp**3
     >             -0.174*exp((1.0-mp)/0.010)
         f3M = -5.866 +11.57*mp -6.665*mp**2 +1.312*mp**3
     >             -0.350*exp((1.0-mp)/0.012)
         lrep = log(re)
         factor = f1M*(lrep-C2)*(lrep-C3)/((C1-C2)*(C1-C3))
     >              +f2M*(lrep-C1)*(lrep-C3)/((C2-C1)*(C2-C3))
     >              +f3M*(lrep-C1)*(lrep-C2)/((C3-C1)*(C3-C2)) 
         rcd1 = rcd_M1 + (rcd_M2-rcd_M1)*factor
       else
         rcd1 = (1.0+0.107*re**0.867) +
     >                (re/24.0)*0.646/(1.0+861.0/re**.634)
       end if ! mp
      endif    ! re

      ! Osnes's volume fraction correction
      b1 = 5.81*phi/((1.0-phi)**2) + 
     >     0.48*(phi**(1.d0/3.d0))/((1.0-phi)**3)

      b2 = ((1.0-phi)**2)*(phi**3)*
     >     re*(0.95+0.61*(phi**3)/((1.0-phi)*2))

      b3 = dmin1(sqrt(20.0d0*mp),1.0d0)*
     >     (5.65*phi-22.0*(phi**2)+23.4*(phi**3))*
     >     (1+tanh((mp-(0.65-0.24*phi))/0.35))


      cd = (24.0/re)*rcd1

      cd = cd/(1.0-phi) + b3 + (24.0/re)*(1.0-phi)*(b1+b2)

      beta = 3.0*rpi*rmu*dp*(re/24.0)*cd

      return
      end
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
! Modified March 5, 2024
!
! Subroutine for quasi-steady force
!
! QS Force calculated as a function of Re, Ma and phi
!
! Use Osnes etal (2023) correlations
! A.N. Osnes, M. Vartdal, M. Khalloufi, 
!    J. Capecelatro, and S. Balachandar.
! Comprehensive quasi-steady force correlations for compressible flow
!    through random particle suspensions.
! International Journal of Multiphase Flow, Vol. 165, 104485, (2023).
! doi: https://doi.org/10.1016/j.imultiphaseflow.2023.104485.
!
! E. Loth, J.T. Daspit, M. Jeong, T. Nagata, and T. Nonomura.
! Supersonic and hypersonic drag coefficients for a sphere.
! AIAA Journal, Vol. 59(8), pp. 3261-3274, (2021).
! doi: https://doi.org/10.2514/1.J060153.
!
! NOTE: Re<45 Rarified formula of Loth et al has been redefined by Balachandar
! to avoid singularity as Ma -> 0.
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_Osnes(i,beta)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 gamma,Knp,fKn,CD1,s,JM,CD2,
     >   cd_loth,CM,GM,HM,b1,b2,b3,cd,beta
      real*8 sgby2, JMt

!
! Code:
!
      gamma = 1.4d0
      mp  = dmax1(rmachp,0.01d0)
      phi = dmax1(rphip,0.0001d0)
      re  = dmax1(rep,0.1d0)

      ! Loth's correlation
      if (re .le. 45.0) then
         ! Rarefied-dominated regime
         Knp = sqrt(0.5d0 * rpi * gamma) * mp / re
         if (Knp > 0.01) then
            fKn = 1.0d0 / (1.0d0 
     >            + Knp*(2.514d0 + 0.8d0*exp(-0.55d0/Knp)))
         else
            fKn = 1.0d0 / (1.0d0 
     >            + Knp*(2.514d0 + 0.8d0*exp(-0.55d0/0.01)))

         end if
         CD1 = (24.0/re)*(1.0d0 + 0.15d0* re**(0.687d0)) * fKn
         s = mp * sqrt(0.5d0 * gamma)
         sgby2 = sqrt(0.5d0 * gamma)
         if (mp <= 1) then
            !JMt = 2.26d0*(mp**4) - 0.1d0*(mp**3) + 0.14d0*mp
            JMt = 2.26d0*(mp**4) + 0.14d0*mp
         else
            JMt = 1.6d0*(mp**4) + 0.25d0*(mp**3) 
     >            + 0.11d0*(mp**2) + 0.44d0*mp
         end if
!
! Reformulated version of Loth et al. to avoid singularity at mp = 0
!
         CD2 = (1.0d0 + 2.0d0*(s**2)) * exp(-s**2) * mp
     >          / ((sgby2**3)*sqrt(rpi)) 
     >          + (4.0d0*(s**4) + 4.0d0*(s**2) - 1.0d0) 
     >          * erf(s) / (2.0d0*(sgby2**4)) 
     >          + (2.0d0*(mp**3) / (3.0d0 * sgby2)) * sqrt(rpi)

         CD2 = CD2 / (1.0d0 + (((CD2/JMt) - 1.0d0) * sqrt(re/45.0d0)))
         cd_loth = CD1 / (1.0d0 + (mp**4)) 
     >          +  CD2 / (1.0d0 + (mp**4))
      else
         ! Compression-dominated regime
         ! TLJ: coefficients tweaked to get continuous values
         !      on the two branches at the critical points
         if (mp < 1.5d0) then
            CM = 1.65d0 + 0.65d0*tanh(4d0*mp - 3.4d0)
         else
            !CM = 2.18d0 - 0.13d0*tanh(0.9d0*mp - 2.7d0)
            CM = 2.18d0 - 0.12913149918318745d0*tanh(0.9d0*mp - 2.7d0)
         end if
         if (mp < 0.8) then
            GM = 166.0d0*(mp**3) + 3.29d0*(mp**2) - 10.9d0*mp + 20.d0
         else
            !GM = 5.0d0 + 40.d0*(mp**(-3))
            GM = 5.0d0 + 47.809331200000017d0*(mp**(-3))
         end if
         if (mp < 1) then
            HM = 0.0239d0*(mp**3) + 0.212d0*(mp**2) 
     >           - 0.074d0*mp + 1.d0
         else
            !HM =   0.93d0 + 1.0d0 / (3.5d0 + (mp**5))
            HM =   0.93967777777777772d0 + 1.0d0 / (3.5d0 + (mp**5))
         end if

         cd_loth = (24.0/re)*(1 + 0.15 * (re**(0.687)))*HM + 
     >      0.42*CM/(1+42500/re**(1.16*CM) + GM/sqrt(re))

      end if

      b1 = 5.81*phi/((1.0-phi)**2) + 
     >     0.48*(phi**(1.d0/3.d0))/((1.0-phi)**3)

      b2 = ((1.0-phi)**2)*(phi**3)*
     >     re*(0.95+0.61*(phi**3)/((1.0-phi)*2))

      b3 = dmin1(sqrt(20.0d0*mp),1.0d0)*
     >     (5.65*phi-22.0*(phi**2)+23.4*(phi**3))*
     >     (1+tanh((mp-(0.65-0.24*phi))/0.35))

      cd = cd_loth/(1.0-phi) + b3 + (24.0/re)*(1.0-phi)*(b1+b2)

      beta = 3.0*rpi*rmu*dp*(re/24.0)*cd


      return
      end
!-----------------------------------------------------------------------
!
! Created Sep. 29, 2025
!      
! Subroutine for Quasi-Steady Drag Model of Gidaspow 
!
! D. Gidaspow, Multiphase Flow and Fluidization (Academic Press, 1994)
!
! Note: Model is provided per cell volume. We convert that to per
! particle using the particle volume fraction and volume
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_Gidaspow(i,beta,cd)
!
      implicit none
!
      include "PPICLF"
!
! Internal:

! Internal variables
      integer*4 i
      real*8 cd, beta, phifRep, phif
!
! Code:

      phi  = dmax1(rphip,0.0001d0)
      phif = dmax1(rphif,0.0001d0)
      re  = dmax1(rep,0.1d0)     

      phifRep = phif*re

      if(phifRep .lt. 1000.0) then
        cd = 24.0/phifRep *(1.0+0.15*(phifRep)**0.687)
      else 
        cd = 0.44
      endif

      if(phif .lt. 0.8) then
        beta = 150.0*((phi**2)*rmu)/(phif * dp**2)
     >          + 1.75*(rhof*phi*vmag/dp)
      else 
        beta = 0.75*cd*phi*rhof*vmag/(dp*phif**1.65)
      endif

      beta = beta*(rpi*dp**3)/(6.0*phi)

      return
      end
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Quasi-steady force fluctuations
! Osnes, Vartdal, Khalloufi, Capecelatro (2023)
!   Comprehensive quasi-steady force correlations
!   for compressible flow through random particle
!   suspensions.
!   International Journal of Multiphase Flows,
!   Vo. 165, 104485.
! Lattanzi, Tavanashad, Subramaniam, Capecelatro (2022)
!   Stochastic model for the hydrodynamic force in
!   Euler-Lagrange silumations of particle-laden flows.
!   Physical Review Fluids, Vol. 7, 014301.
! Note: To compute the granular temperature, we assume
!   the velocity fluctuations are uncorrelated.
! Note: The means are computed using a box filter with an
!   adaptive filter width.
! Compute mean using box filter for langevin model - not for fedback
!
! The mean is calcuated according to Lattanzi etal,
!   Physical Review Fluids, 2022.
!
! Sam - TODO: either couple the projection filter to the
! fluctuation filter or completely decouple them.
! Right now they use the same filter width and are assumed to
! be the same. 
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_fluct_Lattanzi(i,iStage,fqs_fluct)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 fqs_fluct(3)

      real*8 aSDE,bq,bSDE,chi,denum,dW1,dW2,dW3,fq,Fs,gkern,
     >   sigF,tF_inv,theta,upflct,vpflct,wpflct,Z1,Z2,Z3
      real*8 TwoPi

!
! Code:
!
      TwoPi = 2.0d0*acos(-1.0d0)


      if (qs_fluct_filter_flag==0) then
         denum = max(dfloat(icpmean),1.d0)  ! for arithmetic mean
      else if (qs_fluct_filter_flag==1) then
         denum = max(phipmean,1.0e-6)  ! for volume mean
      endif
      

      upmean = upmean / denum
      vpmean = vpmean / denum
      wpmean = wpmean / denum
      u2pmean = u2pmean / denum
      v2pmean = v2pmean / denum
      w2pmean = w2pmean / denum

      if (ppiclf_debug==2) then
         if (ppiclf_nid==0 .and. iStage==1) then
            write(6,"(2x,E16.8,i5,16(1x,F13.8))") ppiclf_time,
     >                denum,upmean,vpmean,wpmean
         endif
      endif

      ! Lattenzi is valid only for incompressible flows,
      !   so here we use the compressible correction of Osnes
      !   Equations (9-12) of Osnes paper, with eqn (9) corrected
      !   Note that sigD has units of force; N = Pa-m^2
      fq = 6.52*rphip - 22.56*(rphip**2) + 49.90*(rphip**3)
      Fs = 3.0*rpi*rmu*dp*(1.0+0.15*((rep*(1.0-rphip))**0.687))
     >          *(1.0-rphip)*vmag
      bq = min(sqrt(20.0*rmachp),1.0)*0.55*(rphip**0.7) 
     >          *(1.0+tanh((rmachp-0.5)/0.2))
      sigF = (fq + bq)*Fs

      chi = (1.0+2.50*rphip+4.51*(rphip**2)+4.52*(rphip**3))
     >         /((1.0-(rphip/0.64)**3)**0.68)

      fqs_fluct = 0.0d0

! Particle velocity fluctuation and Granular Temperature
! Need particle velocity mean
! Though the theory assumes granular temperature to be an 
!    average over neighboring particles, here it is approximated 
!    as that of the chosen particle - Comment 3/6/24
!    This is now fixed - Comment 4/12/24
!
      ! Particle velocity fluctuation
      upflct = ppiclf_y(PPICLF_JVX,i) - upmean
      vpflct = ppiclf_y(PPICLF_JVY,i) - vpmean
      wpflct = ppiclf_y(PPICLF_JVZ,i) - wpmean

      ! Granular temperature
      ! theta = (upflct*upflct + vpflct*vpflct + wpflct*wpflct)/3.0
      ! This is averaged over neighboring particles
      theta  = ((u2pmean + v2pmean + w2pmean) - 
     >          (upmean**2 + vpmean**2 + wpmean**2))/3.0d0

      ! 11/21/24 - Thierry - prevent NaN variables
      if(theta.le.1.d-12) then
        theta = 0.0d0
      endif

      tF_inv = (24.0*rphip*chi)/dp * sqrt(theta/rpi)

      aSDE = tF_inv
      bSDE = sigF*sqrt(2.0*tF_inv)

      call RANDOM_NUMBER(UnifRnd)

      Z1 = sqrt(-2.0d0*log(UnifRnd(1)))*cos(TwoPi*UnifRnd(2))
      Z2 = sqrt(-2.0d0*log(UnifRnd(3)))*cos(TwoPi*UnifRnd(4))
      Z3 = sqrt(-2.0d0*log(UnifRnd(5)))*cos(TwoPi*UnifRnd(6))

      dW1 = sqrt(fac)*Z1
      dW2 = sqrt(fac)*Z2
      dW3 = sqrt(fac)*Z3

      fqs_fluct(1) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFX,i)
     >           + bSDE*dW1
      fqs_fluct(2) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFY,i)
     >           + bSDE*dW2
      fqs_fluct(3) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFZ,i)
     >           + bSDE*dW3


      if (ppiclf_debug==2 .and. (iStage==1 .and. ppiclf_nid==0)) then
         if (ppiclf_time.gt.2.d-8) then
         if (i<=10) then
            write(7350+(i-1)*1,*) i,ppiclf_time,             ! 0-1
     >         rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,  ! 2-10
     >         rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac, ! 11-18
     >         vx,vy,vz,ppiclf_dt,                           ! 19-22
     >         ppiclf_npart,ppiclf_n_bins(1:3),              ! 23-26
     >         ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3), ! 27
     >         ppiclf_binb(1:6),                             ! 28-33
     >         upmean,vpmean,wpmean,phipmean,                ! 34-37
     >         ppiclf_y(PPICLF_JVX:PPICLF_JVZ,i),            ! 38-40
     >         upflct,vpflct,wpflct,icpmean,                 ! 41-44
     >         fq,Fs,bq,theta,chi,tF_inv,                    ! 45-50
     >         aSDE,bSDE,sigF,                               ! 51-53
     >         fqs_fluct(1:3),                               ! 54-56
     >         Z1,Z2,Z3,dW1,dW2,dW3,                         ! 57-62
     >         ppiclf_np
         endif
         endif
      endif


      return
      end
!
!
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024 - T.L. Jackson
! Modified 3/6/24 - Balachandar
!
! Quasi-steady force fluctuations
! Osnes, Vartdal, Khalloufi, Capecelatro (2023)
!   Comprehensive quasi-steady force correlations
!   for compressible flow through random particle
!   suspensions.
!   International Journal of Multiphase Flows,
!   Vo. 165, 104485.
! Lattanzi, Tavanashad, Subramaniam, Capecelatro (2022)
!   Stochastic model for the hydrodynamic force in
!   Euler-Lagrange silumations of particle-laden flows.
!   Physical Review Fluids, Vol. 7, 014301.
! Note: To compute the granular temperature, we assume
!   the velocity fluctuations are uncorrelated.
! Note: The means are computed using a box filter with an
!   adaptive filter width.
! Compute mean using box filter for langevin model - not for fedback
!
! The mean is calcuated according to Lattanzi etal,
!   Physical Review Fluids, 2022.
!
! Sam - TODO: either couple the projection filter to the
! fluctuation filter or completely decouple them.
! Right now they use the same filter width and are assumed to
! be the same. 
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_fluct_Osnes(i,iStage,fqs_fluct,
     >                                      xi_par,xi_perp,xi_T,
     >                                      fqsx,fqsy,fqsz)
!                                                    
      implicit none
!
      include "PPICLF"
!
! Input:
      integer*4 i, iStage
      real*8 fqsx, fqsy, fqsz
!
! Output:
      real*8 xi_par, xi_perp, xi_T
      real*8 fqs_fluct(3)
!
! Internal:
!
      real*8 aSDE,bq,chi,denum,dW1,dW2,dW3,fq,Fs,gkern,
     >   sigD,tF_inv,theta,upflct,vpflct,wpflct,Z1,Z2,Z3
      real*8 TwoPi
      real*8 bSDE_CD, bSDE_CL, bSDE_CT, CD_frac, CD_prime
      real*8 sigT,sigCT
      real*8 sigmoid_cf, f_CF
      real*8 avec(3)
      real*8 bvec(3)
      real*8 cvec(3)
      real*8 dvec(3)
      real*8 cosrand,sinrand
      real*8 eunit(3)

      integer*4 m
      real*8 s_par, s_perp, s_T, Rmean_par, Rmean_perp, R_par, R_perp
      real*8 R(3,3), Q(3,3), Qt(3,3), Tmean_par(3)
      real*8 CD_average
      real*8 k_tilde, k_Mach, b_tilde, b_Mach, b_par, b_perp,
     >       k_Osnes, b_Osnes
      real*8 C1, C2, C3, C4, C5,
     >       D1, D2, D3, D4, D5, D6, D7, D8,
     >       E1, E2, E3, E4
      real*8 F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11,
     >       G1, G2, G3, G4, G5, G6, G7, G8,
     >       H1, H2, H3, H4, H5, H6, H7, H8,
     >       A1, A2, A3, A4
      real*8 D9, D10, D11
      real*8 fit_func
!
! Code:
!
      TwoPi = 2.0d0*acos(-1.0d0)

      if (qs_fluct_filter_flag==0) then
         denum = max(dfloat(icpmean),1.d0)  ! for arithmetic mean
      else if (qs_fluct_filter_flag==1) then
         denum = max(phipmean,1.0e-6)  ! for volume mean
      endif
      
      upmean = upmean / denum
      vpmean = vpmean / denum
      wpmean = wpmean / denum
      u2pmean = u2pmean / denum
      v2pmean = v2pmean / denum
      w2pmean = w2pmean / denum

      if (ppiclf_debug==2) then
         if (ppiclf_nid==0 .and. iStage==1) then
            write(6,*) 'FLUC1 ',
     >          ppiclf_time,denum,upmean,vpmean,wpmean,
     >          abs(upmean),abs(vpmean),abs(wpmean),
     >          u2pmean,v2pmean,w2pmean
         endif
      endif

      ! Equations (9-12) of Osnes paper, with eqn (9) corrected
      ! Note that sigD has units of force; N = Pa-m^2
      fq = 6.52*rphip - 22.56*(rphip**2) + 49.90*(rphip**3)
      Fs = 3.0*rpi*rmu*dp*(1.0+0.15*((rep*(1.0-rphip))**0.687))
     >          *(1.0-rphip)*vmag
      bq = min(sqrt(20.0*rmachp),1.0)*0.55*(rphip**0.7) 
     >          *(1.0+tanh((rmachp-0.5)/0.2))
      sigD = (fq + bq)*Fs

! Particle velocity fluctuation and Granular Temperature
! Need particle velocity mean
! Though the theory assumes granular temperature to be an 
!    average over neighboring particles, here it is approximated 
!    as that of the chosen particle - Comment 3/6/24
!    This is now fixed - Comment 4/12/24
!
      upflct = ppiclf_y(PPICLF_JVX,i) - upmean
      vpflct = ppiclf_y(PPICLF_JVY,i) - vpmean
      wpflct = ppiclf_y(PPICLF_JVZ,i) - wpmean

      ! Granular temperature
      ! theta = (upflct*upflct + vpflct*vpflct + wpflct*wpflct)/3.0
      ! This is averaged over neighboring particles
      theta  = ((u2pmean + v2pmean + w2pmean) - 
     >          (upmean**2 + vpmean**2 + wpmean**2))/3.0d0

      ! 11/21/24 - Thierry - prevent NaN variables
      if(theta.le.1.d-12) then
        theta = 0.0d0
      endif

      chi = (1.0 + 2.50*rphip + 4.51*(rphip**2) + 4.52*(rphip**3))
     >         /((1.0-(rphip/0.64)**3)**0.68)

      tF_inv = (24.0*rphip*chi/dp) * sqrt(theta/rpi)

      aSDE = tF_inv
      bSDE_CD = sigD*sqrt(2.0*tF_inv)  ! Modified 3/6/24

! Fluctuating perpendicular force
! Compare CD' (units of force) against sigma_CD (units of force) 
!    to determine which of 5 bins to use for sigCT
!
! Added 3/6/24 
! Modified 3/14/24 
!
      ! 03/13/2025 - Thierry - if velocity is very small, don't impose fluctuations
      if(vmag > 1.d-8) then
        avec = [vx,vy,vz]/vmag

        CD_prime = ppiclf_rprop(PPICLF_R_FLUCTFX,i)*avec(1) +
     >             ppiclf_rprop(PPICLF_R_FLUCTFY,i)*avec(2) +
     >             ppiclf_rprop(PPICLF_R_FLUCTFZ,i)*avec(3)
        CD_frac  = CD_prime/sigD

      else
        avec     = [1.0, 0.0, 0.0]
        CD_prime = 0.0
        sigD     = 0.0
        CD_frac  = 0.0
      endif

      ! Thierry Daoud - Updated June 2, 2024
      sigmoid_cf = 1.0 / (1.0 + exp(-CD_frac))
      f_CF = 0.39356905*sigmoid_cf + 0.43758848
      sigT  = f_CF*sigD
      bSDE_CL = sigT*sqrt(2.0*tF_inv)

      if (ppiclf_debug==2) then
      if (i<=4) then
         if (ppiclf_nid==0 .and. iStage==1) then
            write(6,*) 'FLUC1 ',i,
     >        CD_prime,CD_frac,sigD,theta,bSDE_CD,bSDE_CL
         endif
      endif
      endif


! Calculate the three orthogonal unit vectors
! The first vector (avec) is vx/vmag, vy/vmag, and vz/vmag
! The second (bvec) is constructued by taking cross-product with eunit
!   Note: if avec is in dir. of e_x=(1,0,0), use e_y=(0,1,0) to get e_z
!       : if avec is in dir. of e_y=(0,1,0), use e_z=(0,0,1) to get e_x
!       : if avec is in dir. of e_z=(0,0,1), use e_x=(1,0,0) to get e_y
! The third (cvec) is cross product of the first two
! written 3/6/24
!
! avec : unit vector in main direction
! bvec, cvec: two orthogonal vectors to avec
      eunit = [1,0,0]
      if (abs(avec(2))+abs(avec(3)) <= 1.d-8) then
         eunit = [0,1,0]
      elseif (abs(avec(1))+abs(avec(3)) <= 1.d-8) then
         eunit = [0,0,1]
      endif

      bvec(1) = avec(2)*eunit(3) - avec(3)*eunit(2)
      bvec(2) = avec(3)*eunit(1) - avec(1)*eunit(3)
      bvec(3) = avec(1)*eunit(2) - avec(2)*eunit(1)
      denum   = max(1.d-8,sqrt(bvec(1)**2 + bvec(2)**2 + bvec(3)**2))
      bvec    = bvec / denum

      cvec(1) = avec(2)*bvec(3) - avec(3)*bvec(2)
      cvec(2) = avec(3)*bvec(1) - avec(1)*bvec(3)
      cvec(3) = avec(1)*bvec(2) - avec(2)*bvec(1)
      denum   = max(1.d-8,sqrt(cvec(1)**2 + cvec(2)**2 + cvec(3)**2))
      cvec    = cvec / denum

      ! Generate  Gaussian Random Values
      call RANDOM_NUMBER(UnifRnd)

! Box-Muller transform for generating two independent standard normal
! (Gaussian) random variables
! Z1 & Z2 are standard normal random variables
      Z1 = sqrt(-2.0d0*log(UnifRnd(1)))*cos(TwoPi*UnifRnd(2))
      Z2 = sqrt(-2.0d0*log(UnifRnd(3)))*cos(TwoPi*UnifRnd(4))

! dW1 & dW2 are scaled stochastic amplitudes       
      dW1 = sqrt(fac)*Z1
      dW2 = sqrt(fac)*Z2

! Random mixture of bvec and cvec - make sure the new one is a unit vector
! Added 3/6/24
! dvec : Random unit direction in perpendicular plane 
      cosrand = cos(TwoPi*UnifRnd(5))
      sinrand = sin(TwoPi*UnifRnd(5)) 
      dvec(1) = bvec(1)*cosrand + cvec(1)*sinrand
      dvec(2) = bvec(2)*cosrand + cvec(2)*sinrand
      dvec(3) = bvec(3)*cosrand + cvec(3)*sinrand
      denum   = max(1.d-8,sqrt(dvec(1)**2 + dvec(2)**2 + dvec(3)**2))
      dvec    = dvec/denum

      ! Store and project back fqs_fluct values only if QSFLUCT = 2
      ! Otherwise, need to call subroutine for PseudoTurbulence
      if(qs_fluct_flag .eq. 0) then
        fqs_fluct = 0.0d0
      else
      fqs_fluct(1) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFX,i)
     >           + bSDE_CD*dW1*avec(1) + bSDE_CL*dW2*dvec(1)
      fqs_fluct(2) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFY,i)
     >           + bSDE_CD*dW1*avec(2) + bSDE_CL*dW2*dvec(2)
      fqs_fluct(3) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFZ,i)
     >           + bSDE_CD*dW1*avec(3) + bSDE_CL*dW2*dvec(3)
      endif

      if (ppiclf_debug==2 .and. (iStage==1 .and. ppiclf_nid==0)) then
         if (ppiclf_time.gt.2.d-8) then
         if (i<=10) then
            write(7350+(i-1)*1,*) i,ppiclf_time,             ! 0-1
     >         rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,  ! 2-10
     >         rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac, ! 11-18
     >         vx,vy,vz,ppiclf_dt,                           ! 19-22
     >         ppiclf_npart,ppiclf_n_bins(1:3),              ! 23-26
     >         ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3), ! 27
     >         ppiclf_binb(1:6),                             ! 28-33
     >         upmean,vpmean,wpmean,phipmean,                ! 34-37
     >         ppiclf_y(PPICLF_JVX:PPICLF_JVZ,i),            ! 38-40
     >         upflct,vpflct,wpflct,icpmean,                 ! 41-44
     >         fq,Fs,bq,theta,chi,tF_inv,                    ! 45-50
     >         aSDE,bSDE_CD,bSDE_CL,                         ! 51-53
     >         sigD,sigT,                                    ! 54-55
     >         CD_prime,CD_frac,sigmoid_cf,f_CF,             ! 56-59
     >         eunit,avec,bvec,                              ! 60-68
     >         cvec,dvec,rpi,                                ! 69-75
     >         fqs_fluct(1:3),                               ! 76-78
     >         Z1,Z2,dW1,dW2,                                ! 79-82
     >         ppiclf_np,                                    ! 83
     >         ppiclf_y(PPICLF_JX:PPICLF_JZ,i)               ! 84-86
         endif
         endif
      endif

!---------------------------------------------------------------------------
! Pseudo-Turbulence Calculations starts here 
      if(pseudoTurb_flag==1) then
        !phi = max(0.05d0, min(0.3d0, rphip))
        !mp = max(0.0d0, min(0.87d0, rmachp))
        !re = max(0.0d0, min(266.0d0, rep))
        phi = rphip
        mp = rmachp
        re = rep
        rem = (1.0-phi)*re

      ! Constants taken from Osnes PseudoTurbulent paper, Table 1
        C1 = -1.2152; D1 = -0.0462; E1 = -0.2906;
        C2 = -7.6314; D2 = -0.1068; E2 =  1.1899; 
        C3 =  0.2889; D3 =  0.6793; E3 =  0.5218; 
        C4 =  0.6143; D4 =  1.1461; E4 =  0.0699;
        C5 =  0.3082; D5 = -2.6886; 
                      D6 = -2.1376;
                      D7 =  0.4873;
                      D8 =  0.2395;
                    ! Dr. Bala's terms
                      D9  =  0.716
                      D10 = -2.14
                      D11 =  1.6

      ! Constants taken from Osnes PseudoTurbulent paper, Table 2
        F1  = -0.0022; G1 = -0.2867; H1 =  0.4992
        F2  = -0.0219; G2 =  0.2176; H2 = -1.3528
        F3  =  0.0932; G3 =  0.2826; H3 = -0.1358
        F4  = -0.0135; G4 = -0.0644; H4 = -0.1463
        F5  =  0.0361; G5 =  0.0466; H5 =  0.2583
        F6  =  0.0403; G6 =  0.0973; H6 = -0.3339
        F7  = -0.0761; G7 = -0.0081; H7 = -0.0407
        F8  =  0.0599; G8 = -0.0235; H8 = -0.0806
        F9  =  0.0164;
        F10 =  0.0453;
        F11 = -0.0265;
     
        ! zero out variables  at first
        Rmean_par = 0.0d0 ; Rmean_perp = 0.0d0
        R = 0.0d0; Rsg = 0.0d0
        Tmean_par = 0.0d0; T_par = 0.0d0
        
        ! CD_average is zero at early time steps

        avec = [vx,vy,vz]/vmag

        CD_average = fqsx*avec(1) +
     >               fqsy*avec(2) +
     >               fqsz*avec(3)

        ! avoiding singularity
        if(CD_average .lt. 1.d-8) return

        ! Reynolds Subgrid Stress Tensor - Eulerian Mean Model 
                                                                       
        ! Reynolds number and vol fraction dependent k^tilde and b_par 
        ! Mehrabadi's terms
        k_tilde = 2.0*phi + 2.5*phi*((1.0-phi)**3) * 
     >         exp(-phi*(rem**0.5))                                 
                                                                    
!        b_par = 0.523/(1.0+0.305*exp(-0.114*rem)) *
!     >          exp(-3.511*phi/(1.0+1.801*exp(-0.005*rem)))

        ! 08/25/2025 - Thierry - Fitted phip function to better match
        ! formulation with Osnes's low Mach number data
        fit_func = -10.18530152*phi**3 + 10.94163073*phi**2
     >              -7.07374862*phi +  0.38424203
        b_par = 0.523/(1.0+0.305*exp(-0.114*rem)) *
     >          exp(fit_func/(1.0+1.801*exp(-0.005*rem)))
                                                                       
        ! Mach number correction provided by Osnes                            
!        k_Mach = phi*(C1 + C2*phi + re**C3) * 
!     >        (tanh(C4/C5) + tanh((mp - C4)/C5))

        ! Mach number correction at Re=100, coeff taken from Osnes
        ! cap vol fraction here at 0.3
        k_Mach = min(phi,0.3)*(-6.918*min(phi,0.3) + 2.238) *
     >        (tanh(C4/C5) + tanh((mp - C4)/C5))

!        b_Mach = (D1 + (re/300.0)*(D2 + D3*re/300.0) +
!     >         phi*(D4 + D5*(re**2/300.0**2) + D6*phi)) *
!     >         (tanh(-D7/D8) - tanh((mp-D7)/D8))

        ! First term was corrected by Dr. Bala (D9, D10, D11)
        b_Mach = phi*(D9 + D10*phi + D11*phi**2) *
     >         (tanh(-D7/D8) - tanh((mp-D7)/D8))

        ! Corrected k^tilde and b_par components                       
        k_Osnes =  k_tilde*(1.0d0 + k_Mach)
        b_Osnes =  b_par  *(1.0d0 + b_Mach)
        b_perp  = -b_Osnes/2.0d0
                                                                       
        ! Mean Eulerian Reynolds Subgrid Stress - Parallel Component   
        Rmean_par  = 2.0d0*k_Osnes*(b_Osnes  + 1.0d0/3.0d0)
                                                                       
        ! Mean Eulerian Reynolds Subgrid Stress - Perpendicular Component
        Rmean_perp = 2.0d0*k_Osnes*(b_perp + 1.0d0/3.0d0)

        ! Mean Pseudo Turbulent Kinetic Energy Model
        Tmean_par = E1 + E2*phi/(E3 + re/300.0) + E4*mp
        
c--  Multiply by the mean relative flow kinetic energy to dimentionalize      
        Rmean_par  = Rmean_par  * 0.5d0 * vmag**2
        Rmean_perp = Rmean_perp * 0.5d0 * vmag**2

c--  Multiply by the mean relative velocity & flow kinetic energy to dimentionalize      
        Tmean_par(1)  = Tmean_par(1) * vx * k_Osnes * 0.5d0 * vmag**2
        Tmean_par(2)  = Tmean_par(2) * vy * k_Osnes * 0.5d0 * vmag**2
        Tmean_par(3)  = Tmean_par(3) * vz * k_Osnes * 0.5d0 * vmag**2

c------ Lagrangian Model
  
        ! 08/15/2025 - Ditch A1 per Dr. Bala, and set A2 to this constant per Osnes
        A1 = 0.0d0
        A2 = 0.064

        A3 = G1 + G2/(min(phi,0.3) + G3) + G4 * mp
        
        ! 09/02/2025 - Cap according to Osnes model range
        A4 = H1 + H2*max(0.0d0, min(0.3d0, phi))
     >          + H3*max(0.0d0, min(0.87d0, mp)) 
     >          + H4*max(30.0d0, min(266.0d0, re))/300.0d0
  
        s_par = F8 + F9/(phi + F10) + F11 * mp
        !s_perp = G5/(phi + G6) + (G7*re)/(300.0*phi) + G8
c---     We ditch Osnes's expression for s_perp and assume it as big as s_par
        s_perp = s_par
        
        ! 09/02/2025 - Cap according to Osnes model range
        s_T = H5 + H6*max(0.0d0, min(0.3d0, phi))
     >           + H7*max(0.0d0, min(0.87d0, mp)) 
     >           + H8*max(30.0d0, min(266.0d0, re))/300.0d0

        tF_inv = (24.0*phi*chi/dp) * sqrt(theta/rpi)
        aSDE = tF_inv
        bSDE_CD = s_par *sqrt(2.0*tF_inv)
        bSDE_CL = s_perp*sqrt(2.0*tF_inv)
        bSDE_CT = s_T *sqrt(2.0*tF_inv)
  
        call RANDOM_NUMBER(UnifRnd)
        
        ! Box-Muller transform for generating two independent standard normal
        ! (Gaussian) random variables
        ! Z1 & Z2 are standard normal random variables
        Z1 = sqrt(-2.0d0*log(UnifRnd(1))) * cos(TwoPi*UnifRnd(2))
        Z2 = sqrt(-2.0d0*log(UnifRnd(3))) * sin(TwoPi*UnifRnd(4))
        Z3 = sqrt(-2.0d0*log(UnifRnd(5))) * cos(TwoPi*UnifRnd(6))

        ! dW1 & dW2 are scaled stochastic amplitudes       
        dW1 = sqrt(fac)*Z1
        dW2 = sqrt(fac)*Z2
        dW3 = sqrt(fac)*Z3
  
        ! Langevin Model implemented for xi_par, xi_perp, xi_T
        xi_par = (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_XIPAR,i)
     >            + bSDE_CD*dW1
        xi_perp = (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_XIPERP,i)
     >            + bSDE_CL*dW2
        xi_T = (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_XIT,i)
     >            + bSDE_CT*dW3

        ! CD_prime has unit of Force
        ! CD_average has unit of Force

        ! Lagrangian Reynolds Subgrid Stress - Parallel Component
        R_par = 1.0 + A1 + A2 * CD_prime / CD_average + xi_par
  
        ! Lagrangian Reynolds Subgrid Stress - Perpendicular Component
        R_perp = 1.0 + A3 * CD_prime / CD_average + xi_perp

c--  Multiply Lagrangian Model by the Eulerian Mean Model
        R_par  = R_par  * Rmean_par 
        R_perp = R_perp * Rmean_perp

c---  Q = [avec | bvec | cvec], 3x3 matrix
c---  avec : unit vector in main direction
c---  bvec, cvec: two orthogonal vectors to avec

        ! 08/15/2025 - Thierry - still need to finalize how to do the
        ! rotation of the R tensor 
        do m=1,3
          Q(m,1) = avec(m)
          Q(m,2) = bvec(m)
          Q(m,3) = cvec(m)
        enddo
  
        Qt = transpose(Q)
  
c--- R = |R_par,   0   ,   0   |
c---     | 0   , R_perp,   0   |
c---     | 0       0   , R_perp|
  
c---  R matrix only has diagonal components
        R(1,1) = R_par
        R(2,2) = R_perp
        R(3,3) = R_perp
  
c--- Now Rotate the matrix, Rsg = Q . R . Q^T
  
       Rsg = matmul(Q, matmul(R,Qt))

c--- Osnes Formulation for PTKE

      T_par = A4 * CD_prime/CD_average + xi_T

c--  Multiply by the mean relative velocity & flow kinetic energy to dimentionalize      
c--  then add mean PTKE
       T_par(1) = T_par(1) * vx * k_Osnes * 0.5d0 * vmag**2 
     >            + Tmean_par(1)

       T_par(2) = T_par(2) * vy * k_Osnes * 0.5d0 * vmag**2 
     >            + Tmean_par(2)

       T_par(3) = T_par(3) * vz * k_Osnes * 0.5d0 * vmag**2 
     >            + Tmean_par(3)

      endif ! pseudoTurb_flag

      return
      end
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for added mass
!   also called the quasi-unsteady force,
!   or the inviscid unsteady force in case of the Euler equations
!      
! ADDED from rocflu
! Added mass force (phi corrections):
!   On the dispersed two-phase flow in the laminar flow regime
!   - Zuber (1964), Chem. Engng Sci.
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_AM_Parmar(i,iStage,
     >                 famx,famy,famz,rmass_add)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz, rmass_add
      real*8 rcd_am
      real*8 SDrho
      real*8 ug,vg,wg,vgradrho
      real*8 famx_Ling
      real*8 famx_Brad

!
! Code:
!
      ! Mach number correction
      if (rmachp .gt. 0.6) then        
         rcd_am = 1.0 + 1.8*((0.6)**2)  + 7.6*((0.6)**4)
      else
         rcd_am = 1.0 + 1.8*(rmachp**2) + 7.6*(rmachp**4)
      endif
      rcd_am = rcd_am * 0.5

      ! Sangani's volume fraction correction for dilute random arrays
      ! Capping volume fraction at 0.5 
      ! 09/20/2025 - Thierry - Sangani's volume fraction correction
      ! overshoots Unary added-mass term A LOT. 
      !rcd_am = rcd_am*(1.0+3.32*min(rphip,0.5))

      ! Adopting the volume fraction correction from Beguin & Etienne
      ! (2016).
      rcd_am = rcd_am*(1.0+0.68*rphip**2)
      rmass_add = rhof*ppiclf_rprop(PPICLF_R_JVOLP,i)*rcd_am

      !NEW Added mass, using how rocflu does it
      !1st Derivative, substantial how rocflu does it
      SDrho = ppiclf_rprop(PPICLF_R_JRHSR,i) 
     >      + ppiclf_y(PPICLF_JVX,i) * ppiclf_rprop(PPICLF_R_JPGCX,i)
     >      + ppiclf_y(PPICLF_JVY,i) * ppiclf_rprop(PPICLF_R_JPGCY,i)
     >      + ppiclf_y(PPICLF_JVZ,i) * ppiclf_rprop(PPICLF_R_JPGCZ,i)

      ! 03/11/2025 - Thierry - substantial derivative from Rocflu is 
      !              weighted by \phi^g.
      ! d(rho^g phi^g)/dt = rho^g * d(phi^g)/dt + phi^g * d(rho^g)/dt
      !                   = phi^g * d(rho^g)/dt
      !  
      !     d(rho^g)/dt   = SDrho = d(rho phi^g)/dt / phi^g
      SDrho = SDrho / (rphif) 

      ! 03/23/2025 - TLJ - added extra term involving grad(rhog)
      vgradrho = vx*ppiclf_rprop(PPICLF_R_JRHOGX,i) +
     >           vy*ppiclf_rprop(PPICLF_R_JRHOGY,i) +
     >           vz*ppiclf_rprop(PPICLF_R_JRHOGZ,i)

      ug = ppiclf_rprop(PPICLF_R_JUX,i)
      vg = ppiclf_rprop(PPICLF_R_JUY,i)
      wg = ppiclf_rprop(PPICLF_R_JUZ,i)

      famx = rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i) *
     >   (vx*SDrho + rhof*ppiclf_rprop(PPICLF_R_JSDRX,i) + ug*vgradrho)

      famy = rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i) *
     >   (vy*SDrho + rhof*ppiclf_rprop(PPICLF_R_JSDRY,i) + vg*vgradrho)

      famz = rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i) *
     >   (vz*SDrho + rhof*ppiclf_rprop(PPICLF_R_JSDRZ,i) + wg*vgradrho)


      if (1==2) then
      ! This is Ling's 2012 formulation where he replaced
      !   D(rhog*ug) with -grad(pg)
      famx_Ling = rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i)*
     > (-ppiclf_rprop(PPICLF_R_JDPDX,i) - ppiclf_y(PPICLF_JVX,i)*SDrho)

      ! Original version
      ! /home/tlj/Codes_Rocflu/Rocflu_picl_tlj/ppiclf/source/ppiclf_user.f
      ! In the original version JSDRX was assumed to be D(rhog*ug)/Dt,
      !    but this was before we realized that the conserved Rocflu
      !    variable was phig*rhog*ug
      famx_Brad = rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i) *
     >   (ppiclf_rprop(PPICLF_R_JSDRX,i) - ppiclf_y(PPICLF_JVX,i)*SDrho)

      ! writing data only for the median particle
      if((ppiclf_iprop(5,i).eq.29.0) .and. (ppiclf_iprop(6,i).eq.0.0)
     >   .and. (ppiclf_iprop(7,i).eq.151.0)) then
      
      open(unit=20,file='fort.20',position='append') 
      write(20,*) ppiclf_time,famx-rmass_add*ppiclf_ydot(PPICLF_JVX,i),
     >            rcd_am, ppiclf_rprop(PPICLF_R_JVOLP,i),
     >            rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i),
     >            vx, SDrho, vx*SDrho,
     >            rhof, ppiclf_rprop(PPICLF_R_JSDRX,i),
     >            rhof*ppiclf_rprop(PPICLF_R_JSDRX,i),
     >            ug, vgradrho,
     >            ug*vgradrho,
     >            famx_Ling-rmass_add*ppiclf_ydot(PPICLF_JVX,i),
     >            famx_Brad-rmass_add*ppiclf_ydot(PPICLF_JVX,i)
      flush(20)
      endif
      endif


      return
      end
!
!
!-----------------------------------------------------------------------
!
! Created May 20, 2024
!
! Subroutine for added mass
!   also called the quasi-unsteady force,
!   or the inviscid unsteady force in case of the Euler equations
!      
! Implementing Added Mass Algorithm from S.Briney (2025)
!  
! n       = number of points
! alpha   = volume fraction
! rad     = particle radius
! d       = center-to-center distance
! rmax    = center-to-center max neighbor distance
! R       = resistance matrix (output)
! dr_max  = max interaction distance between particles considered 
! poins   = 3xn array of points x, y, z. Initialized as points(3,n)
!
! correction_analytical_always = 
!    if true, always use the analytical distant neighbor correction
!    if false, use numerical when dr_max/rad < 3.49
! 
! This subroutine corresponds to the flag:  am_flag = 2 
!
!
! IMPORTANT NOTE: THIS SUBROUTINE IS CURRENTLY ONLY 
!    VALID FOR MONODISPERSE PARTICLE BEDS 
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!                          Unary part 
!-----------------------------------------------------------------------
!
!
!
      subroutine ppiclf_user_AM_Briney_Unary(i,iStage,
     >                 famx,famy,famz,rmass_add)
!
      implicit none
!
      include "PPICLF"
!
      integer i, j, k, l, n, jj
      integer*4 iStage
      real*8 rad
      real*8 famx, famy, famz, rmass_add
      real*8 rcd_am
      real*8 SDrho
      real*8 ug,vg,wg,vgradrho

!
! Code:
!
      ! Thierry - need to decide if we want to keep 
      !    1) Mach no. correction  - yes, according to Bala
      !    2) Vol frac. correction - no, this is in binary term
      !

      ! Mach number correction
      if (rmachp .gt. 0.6) then
         rcd_am = 1.0 + 1.8*((0.6)**2)  + 7.6*((0.6)**4)
      else
         rcd_am = 1.0 + 1.8*(rmachp**2) + 7.6*(rmachp**4)
      endif
      rcd_am = rcd_am * 0.5

      ! Volume fraction correction - In the Briney model the
      !    volume fraction correction is in the binary term,
      !    not the unary term
      !rcd_am = rcd_am*(1.0+2.0*rphip)

      rmass_add = rhof*ppiclf_rprop(PPICLF_R_JVOLP,i)*rcd_am

      !NEW Added mass, using how rocflu does it
      !1st Derivative, substantial how rocflu does it
      SDrho = ppiclf_rprop(PPICLF_R_JRHSR,i)
     >      + ppiclf_y(PPICLF_JVX,i) * ppiclf_rprop(PPICLF_R_JPGCX,i)
     >      + ppiclf_y(PPICLF_JVY,i) * ppiclf_rprop(PPICLF_R_JPGCY,i)
     >      + ppiclf_y(PPICLF_JVZ,i) * ppiclf_rprop(PPICLF_R_JPGCZ,i)
      ! material derivative is phi weighted in Rocflu
      ! drho/dt
      SDrho = SDrho / (rphif) 

      ! 03/23/2025 - TLJ - added extra term involving grad(rhog)
      vgradrho = vx*ppiclf_rprop(PPICLF_R_JRHOGX,i) +
     >           vy*ppiclf_rprop(PPICLF_R_JRHOGY,i) +
     >           vz*ppiclf_rprop(PPICLF_R_JRHOGZ,i)

      ug = ppiclf_rprop(PPICLF_R_JUX,i)
      vg = ppiclf_rprop(PPICLF_R_JUY,i)
      wg = ppiclf_rprop(PPICLF_R_JUZ,i)

      ! Take care of volume in Binary subroutine
      famx = rcd_am*
     >   (vx*SDrho + rhof*ppiclf_rprop(PPICLF_R_JSDRX,i) + ug*vgradrho)

      famy = rcd_am*
     >   (vy*SDrho + rhof*ppiclf_rprop(PPICLF_R_JSDRY,i) + vg*vgradrho)

      famz = rcd_am*
     >   (vz*SDrho + rhof*ppiclf_rprop(PPICLF_R_JSDRZ,i) + wg*vgradrho)

      ! Multiply by neighbors here for storing
      FamUnary(1) = famx*ppiclf_rprop(PPICLF_R_JVOLP,i)
      FamUnary(2) = famy*ppiclf_rprop(PPICLF_R_JVOLP,i)
      FamUnary(3) = famz*ppiclf_rprop(PPICLF_R_JVOLP,i)

      ! Do not multiply by volume for Fam, as this is done
      ! in user file (if nneighbors=0) or Binary subroutine (if nneighbors>0)
      Fam(1) = famx
      Fam(2) = famy
      Fam(3) = famz

      if (ppiclf_debug==2) then
      if (ppiclf_nid .eq. 0 .or. ppiclf_np == 1) then
      if (i<=5 .and. iStage==1) then  
         open(unit=7051,file='fort.7051',position='append')
         open(unit=7052,file='fort.7052',position='append')
         open(unit=7053,file='fort.7053',position='append')
         open(unit=7054,file='fort.7054',position='append')
         open(unit=7055,file='fort.7055',position='append')
         write(7050+i,*) i, ppiclf_nid, ppiclf_np, ppiclf_time, ! 0-3
     >      ppiclf_rprop(PPICLF_R_JRHSR,i),        ! 4
     >      ppiclf_rprop(PPICLF_R_JPGCX,i),        ! 5
     >      ppiclf_rprop(PPICLF_R_JPGCY,i),        ! 6
     >      ppiclf_rprop(PPICLF_R_JPGCZ,i),        ! 7
     >      SDrho,                                 ! 8
     >      rhof, rphip, rmachp, SDrho,            ! 9-12
     >      ppiclf_rprop(PPICLF_R_WDOTX,i),        ! 13
     >      ppiclf_rprop(PPICLF_R_WDOTY,i),        ! 14
     >      ppiclf_rprop(PPICLF_R_WDOTZ,i),        ! 15
     >      Wdot_neighbor_mean(1:3), nneighbors,   ! 16-19
     >      famx, famy, famz,                      ! 20-22
     >      ppiclf_rprop(PPICLF_R_JVOLP,i)*Fam(1), ! 23
     >      ppiclf_rprop(PPICLF_R_JVOLP,i)*Fam(2), ! 24
     >      ppiclf_rprop(PPICLF_R_JVOLP,i)*Fam(3)  ! 25
         flush(7051)
         flush(7052)
         flush(7053)
         flush(7054)
         flush(7055)
      end if  
      end if  
      end if  

      return
      end
!
!-----------------------------------------------------------------------
!
!
! IMPORTANT NOTE: THIS SUBROUTINE IS CURRENTLY ONLY 
!    VALID FOR MONODISPERSE PARTICLE BEDS 
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!                          Binary part 
!-----------------------------------------------------------------------
!
!
      subroutine ppiclf_user_AM_Briney_Binary(i,iStage,
     >                 famx,famy,famz,rmass_add)
!
      implicit none
!
      include "PPICLF"
!
      integer i, j, k, l, n, jj
      integer*4 iStage
      real*8 rad
      real*8 dr_max
      real*8 IA, II
      real*8 famx, famy, famz, rmass_add
      real*8 alpha

! Declare functions
      real*8 IA_analytical, IA_numerical, II_analytical, II_numerical
!
! Code:
!
      ! particle radius
      rad = ppiclf_rprop(PPICLF_R_JDP,i) * 0.5d0
      
      ! ppiclf_nndist is neighbor width - MIN(user defined,4*P_dia)
      dr_max = ppiclf_nndist 

      ! In the example program binary_model.f90, the nearest
      ! neighbor calculations are done here at this point to
      ! calculate wdot. However, these have been previously 
      ! calculated in ppiclf_user.f in a separate do loop
      ! over 1 <= i <= ppiclf_npart

      ! Model only valid for local volume fraction
      ! less than 0.4, so we limit it here without
      ! over riding rphip
      alpha = min(0.4, rphip) ! limit alpha to mitigate misuse

      ! we need to make the numerical integration more efficient
      ! do it in a table and save it 
      if (dr_max / rad >= 3.49) then
         IA = IA_analytical(dr_max, rad, alpha) ! self acceleration
         II = II_analytical(dr_max, rad, alpha) ! neighbor acceleration (induced)
      else
         if (ppiclf_nid==0 .and. iStage==1) then
             print*, "***WARNING*** - NUMERICAL 
     >               FUNCTIONS USED IN ADDED MASS"
         endif
         IA = IA_numerical(dr_max, rad, alpha)
         II = II_numerical(dr_max, rad, alpha)
      end if

      do j=1,3
         jj = PPICLF_R_WDOTX + (j-1)
         Fam(j) = Fam(j) + IA*ppiclf_rprop(jj, i) ! added mass
         Fam(j) = Fam(j) + II*Wdot_neighbor_mean(j) 
     >                               / nneighbors ! induced added mass
      end do

      ! multiply by volume before adding unary term
      ! doing so here implies that the particle volume is
      ! the same for all particles; i.e., monodisperse packs
      do j=1,3
         Fam(j) = ppiclf_rprop(PPICLF_R_JVOLP,i)*Fam(j) 
      end do

      famx = Fam(1)
      famy = Fam(2)
      famz = Fam(3)
        
      if (ppiclf_debug==2) then
      if (ppiclf_nid .eq. 0 .or. ppiclf_np == 1) then
      if (i<=5 .and. iStage==1) then
         open(unit=7061,file='fort.7061',position='append')
         open(unit=7062,file='fort.7062',position='append')
         open(unit=7063,file='fort.7063',position='append')
         open(unit=7064,file='fort.7064',position='append')
         open(unit=7065,file='fort.7065',position='append')
         write(7060+i,*) i,iStage, ppiclf_time,     ! 0-2
     >      rphip, rmachp,                          ! 3-4
     >      IA, II,                                 ! 5-6
     >      ppiclf_rprop(PPICLF_R_WDOTX,i),         ! 7
     >      ppiclf_rprop(PPICLF_R_WDOTY,i),         ! 8
     >      ppiclf_rprop(PPICLF_R_WDOTZ,i),         ! 9
     >      Wdot_neighbor_mean(1:3), nneighbors,    ! 10-13
     >      famx, famy, famz,                       ! 14-16
     >      Fam(1), Fam(2), Fam(3)                  ! 17-19
         flush(7061)
         flush(7062)
         flush(7063)
         flush(7064)
         flush(7065)
      end if
      end if
      end if
        

      return
      end
!------------------------------------------------------------------------
!
! Created May 20, 2024
!
! Subroutine for added mass
!   also called the quasi-unsteady force,
!   or the inviscid unsteady force in case of the Euler equations
!
! Contains functions:
!    real*8 function B11_11(d, alpha)
!    real*8 function B11_22(d, alpha)
!    real*8 function B12_11(d, alpha)
!    real*8 function B12_22(d, alpha)
!    real*8 function IA_analytical(rmax, rad, alpha)
!    real*8 function II_analytical(rmax, rad, alpha)
!    real*8 function rdf_analytical(r, alpha)
!    real*8 function IA_numerical_integrand(d, rad, alpha)
!    real*8 function II_numerical_integrand(d, rad, alpha)
!    real*8 function IA_numerical(rmax, rad, alpha)
!    real*8 function II_numerical(rmax, rad, alpha)
!    
! Implementing Added Mass Algorithm from S.Briney (2024)
!  
! n       = number of points
! alpha   = volume fraction
! rad     = particle radius
! d       = center-to-center distance
! rmax    = center-to-center max neighbor distance
! R       = resistance matrix (output)
! x       = x_2 - x_1
! y       = y_2 - y_1
! z       = z_2 - z_1
! dr_max  = max interaction distance between particles considered 
! poins   = 3xn array of points x, y, z. Initialized as points(3,n)
!
! correction_analytical_always 
!    = if true, always use the analytical distant neighbor correction
!    > if false, use numerical when dr_max/rad < 3.49
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Functions for calculating the four curves that define the binary model
! These the the vol fraction-corrected binary added mass terms
!   taken from Appendix C of Briney et at (2024).
!
! parallel added mass
      real*8 function B11_11(d, alpha)
      implicit none
       
      ! input
      real*8 d, alpha
      
      ! local vars
      ! asymptotic solution (Beguin et al. 2016)
      ! empirical higher order terms
      ! empirical volume fraction correction
      real*8 asym, hot, alpha_corr
      
      ! Beguin et al. 2016
      asym = 0.5*(3.0/64.0 * (2.0/d)**6 + 9.0/256.0*(2.0/d)**8) 

      hot = 0.059/((d-1.9098)*(d-0.4782)**2 * (d**3))
      alpha_corr = (alpha*alpha - 0.0902*alpha) 
     >                 * 70.7731/((d+2.0936)**6)
      
      B11_11 = asym + hot + alpha_corr
      
      return
      end

! perpendicular added mass
      real*8 function B11_22(d, alpha)
      implicit none
      
      ! input
      real*8 d, alpha
      
      ! local vars
      ! hot combined with alpha correction in this case
      ! See comment in calc_B11_11
      real*8 asym, hot
      real*8 A, B
          
      ! Beguin et al. 2016
      asym = 0.5*(3.0/256.0 * (2.0/d)**6 + 3.0/256.0 
     >               *(2.0/d)**8) 
      
      A = 0.0003 + 0.0262*alpha*alpha - 0.0012*alpha
      B = 1.3127 + 1.0401*alpha*alpha - 1.2519*alpha
      hot = A / ((d-B)**6)
      
      B11_22 = asym + hot
      
      return
      end
      
! parallel induced added mass
      real*8 function B12_11(d, alpha)
      implicit none
         
      ! input
      real*8 d, alpha
      
      ! local vars
      real*8 asym, hot, alpha_corr
      
      ! Beguin et al. 2016
      asym = 0.5*(-3.0/8.0*(2.0/d)**3 - 3.0/512.0
     >               *(2.0/d)**9) 

      hot = -0.0006/((d-1.5428)**5)
      alpha_corr = -0.7913*alpha
     >               *exp(-(0.9801 - 0.1075*alpha)*d)
      
      B12_11 = asym + hot + alpha_corr

      return
      end
      
! perpendicular induced added mass
      real*8 function B12_22(d, alpha)
      implicit none
          
      ! input
      real*8 d, alpha
      
      ! local vars
      real*8 asym, hot
          
      ! Beguin et al. 2016
      asym = 0.5*(3.0/16.0 * (2.0/d)**3 + 3.0/4096.0 * (2.0/d)**9) 
          
      ! includes alpha correction
      hot = (-0.1985 + 16.7372*alpha)/(d**6)
     >           + (1.3907 - 48.2604*alpha)/(d**8) 
      
      B12_22 = asym + hot
      
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Analytical correction functions for distant neighbors. 
!     Assumes g(r) = 1. Good for maxr >= 3.5.
!      
!     These can be found in Appendix D of Briney etal (2024).
!
!     added mass
!     rmax = center-to-center maximum neighbor distance
!     rad = particle radius
!     alpha = volume fraction
      real*8 function IA_analytical(rmax, rad, alpha)
      implicit none
      
      ! input
      real*8 rmax, rad, alpha
      
      ! local vars
      real*8 term1, term2, term3, c, numerator, B11, B22, maxr
      
      maxr = rmax / rad
      
      ! B11
      term1 = (5.0*maxr**2 + 9.0)/(10.0*maxr**5)
          
      term2 = ((0.00720826 - 0.0150737 * maxr) * log(maxr - 1.9098)
     >         - 0.120023 * maxr * log(maxr - 0.4782) 
     >         + 0.057395 * log(maxr - 0.4782) 
     >         + (0.135097 * maxr - 0.0646033) * log(maxr) - 0.0861828)
     >                                           /(maxr - 0.4782)
      
      term3 = (alpha**2 - 0.0902*alpha) 
     >     * (0.333333 * maxr**2 + 0.348933 * maxr + 0.146105)
     >     /(maxr + 2.0936)**5
      
      B11 = term1 + term2 + term3
      
      ! B22
      term1 = (5.0*maxr**2 + 12.0)/(40.0*maxr**5)
      
      numerator = 0.0262*alpha**2 - 0.0012*alpha + 0.0003
      c = -1.0401*alpha**2 + 1.2519*alpha - 1.3127
      term2 = numerator * (10.0*maxr**2 + 5.0*maxr*c + c**2) 
     >             / (30.0*(maxr+c)**5)
      
      B22 = term1 + term2
      
      IA_analytical = alpha*(B11 + 2.0*B22)
          
!     write(1,*) IA_analytical, rmax, rad, alpha,
!    >                maxr, term1, term2, term3, B11, B22

      return
      end

      ! induced added mass
      real*8 function II_analytical(rmax, rad, alpha)
      implicit none
      
      ! input
      real*8 rmax, rad, alpha
      
      ! local vars
      real*8 term1, term2, term3, B11, B22, c, maxr, prefactor
      
      maxr = rmax / rad ! scale the filter width
      
      ! B11
      term1 = -1.0/(4.0*maxr**6) 
      
      term2 = -6.0e-4 * (0.5*maxr**2 - 0.516067*maxr + 0.199744)
     >             /(1.5482 - maxr)**4
      
      prefactor = -0.7913*alpha
      c = 0.9801 - 0.1075*alpha
      term3 = prefactor * (maxr *c *(maxr *c + 2.0) + 2.0) 
     >                       * exp((-maxr *c))/c**3
      
      B11 = term1 + term2 + term3
      
      ! B22
      term1 = 1.0 / (32.0 * maxr**6)
      
      term2 = (-0.1985 + 16.7372*alpha) / (3.0*maxr**3)
      
      term3 = (1.3907 - 48.2604*alpha) / (5.0*maxr**5)
      
      B22 = term1 + term2 + term3
      
      II_analytical = alpha*(B11 + 2.0*B22)
          
!          write(2,*) II_analytical, rmax, rad, alpha,
!     >                maxr, term1, term2, term3, B11, B22
      
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Numerical distant neighbor corrections
!      
!     Uses the rdf function
!      
!     r = dist between particles / diameter 
!          - normalization different than my convention!
!     phi = particle volume fraction
!     Trokhymchuk et al. 2005. The Journal of Chemical Physics.
!     https://doi.org/10.1063/1.1979488
!     erratum: https://doi.org/10.1063/1.2188941

      real*8 function rdf_analytical(r, alpha)
      implicit none
      
      ! inputs
      real*8 r, alpha
      
      ! local vars
      real*8  eta, d, mu, alpha0, beta0, denom_gamma, 
     >              gamma, omega, kappa
      real*8  alpha1, beta, rstar, gm, gsig, Brad, Arad, delta, 
     >               Crad, g2
      
      eta = alpha
      
      d   = (2.0*eta*(eta*eta - 3.0*eta - 3.0
     >     + sqrt(3.0*(eta**4-2.0*(eta**3)
     >     +eta**2+6.0*eta+3.0))))**(1.0/3.0)
      
      !mu    = (2.0*eta/(1.0-eta))*(-1.0-(d/(2.0*eta))-(eta/d))
      mu    = (2.0*eta/(1.0-eta))*(-1.0-(d/(2.0*eta))+(eta/d)) ! see erratum
      
      alpha0 = (2.0*eta/(1.0-eta))
     >         *(-1.0+(d/(4.0*eta))-(eta/(2.0*d)))
      
      beta0 = (2.0*eta/(1.0-eta))*sqrt(3.0)
     >         *(-(d/(4.0*eta))-(eta/(2.0*d)))
      
      denom_gamma = (alpha0**2 + beta0**2 - 2.0*mu*alpha0)
     >         *(1.0 + 0.5*eta) - mu*(1.0 + 2.0*eta) ! erratum
      
      gamma = atan(-(1.0/beta0)*((alpha0*(alpha0**2+beta0**2)
     >         -mu*(alpha0**2-beta0**2))*(1.0+0.5*eta)
     >         +(alpha0**2+beta0**2-mu*alpha0)*(1.0+2.0*eta))
     >         /(denom_gamma)) ! erratum
      
      !gamma = atan(-(1.0/beta0)*((alpha0*(alpha0**2+beta0**2)
      ! > -mu*(alpha0**2+beta0**2))*(1.0+0.5*eta)+
      ! > (alpha0**2+beta0**2-mu*alpha0)*(1.0+2.0*eta)));
      
      omega = -0.682*exp(-24.697*eta)+4.72+4.45*eta
      
      kappa = 4.674*exp(-3.935*eta)+3.536*exp(-56.27*eta)
      
      alpha1 = 44.554+79.868*eta+116.432*eta*eta-44.652*exp(2.0*eta)
      
      beta  = -5.022+5.857*eta+5.089*exp(-4.0*eta)
      
      rstar = 2.0116-1.0647*eta+0.0538*eta*eta
      
      gm    = 1.0286-0.6095*eta+3.5781*(eta**2)-21.3651*(eta**3)
     >             +42.6344*(eta**4)-33.8485*(eta**5)
      
      gsig  = (1.0/(4.0*eta))*(((1.0+eta+(eta**2)-(2.0/3.0)*(eta**3)
     >             -(2.0/3.0)*(eta**4))/((1.0-eta)**3))-1.0)
      
      Brad = (gm-(gsig/rstar)*exp(mu*(rstar - 1.0)))
     >           *rstar/(cos(beta*(rstar-1.0)+gamma)
     >           *exp(alpha1*(rstar-1.0))
     >           -cos(gamma)*exp(mu*(rstar-1.0)))
      
      Arad = gsig - Brad*cos(gamma)
      
      delta = -omega*rstar - atan((kappa*rstar+1.0)/(omega*rstar))
      
      Crad = rstar*(gm-1.0)*exp(kappa*rstar)/cos(omega*rstar+delta)
      
      g2 = (Arad/r)*exp(mu*(r-1.0)) + (Brad/r)
     >         *cos(beta*(r-1.0)+gamma)*exp(alpha1*(r-1.0))
      
      if (r > rstar) then
         g2 = 1.0 +(Crad/r)*cos(omega*r+delta)*exp(-kappa*r)
      else if (r < 1) then
         g2 = 0.0
      end if
      
      rdf_analytical = g2
      
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     d = center to center distance
!     rad = radius
!     alpha = volume fraction

      real*8 function IA_numerical_integrand(d, rad, alpha)
      implicit none
      
      ! input
      real*8 d, rad, alpha
      
      ! local vars
      real*8 g, f11, f22, rho, pi, w

      ! declare funcitons used
      real*8 rdf_analytical, B11_11, B11_22
      
      pi = 4.0*ATAN(1.0)
      
      g = rdf_analytical((d/(2.0*rad)), alpha)
      
      f11 = B11_11(d/rad, alpha)
      f22 = B11_22(d/rad, alpha)
      
      rho = alpha / (4.0/3.0*pi) ! number density
      w = 1.0/3.0 * rho * 4.0*pi*(d/rad)**2
      
      IA_numerical_integrand = w*g*(f11 + 2.0*f22)

      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     d = center to center distance
!     rad = radius
!     alpha = volume fraction

      real*8 function II_numerical_integrand(d, rad, alpha)
      implicit none
      
      ! input
      real*8 d, rad, alpha
      
      ! local vars
      real*8 g, f11, f22, rho, pi, w
          
      ! declare funcitons used
      real*8 rdf_analytical, B12_11, B12_22
      
      pi = 4.0*ATAN(1.0)
      
      g = rdf_analytical((d/(2.0*rad)), alpha)
      
      f11 = B12_11(d/rad, alpha)
      f22 = B12_22(d/rad, alpha)
      
      rho = alpha / (4.0/3.0*pi) ! number density
      w = 1.0/3.0 * rho * 4.0*pi*(d/rad)**2
      
      II_numerical_integrand = w*g*(f11 + 2.0*f22)

      return
      end
!
!-----------------------------------------------------------------------
!
      
      real*8 function IA_numerical(rmax, rad, alpha)
      implicit none
      
      ! input
      real*8 rmax, rad, alpha
      
      ! local vars
      real*8 maxr, dr, r, integral, f, coef(2)
      integer npts, i, j
          
      ! declare funcitons used
      real*8 IA_numerical_integrand
      
      npts = 1001 ! must be odd
      
      maxr = rad*20.0 ! essentially infinity
      
      coef(1) = 4.0
      coef(2) = 2.0
      
      if (maxr > rmax) then
         dr = (maxr - rmax) / (npts + 1)
     
         integral = 0.0
         r = rmax
      
         ! Simpson's rule
         ! endpoints first
         f = IA_numerical_integrand(r, rad, alpha)
         integral = integral + f
      
         f = IA_numerical_integrand(maxr, rad, alpha)
         integral = integral + f
      
         ! middle points
         do i=2,npts-1
            r = r + dr
            f = IA_numerical_integrand(r, rad, alpha)
      
            j = mod(i, 2) + 1
            integral = integral + coef(j)*f ! 4, 2, 4, 2, ...
         end do
      
         IA_numerical = dr * integral / 3.0
      
      else
         IA_numerical = 0.0
      end if
      
      return
      end
      
!
!-----------------------------------------------------------------------
!
      
      real*8 function II_numerical(rmax, rad, alpha)
      implicit none
      
      ! input
      real*8 rmax, rad, alpha
      
      ! local vars
      real*8 maxr, dr, r, integral, f, coef(2)
      integer npts, i, j
      
      ! declare funcitons used
      real*8 II_numerical_integrand
          
      npts = 1001 ! must be odd
      
      maxr = rad*20.0 ! essentially infinity
      
      coef(1) = 4.0
      coef(2) = 2.0
      
      if (maxr > rmax) then
         dr = (maxr - rmax) / (npts + 1)
      
         integral = 0.0
         r = rmax
      
         ! Simpson's rule
         ! endpoints first
         f = II_numerical_integrand(r, rad, alpha)
         integral = integral + f
      
         f = II_numerical_integrand(maxr, rad, alpha)
         integral = integral + f
      
         ! middle points
         do i=2,npts-1
            r = r + dr
            f = II_numerical_integrand(r, rad, alpha)
      
            j = mod(i, 2) + 1
            integral = integral + coef(j)*f ! 4, 2, 4, 2, ...
         end do
      
         II_numerical = dr * integral / 3.0
      
      else
         II_numerical = 0.0
      end if
      
      return 
      end
!------------------------------------------------------------------------
!
! Created May 20, 2024
!
! Subroutine for added mass
!   also called the quasi-unsteady force,
!   or the inviscid unsteady force in case of the Euler equations
!
! Contains subroutines:
!     rotation_matrix(x, y, z, Q)
!     resistance_pair(x, y, z, alpha, rad, R)
!
! Implementing Added Mass Algorithm from S.Briney (2024)
!  
! n       = number of points
! alpha   = volume fraction
! rad     = particle radius
! d       = center-to-center distance
! rmax    = center-to-center max neighbor distance
! R       = resistance matrix (output)
! x       = x_2 - x_1
! y       = y_2 - y_1
! z       = z_2 - z_1
! dr_max  = max interaction distance between particles considered 
! poins   = 3xn array of points x, y, z. Initialized as points(3,n)
!
! correction_analytical_always 
!    = if true, always use the analytical distant neighbor correction
!    > if false, use numerical when dr_max/rad < 3.49
!
!
!-----------------------------------------------------------------------
!
! Calculates the rotation matrix Q to align the coordinate 
!    system such that the point (x, y, z) is on the x-axis

      subroutine rotation_matrix(x, y, z, Q)
       
      ! inputs
      real*8 x, y, z
       
      ! output rotation matrix
      real*8 Q(3, 3)
       
      ! local vars
      real*8 gamma, beta
           
      gamma =  atan2(y, x)
      beta  = -atan2(z, sqrt(x*x + y*y))
       
      Q(1, 1) = cos(beta)*cos(gamma)
      Q(1, 2) = cos(beta)*sin(gamma)
      Q(1, 3) = -sin(beta)
       
      Q(2, 1) = -sin(gamma)
      Q(2, 2) = cos(gamma)
      Q(2, 3) = 0.0
       
      Q(3, 1) = sin(beta)*cos(gamma)
      Q(3, 2) = sin(beta)*sin(gamma)
      Q(3, 3) = cos(beta)
       
      return
      end
!
!-----------------------------------------------------------------------
!
      ! This is the one of the functions the user should typically call.
      ! Returns the resistance matrix for a 2 particle system 
      !    of arbitrary orientation
      ! x = x2 - x1, y = y2 - y1, z = z2 - z1 (relative position)
      ! rad = particle radius
      ! R = resistance matrix (output)

      subroutine resistance_pair(x, y, z, alpha, rad, R)
      implicit none
       
      ! input: x, y, z of second particle relative to the second particle
      ! x = x2 - x1, etc.
      ! rad = particle radius (monodisperse)
      real*8 x, y, z, alpha, rad
       
      ! output: resistance matrix
      real*8 R(6, 6)
       
      ! local vars
      real*8, dimension(3, 3) :: Q, Qt, B11, B12
       
      real*8 dist
       
      ! declare functions
      real*8 B11_11, B11_22, B12_11, B12_22
       
      ! get rotation matrix Q
      call rotation_matrix(x, y, z, Q)
       
      ! normalize the distance by the particle radius
      dist = sqrt(x*x + y*y + z*z) / rad 
       
      ! initially set coefficients to 0
      B11(1:3, 1:3) = 0.0
      B12(1:3, 1:3) = 0.0
      
      ! self acceleration
      B11(1, 1) = B11_11(dist, alpha)
      B11(2, 2) = B11_22(dist, alpha)
      B11(3, 3) = B11(2, 2)
       
      ! neighbor acceleration
      B12(1, 1) = B12_11(dist, alpha)
      B12(2, 2) = B12_22(dist, alpha)
      B12(3, 3) = B12(2, 2)
       
      ! rotate
      Qt  = transpose(Q)
      B11 = matmul(Qt, matmul(B11, Q))
      B12 = matmul(Qt, matmul(B12, Q))
            
      ! store resitance matrix
      R(1:3, 1:3) = B11
      R(1:3, 4:6) = B12
       
      R(4:6, 1:3) = B12
      R(4:6, 4:6) = B11
          
      ! write(7059,*) x, y, z, alpha, rad, dist 
      !
      ! write(7060,*) B11(1,1:3), B11(2,1:3), B11(3,1:3)
      ! write(7061,*) B12(1,1:3), B12(2,1:3), B12(3,1:3)
      ! write(7062,*) C11(1,1:3), C11(2,1:3), C11(3,1:3)
      ! write(7063,*) C12(1,1:3), C12(2,1:3), C12(3,1:3)
      ! write(7064,*) Qt(1,1:3), Qt(2,1:3), Qt(3,1:3)
      !
      ! write(7069,*) R(1,1:6)
      ! write(7069,*) R(2,1:6)
      ! write(7069,*) R(3,1:6)
      ! write(7069,*) R(4,1:6)
      ! write(7069,*) R(5,1:6)
      ! write(7069,*) R(6,1:6)
      ! write(7069,*) "------------------------------------"

      return
      end
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer models for non-burning
!    particles
! May update if we include particle combustion
!
! if heattransfer_flag = 0  ignore heat transfer
!                      = 1  Stokes
!                      = 2  Ranz-Marshall (1952)
!                      = 3  Gunn (1977)
!                      = 4  Fox (1978)
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_HT_driver(i,qq)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 qq, Nuss, Q_conv

!
! Code:
!
      Q_conv = rpi*rkappa*dp*(ppiclf_rprop(PPICLF_R_JT,i) -
     >                          ppiclf_y(PPICLF_JT,i) )

      Nuss = 0.0d0
      if (heattransfer_flag == 1) then
         call HTModel_Stokes(i,Nuss)
      elseif (heattransfer_flag == 2) then
         call HTModel_RM(i,Nuss)
      elseif (heattransfer_flag == 3) then
         call HTModel_Gunn(i,Nuss)
      elseif (heattransfer_flag == 4) then
         call HTModel_Fox(i,Nuss)
      else
         call ppiclf_exittr('Unknown heat transfer model$', 0.0d0, 0)
      endif

      qq = qq + Q_conv*Nuss


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: Stokes limit
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_Stokes(i,Nuss)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 Nuss
!
! Code:
!
      Nuss = 2.0d0

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: Ranz-Marshall
!    define Nusselt number Nu = Nu(Pr,Re)
!
! (1) Ling et al (2012) 
!     "Interaction of a planar shock wave with a dense particle
!     curtain: Modeling and experiments"
!     Physics of Fluids, Vol. 24, 113301
! (2) Durant et al (2022) 
!     "Explosive dispersal of particles in high speed environments"
!     Journal of Applied Physics, Vol. 132, 184902
!
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_RM(i,Nuss)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 Nuss
!
! Code:
!
      Nuss = 2.0d0+0.6d0*(rep**0.5d0)*(rpr**OneThird)

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: Gunn
!    define Nusselt number Nu = Nu(Pr,Re,phi)
!
! (1)  Gunn (1978)
!        "Transfer of heat or mass to particles in
!        fixed and fluidised beds"
!        Int. J. Heat Mass Transfer, Vol. 21, pp. 467-476
! (2)  Houim and Oran (2016)
!        "A multiphase model for compressible granular-gaseous
!        flows: formulation and initial tests"
!        J. Fluid Mech., Vol. 789, pp. 166-220
! (3) Boniou and Fox (2023)
!        "Shock-particle-curtain-interaction study with
!        a hyperbolic two-fluid model: Effect of particle
!        force models"
!        Int. Journal of Mutiphase Flow, Vol. 169, 104591
!
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_Gunn(i,Nuss)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 vg
      real*8 Nuss
!
! Code:
!
      ! define bed voidage = ratio of free volume avaliable
      ! for flow to the total volume of bed; aka, volume
      ! fraction of the gas phase
      ! vg = 1.0 - rphip == rphif
      vg = rphif

      ! define Nusselt number Nu = Nu(Pr,Re,phi)
      Nuss = (7.0d0-10.0d0*vg+5.0d0*vg*vg)
     >           *(1.0d0+0.7d0*(rep**0.2d0)*(rpr**OneThird))
     >     + (1.33d0-2.4d0*vg+1.2d0*vg*vg)*(rep**0.7d0)*(rpr**OneThird)

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: Fox
!    define Nusselt number Nu = Nu(Pr,Re,M)
!
! (1)  Fox, Rackett, and Nicholls (1978)
!        "Shock wave ignition of magnesium powders"
!        in Proc. 11th Int. Symp. Shock Tubes and Waves.
!        University of Washington Press, Seattle, WA, pp. 262-268
! (2)  Ling, Haselbacher, and Balachandar (2011)
!        "Importance of unsteady contributions to force and heating
!        for particles in compressible flows. Part I: Modeling and
!        anaylsis for shock-particle interaction"
!        Int. Journal of Multiphase Flow, Vol. 37, pp. 1026-1044
!
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_Fox(i,Nuss)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 Nuss
!
! Code:
!
      ! define Nusselt number Nu = Nu(Pr,Re,M)
      Nuss = 2.0d0*exp(-rmachp)/(1.0d0+17.0d0*rmachp/rep)
     >     + 0.495d0*(rpr**OneThird)*(rep**0.55d0)*
     >       ((1.0d0+0.5d0*exp(-17.0d0*rmachp/rep))/1.5d0)

      return
      end
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for computing the torque terms on the
!    RHS of the angular velocity equations
!
!
! if collisional_flag = 1  F = Fn
!                     = 2  F = Fn + Ft + Tt
!                     = 3  F = Fn + Ft + Tt + Th + Tr
! where Tt = collisional torque
!       Th = hydrodynamic torque
!       Tr = rolling torque
!
! Note that the collisional and rolling torques are due to
!     particle-particle interactions and thus are evaulated
!     in ppiclf_user_EvalNearestNeighbor. Therefore, only the
!     hydrodynamic torque is left to be calculated.
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_Torque_driver(i,iStage,taux,tauy,tauz,
     >                                 taux_hydro,tauy_hydro,tauz_hydro)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i,iStage
      real*8 taux, tauy, tauz
      real*8 taux_hydro, tauy_hydro, tauz_hydro
      real*8 taux_undist, tauy_undist, tauz_undist
      real*8 rmass_local

!
! Code:
!
      taux_hydro = 0.0d0
      tauy_hydro = 0.0d0
      tauz_hydro = 0.0d0
      taux_undist = 0.0d0
      tauy_undist = 0.0d0
      tauz_undist = 0.0d0

      if (collisional_flag >= 3) then
         call Torque_Hydro(i,taux_hydro,tauy_hydro,tauz_hydro)
         call Torque_Undisturbed(i,taux_undist,tauy_undist,tauz_undist)
      endif

      taux = taux + taux_hydro + taux_undist
      tauy = tauy + tauy_hydro + tauy_undist
      tauz = tauz + tauz_hydro + tauz_undist

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for hydrodynamic torque
!
!
!-----------------------------------------------------------------------
!
      subroutine Torque_Hydro(i,taux_hydro,tauy_hydro,tauz_hydro)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 taux_hydro, tauy_hydro, tauz_hydro
      real*8 omgrx, omgry, omgrz, omgr_mag
      real*8 Ct1, Ct2, Ct3, Ct
      real*8 reyr, beta, rIp, factor

!
! Code:
!
      ! Compute relative angular velocity components
      !    and magnitude
      omgrx = 0.5d0*ppiclf_rprop(PPICLF_R_JXVOR,i)
     >          - ppiclf_y(PPICLF_JOX,i)
      omgry = 0.5d0*ppiclf_rprop(PPICLF_R_JYVOR,i)
     >          - ppiclf_y(PPICLF_JOY,i)
      omgrz = 0.5d0*ppiclf_rprop(PPICLF_R_JZVOR,i)
     >          - ppiclf_y(PPICLF_JOZ,i)
      omgr_mag = sqrt(omgrx*omgrx + omgry*omgry + omgrz*omgrz)

      ! Particle rotational Reynolds number
      reyr = rhof*dp*dp*omgr_mag/(4.0d0*rmu)
      reyr = max(reyr,0.001d0)

      ! Compute the hydrodynamic torque parameter Ct=Ct(Re_r)
      if (reyr < 1) then
         Ct1 = 0.0d0
         Ct2 = 16.0d0*rpi
         Ct3 = 0.0d0
      elseif (reyr < 10) then
         Ct1 = 0.0d0
         Ct2 = 16.0d0*rpi
         Ct3 = 0.0418d0
      elseif (reyr < 20) then
         Ct1 = 5.32d0
         Ct2 = 37.2d0
         Ct3 = 0.0d0
      elseif (reyr < 50) then
         Ct1 = 6.44d0
         Ct2 = 32.2d0
         Ct3 = 0.0d0
      elseif (reyr < 100) then
         Ct1 = 6.45d0
         Ct2 = 32.1d0
         Ct3 = 0.0d0
      else
         call ppiclf_exittr('Re rotational too large$', reyr, 0)
      endif

      Ct = Ct1/sqrt(reyr) + Ct2/Reyr + Ct3*reyr

      ! Now compute hydrodynamic torque components
      beta = rhop/rhof
      rIp  = rmass*dp*dp/10.0d0
      factor = rIp*60.0d0*Ct*omgr_mag/(64.0d0*rpi*beta)

      taux_hydro = factor*omgrx
      tauy_hydro = factor*omgry
      tauz_hydro = factor*omgrz


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created April 01, 2025
!
! Subroutine for undisturbed torque
!
!
!-----------------------------------------------------------------------
!
      subroutine Torque_Undisturbed(i, 
     >           taux_undist,tauy_undist,tauz_undist)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 taux_undist, tauy_undist, tauz_undist
      real*8 rIf

!
! Code:
!

      ! Moment of interia with respect to gas
      rIf = rhof*dp*dp*ppiclf_rprop(PPICLF_R_JVOLP,i)/10.0d0

      ! Undisturbed torque component
      ! Written using angular velocity = 0.5*vorticity
      taux_undist = 0.5d0*rIf*ppiclf_rprop(PPICLF_R_JSDOX,i)
      tauy_undist = 0.5d0*rIf*ppiclf_rprop(PPICLF_R_JSDOY,i)
      tauz_undist = 0.5d0*rIf*ppiclf_rprop(PPICLF_R_JSDOZ,i)


      return
      end
!-----------------------------------------------------------------------
!
! Created April 01, 2025
!
! Subroutine for computing the lift terms
!    Lift components first requires computing gas-phase vorticity
!    and particle angular velocity
!
!
! if collisional_flag = 4  Add Saffman and Magnus lift
!
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_Lift_driver(i,iStage,liftx,lifty,liftz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i,iStage
      real*8 liftx, lifty, liftz

!
! Code:
!
      liftx = 0.0d0
      lifty = 0.0d0
      liftz = 0.0d0

      if (collisional_flag >= 4) then
         call Lift_Saffman(i,liftx,lifty,liftz)
         call Lift_Magnus (i,liftx,lifty,liftz)
      endif


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created April 01, 2025
!
! Subroutine for Saffman lift - shear-induced lift
!
! Requires gas-phase vorticity to be computed
! Valid for Rep < 50 and omg* < 0.8 (see Loft, "Lift of a spherical
!    particle subject to vorticity and/or spin", AIAA J., 
!    Vol. 46,  pp. 801-809, 2008)
!      
! References:
! 1) Fundamentals of Dispersed Multiphase Flows (S.Balachandar), Chap.5
!
!-----------------------------------------------------------------------
!
      subroutine Lift_Saffman(i,liftx,lifty,liftz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 liftx, lifty, liftz
      real*8 omgx, omgy, omgz, omg_mag, omg_star
      real*8 epi, Jepi
      real*8 d1, d2, d3
      real*8 factor
      real*8 elx, ely, elz, elm, ielm

!
! Code:
!
      if (vmag .lt. 1.d-8) return

      ! Compute gas-phase vorticity components and magnitude
      omgx = ppiclf_rprop(PPICLF_R_JXVOR,i)
      omgy = ppiclf_rprop(PPICLF_R_JYVOR,i)
      omgz = ppiclf_rprop(PPICLF_R_JZVOR,i)
      omg_mag = sqrt(omgx*omgx + omgy*omgy + omgz*omgz)

      ! Compute Mei correction
      omg_star = omg_mag*dp/vmag
      epi = sqrt(omg_star/rep)

      d1 = 1.0d0 + tanh(2.5d0*(log10(epi)+0.191d0))
      d2 = 0.667d0 + tanh(6.0d0*epi-1.92d0)
      Jepi = 0.3d0*d1*d2

      factor = 1.615d0*rmu*(dp*dp)*vmag*sqrt(omg_mag/rnu)

      ! Compute lift components
      elx = vy*omgz - vz*omgy
      ely = vz*omgx - vx*omgz
      elz = vx*omgy - vy*omgx
      elm = sqrt(elx*elx + ely*ely +elz*elz)
      elm = max(1.0d-20,elm)
      ielm = 1.0d0/elm

      liftx = liftx + factor*Jepi*elx*ielm
      lifty = lifty + factor*Jepi*ely*ielm
      liftz = liftz + factor*Jepi*elz*ielm


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created April 01, 2025
!
! Subroutine for Magnus lift - lift induced by particle rotation
!
! Requires particle angular velocity to be calculated
!      
! References:
! 1) Fundamentals of Dispersed Multiphase Flows (S.Balachandar), Chap.5
! 2) Loth and Drogan, "An equation of motion for particles of finite
! Reynolds number and size", (2009). 
!      
!-----------------------------------------------------------------------
!
      subroutine Lift_Magnus(i,liftx,lifty,liftz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 liftx, lifty, liftz
      real*8 omgx, omgy, omgz, omg_mag, omg_star
      real*8 epi, CL
      real*8 d1
      real*8 factor
      real*8 elx, ely, elz

!
! Code:
!
      if (vmag .lt. 1.d-8) return

      ! Compute particle angular velocity
      omgx = ppiclf_y(PPICLF_JOX,i)
      omgy = ppiclf_y(PPICLF_JOY,i)
      omgz = ppiclf_y(PPICLF_JOZ,i)
      omg_mag = sqrt(omgx*omgx + omgy*omgy + omgz*omgz)

      ! Correction to lift
      omg_star = omg_mag*dp/vmag
      epi = omg_star
      d1 = 0.675d0 + 0.15d0*(1.0d0 + tanh(0.28d0*(epi-2.0d0)))
      CL = 1.0d0 - d1*tanh(0.18*sqrt(rep))

      factor = 0.125d0*dp*dp*dp*rhof

      ! Compute lift components
      elx = vy*omgz - vz*omgy
      ely = vz*omgx - vx*omgz
      elz = vx*omgy - vy*omgx

      liftx = liftx + factor*CL*elx
      lifty = lifty + factor*CL*ely
      liftz = liftz + factor*CL*elz


      return
      end
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for viscous unsteady force with history kernel
!
! Mei-Adrian history kernel
!
! Copied from either files in rocintereact/
!   INRT_CalcDragUnsteady_AMImplicit.F90
!   INRT_CalcDragUnsteady_AMExplicit.F90
!
! The number of time steps kept for the history
!   kernel is set in libpicl/ppiclF/source/PPICLF_STD.h
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_VU_Rocflu(i,iStage,fvux,fvuy,fvuz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage, iT
      real*8 fvux,fvuy,fvuz
      real*8 time,fH,factor,A,B,kernelVU

!
! Code:
!
      fvux = 0.0d0
      fvuy = 0.0d0
      fvuz = 0.0d0
      iT   = 1
      time = 0.0d0

      fH     = 0.75d0 + .105d0*reyL
      ! Sangani's volume fraction correction for dilute random arrays
      ! Capping volume fraction at 0.5 
      factor = 3.0d0*rpi*rnu*dp*ppiclf_dt*(1.0+2.28*min(rphip,0.5))

      if (ppiclf_nTimeBH > 1) then
         do iT = 2,ppiclf_nTimeBH-1
            time = ppiclf_timeBH(iT)

            A  = (4.0d0*rpi*time*rnu/(dp**2))**(.25d0)
            B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >                 (0.5d0*dp*rnu*(fH**3)))**(.5d0)

            kernelVU = factor*(A+B)**(-2)

            fvux = fvux + kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JX,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JX,iT,i) )
            fvuy = fvuy + kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JY,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JY,iT,i) )
            fvuz = fvuz + kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JZ,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JZ,iT,i) )
         enddo

         iT = ppiclf_nTimeBH
         time = ppiclf_timeBH(iT)

         A  = (4.0d0*rpi*time*rnu/(dp**2))**(.25d0)
         B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >                 (0.5d0*dp*rnu*(fH**3)))**(.5d0)

         kernelVU = 0.5d0*factor*(A+B)**(-2)

         fvux = fvux + kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JX,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JX,iT,i) )
         fvuy = fvuy + kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JY,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JY,iT,i) )
         fvuz = fvuz + kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JZ,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JZ,iT,i) )
      endif


      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Shift arrays for Viscous Unsteady Force
!
! See rocpart/PLAG_RFLU_ShiftUnsteadyData.F90
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_ShiftUnsteadyData
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iT

!
! Code:
!
      do i=1,ppiclf_npart
         do iT = ppiclf_nUnsteadyData,2,-1
            ppiclf_drudtMixt(PPICLF_JX,iT,i) = 
     >                      ppiclf_drudtMixt(PPICLF_JX,iT-1,i)
            ppiclf_drudtMixt(PPICLF_JY,iT,i) = 
     >                      ppiclf_drudtMixt(PPICLF_JY,iT-1,i)
            ppiclf_drudtMixt(PPICLF_JZ,iT,i) = 
     >                      ppiclf_drudtMixt(PPICLF_JZ,iT-1,i)

            ppiclf_drudtPlag(PPICLF_JX,iT,i) = 
     >                      ppiclf_drudtPlag(PPICLF_JX,iT-1,i)
            ppiclf_drudtPlag(PPICLF_JY,iT,i) = 
     >                      ppiclf_drudtPlag(PPICLF_JY,iT-1,i)
            ppiclf_drudtPlag(PPICLF_JZ,iT,i) = 
     >                      ppiclf_drudtPlag(PPICLF_JZ,iT-1,i)
         enddo
      enddo


      if (ppiclf_nTimeBH < ppiclf_nUnsteadyData) then
            ppiclf_nTimeBH = ppiclf_nTimeBH + 1
      endif

      do iT = ppiclf_nTimeBH,2,-1
            ppiclf_timeBH(it) = ppiclf_timeBH(iT-1) + ppiclf_dt
      enddo


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Update arrays for Viscous Unsteady Force for JT=1 (current time step)
!
! See libpicl/user_files/ppiclf_user_AddedMass.f
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_UpdatePlag(i)
!
      implicit none
!
      include "PPICLF"
!
      integer*4 i
      real*8 SDrho
      real*8 ug,vg,wg
      real*8 up,vp,wp
      real*8 vgradrho

!
! Code:
!
      SDrho = ppiclf_rprop(PPICLF_R_JRHSR,i)
     >         + ppiclf_y(PPICLF_JVX,i) * ppiclf_rprop(PPICLF_R_JPGCX,i)
     >         + ppiclf_y(PPICLF_JVY,i) * ppiclf_rprop(PPICLF_R_JPGCY,i)
     >         + ppiclf_y(PPICLF_JVZ,i) * ppiclf_rprop(PPICLF_R_JPGCZ,i)

      ! 03/11/2025 - Thierry - substantial derivative from Rocflu is
      !              weighted by \phi^g.
      ! d(rho^g phi^g)/dt = rho^g * d(phi^g)/dt + phi^g * d(rho^g)/dt
      !                   = phi^g * d(rho^g)/dt
      !
      !     d(rho^g)/dt   = SDrho = d(rho phi^g)/dt / phi^g
      SDrho = SDrho / (rphif)

      ! 03/23/2025 - TLJ - added extra term involving grad(rhog)
      vgradrho = vx*ppiclf_rprop(PPICLF_R_JRHOGX,i) +
     >           vy*ppiclf_rprop(PPICLF_R_JRHOGY,i) +
     >           vz*ppiclf_rprop(PPICLF_R_JRHOGZ,i)

      ug = ppiclf_rprop(PPICLF_R_JUX,i)
      vg = ppiclf_rprop(PPICLF_R_JUY,i)
      wg = ppiclf_rprop(PPICLF_R_JUZ,i)
      up = ppiclf_y(PPICLF_JVX,i)
      vp = ppiclf_y(PPICLF_JVY,i)
      wp = ppiclf_y(PPICLF_JVZ,i)

      ! D(rhog*ug)/Dt
      ppiclf_drudtMixt(PPICLF_JX,1,i) =
     >   ug*(SDrho+vgradrho) + rhof*ppiclf_rprop(PPICLF_R_JSDRX,i)
      ppiclf_drudtMixt(PPICLF_JY,1,i) =
     >   vg*(SDrho+vgradrho) + rhof*ppiclf_rprop(PPICLF_R_JSDRY,i)
      ppiclf_drudtMixt(PPICLF_JZ,1,i) =
     >   wg*(SDrho+vgradrho) + rhof*ppiclf_rprop(PPICLF_R_JSDRZ,i)

      ! d(rhog*up)/dt
      ppiclf_drudtPlag(PPICLF_JX,1,i) =
     >   up*SDrho + rhof*ppiclf_ydot(PPICLF_JVX,i)
      ppiclf_drudtPlag(PPICLF_JY,1,i) =
     >   vp*SDrho + rhof*ppiclf_ydot(PPICLF_JVY,i)
      ppiclf_drudtPlag(PPICLF_JZ,1,i) =
     >   wp*SDrho + rhof*ppiclf_ydot(PPICLF_JVZ,i)


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Sets drudtMixt and drudtPlag from rprop3
! Needed for proper particle tracking
! Load communication buffers rprop3 into particle data
! See rocpart/PLAG_RFLU_ModComm.F90:
!     SUBROUTINE PLAG_RFLU_LoadBuffersSend(pRegion)
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_prop2plag
!
      implicit none
!
      include "PPICLF"
!
      integer*4 i,k,ic,iT
!
! Code:
!
      do i=1,ppiclf_npart
         k = 0
         do ic = 1,3
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_drudtMixt(ic,iT,i) = ppiclf_rprop3(k,i)
         enddo
         enddo
         do ic = 1,3
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_drudtPlag(ic,iT,i) = ppiclf_rprop3(k,i)
         enddo
         enddo
      enddo


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Sets rprop3 from drudtMixt and drudtPlag
! Needed for proper particle tracking
! Load particle data into communication buffers rprop3
! See rocpart/PLAG_RFLU_ModComm.F90:
!     SUBROUTINE PLAG_RFLU_UnloadBuffersRecv(pRegion)
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_plag2prop
!
      implicit none
!
      include "PPICLF"
!
      integer*4 i,k,ic,iT
!
! Code:
!
      do i=1,ppiclf_npart
         k = 0
         do ic = 1,3
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_drudtMixt(ic,iT,i)
         enddo
         enddo
         do ic = 1,3
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_drudtPlag(ic,iT,i)
         enddo
         enddo
      enddo


      return
      end
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
      subroutine ppiclf_user_EvalNearestNeighbor
     >    (i,j,yi,rpropi,yj,rpropj)
!
      implicit none
!
      include "PPICLF"
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

      integer*4 i
      integer*4 j
      real*8 yi    (PPICLF_LRS)    
      real*8 rpropi(PPICLF_LRP)
      real*8 yj    (PPICLF_LRS)    
      real*8 rpropj(PPICLF_LRP)
!
! Internal:
!
      real*8 pi2, rthresh, rxdiff, rydiff, rzdiff, rdiff, rm1, rm2,
     >       rmult, eta_n, rbot, rn_12x, rn_12y, rn_12z, rdelta12,
     >       rv12_mag, rv12_mage, rksp_max, rnmag, rksp_wall, rextra,
     >       JDP_i,JDP_j
      real*8 eps

      ! Sam - from TLJ for box filter
      integer*4 ifilt
      real*8 adptfilter, dpl, phip, dist2, xdist2, ydist2, zdist2
      real*8 dist, rsig
      real*8 sig2, gkern, pi

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
      real*8 rp, yp, zp, vp, wp, 
     >       rp_new, yp_new, zp_new, vp_new, wp_new,
     >       rbound, urp, thetap

!
! Code:
!
      pi2  = rpi*rpi

      ! other particles
      if (j .ne. 0) then
         !Added spload and radius factor

         ! Compute mean particle diameter between i and j; delta_{ij}
         rthresh  = 0.5d0*(rpropi(PPICLF_R_JDP) + rpropj(PPICLF_R_JDP))

         ! Compute vector components and distance between 
         !    centers of particles i and j; D_{ij}
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
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
            rad = 0.5d0*rpropi(PPICLF_R_JDP)

            call resistance_pair(rxdiff1, rydiff1, rzdiff1, 
     >           alpha_local, rad, R_pair)
            
            ! accumulate number of neighbors
            nneighbors = nneighbors + 1
            
            do k=1,3
               do l=1,3
                 ll = PPICLF_R_WDOTX + (l-1)
                 ! added mass
                 Fam(k) = Fam(k) + R_pair(k,l)   * rpropi(ll)
                 ! induced added mass
                 Fam(k) = Fam(k) + R_pair(k,l+3) * rpropj(ll)
               end do ! l-loop

               ! accumulate neighbor acceleration
               kk = PPICLF_R_WDOTX + (k-1)
               Wdot_neighbor_mean(k) = Wdot_neighbor_mean(k)
     >                               + rpropj(kk)
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
            r1 = 0.5d0*rpropi(PPICLF_R_JDP)
            r2 = 0.5d0*rpropj(PPICLF_R_JDP)
            Rstar = r1*r2/(r1+r2)
            ksp2 = (4.0d0/3.0d0)*Estar*sqrt(Rstar)
            ksp2 = ksp2*sqrt(abs(rdiff-rthresh))
            ! kn = min(k1,k2)
            ksp_min = min(ksp1,ksp2)

            rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
            rm2 = rpropj(PPICLF_R_JRHOP)*rpropj(PPICLF_R_JVOLP)
         
            rmult = (rm1*rm2)/(rm1+rm2)
            eta_n = -2.0d0*sqrt(ksp_min)*log(erest)
     >              /sqrt(log(erest)**2+pi2)
     >              *sqrt(rmult)

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
            u12x = yi(PPICLF_JVX)-yj(PPICLF_JVX)
            u12y = yi(PPICLF_JVY)-yj(PPICLF_JVY)
            u12z = yi(PPICLF_JVZ)-yj(PPICLF_JVZ)

            if (collisional_flag>=2) then
               ! Add contribution from angular velocity
               rad1 = 0.5d0*rpropi(PPICLF_R_JDP)
               rad2 = 0.5d0*rpropj(PPICLF_R_JDP)
               A12x = rad1*yi(PPICLF_JOX) + rad2*yj(PPICLF_JOX)
               A12y = rad1*yi(PPICLF_JOY) + rad2*yj(PPICLF_JOY)
               A12z = rad1*yi(PPICLF_JOZ) + rad2*yj(PPICLF_JOZ)

               u12x = u12x + (A12y*rn_12z - A12z*rn_12y)
               u12y = u12y + (A12z*rn_12x - A12x*rn_12z)
               u12z = u12z + (A12x*rn_12y - A12y*rn_12x)
            endif

            ! Compute (u_ij \cdot n_ij)
            rv12_mag = u12x*rn_12x
     >               + u12y*rn_12y
     >               + u12z*rn_12z
         
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
               rad1 = 0.5d0*rpropi(PPICLF_R_JDP)
               tcx = rad1*(rn_12y*Ftz - rn_12z*Fty)
               tcy = rad1*(rn_12z*Ftx - rn_12x*Ftz)
               tcz = rad1*(rn_12x*Fty - rn_12y*Ftx)

               if (collisional_flag>=3) then
                  ! Add Rolling torque contribution
                  thetar = 0.06  ! Needs to be calibrated
                  dp1 = rpropi(PPICLF_R_JDP)
                  dp2 = rpropj(PPICLF_R_JDP)
                  r12 = 0.5d0*(dp1*dp2)/(dp1+dp2)
                  omgrx = yi(PPICLF_JOX) - yj(PPICLF_JOX)
                  omgry = yi(PPICLF_JOY) - yj(PPICLF_JOY)
                  omgrz = yi(PPICLF_JOZ) - yj(PPICLF_JOZ)
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
            ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_ydotc(PPICLF_JVX,i)
     >                                 + rnmag*rn_12x
     >                                 + Ftmin*rt_12x
            ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_ydotc(PPICLF_JVY,i)
     >                                 + rnmag*rn_12y
     >                                 + Ftmin*rt_12y
            ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_ydotc(PPICLF_JVZ,i)
     >                                 + rnmag*rn_12z
     >                                 + Ftmin*rt_12z

            ! Particle angular velocities
            ppiclf_ydotc(PPICLF_JOX,i) = ppiclf_ydotc(PPICLF_JOX,i)
     >                                 + tcx + trx
            ppiclf_ydotc(PPICLF_JOY,i) = ppiclf_ydotc(PPICLF_JOY,i)
     >                                 + tcy + try
            ppiclf_ydotc(PPICLF_JOZ,i) = ppiclf_ydotc(PPICLF_JOZ,i)
     >                                 + tcz + trz

         end if ! rdiff lt rthresh

!-----------------------------------------------------------------------
!
         ! Feedback fluctuation mean

         dist2 = MAX(ppiclf_filter(1),ppiclf_filter(2),ppiclf_filter(3))

         ! Box filter half-width dist2
         IF(qs_fluct_filter_adapt_flag.NE.0) THEN
            ! Adaptive filter defined wrt particle i
            ! Used for adaptive box or gaussian
            dpl = rpropi(PPICLF_R_JDP)
            phip = rpropi(PPICLF_R_JPHIP)
            adptfilter = ( 10.*(dpl**3)/max(1.e-4,phip) )**(1./3.)
            adptfilter = adptfilter/2.0
            IF(adptfilter .GT. dist2) dist2 = adptfilter
         END IF

         ! Check if particle lies inside box or gaussian filter
         xdist2 = abs(yi(PPICLF_JX)-yj(PPICLF_JX))
         if (xdist2 .gt. dist2) return

         ydist2 = abs(yi(PPICLF_JY)-yj(PPICLF_JY))
         if (ydist2 .gt. dist2) return

         if (ppiclf_ndim .eq. 3) then
           zdist2 = abs(yi(PPICLF_JZ)-yj(PPICLF_JZ))
           if (zdist2 .gt. dist2) return
         endif

         !
         ! The mean is calcuated according to Lattanzi etal,
         !   Physical Review Fluids, 2022.
         !
         if (j.ne.0) then
         if (qs_fluct_filter_flag==0) then
            upmean   = upmean + yj(PPICLF_JVX)
            vpmean   = vpmean + yj(PPICLF_JVY)
            wpmean   = wpmean + yj(PPICLF_JVZ)
            u2pmean  = u2pmean + yj(PPICLF_JVX)**2
            v2pmean  = v2pmean + yj(PPICLF_JVY)**2
            w2pmean  = w2pmean + yj(PPICLF_JVZ)**2
            icpmean  = icpmean + 1
         else if (qs_fluct_filter_flag==1) then
            ! See https://dpzwick.github.io/ppiclF-doc/algorithms/overlap_mesh.html
            dist = sqrt(xdist2**2 + ydist2**2 + zdist2**2)
            gkern = sqrt(pi*dist2**2/
     >              (4.0d0*log(2.0d0)))**(-ppiclf_ndim) * 
     >              exp(-dist**2/(dist2**2/(4.0d0*log(2.0d0))))

            phipmean = phipmean + gkern*rpropj(PPICLF_R_JVOLP)
            upmean   = upmean +
     >                 gkern*yj(PPICLF_JVX)*rpropj(PPICLF_R_JVOLP)
            vpmean   = vpmean +
     >                 gkern*yj(PPICLF_JVY)*rpropj(PPICLF_R_JVOLP)
            wpmean   = wpmean +
     >                 gkern*yj(PPICLF_JVZ)*rpropj(PPICLF_R_JVOLP)
            u2pmean  = u2pmean +
     >                gkern*(yj(PPICLF_JVX)**2)*rpropj(PPICLF_R_JVOLP)
            v2pmean  = v2pmean +
     >                gkern*(yj(PPICLF_JVY)**2)*rpropj(PPICLF_R_JVOLP)
            w2pmean  = w2pmean +
     >                gkern*(yj(PPICLF_JVZ)**2)*rpropj(PPICLF_R_JVOLP)
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
         rthresh  = (0.5d0+rextra)*rpropi(PPICLF_R_JDP)
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = rxdiff**2 + rydiff**2 + rzdiff**2
         rdiff = sqrt(rdiff)
         
         if (rdiff .gt. rthresh) return

         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)

         ! Compute spring stiffness constant dynamically, 
         !   which overrides the user defined value
         ! Need to make sure this formula is valid for a wall
         ! k1 = k_{n,limit}
         ksp1 = rm1*rpi*rpi/((ksp*ppiclf_dt)**2)
         ! k2 = k_{hertzian}
         E1  = 1.0d9  ! Assumed value for Young's modulus
         nu1 = 0.35d0 ! Assumed value for Poisson's ratio
         Estar = E1/(1.0d0-nu1*nu1)
         r1 = 0.5d0*rpropi(PPICLF_R_JDP)
         r2 = r1
         Rstar = r1*r2/(r1+r2)
         ksp2 = (2.0d0/3.0d0)*Estar*sqrt(Rstar)
         ksp2 = ksp2*sqrt(abs(rdiff-rthresh))
         ! kn = min(k1,k2)
         rksp_wall = min(ksp1,ksp2)
         
         rmult = sqrt(rm1)
         eta_n = 2.0d0*sqrt(rksp_wall)*log(erest)
     >           /sqrt(log(erest)**2+pi2)*rmult
         
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
        
         rdelta12 = rthresh - rdiff
         
         rv12_mag = -yi(PPICLF_JVX)*rn_12x
     >              -yi(PPICLF_JVY)*rn_12y
     >              -yi(PPICLF_JVZ)*rn_12z

         rv12_mage = rv12_mag*eta_n
         rksp_max  = rksp_wall*rdelta12
         rnmag     = -rksp_max - rv12_mage

         
!         ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_ydotc(PPICLF_JVX,i)
!     >                              + rnmag*rn_12x
!         ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_ydotc(PPICLF_JVY,i)
!     >                              + rnmag*rn_12y
!         ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_ydotc(PPICLF_JVZ,i)
!     >                              + rnmag*rn_12z
!


         ! Particles leavind the domain with wall collisions
         ! Simple fix for a conical geometry
         yp = yi(PPICLF_JY)
         zp = yi(PPICLF_JZ)
         vp = yi(PPICLF_JVY)
         wp = yi(PPICLF_JVZ)

         rbound = sqrt(yj(PPICLF_JY)**2 + yj(PPICLF_JZ)**2)
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

           ppiclf_y(PPICLF_JY,i) = yp_new
           ppiclf_y(PPICLF_JZ,i) = zp_new

           ppiclf_y(PPICLF_JVY,i) = vp_new
           ppiclf_y(PPICLF_JVZ,i) = wp_new
         endif

         !write(*,*) "Wall NEAR",i,ppiclf_ydotc(PPICLF_JVY,i)  
      endif


      return
      end
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine to set user-defined values at time t=0
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_InitZero
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
!
! Code:
!
      ppiclf_TimeBH = 0.0d0

      ppiclf_drudtMixt = 0.0d0
      ppiclf_drudtPlag = 0.0d0


      return
      end

!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for output if ppiclf_debug=1
! fort.72## is reserved for debug
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_debug
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Internal:
!
      integer*4 i, n, ic, k, iStage

! Needed for allreduce
      integer*4 ngop
      parameter(ngop = 32)
      real*8 xin(ngop),wout(ngop)

! Needed for viscous unsteady
      integer*4 iT,ii
      real*8 time, factor, A, B, fH, kernelVU
      real*8 FVUoutput
      real*8 ppiclf_npart_sum
      integer*4 npart_tot, npart_max, npart_min
      integer*4 ppiclf_iglsum,ppiclf_iglmax,ppiclf_iglmin
      external  ppiclf_iglsum,ppiclf_iglmax,ppiclf_iglmin

      iStage = 1 ! for internal use only

!
! Code:
!
      ! Use ppiclf ALLREDUCE to compute values across processors
      ! Note that ALLREDUCE uses MPI_BARRIER, which is cpu expensive
      ! Print out every 10th iStage=1 counts

      !xin(1) = dfloat(ppiclf_npart)
      !call ppiclf_gop(xin, wout, '+  ', 1)
      !ppiclf_npart_sum = wout(1)
      npart_tot = ppiclf_iglsum(PPICLF_NPART,1)
      npart_max = ppiclf_iglmax(PPICLF_NPART,1)
      npart_min = ppiclf_iglmin(PPICLF_NPART,1)


      xin=(/phimax,
     >         fqsx_max,fqsy_max,fqsz_max,
     >         famx_max,famy_max,famz_max, 
     >         fdpdx_max,fdpdy_max,fdpdz_max, 
     >         fcx_max,fcy_max,fcz_max,
     >         umean_max,vmean_max,wmean_max,
     >         fqs_mag,fam_mag,fdp_mag,fc_mag,
     >         fqsx_fluct_max,fqsy_fluct_max,fqsz_fluct_max,
     >         fqsx_total_max,fqsy_total_max,fqsz_total_max,
     >         fvux_max,fvuy_max,fvuz_max,
     >         qq_max,tau_max,lift_max/)
      call ppiclf_gop(xin, wout, 'M  ', ngop)
      phimax     = wout(1)
      fqsx_max   = wout(2)
      fqsy_max   = wout(3)
      fqsz_max   = wout(4)
      famx_max   = wout(5)
      famy_max   = wout(6)
      famz_max   = wout(7)
      fdpdx_max  = wout(8)
      fdpdy_max  = wout(9)
      fdpdz_max  = wout(10)
      fcx_max    = wout(11)
      fcy_max    = wout(12)
      fcz_max    = wout(13)
      umean_max  = wout(14)
      vmean_max  = wout(15)
      wmean_max  = wout(16)
      fqs_mag    = wout(17)
      fam_mag    = wout(18)
      fdp_mag    = wout(19)
      fc_mag     = wout(20)
      fqsx_fluct_max = wout(21)
      fqsy_fluct_max = wout(22)
      fqsz_fluct_max = wout(23)
      fqsx_total_max = wout(24)
      fqsy_total_max = wout(25)
      fqsz_total_max = wout(26)
      fvux_max   = wout(27)
      fvuy_max   = wout(28)
      fvuz_max   = wout(29)
      qq_max     = wout(30)
      tau_max    = wout(31)
      lift_max   = wout(32)

      ! Sam - logging for debugging purposes
      ! TLJ - below is a mess I created, need to clean up
      if (ppiclf_nid.eq.0) then

         if (ViscousUnsteady_flag>=1) then
            fH     = 0.75d0 + .105d0*reyL
            factor = 3.0d0*rpi*rnu*dp*fac
            FVUoutput = 0.0
            if (ppiclf_nTimeBH>1) then
               do iT = 2,ppiclf_nTimeBH-1
                  time = ppiclf_timeBH(iT)
                  A  = (4.0d0*rpi*time*rnu/dp**2)**(.25d0)
                  B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >              (0.5d0*dp*rnu*fH**3))**(.5d0)
                  kernelVU = factor*(A+B)**(-2)
                  FVUoutput = FVUoutput + kernelVU*
     >               (ppiclf_drudtMixt(1,iT,1)-ppiclf_drudtPlag(1,iT,1))
                  if (abs(FVUoutput) < 1.d-20) FVUoutput = 0.0d0
               enddo
               iT = ppiclf_nTimeBH
               time = ppiclf_timeBH(iT)
               A  = (4.0d0*rpi*time*rnu/dp**2)**(.25d0)
               B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >              (0.5d0*dp*rnu*fH**3))**(.5d0)
               kernelVU = 0.5*factor*(A+B)**(-2)
               FVUoutput = FVUoutput + kernelVU*
     >             (ppiclf_drudtMixt(1,iT,1)-ppiclf_drudtPlag(1,iT,1))
            endif

            WRITE(7225,"(700(1x,E14.6))") ppiclf_time, 
     >        ((ppiclf_drudtMixt(1,i,1)-ppiclf_drudtPlag(1,i,1))
     >        ,i=1,6)
            WRITE(7227,"(2i5,2x,700(2x,E14.6))") iT,iStage,ppiclf_time,
     >        time,FVUoutput,A,B,kernelVU
            WRITE(7229,"(i5,2x,i5,2x,800(1x,E14.6))") ppiclf_nTimeBH,
     >          ppiclf_nUnsteadyData,ppiclf_dt,
     >          ppiclf_time,ppiclf_TimeBH(1:6)

         endif

         WRITE(7226,"(700(1x,E14.6))") ppiclf_time,
     >        ((ppiclf_drudtMixt(1,i,1)-ppiclf_drudtPlag(1,i,1))
     >        ,i=1,ppiclf_nUnsteadyData)
         WRITE(7228,"(70(1x,E14.6))") ppiclf_time,
     >        ((ppiclf_drudtMixt(3,i,1)-ppiclf_drudtPlag(3,i,1))
     >        ,i=1,ppiclf_nUnsteadyData)

         WRITE(7230,"(27(1x,E23.16))") ppiclf_time, ppiclf_y(1:12, 1)

         WRITE(7231,"(28(1x,E23.16))") ppiclf_time,phimax,
     >             fqsx_max, fqsy_max, fqsz_max,
     >             famx_max, famy_max, famz_max,
     >             fdpdx_max, fdpdy_max, fdpdz_max,
     >             fcx_max, fcy_max, fcz_max,
     >             qq_max,tau_max,lift_max,
     >             fqsx_total_max,fqsy_total_max,fqsz_total_max,
     >             fvux_max, fvuy_max, fvuz_max
         WRITE(7232,"(26(1x,F13.8))") ppiclf_time,
     >             umean_max,vmean_max,wmean_max
         WRITE(7233,"(i5,2x,28(1x,E23.16))")
     >             ppiclf_nid,ppiclf_dt,ppiclf_time,
     >             fac, phimax,
     >             fqsx_fluct_max, fqsy_fluct_max, fqsz_fluct_max
         WRITE(7234,*) ppiclf_nid,istage,PPICLF_LRS ,PPICLF_LPART,
     >             PPICLF_NPART,ppiclf_time,
     >             ppiclf_rprop(PPICLF_R_FLUCTFX:PPICLF_R_FLUCTFZ,1),
     >             ppiclf_ydotc(PPICLF_JVX:PPICLF_JT,1)
         WRITE(7235,"(26(1x,F13.8))") ppiclf_time,
     >             fqs_mag,fam_mag,fdp_mag,fc_mag
         WRITE(7236,"(26(1x,F13.8))") ppiclf_time,
     >             fcx_max, fcy_max, fcz_max
         WRITE(7237,"(26(1x,F13.8))") ppiclf_time,UnifRnd
         WRITE(7240,"(26(1x,F13.8))") ppiclf_time,
     >             fqsx_max, fqsy_max, fqsz_max,
     >             fqsx_fluct_max, fqsy_fluct_max, fqsz_fluct_max,
     >             fqsx_total_max,fqsy_total_max,fqsz_total_max
         WRITE(7241,"(5(1x,E23.16))") ppiclf_time,
     >             phipmean, upmean, vpmean, wpmean
         WRITE(7243,"(1x,E23.16,5(2x,I8))") ppiclf_time,
     >             npart_tot,npart_max,npart_min

         do i = 1,4
            write(7250+i,*) ppiclf_time,
     >              ppiclf_rprop(PPICLF_R_JSDRX:PPICLF_R_JSDRZ,i), ! Du/Dt
     >              ppiclf_rprop(PPICLF_R_JSDOX:PPICLF_R_JSDOZ,i)  ! DOmega/Dt
         enddo


      endif


      return
      end
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for burn rate models for reactive particles
!
! if heattransfer_flag = 0  ignore heat transfer
!                      = 1  Stokes
!                      = 2  Ranz-Marshall (1952)
!                      = 3  Gunn (1977)
!                      = 4  Fox (1978)
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_BR_driver(i,iStage,burnrate_model,
     >   qq,mdot_me,mdot_ox)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage, burnrate_model
      real*8 qq
      real*8 mdot_me, mdot_ox

!
! Code:
!
      if (burnrate_model == 1) then
         call AL_CombModel(i,iStage,qq,mdot_me,mdot_ox)
      elseif (burnrate_model == 2) then
         !call Carbon_CombModel(i,iStage,qq,mdot_me,mdot_ox)
      elseif (burnrate_model == 3) then
         !call Mg_CombModel(i,iStage,qq,mdot_me,mdot_ox)
      else
         call ppiclf_exittr('Unknown combustion model$', 0.0d0, 0)
      endif


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
!-----------------------------------------------------------------------
!
      subroutine AL_CombModel(i,iStage,qq,mdot_me,mdot_ox)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i,iStage
      real*8 qq
      real*8 mdot_me, mdot_ox

      REAL*8 emi, A, E_a, hr, Ks, r_gas
      REAL*8 Q_HSR, Q_comb, T_part, T_ign, T_evap
      REAL*8 Dp_me, m_me, h_comb, h_evap, rho_me

      !
      real*8 rho_ox, m_ox
      real*8 MW_me, MW_o, MW_me2o3
      real*8 oxide_t, V_me, V_ox
      real*8 vol_avg, rho_p, Dia, D0
      real*8 Cp_me, Cp_ox, t
      real*8 rmass_therm

      REAL*8 T_flame
      REAL*8 C, Cs
      REAL*8 xi_co2, xi_h2, xi_h2o, xi_o2, xi_eff
      REAL*8 psi_me, V_p
      REAL*8 Pres
      INTEGER*4 IVALUE, OxideFilm


!===============================================================
!     PARTICLE PROPERTIES
!===============================================================

      T_part = PPICLF_y(PPICLF_JT,i)
      Pres   = PPICLF_RPROP(PPICLF_R_JP,i)

      m_me = PPICLF_y(PPICLF_JMETAL,i)
      m_ox = PPICLF_y(PPICLF_JOXIDE,i)

      Dia    = PPICLF_RPROP(PPICLF_R_JDP,i)
      V_p    = PPICLF_RPROP(PPICLF_R_JVOLP,i)

      D0 = PPICLF_RPROP(PPICLF_R_JIDP,i)

      V_me = m_me / rho_me
      psi_me = V_me / V_p
      Dp_me = ( 6.0d0 * m_me / rpi / rho_me )**(1.0d0/3.0d0)


!-----Al2O3 ignition as a function of diameter, Sundaram et al.--

      if (D0 .le. 100.0d-6) then
         T_ign = 368.0d0 * (D0*1.0d6)**0.268d0 + 780.0d0
      elseif (D0 .le. 1000.0d-6 .and. D0 .gt. 100.0d-6) then
         T_ign = 0.1617d0 * (D0*1.0d6) + 2040.0d0
      else
         T_ign = 3000.0d0
      endif

!-----Updating density--------------------------------------------

      !metal oxide
      if (T_part .le. 3250.0d0 .and. T_part .gt. 2345.0d0) then
         rho_ox = 5632.0d0 - 1.127d0*T_part
      elseif (T_part .le. 2345.0d0) then
         rho_ox = 3970.0d0 * (1.0d0 - 8.0d-6 * (T_part - 300.0d0))
      else
         rho_ox = 5632.0d0 - 1.127d0*3250.0d0
      endif

      !metal
      if (T_part .le. 2743.0d0 .and. T_part .gt. 933.0d0) then
         rho_me = 2385.0d0 - 0.280d0 * (T_part - 933.0d0)
      elseif (T_part .le. 933.0d0) then
         rho_me = 2700.0d0 * (1.0d0 - 23.1d-6 * (T_part - 300.0d0))
      else
         rho_me = 2385.0d0 - 0.280d0 * (2743.0d0 - 933.0d0)
      endif

      !----Updating particle properties--------------------------------

      V_me = m_me/rho_me
      V_ox = m_ox/rho_ox
      vol_avg  = V_me + V_ox

!      phi_me = m_me / (m_me + m_ox)
!      phi_ox = 1.0d0 - phi_me
!      psi_me = ( m_me/rho_me ) / vol_avg
!      psi_ox = 1.0d0 - psi_me

      rho_p = ( m_me + m_ox ) / vol_avg
      Dia = ( 6.0d0 * vol_avg / rpi )**( 1.0d0 / 3.0d0 )

      !-----Molar Weight------------------------------------------------

      MW_me = 0.0269815384d0       ! kg/mol
      MW_o = 0.015999d0            ! kg/mol
      MW_me2o3 = MW_me * 2.0d0 + MW_o * 3.0d0 ! kg/mol

!----------------------------------------------------------------
!     Heat capacity for metal and oxide
!----------------------------------------------------------------

      t = T_part / 1000.0d0

      !Metal heat capacity
      if (T_part .le. 933.0d0) then
         !solid phase heat capacity (Shomate Equation)
         Cp_me = 28.08920d0 - 5.414849d0 * t + 8.560423d0 * t**2.0
     >      + 3.427370d0 * t**3.0d0 - 0.277375d0 / t**2.0

         Cp_me = Cp_me / MW_me ! J/mol*K to J/kg*K
      else
         !liquid phase heat capacity (Shomate Equation)
         Cp_me = 31.75104d0 + 3.935826d-8 * t - 1.786515d-8 * t**2.0
     >       + 2.694171d-9 * t**3.0d0 + 5.480037d-9 / t**2.0

         Cp_me = Cp_me / MW_me ! J/mol*K to J/kg*K
      endif


      !Metal oxide heat capacity
      if (T_part .le. 2327.0d0) then
         !solid phase heat capacity (Shomate Equation)
         Cp_ox = 102.4290d0 + 38.7498d0 * t - 15.91090d0 * t**2.0
     >      + 2.628181d0 * t**3.0d0 - 3.007551d0 / t**2.0

         Cp_ox = Cp_ox / MW_me2o3 ! J/mol*K to J/kg*K
      else
         !liquid phase heat capacity (Shomate Equation)
         Cp_ox = 192.4640d0 + 9.519856d-8 * t - 2.858928d-8 * t**2.0
     >      + 2.929147d-9 * t**3.0d0 + 5.599405d-8 / t**2.0

         Cp_ox = Cp_ox / MW_me2o3 ! J/mol*K to J/kg*K
      endif

      rmass_therm = m_me * Cp_me + m_ox * Cp_ox

!----Updating PPICLF_RPOP values---------------------------------

      ! TEMP FIX
      if (Dia.gt.PPICLF_RPROP(PPICLF_R_JIDP,i)) then
!         print*,'Warning Dia too big',
!     >    i, PPICLF_RPROP(PPICLF_R_JIDP,i),Dia
         Dia = PPICLF_RPROP(PPICLF_R_JIDP,i)
      endif
      PPICLF_RPROP(PPICLF_R_JDP,i)   = Dia
      PPICLF_RPROP(PPICLF_R_JRHOP,i) = rho_p
      PPICLF_RPROP(PPICLF_R_JVOLP,i) = vol_avg

!===============================================================
!     COMBUSTION HEAT TRANSFER
!===============================================================

      !constants - Al
      emi=0.1d0
      A=200.0d0
      E_a=95395.0d0
      hr=3.1d7
      h_evap=1.183d7
      T_evap=2743.d0

      r_gas=8.314d0 ! Gas Constant

!----------------------------------------------------------------
!     Source terms for heat transfer
!----------------------------------------------------------------

      !initialize source terms to 0

      Q_HSR = 0.0d0
      Q_comb = 0.0d0

      !preheat - HSR
      if (T_part .ge. 933.0d0 .and. T_part .lt. T_ign) then
         Ks = A * exp(-E_a / (r_gas * T_part))
         Q_HSR = hr * Ks * Dp_me**2.0d0 * rpi
      endif

      !combustion - comb
      if (T_part .ge. T_ign) then
            if (T_part .le. 3300.d0) then
               h_comb = -4.3334d7
            else
               h_comb = 3.3d6
            endif
         Q_comb = -mdot_me * (h_comb + h_evap)
      endif

      qq = qq + Q_HSR + Q_comb


!===============================================================
!     COMBUSTION MODEL
!===============================================================

!----------------------------------------------------------------
!     Setup information
!----------------------------------------------------------------

      !particle
      IVALUE = 0                   ! 0 for B&W; 1 for Maggi

      !constants
      C = 5.9427624643167829d-7    ! (m**1.8,s**-1,K**-0.2,Pa**-0.1)
      Cs = 0.0d0                   ! kg/m^3

!----------------------------------------------------------------
!         WIDENER and Beckstead, AIAA-98-382
!         With Fadi Najjar's correction, J. Spacecraft & Rockets,
!           43(6), 2006, pp.1258--1270
!----------------------------------------------------------------

!-------------MOLAR FRACTIONS IN THE GAS PHASE-------------------

      if (IVALUE==0) then  ! Widener and Beckstead (1998)
          xi_o2  = 0.02d0
          xi_h2o = 0.2d0
          xi_co2 = 0.2d0
          xi_h2  = 0.2d0
      endif
      if (IVALUE==1) then  ! Filippo Maggi (CSAR, 2007)
          xi_o2  = 0.00002d0
          xi_h2o = 0.04809d0
          xi_co2 = 0.00394d0
          xi_h2  = 0.33582d0
      endif

!-----Relative oxidizer concentration----------------------------

      xi_eff = 0.0291d0
!      xi_eff = xi_o2 + 0.58d0*xi_h2o + 0.22d0*xi_co2



!----------------------------------------------------------------
!             Aluminum mass burning rate (UNITS: SI)
!             Weighted by aluminum volume fraction
!             this formula takes diameter in meters
!----------------------------------------------------------------

!-------------Check if particle is burning-----------------------

              if ((T_part .ge. T_ign .or. PPICLF_RPROP(PPICLF_R_JBRNT,i)
     >           .gt. 0.0) .and. Dp_me .gt. 5.0d-6) then

                 !Burn law
                 mdot_me = C*rho_me*(Dia**1.2d0) *
     >              xi_eff*(T_part**0.2d0) *
     >              (pres**0.1d0) * psi_me

                 mdot_ox = 0.25d0 * rpi * Dia**2.0d0 *
     >              abs(vmag) * 0.25d0 * Cs

                 !update total burn time
                 PPICLF_RPROP(PPICLF_R_JBRNT,i) =
     >              PPICLF_RPROP(PPICLF_R_JBRNT,i) + ppiclf_dt
              else
                 mdot_me = 0.0d0
                 mdot_ox = 0.0d0
              endif

      return
      end
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
      IF (PPICLF_LINIT .OR. PPICLF_OVERLAP)
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
      PPICLF_LCOMM = .TRUE.

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_InitCrystal
!
!     This subroutine is called form ppiclf_comm_InitMPI
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"

!
! Code:
!
      ! GSlib call to setup crystal router for communication across
      ! processors
      CALL pfgslib_crystal_setup(ppiclf_cr_hndl,ppiclf_comm,ppiclf_np)

      RETURN
      END
!-----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_comm_InitOverlapGrid(ncell,fluidGrid)
     > bind(C, name="ppiclc_comm_InitOverlapGrid")
#else
      SUBROUTINE ppiclf_comm_InitOverlapGrid(ncell,fluidGrid)
#endif
!
! This subroutine is called from rocpicl/PICL_TEMP_InitSolver.F90
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4 ncell
      REAL*8    fluidGrid(7,ncell)
      ! Expected size: fluidGrid(7,ncell)
      ! Indicies 1-3: Centroid x,y,z position
      ! Indicies 4-6: Max cell dx,dy,dz based on any vertex combination
      ! Index      7: Cell Volume 
!
! External:
!
      INTEGER*4 ie, i, ierr
!
      ppiclf_overlap = .TRUE.

      IF(.NOT. PPICLF_LCOMM)
     > CALL ppiclf_exittr('InitMPI must be before InitOverlap$',0.0d0,0)
      IF(.NOT. PPICLF_LINIT)
     > CALL ppiclf_exittr('InitParticle must be before InitOverlap$'
     >                  ,0.0d0,0)

      IF(ncell .GT. PPICLF_LEE .OR. ncell .LT. 0) THEN
        PRINT*, '***ERROR*** PPICLF_LEE', PPICLF_LEE, 'in', 
     >   'InitMapOverlapGrid must be greater than', ncell
        CALL MPI_BARRIER(ppiclf_comm,ierr) 
        CALL ppiclf_exittr('Increase LEE in InitOverlap$',0.0d0,ncell)
      END IF

      ! Number of finite volume cells from fluid solver
      ppiclf_nFVCells = ncell

      DO ie=1,ppiclf_nFVCells
        DO i = 1,7
          ppiclf_fluid_grid(i,ie) = fluidGrid(i,ie)
        END DO
      END DO

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_CreateBin
!
! This subroutine is called from ppiclf_solve_InitSolve
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal:
!
      INTEGER*4 ix, iy, iz, npt_total, i, j, k, idum, jdum, kdum, 
     >          total_bin, targetTotBin, idealBin(3), Temp_iBin(3),
     >          iBin(3), iBinTot, tempi, ideal_bin_index(3), largeBin, 
     >          medBin, smallBin,
     >          ppiclf_iglsum, NBMax, ierr, MaxPotentialBins(3), 
     >          maxbincount1, maxbincount2, BinCheck, minbin(3),
     >          ppiclf_iglmax, binNegBound, binPosBound, binIterations
     >          
      REAL*8    xmin, ymin, zmin, xmax, ymax, zmax, temp1, temp2,
     >          binb_length(3), BinMinLen(3), ppiclf_glmin, 
     >          ppiclf_glmax, ppiclf_glsum, periodicDistCheck,
     >          BinBuffer(3), binsReal(3), binError, minBinError,
     >          increaseRatio
      EXTERNAL  ppiclf_iglsum, ppiclf_glmin, ppiclf_glmax, ppiclf_glsum,
     >          ppiclf_iglmax
      LOGICAL   MaxBinsAchieved(3), TwoSmallBins
!

      ix = 1
      iy = 2
      iz = 3
        
      ! iglsum is integer global sum across MPI ranks.
      npt_total = ppiclf_iglsum(ppiclf_npart,1)

      ! Bin must be larger than nearest neighbor search distance
      ! and the ppiclf_filter(1:3).  This makes a buffer around the bin
      ! domain. Increase if you desire to bin less frequently.
      ! ppiclf_filter is 2x cell length in each direction
      DO i = 1,3
        BinMinLen(i) = MAX(ppiclf_filter(i),ppiclf_nndist)
        ! Need ppiclf_filter to make sure you have 1 layer
        ! of outer fluid cells
        ! Need ppiclf_nndist/2 to ensure BinMinLen is never violated
        BinBuffer(i) = MAX(ppiclf_filter(i)/2,ppiclf_nndist/2)
        MaxBinsAchieved(i) = .FALSE.
      END DO

      xmin =  1D10
      ymin =  1D10
      zmin =  1D10
      xmax = -1D10
      ymax = -1D10
      zmax = -1D10

      ! Looping through particles on this processor
      ! to find bin boundary locations
      DO i=1,ppiclf_npart
         ! Finding min/max particle extremes.
         ! Add buffer so that layers of outer cells 
         ! are available for interpolation/projection.
         temp1 = ppiclf_y(ix,i) - BinBuffer(ix)
         temp2 = ppiclf_y(ix,i) + BinBuffer(ix)
         IF(temp1 .LT. xmin) xmin = temp1
         IF(temp2 .GT. xmax) xmax = temp2

         temp1 = ppiclf_y(iy,i) - BinBuffer(iy)
         temp2 = ppiclf_y(iy,i) + BinBuffer(iy)
         IF(temp1 .LT. ymin) ymin = temp1
         IF(temp2 .GT. ymax) ymax = temp2

         temp1 = ppiclf_y(iz,i) - BinBuffer(iz)
         temp2 = ppiclf_y(iz,i) + BinBuffer(iz)
         IF(temp1 .LT. zmin) zmin = temp1
         IF(temp2 .GT. zmax) zmax = temp2
      END DO

      ! Finds global bin domain boundaries across MPI ranks
      ppiclf_binb(1) = ppiclf_glmin(xmin,1)
      ppiclf_binb(2) = ppiclf_glmax(xmax,1)
      ppiclf_binb(3) = ppiclf_glmin(ymin,1)
      ppiclf_binb(4) = ppiclf_glmax(ymax,1)
      ppiclf_binb(5) = ppiclf_glmin(zmin,1)
      ppiclf_binb(6) = ppiclf_glmax(zmax,1)

      CALL MPI_BARRIER(ppiclf_comm,ierr)

      ! If all particles within last RK Stage binbound, do not calculate
      ! bins again and do not remap overlap grid.
      BinCheck = 0
      DO i = 1,3
        IF((ppiclf_binb(2*i-1) + BinBuffer(i)) .LT.
     >                             ppiclf_previousbinb(2*i-1)) THEN
          BinCheck = 1
          EXIT
        END IF
        IF((ppiclf_binb(2*i)   - BinBuffer(i)) .GT.
     >                             ppiclf_previousbinb(2*i)) THEN
          BinCheck = 1
          EXIT
        END IF
      END DO

      CALL MPI_BARRIER(ppiclf_comm,ierr)

      BinCheck = ppiclf_iglmax(BinCheck,1)

#ifdef TEST
      ! So that CreateBin can be called repeatedly
      ! for different number of processors 
      BinCheck = 1
#endif

      IF(BinCheck .EQ. 0) THEN
        DO i = 1,3
          ppiclf_binb(2*i-1) = ppiclf_previousbinb(2*i-1)
          ppiclf_binb(2*i)   = ppiclf_previousbinb(2*i)
        END DO
        ppiclf_binchanged = .FALSE.
        RETURN
      ELSE
        ppiclf_binchanged  = .TRUE.
        ppiclf_printbinvtu = .TRUE.
      END IF

      ! Ensuring ppiclf_binb not greater than 
      ! cartesian fluid domain extremes.
      ! If dist within ppiclf_nndist, set ppiclf_binb
      ! equal to fluid domain for periodic ghost particles.
      ! Needed to know when to use linear periodic
      
      ppiclf_EqualDomain(1) = .FALSE.
      ppiclf_EqualDomain(2) = .FALSE.
      ppiclf_EqualDomain(3) = .FALSE.

      DO i = 1,3
        ! Check bin min domain
        periodicDistCheck = MAX(ppiclf_nndist,ppiclf_filter(i))
        IF(ppiclf_binb(i*2-1) - periodicDistCheck .LE. 
     >                          ppiclf_xdrange(1,i)) THEN
          ppiclf_binb(i*2-1) = ppiclf_xdrange(1,i)
          ppiclf_EqualDomain(i) = .TRUE.
        END IF
        ! Check bin max domain
        IF(ppiclf_binb(i*2)+periodicDistCheck .GE. 
     >                          ppiclf_xdrange(2,i)) THEN
          ppiclf_binb(i*2) = ppiclf_xdrange(2,i)
        ELSE
          ppiclf_EqualDomain(i) = .FALSE.
        END IF
      END DO

      ! Set previous bin bound for next RK Stage check
      DO i = 1,6
        ppiclf_previousbinb(i) = ppiclf_binb(i)
      END DO

      ! Return from subroutine if no particles present      
      IF(npt_total .LT. 1) RETURN

      ! Find ppiclf bin domain lengths
      DO i = 1,3
        binb_length(i) = ppiclf_binb(2*i) -
     >                         ppiclf_binb(2*i-1)
      END DO

!*** Start active bin iteration loop here      
      ! Update with targetTotBin based on active/inactive
      targetTotBin = ppiclf_np
      ! This sets the number of bin permutations to try to find optimal
      ! balance
      binPosBound = 1
      binNegBound = 1
      binIterations = binPosBound + binNegBound 
      DO i = 1,3
        MaxPotentialBins(i) = FLOOR(binb_length(i)/BinMinLen(i))
        IF(MaxPotentialBins(i) .LT. 1) THEN
          CALL ppiclf_exittr('BinMinLen() criteria violated.',0.0D0,0)
        END IF
      END DO

      ! Number of bins calculated based on bin surface
      ! area minimization and bin aspect ratio close to 1
      binsReal(1) = (targetTotBin**(1.0D0/3.0D0))*
     >              (binb_length(1)**(2.0D0/3.0D0))/ 
     >              ((binb_length(2)**(1.0D0/3.0D0))*
     >              (binb_length(3))**(1.0D0/3.0D0)) 
      
      binsReal(2) = (targetTotBin**(1.0D0/3.0D0))*
     >              (binb_length(2)**(2.0D0/3.0D0))/ 
     >              ((binb_length(1)**(1.0D0/3.0D0))*
     >              (binb_length(3))**(1.0D0/3.0D0)) 
     
      binsReal(3) = (targetTotBin**(1.0D0/3.0D0))*
     >              (binb_length(3)**(2.0D0/3.0D0))/ 
     >              ((binb_length(2)**(1.0D0/3.0D0))*
     >              (binb_length(1))**(1.0D0/3.0D0)) 

      ! The loop below ensures bins don't exceed number of
      ! processors or minimum bin length criteria, while
      ! maximizing the other dimensions.
      maxbincount2 = 0
      tempi = 0
      DO
        tempi = tempi + 1
        IF(tempi .GT. 100) THEN
          PRINT*, 'CreateBin stuck in infinite loop'
          CALL ppiclf_exittr('Error in Createbins',0.0,0)
        END IF
        ! Ensures most bins don't exceed number of 
        ! processors in case where one direction is
        ! much longer than the others.
        ! Ranking the bins Real numbers
        largeBin = 1 
        smallBin = 1
        medBin   = 1
        DO i = 1,3
          ! INT(x+0.5) is equivalent to ROUND(x). 
          ! Round isn't a built in fortran function.
          IF(binsReal(i) .LT. 0.50D0) binsReal(i) = 0.51D0
          ppiclf_n_bins(i) = INT(binsReal(i)+0.5D0)
          minBin(i) = ppiclf_n_bins(i) - binNegBound
          IF(minBin(i) .LT. 1) minBin(i) = 1
          IF(binsReal(i) .GT. binsReal(largeBin)) largeBin = i
          IF(binsReal(i) .LT. binsReal(smallBin)) smallBin = i
        END DO
        DO i = 1,3
          IF(i .EQ. largeBin) CYCLE
          IF(i .EQ. smallBin) CYCLE
          medBin = i
        END DO
        ! This checks to make sure at least one combination is within
        ! criteria of the largest possible total bins
        IF(minBin(1)*minBin(2)*minBin(3) .GT. targetTotBin) THEN
          IF(binsReal(medBin) .LT. 1) THEN
            ! Make largest dimension equal to number of processors.
            ! Other two dimensions are small
            DO i = 1,3
              IF(i .EQ. largeBin) THEN
                ppiclf_n_bins(largeBin) = targetTotBin 
                binsReal(largeBin) = REAL(ppiclf_n_bins(largeBin))
              ELSE
                ppiclf_n_bins(i) = 1
                binsReal(i) = 0.51D0
              END IF
            END DO
          ELSE
            ! Set bins of two large dimensions proportionally by length
            ! Set small dimension with 1 bins.
            binsReal(smallBin) = 0.51D0
            ppiclf_n_bins(smallBin) = 1

            binsReal(medBin) = SQRT(REAL(targetTotBin)) *
     >                    binb_length(medBin)/binb_length(largeBin)
            ppiclf_n_bins(medBin) = INT(binsReal(medBin)+0.5D0)

            binsReal(largeBin) = SQRT(REAL(targetTotBin)) *
     >                     binb_length(largeBin)/binb_length(medBin)
            ppiclf_n_bins(largeBin) = INT(binsReal(largeBin)+0.5D0)
          END IF 
        END IF
        ! Ensures bins don't violate minimum length criteria
        maxbincount1  = 0
        increaseRatio = 1.0D0
        DO i = 1,3
          IF(ppiclf_n_bins(i) .GE. MaxPotentialBins(i)) THEN
            ppiclf_n_bins(i) = MaxPotentialBins(i)
            ! This will be used to increase non-maximized bins
            increaseRatio  = (binsReal(i)/REAL(MaxPotentialBins(i)))* 
     >                        increaseRatio
            binsReal(i) = REAL(MaxPotentialBins(i))
            MaxBinsAchieved(i) = .TRUE.
            maxbincount1 = maxbincount1 + 1
          END IF
        END DO
        ! Either they both equal 0 and exit on the first pass,
        ! or they iterate and are equal on second or third pass.
        IF(maxbincount1 .EQ. maxbincount2) EXIT

        ! Increases other bins if one or two dimensions 
        ! are at maximum size based on minimum length criteria
        IF(MaxBinsAchieved(1) .OR. MaxBinsAchieved(2) .OR.
     >                             MaxBinsAchieved(3)) THEN
          maxbincount2 = 0
          DO i = 1,3
            IF(MaxBinsAchieved(i)) THEN
              maxbincount2 = maxbincount2 + 1
            END IF
          END DO
          DO i = 1,3
            IF(MaxBinsAchieved(i)) CYCLE
            IF(maxbincount2 .EQ. 1) THEN
              binsReal(i) = SQRT(increaseRatio)*binsReal(i)
            ELSE
              binsReal(i) = binsReal(i)*increaseRatio
            END IF
          END DO
        END IF
      END DO

      ! Since bin must be an integer, check -1:+1 number of bins
      ! for each bin dimension. Ideal number of bins will be max value
      ! while less than number of total target of bins.
      ! Minimize total error of sum of real bin calc - bin integer
      ! This will ensure larger dimensions get more bins.
      total_bin = 0
      minBinError = 1.0D9
      ideal_bin_index = 0
      DO ix = 0,binIterations
        iBin(1) = ppiclf_n_bins(1) + (ix-binNegBound)
        ppiclf_bins_dx(1) = binb_length(1)/iBin(1)
        IF(ppiclf_bins_dx(1) .LT. BinMinLen(1) .OR.
     >                           iBin(1) .LT. 1) CYCLE
        DO iy = 0,binIterations
          iBin(2) = ppiclf_n_bins(2) + (iy-binNegBound)
          ppiclf_bins_dx(2) = binb_length(2)/iBin(2)
          IF(ppiclf_bins_dx(2) .LT. BinMinLen(2) .OR.
     >                             iBin(2) .LT. 1) CYCLE
          DO iz = 0,binIterations
            iBin(3) = ppiclf_n_bins(3) + (iz-binNegBound)
            ppiclf_bins_dx(3) = binb_length(3)/iBin(3)
            IF(ppiclf_bins_dx(3) .LT. BinMinLen(3) .OR.
     >                               iBin(3) .LT. 1) CYCLE
            iBinTot = iBin(1)*iBin(2)*iBin(3)
            IF(iBinTot .GE. total_bin .AND. iBinTot .GT. 0 .AND.
     >                     iBinTot .LE. targetTotBin) THEN
              ! This resets the error min when a larger amount of
              ! total bins are found.
              IF(iBinTot .GT. total_bin) THEN
                minBinError = 1.0D9
              END IF
              ! If this is equal maximum # bins, find the best
              ! combination of bins per dimension
              binError = 0.0D0
              DO i = 1,3
                binError = binError + ABS(binsReal(i) - REAL(iBin(i))) 
              END DO
              IF(binError .LT. minBinError) THEN
                  minBinError = binError
                  ideal_bin_index(1) = ix
                  ideal_bin_index(2) = iy
                  ideal_bin_index(3) = iz
                  total_bin = iBinTot
              END IF            
            END IF
          END DO !iz
        END DO !iy
      END DO !ix
      IF(total_bin .EQ. 0 .AND. ppiclf_nid .EQ. 0) THEN
        PRINT*, 'correct bin combination not found'
      END IF

      tempi = 0
      total_bin = 1      
      DO i = 1,3
        ! Set common ppiclf bin arrays based on above calculation
        ! same number of bin equation as in above loops with best
        ! indices
        ppiclf_n_bins(i) = INT(binsReal(i)+0.5D0)
     >                     + (ideal_bin_index(i) - binNegBound)
        ppiclf_bins_dx(i) = binb_length(i)/ppiclf_n_bins(i)
        total_bin = total_bin*ppiclf_n_bins(i)
        IF(total_bin .GT. ppiclf_np) THEN
          IF(ppiclf_nid .EQ. 0) THEN
            PRINT*, 'binERROR: Num Bins > NumProcessors',total_bin,
     >            ppiclf_np, targetTotBin
          END IF
          CALL ppiclf_exittr('Error in Createbins',0.0,0)
        END IF
        IF(ppiclf_n_bins(i) .LT. 1) THEN
          PRINT*, 'binERROR: this dimension has negative bins', i,
     >            'Bins per dimension:', ppiclf_bins_dx(1),
     >            ppiclf_bins_dx(2), ppiclf_bins_dx(3)
        END IF
        IF(ppiclf_n_bins(i) .GT. tempi) THEN
          NBMax = i ! Dimension with max number of bins
          tempi = ppiclf_n_bins(i)
        END IF
      END DO
     
      ! Loop to see if we can add one to dimension with largest number of bins
      ! Choose this dimension because it is smallest incremental increase to total bins 
      DO
        IF((total_bin/ppiclf_n_bins(NBMax))*
     >      (ppiclf_n_bins(NBMax)+1) .LT. targetTotBin) THEN
          ! Add a bin and set new bin dx length
          ppiclf_n_bins(NBMax) = ppiclf_n_bins(NBMax)+1
          ppiclf_bins_dx(NBMax) = binb_length(NBMax)/
     >                              ppiclf_n_bins(NBMax)
          IF(ppiclf_bins_dx(NBMax) .LT. BinMinLen(NBMax)) THEN
            ! If BinMinLen criteria violated, return to previous bin configuration
            ppiclf_n_bins(NBMax) = ppiclf_n_bins(NBMax)-1
            ppiclf_bins_dx(NBMax) = binb_length(NBMax)/
     >                                ppiclf_n_bins(NBMax)
            EXIT
          END IF
          total_bin = 1
          DO i = 1,3
            total_bin = total_bin*ppiclf_n_bins(i)
          END DO
        ELSE
          EXIT
        END IF
      END DO

! CALL ActiveBinCounter
! End active bin check loop
! CALL BinToProcessorMap
      ! *** AVERY
      ! Find this processor's x,y,z bin indicies
      idum = modulo(ppiclf_nid,ppiclf_n_bins(1))
      jdum = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
      kdum = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))

      ! Calculate this processor's bin min/max position in each dimension.
      ! Note that this value stays with this MPI rank.
      ppiclf_bin_pos(1,1) = ppiclf_binb(1) + idum    *ppiclf_bins_dx(1)
      ppiclf_bin_pos(2,1) = ppiclf_binb(1) + (idum+1)*ppiclf_bins_dx(1)
      ppiclf_bin_pos(1,2) = ppiclf_binb(3) + jdum    *ppiclf_bins_dx(2)
      ppiclf_bin_pos(2,2) = ppiclf_binb(3) + (jdum+1)*ppiclf_bins_dx(2)
      ppiclf_bin_pos(1,3) = ppiclf_binb(5) + kdum    *ppiclf_bins_dx(3)
      ppiclf_bin_pos(2,3) = ppiclf_binb(5) + (kdum+1)*ppiclf_bins_dx(3)

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_FindParticle
!
! This subroutine is called from ppiclf_solve_InitSolve
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal:
!
      INTEGER*4  i, ii, jj, kk, nrank, ierr, partcheck
      EXTERNAL   ppiclf_iglmax
      INTEGER*4  ppiclf_iglmax
!
      partcheck = 0
      DO i=1,ppiclf_npart
         ! Calculates particle's bin index
         ii  = FLOOR((ppiclf_y(1,i)-ppiclf_binb(1))/ppiclf_bins_dx(1))
         jj  = FLOOR((ppiclf_y(2,i)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
         kk  = FLOOR((ppiclf_y(3,i)-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
        
         ! Calculates particle's bin
         nrank  = ii + ppiclf_n_bins(1)*jj + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
         IF(nrank .NE. ppiclf_iprop(4,i)) partcheck = 1

         ! Maps particle to correct processor based on active bin number
         ! ***Use BinToProcMap for active/inactive bin***
         ppiclf_iprop(4,i) = nrank ! Processor to send to
         ppiclf_iprop(5,i) = ii    ! x bin #
         ppiclf_iprop(6,i) = jj    ! y bin #
         ppiclf_iprop(7,i) = kk    ! z bin #
         ppiclf_iprop(8,i) = nrank ! total bin number
      END DO
      ppiclf_particleMoved = ppiclf_iglmax(partcheck,1)
      CALL mpi_barrier(ppiclf_comm,ierr)

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_MoveParticle
!
! This subroutine is called from ppiclf_solve_InitSolve
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal:
!
      LOGICAL   partl ! dummy variable    
      INTEGER*4 rtempLim
      PARAMETER(rtempLim = PPICLF_LRS*4 + PPICLF_LRP + PPICLF_LRP2
     >       + PPICLF_LRP3 + PPICLF_LRP4 + PPICLF_LRP5 + PPICLF_LRP_PRO)
      REAL*8    rtemp(rtempLim,PPICLF_LPART)
      INTEGER*4 i, icount, j0
!
      ! copy particle y, rprop, rprop2, rprop3 arrays into rtemp
      ! array for communication
      DO i=1,ppiclf_npart
        icount = 1
        CALL ppiclf_copy(rtemp(icount,i),ppiclf_y(1,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(rtemp(icount,i),ppiclf_y1(1,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(rtemp(icount,i),ppiclf_ydot(1,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(rtemp(icount,i),ppiclf_ydotc(1,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(rtemp(icount,i),ppiclf_rprop(1,i),PPICLF_LRP)
        icount = icount + PPICLF_LRP
        IF(PPICLF_LRP2 .GT. 1) THEN
          CALL ppiclf_copy(rtemp(icount,i),
     >                     ppiclf_rprop2(1,i),PPICLF_LRP2)
          icount = icount + PPICLF_LRP2
        END IF
        IF(PPICLF_LRP3 .GT. 1) THEN
          CALL ppiclf_copy(rtemp(icount,i),
     >                     ppiclf_rprop3(1,i),PPICLF_LRP3)
          icount = icount + PPICLF_LRP3
        END IF
        IF(PPICLF_LRP4 .GT. 1) THEN
          CALL ppiclf_copy(rtemp(icount,i),
     >                     ppiclf_rprop4(1,i),PPICLF_LRP4)
          icount = icount + PPICLF_LRP4
        END IF
        IF(PPICLF_LRP5 .GT. 1) THEN
          CALL ppiclf_copy(rtemp(icount,i),
     >                     ppiclf_rprop5(1,i),PPICLF_LRP5)
          icount = icount + PPICLF_LRP5
        END IF
        CALL ppiclf_copy(rtemp(icount,i),
     >                   ppiclf_feedbk(1,i),PPICLF_LRP_PRO)
      END DO
      
      j0 = 4 ! index of ppiclf_iprop that contains rank to send to
      CALL pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >             ,ppiclf_npart,PPICLF_LPART ! Setup
     >             ,ppiclf_iprop,PPICLF_LIP   ! Integer Comm
     >             ,partl,0                   ! Logical Comm
     >             ,rtemp,rtempLim            ! Real Comm
     >             ,j0)                       ! Receiver processor index

      IF(ppiclf_npart .GT. PPICLF_LPART .OR. ppiclf_npart .LT. 0) THEN
        PRINT*,'Increase LPART. Processor:',ppiclf_nid,
     >   'LPART should be greater than:',ppiclf_npart
        CALL ppiclf_exittr('Increase LPART$',0.0d0,ppiclf_npart)
      END IF
 
      ! Update processor particle values with newly transfered rtemp
      ! array from communication
      DO i=1,ppiclf_npart
        icount = 1
        CALL ppiclf_copy(ppiclf_y(1,i),rtemp(icount,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(ppiclf_y1(1,i),rtemp(icount,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(ppiclf_ydot(1,i),rtemp(icount,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(ppiclf_ydotc(1,i),rtemp(icount,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(ppiclf_rprop(1,i),rtemp(icount,i),PPICLF_LRP)
        icount = icount + PPICLF_LRP
        IF(PPICLF_LRP2 .GT. 1) THEN
        CALL ppiclf_copy(ppiclf_rprop2(1,i),rtemp(icount,i),
     >                   PPICLF_LRP2)
        icount = icount + PPICLF_LRP2
        END IF
        IF(PPICLF_LRP3 .GT. 1) THEN
          CALL ppiclf_copy(ppiclf_rprop3(1,i),rtemp(icount,i),
     >                     PPICLF_LRP3)
          icount = icount + PPICLF_LRP3
        END IF
        IF(PPICLF_LRP4 .GT. 1) THEN
          CALL ppiclf_copy(ppiclf_rprop4(1,i),rtemp(icount,i),
     >                     PPICLF_LRP4)
          icount = icount + PPICLF_LRP4
        END IF
        IF(PPICLF_LRP5 .GT. 1) THEN
          CALL ppiclf_copy(ppiclf_rprop5(1,i),rtemp(icount,i),
     >                     PPICLF_LRP5)
          icount = icount + PPICLF_LRP5
        END IF
        CALL ppiclf_copy(ppiclf_feedbk(1,i),rtemp(icount,i),
     >           PPICLF_LRP_PRO)
      END DO
        
      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_MapOverlapGrid
!
! This subroutine is called from ppiclf_solve_InitSolve
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
      INTEGER*4 nkey(2), i, j, k, l, ie, iee, ii, jj, kk, nrank,
     >          nl, nii, njj, nrr, iic, jjc, kkc il, ierr
      INTEGER*4 ix, iy, iz, ixLow, ixHigh, iyLow,
     >          iyHigh, izLow, izHigh 
      REAL*8    rxval, ryval, rzval, EleSizei(3), MaxPoint(3),
     >          MinPoint(3), ppiclf_vlmin, ppiclf_vlmax,
     >          centeri(3), exchCellMultiplier
      LOGICAL   partl, ErrorFound
      EXTERNAL  ppiclf_vlmin, ppiclf_vlmax
!
! Code Start:
!
      ! Number of fluid finite volume cells that map to particle bins
      ppiclf_nCells_FV2PICL = 0 
      
      ! Multiplies by x.4999 the cell length to ensure that x layers of cells
      ! outside of the ppiclf bin are mapped for interpolation.
      ! This could be changed based on the desired frequency of
      ! ppiclf bin creation and mapping.
      exchCellMultiplier  = 2.499D0

      ! Loops through number of fluid FV cells on this processor
      DO ie=1,ppiclf_nFVCells  
        DO l=1,3
         !indicies 1:3
         centeri(l) = ppiclf_fluid_grid(l,ie)
         !indicies 4:6
         EleSizei(l) =  exchCellMultiplier*ppiclf_fluid_grid(3+l,ie)
         IF(EleSizei(l) .GT. ppiclf_bins_dx(l)/2) THEN
           EleSizei(l) =  1.499D0*ppiclf_fluid_grid(3+l,ie)
         END IF
        END DO !l 

        ! Fluid Cell vertex position without additional length
        rxval = centeri(1)
        ryval = centeri(2)
        rzval = centeri(3)
      
        ! Exits if fluid cell center is outside of any bin boundaries 
        IF (rxval .GT. ppiclf_binb(2)) CYCLE
        IF (rxval .LT. ppiclf_binb(1)) CYCLE
        IF (ryval .GT. ppiclf_binb(4)) CYCLE
        IF (ryval .LT. ppiclf_binb(3)) CYCLE
        IF (rzval .GT. ppiclf_binb(6)) CYCLE
        IF (rzval .LT. ppiclf_binb(5)) CYCLE
 
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
        
        IF(FLOOR((rxval + EleSizei(1) - ppiclf_binb(1))
     >       /ppiclf_bins_dx(1)) .NE. ii)  ixHigh = 3

        IF(FLOOR((rxval - EleSizei(1) - ppiclf_binb(1))
     >       /ppiclf_bins_dx(1)) .NE. ii)  ixLow = 1
        IF(FLOOR((ryval + EleSizei(2) - ppiclf_binb(3))
     >       /ppiclf_bins_dx(2)) .NE. jj)  iyHigh = 3

        IF(FLOOR((ryval - EleSizei(2) - ppiclf_binb(3))
     >       /ppiclf_bins_dx(2)) .NE. jj)  iyLow = 1

        IF(ppiclf_ndim .GT. 2 .AND. FLOOR((rzval + EleSizei(3)
     >    - ppiclf_binb(5))/ppiclf_bins_dx(3)) .NE. kk)  izHigh = 3

        IF(ppiclf_ndim .GT. 2 .AND. FLOOR((rzval - EleSizei(3)
     >    - ppiclf_binb(5))/ppiclf_bins_dx(3)) .NE. kk)  izLow = 1

        DO ix=ixLow,ixHigh
          DO iy=iyLow,iyHigh
            DO iz=izLow,izHigh
              ! Change cell position by EleSizei if ix,iy,or iz NE 2
              rxval = centeri(1) + (ix-2)*EleSizei(1)
              ryval = centeri(2) + (iy-2)*EleSizei(2)
              rzval = centeri(3) + (iz-2)*EleSizei(3)
              ! Find bin for adjusted rval
              ii    = FLOOR((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
              jj    = FLOOR((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
              kk    = FLOOR((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
              

              ! This covers ghost exchanged cells for linear periodicity
              ! Maps cells greater than ppiclf bin domain to first bin
              ! Maps cells less than ppiclf bin domain to last bin
              IF(ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)
     >                               .AND. ppiclf_n_bins(1) .GT. 1) THEN
                IF(ii .EQ. ppiclf_n_bins(1)) ii = 0
                IF(ii .EQ. -1) ii = ppiclf_n_bins(1) - 1
              END IF
              IF(ppiclf_linperiodic(2) .AND. ppiclf_EqualDomain(2)
     >                               .AND. ppiclf_n_bins(2) .GT. 1) THEN
                IF(jj .EQ. ppiclf_n_bins(2)) jj = 0
                IF(jj .EQ. -1) jj = ppiclf_n_bins(2) - 1
              END IF
              IF(ppiclf_linperiodic(3) .AND. ppiclf_EqualDomain(3)
     >                               .AND. ppiclf_n_bins(3) .GT. 1) THEN
                IF(kk .EQ. ppiclf_n_bins(3)) kk = 0
                IF(kk .EQ. -1) kk = ppiclf_n_bins(3) - 1
              END IF
              
              ! Ensures duplicate cells don't get sent to same processor
              IF (ii .LT. 0 .OR. ii .GT. ppiclf_n_bins(1)-1) CYCLE
              IF (jj .LT. 0 .OR. jj .GT. ppiclf_n_bins(2)-1) CYCLE
              IF (kk .LT. 0 .OR. kk .GT. ppiclf_n_bins(3)-1) CYCLE


              ! Calculates processor rank
              nrank  = ii + ppiclf_n_bins(1)*jj + 
     >                     ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
              ppiclf_nCells_FV2PICL = ppiclf_nCells_FV2PICL + 1
              IF(ppiclf_nCells_FV2PICL .GT. PPICLF_LEE) THEN
                PRINT*, '***ERROR*** PPICLF_LEE',PPICLF_LEE, 'in', 
     >           'MapOverlapGrid must be greater than',
     >            ppiclf_nCells_FV2PICL 
                CALL ppiclf_exittr('Increase PPICLF_LEE$ (MapOverlap)',0.0D0
     >               ,ppiclf_nCells_FV2PICL)
              END IF

              ! make sure it is mapped to active nrank and map rank to
              ! processor. *** FOR ACTIVE BINNING ***

              ! Stores cells to rank mapping.
              ! Fluid solver cell ID
              ppiclf_cell_map(1,ppiclf_nCells_FV2PICL) = ie
              ! Fluid solver cell rank
              ppiclf_cell_map(2,ppiclf_nCells_FV2PICL) = ppiclf_nid
              ! Particle solver cell rank and bin indicies
              ppiclf_cell_map(3,ppiclf_nCells_FV2PICL) = nrank
              ppiclf_cell_map(4,ppiclf_nCells_FV2PICL) = ii
              ppiclf_cell_map(5,ppiclf_nCells_FV2PICL) = jj
              ppiclf_cell_map(6,ppiclf_nCells_FV2PICL) = kk
            END DO !iz
          END DO !iy
        END DO !ix
      END DO !ie

      DO ie=1,ppiclf_nCells_FV2PICL 
        ! These copy all indicies since Fortran is column-major
        iee = ppiclf_cell_map(1,ie)
        CALL ppiclf_copy(ppiclf_picl_grid(1,ie)
     >                 ,ppiclf_fluid_grid(1,iee),7)
 
        ! ppiclf_filter initially set in PICL_TEMP_InitSolver
        ! Want to only consider cells that reside in the particle domain
        ! Update ppiclf_filter for next binning cycle 2.1*dx since next
        ! layer of cells in a growing particle domain may be slightly larger
        ! and we want filter equal a minimum of 2 cells.
        DO l = 1,3
          IF(ppiclf_filter(l) .LT. ppiclf_fluid_grid(3+l,iee)) THEN
            ppiclf_filter(l) = 2.1D0*ppiclf_fluid_grid(3+l,iee)
          END IF
        END DO
      END DO

      ! Copy mapping since it is need to send fluid properties in interp
      ppiclf_nCells_FV2PICL_Orig = ppiclf_nCells_FV2PICL
      DO ie=1,ppiclf_nCells_FV2PICL_Orig
         ! Copies cells to rank mapping (integer copy)
         CALL ppiclf_icopy(ppiclf_cell_map_Orig(1,ie)
     >            ,ppiclf_cell_map(1,ie),PPICLF_LRMAX)
      END DO

      ! GSLIB required info
      ! NumPiclCells - number of columns to transfer
      ! PPICLF_LEE - number of columns declared
      ! nl - partl row size (dummy logical variable)
      nl   = 0
      ! nii - ppiclf_cell_map row size declared
      nii  = PPICLF_LRMAX
      ! njj - Row index of ppiclf_cell_map with receiver processor/rank
      njj  = 3
      ! nrr - ppiclf_rocGrid row size declared
      nrr  = 7
      ! Defines sorting order
      nkey(1) = 2
      nkey(2) = 1

      CALL pfgslib_crystal_tuple_transfer(
     >        ppiclf_cr_hndl,ppiclf_nCells_FV2PICL,PPICLF_LEE !setup
     >        ,ppiclf_cell_map,nii ! Integer Comm
     >        ,partl,nl                 ! Logical Comm
     >        ,ppiclf_picl_grid,nrr      ! Real Comm
     >        ,njj)                      ! Receiver processor index
      CALL pfgslib_crystal_tuple_sort(
     >        ppiclf_cr_hndl,ppiclf_nCells_FV2PICL !setup
     >        ,ppiclf_cell_map,nii !Integer to sort
     >        ,partl,nl                 !Logical to sort
     >        ,ppiclf_picl_grid,nrr      !Real to sort
     >        ,nkey,2)                  !sorting method


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
      SUBROUTINE ppiclf_comm_CreateGhost
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal:
!
      REAL*8     GhostPos(3), PeriodicShift(3), 
     >           distSQ(3), distCheckSQ, buffer
      INTEGER*4  ip, idum, iip, jjp, kkp, iig, jjg, kkg, nrank, 
     >           j, k, l, GhostInc(3), ix, iy, iz

      ! Calculate the linear periodicity shift in each dimension
      DO l = 1,3
        IF(ppiclf_linperiodic(l)) THEN
          PeriodicShift(l) = ppiclf_xdrange(2,l) - ppiclf_xdrange(1,l)
        ELSE
          PeriodicShift(l) = 0.0D0
        END IF
      END DO

      ppiclf_npart_gp = 0
      
      DO ip=1,ppiclf_npart
        idum = 0
        ! Copy particle solution variables
        DO j=1,PPICLF_LRS
           idum = idum + 1
           ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip)
        END DO
        ! Copy particle property variables
        DO j=1,PPICLF_LRP
           idum = idum + 1
           ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip)
        END DO

        ! GP Bin Index
        iip    = ppiclf_iprop(5,ip)
        jjp    = ppiclf_iprop(6,ip)
        kkp    = ppiclf_iprop(7,ip)
   
        ! Found that buffer was needed in unit testing 
        ! due to round-off errors with periodicity
        buffer = 1.02 
        distCheckSQ = (ppiclf_nndist*(buffer**3))**2

        DO ix = 1,3
          distSQ = 0.0D0
          GhostPos(1) = ppiclf_cp_map(1,ip)
          IF(ix .LT. 3) THEN
            CALL ppiclf_comm_GhostDistCheck(ix,GhostPos(1),
     >                     ppiclf_nndist*buffer,GhostInc(1),1,distSQ(1))
            IF(GhostInc(1) .EQ. 0) CYCLE
          ELSE
            GhostInc(1) = 0 !For ghosts in other 2 dimensions only
          END IF
          iig = iip + GhostInc(1)

          ! Angular Periodicity Check
          ! *** Add here ***

          ! If ghost is outside of ppiclf domain:
          IF(iig .LT. 0 .OR. iig .GT. ppiclf_n_bins(1)-1) THEN
            IF(ppiclf_linperiodic(1).AND. ppiclf_EqualDomain(1)) THEN
              CALL ppiclf_comm_LinearPeriodicityGhost
     >                       (iig,1,GhostPos(1),PeriodicShift(1))
            ELSE
              CYCLE
            END IF
          END IF

          DO iy = 1,3
            GhostPos(2) = ppiclf_cp_map(2,ip)
            IF(iy .LT. 3) THEN
              CALL ppiclf_comm_GhostDistCheck(iy,GhostPos(2),
     >                    ppiclf_nndist*buffer,GhostInc(2),2,distSQ(2))
              IF(GhostInc(2) .EQ. 0.) CYCLE
              ! This corner/edge check caused issues when unit testing
              ! Removing this makes extra ghost particles
              !IF(distSQ(1)+distSQ(2) .GT. distCheckSQ) CYCLE !corner/edge check
            ELSE
              GhostInc(2) = 0 !For ghosts in other 2 dimensions only
            END IF
            jjg = jjp + GhostInc(2)

          ! Angular Periodicity Check
          ! *** Add here ***

            ! If ghost is outside of ppiclf domain:
            IF(jjg .LT. 0 .OR. jjg .GT. ppiclf_n_bins(2)-1) THEN
              IF(ppiclf_linperiodic(2) .AND.
     >                     ppiclf_EqualDomain(2)) THEN
                CALL ppiclf_comm_LinearPeriodicityGhost(jjg,2,
     >                           GhostPos(2),PeriodicShift(2))
              ELSE
                CYCLE
              END IF
            END IF

            DO iz = 1,3
              GhostPos(3) = ppiclf_cp_map(3,ip)
              IF(iz .LT. 3) THEN
                CALL ppiclf_comm_GhostDistCheck(iz,GhostPos(3),
     >                    ppiclf_nndist*buffer,GhostInc(3),3,distSQ(3))
                IF(GhostInc(3) .EQ. 0) CYCLE
                ! This corner/edge check caused issues when unit testing
                ! Removing this makes extra ghost particles
                !corner/edge check
                !IF(distSQ(1)+distSQ(2)+distSQ(3) .GT. distCheckSQ) CYCLE
              ELSE
                GhostInc(3) = 0
              END IF
              kkg = kkp + GhostInc(3)

          ! Angular Periodicity Check
          ! *** Add here ***              

              ! If ghost is outside of ppiclf domain:
              IF(kkg .LT. 0 .OR. kkg .GT. ppiclf_n_bins(3)-1) THEN
                IF(ppiclf_linperiodic(3) .AND.
     >                        ppiclf_EqualDomain(3)) THEN
                  CALL ppiclf_comm_LinearPeriodicityGhost
     >                           (kkg,3,GhostPos(3),PeriodicShift(3))
                ELSE
                  CYCLE
                END IF
              END IF

              ! This prevents ghost in same bin.
              IF(GhostInc(1) .EQ. 0 .AND. GhostInc(2) .EQ. 0 .AND.
     >           GhostInc(3) .EQ. 0 ) CYCLE
              
              ! Add ghost particle and map integer and real properties
              nrank = iig + ppiclf_n_bins(1)*jjg 
     >               + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
              !ghostsMade = ghostsMade + 1
              ppiclf_npart_gp = ppiclf_npart_gp + 1
              ! Copy particle ID info
              ppiclf_iprop_gp(1,ppiclf_npart_gp) = ppiclf_iprop(1,ip)
              ppiclf_iprop_gp(2,ppiclf_npart_gp) = ppiclf_iprop(2,ip)
              ppiclf_iprop_gp(3,ppiclf_npart_gp) = ppiclf_iprop(3,ip)
              ppiclf_iprop_gp(4,ppiclf_npart_gp) = nrank !*** change to processor
              ppiclf_iprop_gp(5,ppiclf_npart_gp) = iig
              ppiclf_iprop_gp(6,ppiclf_npart_gp) = jjg
              ppiclf_iprop_gp(7,ppiclf_npart_gp) = kkg
              ppiclf_iprop_gp(8,ppiclf_npart_gp) = nrank

              ppiclf_rprop_gp(1,ppiclf_npart_gp) = GhostPos(1)
              ppiclf_rprop_gp(2,ppiclf_npart_gp) = GhostPos(2)
              ppiclf_rprop_gp(3,ppiclf_npart_gp) = GhostPos(3)

              DO k=4,PPICLF_LRP_GP
                ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
              END DO
            END DO !iz = 1:3
          END DO !iy = 1:3
        END DO !ix = 1:3
      END DO !ip = 1:ppiclf_npart

      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_GhostDistCheck(ix,Pos,distchk,
     >                                      GhostInc,l,dSQ)
      
      IMPLICIT NONE
      
      INCLUDE "PPICLF"

      ! ix: ghostcheck loop counter, GhostInc: bin +/-, l: dimenison
      INTEGER*4 ix, GhostInc, l
      ! Pos: Position of Ghost Particle, distchk: criteria to create
      ! ghost particle
      ! distSQ: used to evaluate distance ghost in edge & corner case
      REAL*8    Pos, distchk, dSQ

      ! ppiclf_bin_pos(1,1) is bin min position in x
      ! ppiclf_bin_pos(2,1) is bin max position in x
      IF(ABS(Pos - ppiclf_bin_pos(ix,l)) 
     >                          .LT. distchk) THEN
        dSQ = (Pos-ppiclf_bin_pos(ix,l))**2
        IF(ix .EQ. 1) GhostInc = -1 ! close to bin min
        IF(ix .EQ. 2) GhostInc =  1 ! clost to bin max
      ELSE
        GhostInc = 0
      END IF
      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_LinearPeriodicityGhost(iig,l,Pos,PerShift)
      
      IMPLICIT NONE
      
      INCLUDE "PPICLF"

      !iig: Ghost bin index, l: dimension Number (1:x,2:y,3:z)
      INTEGER*4 iig, l 
      ! Pos: GhostPos(l), PerShift: PeriodicShift(l)
      REAL*8    Pos, PerShift
      IF(iig .LT. 0) THEN
        iig = ppiclf_n_bins(l) - 1
        Pos = Pos + PerShift
      ELSE IF (iig .GT. ppiclf_n_bins(l) - 1) THEN
        iig = 0
        Pos = Pos - PerShift
      END IF

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_MoveGhost
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Internal:
!
      INTEGER*4 iprop_proc_index
      LOGICAL   partl  ! Dummy variable       
!
      iprop_proc_index = 4 ! since ppiclf_iprop(4,np) contains processor
                           ! that should receive ghost particle
      CALL pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >             ,ppiclf_npart_gp,PPICLF_LPART_GP ! Setup
     >             ,ppiclf_iprop_gp,PPICLF_LIP_GP   ! Integer Comm
     >             ,partl,0                         ! Logical Comm
     >             ,ppiclf_rprop_gp,PPICLF_LRP_GP   ! Real Comm
     >             ,iprop_proc_index)               ! Receiver processor index

      RETURN
      END

!----------------------------------------------------------------------
!      subroutine ppiclf_comm_AngularCreateGhost
!!
!      implicit none
!!
!      include "PPICLF"
!!
!! Internal:
!!
!      real*8 xdlen,ydlen,zdlen,rxdrng(3),rxnew(3), rfac, rxval, ryval,
!     >       rzval, rxl, ryl, rzl, rxr, ryr, rzr, distchk, dist, gFilt
!      integer*4 iadd(3),gpsave(27)
!      real*8 map(PPICLF_LRP_PRO)
!      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
!     >           nfacegp, nedgegp, ncornergp, iperiodicx, iperiodicy,
!     >           iperiodicz, jx, jy, jz, ip, idum, iip, jjp, kkp, ii1,
!     >           jj1, kk1, iig, jjg, kkg, iflgx, iflgy, iflgz,
!     >           isave, iflgsum, ndumn, nrank, ibctype, i, ifc, ist, j,
!     >           k
!      ! 08/27/24 - Thierry - added for angular periodicty starts here
!      real*8 alpha
!      integer*4 xrank, yrank, zrank
!      ! 08/27/24 - Thierry - added for angular periodicty ends here
!      ! 09/26/24 - Thierry - added for angular periodicty starts here
!      real*8 dist1, dist2
!      ! 09/26/24 - Thierry - added for angular periodicty ends here
!!
!
!c     face, edge, and corner number, x,y,z are all inline, so stride=3
!      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
!      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
!     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
!     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
!      el_corner_num = (/ -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
!     >                   -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)
!
!      nfacegp   = 4  ! number of faces
!      nedgegp   = 4  ! number of edges
!      ncornergp = 0  ! number of corners
!
!      if (ppiclf_ndim .gt. 2) then
!         nfacegp   = 6  ! number of faces
!         nedgegp   = 12 ! number of edges
!         ncornergp = 8  ! number of corners
!      endif
!
!      iperiodicx = ppiclf_iperiodic(1)
!      iperiodicy = ppiclf_iperiodic(2)
!      iperiodicz = ppiclf_iperiodic(3)
!
!! ------------------------
!c CREATING GHOST PARTICLES
!! ------------------------
!      jx    = 1
!      jy    = 2
!      jz    = 3
!
!      ! Thierry - we dont use xdlen and ydlen in this algorithm. no need to modify them.
!      xdlen = ppiclf_binb(2) - ppiclf_binb(1)
!      ydlen = ppiclf_binb(4) - ppiclf_binb(3)
!      zdlen = -1.
!      if (ppiclf_ndim .gt. 2) 
!!     >   zdlen = ppiclf_binb(6) - ppiclf_binb(5)
!     >   zdlen = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
!      if (iperiodicx .ne. 0) xdlen = -1
!      if (iperiodicy .ne. 0) ydlen = -1
!      if (iperiodicz .ne. 0) zdlen = -1
!
!      rxdrng(1) = xdlen
!      rxdrng(2) = ydlen
!      rxdrng(3) = zdlen
!
!      ppiclf_npart_gp = 0
!
!      rfac = 1.0d0
!      gFilt = MAX(ppiclf_nndist,ppiclf_filter(1),
!     >            ppiclf_filter(2),ppiclf_filter(3))
!
!
!      do ip=1,ppiclf_npart
!
!         call ppiclf_user_MapProjPart(map,ppiclf_y(1,ip)
!     >         ,ppiclf_ydot(1,ip),ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip))
!
!c        idum = 1
!c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
!c        idum = 2
!c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
!c        idum = 3
!c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
!
!         idum = 0
!         do j=1,PPICLF_LRS
!            idum = idum + 1
!            ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip) ! ppiclf_y(PPICLF_JX/ JY/ JZ/ JVX/ JVY/ JVZ/ JT, ip)
!         enddo
!         idum = PPICLF_LRS
!         do j=1,PPICLF_LRP
!            idum = idum + 1
!            ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip) ! ppiclf_rprop(PPICLF_R_JRHOP/ R_JRHOF/ .../ R_WDOTZ, ip)
!         enddo
!         idum = PPICLF_LRS+PPICLF_LRP
!         do j=1,PPICLF_LRP_PRO
!            idum = idum + 1
!            ppiclf_cp_map(idum,ip) = map(j) ! map(PPICLF_P_JPHIP/ JFX/ .../ JPHIPW) - these are found in ppiclf_user_MapProjPart
!         enddo
!
!         rxval = ppiclf_cp_map(1,ip) ! ppiclf_y(PPICLF_JX,ip)
!         ryval = ppiclf_cp_map(2,ip) ! ppiclf_y(PPICLF_JY,ip)
!         rzval = 0.0d0
!         if (ppiclf_ndim .gt. 2) rzval = ppiclf_cp_map(3,ip) ! ppiclf_y(PPICLF_JZ,ip)
!
!         iip    = ppiclf_iprop(4,ip) ! ith coordinate of bin
!         jjp    = ppiclf_iprop(5,ip) ! jth coordinate of bin
!         kkp    = ppiclf_iprop(6,ip) ! kth coordinate of bin
!
!         rxl = ppiclf_binb(1) + ppiclf_bins_dx(1)*iip ! min x of bin
!         rxr = rxl + ppiclf_bins_dx(1)                ! max x of bin
!         ryl = ppiclf_binb(3) + ppiclf_bins_dx(2)*jjp
!         ryr = ryl + ppiclf_bins_dx(2)
!         rzl = 0.0d0
!         rzr = 0.0d0
!         if (ppiclf_ndim .gt. 2) then
!            rzl = ppiclf_binb(5) + ppiclf_bins_dx(3)*kkp
!            rzr = rzl + ppiclf_bins_dx(3)
!         endif
!
!         isave = 0
!
!         ! faces
!         do ifc=1,nfacegp
!            ist = (ifc-1)*3
!            ii1 = iip + el_face_num(ist+1) 
!            jj1 = jjp + el_face_num(ist+2)
!            kk1 = kkp + el_face_num(ist+3)
!
!            iig = ii1
!            jjg = jj1
!            kkg = kk1
!
!            distchk = 0.0d0
!            dist = 0.0d0
!            if (ii1-iip .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
!               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
!            endif
!            if (jj1-jjp .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
!               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
!            endif
!            if (ppiclf_ndim .gt. 2) then
!            if (kk1-kkp .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
!               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
!            endif
!            endif
!            distchk = sqrt(distchk)
!            dist = sqrt(dist)
!
!            if (ang_case==1) then  ! for wedge geometry
!
!               ! Thierry - I dont think it's code efficient to call this subroutine
!               !           for every particle, every ghost face, at every time step
!               !           I'm wondering if it's better if we make the plane values 
!               !           as global values that are initialized in the beginning 
!            
!               call ppiclf_solve_InitAngularPlane(ip,
!     >                                 ang_per_rin  , ang_per_rout  ,
!     >                                 ang_per_angle, ang_per_xangle,
!     >                                 dist1, dist2)
!               if ((dist .gt. distchk).and.(dist1.gt.distchk)
!     >           .and.(dist2.gt.distchk)) cycle
!            else
!               if (dist .gt. distchk) cycle
!            endif
!
!            iflgx = 0
!            iflgy = 0
!            iflgz = 0
!!-----------------------------------------------------------------------
!            ! 08/27/24 - Thierry - modification for angular periodicty starts here
!
!               ! angle between particle and x-axis
!                alpha = atan2(ppiclf_y(PPICLF_JY,ip), 
!     >                        ppiclf_y(PPICLF_JX,ip))
!                
!
!                call ppiclf_solve_InvokeAngularPeriodic(ip, 
!     >                                                  ang_per_flag,
!     >                                                  alpha,         
!     >                                                  ang_per_angle,  
!     >                                                  ang_per_xangle, 
!     >                                                  0)
!
!              ! Thierry - this is how FindParticle implements it
!              ! need to find a way to make the code deal with negative xrot values
!
!            xrank = iig ; yrank=jjg; zrank = kkg
!            ! Thierry - previously placed before the CheckPeriodicBC call, had to move them for the periodic check
!            iadd(1) = ii1
!            iadd(2) = jj1
!            iadd(3) = kk1
!            rxnew(1) = rxval
!            rxnew(2) = ryval
!            rxnew(3) = rzval ! z-coordinate does not change when angular periodicity is invoked
!            
!            ! Angular periodicity check in x- and y-directions
!            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
!              iflgx = 1
!              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
!              if (iperiodicx .ne. 0) cycle
!              iig = xrank
!              jjg = yrank
!            end if
!            
!            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
!              iflgy = 1
!              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
!              if (iperiodicy .ne. 0) cycle
!              iig = xrank
!              jjg = yrank
!            end if
!            
!            ! Linear periodicity check in z-direction
!            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
!              iflgz = 1
!              kkg =modulo(kkg,ppiclf_n_bins(3))
!              if (iperiodicz .ne. 0) cycle
!              ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
!              ! rxdrng(3) = -1.0  if not periodic in Z
!              if (rxdrng(3) .gt. 0) then 
!                if (iadd(3) .ge. ppiclf_n_bins(3)) then ! particle leaving from max z-face
!                  rxnew(3) = rxnew(3) - rxdrng(3)
!                elseif (iadd(3) .lt. 0) then ! particle leaving from min z-face
!                  rxnew(3) = rxnew(3) + rxdrng(3)
!                end if ! iadd
!              end if ! rxrdrng
!            else ! z-linear periodicity not applicable
!              kkg = zrank
!            end if ! kkg
!            
!            iflgsum = iflgx + iflgy + iflgz
!            ndumn  = iig + ppiclf_n_bins(1)*jjg + 
!     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
!             nrank = ndumn
!
!            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle
!
!            ! 08/27/24 - Thierry - modification for angular periodicty ends here
!!-----------------------------------------------------------------------
!
!            do i=1,isave
!               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 111
!            enddo
!            isave = isave + 1
!            gpsave(isave) = nrank
!
!            ibctype = iflgx+iflgy+iflgz
!            
!            rxnew(1) = xrot(1)
!            rxnew(2) = xrot(2)
!            ppiclf_cp_map(4,ip) = vrot(1)
!            ppiclf_cp_map(5,ip) = vrot(2)
!                 
!            ppiclf_npart_gp = ppiclf_npart_gp + 1
!            ppiclf_iprop_gp(1,ppiclf_npart_gp) = ppiclf_iprop(1,ip)
!            ppiclf_iprop_gp(2,ppiclf_npart_gp) = ppiclf_iprop(2,ip)
!            ppiclf_iprop_gp(3,ppiclf_npart_gp) = nrank
!            ppiclf_iprop_gp(4,ppiclf_npart_gp) = iig
!            ppiclf_iprop_gp(5,ppiclf_npart_gp) = jjg
!            ppiclf_iprop_gp(6,ppiclf_npart_gp) = kkg
!            ppiclf_iprop_gp(7,ppiclf_npart_gp) = nrank
!
!            ! Thierry - we don't need ppiclf_comm_CheckPeriodicBC anymore for the angular periodic ghost algorithm
!            !           as this is now taken care of when anticipating where the particle might be when calling
!            !           ppiclf_comm_InvokeAngularPeriodic
!            !           we only need to assign xr and vr to ppiclf_rprop_gp
!
!            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1) ! ppiclf_y(PPICLF_JX, ip) for the periodic ghost particle
!            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2) ! JY
!            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3) ! JZ
!            
!            do k=4,PPICLF_LRP_GP
!               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
!            enddo
!  111 continue
!         enddo
!
!         ! edges
!         do ifc=1,nedgegp
!            ist = (ifc-1)*3
!            ii1 = iip + el_edge_num(ist+1) 
!            jj1 = jjp + el_edge_num(ist+2)
!            kk1 = kkp + el_edge_num(ist+3)
!
!            iig = ii1
!            jjg = jj1
!            kkg = kk1
!
!            distchk = 0.0d0
!            dist = 0.0d0
!            if (ii1-iip .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
!               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
!            endif
!            if (jj1-jjp .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
!               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
!            endif
!            if (ppiclf_ndim .gt. 2) then
!            if (kk1-kkp .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
!               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
!            endif
!            endif
!            distchk = sqrt(distchk)
!            dist = sqrt(dist)
!
!            if (ang_case==1) then  ! for wedge geometry
!
!               call ppiclf_solve_InitAngularPlane(ip,
!     >                                 ang_per_rin  , ang_per_rout  ,
!     >                                 ang_per_angle, ang_per_xangle,
!     >                                 dist1, dist2)
!               if ((dist .gt. distchk).and.(dist1.gt.distchk)
!     >           .and.(dist2.gt.distchk)) cycle
!            else
!               if (dist .gt. distchk) cycle
!            endif
!
!            iflgx = 0
!            iflgy = 0
!            iflgz = 0
!            ! periodic if out of domain - add some ifsss
!!-----------------------------------------------------------------------
!            ! 08/27/24 - Thierry - modification for angular periodicty starts here
!
!               ! angle between particle and x-axis
!                alpha = atan2(ppiclf_y(PPICLF_JY,ip), 
!     >                        ppiclf_y(PPICLF_JX,ip))
!                
!
!                call ppiclf_solve_InvokeAngularPeriodic(ip, 
!     >                                                  ang_per_flag,
!     >                                                  alpha,         
!     >                                                  ang_per_angle,  
!     >                                                  ang_per_xangle, 
!     >                                                  0)
!
!              ! Thierry - this is how FindParticle implements it
!              ! need to find a way to make the code deal with negative xrot values
!
!            xrank = iig ; yrank=jjg; zrank = kkg
!            ! Thierry - previously placed before the CheckPeriodicBC call, had to move them for the periodic check
!            iadd(1) = ii1
!            iadd(2) = jj1
!            iadd(3) = kk1
!            rxnew(1) = rxval
!            rxnew(2) = ryval
!            rxnew(3) = rzval ! z-coordinate does not change when angular periodicity is invoked
!            
!            ! Angular periodicity check in x- and y-directions
!            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
!              iflgx = 1
!              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
!              if (iperiodicx .ne. 0) cycle
!              iig = xrank
!              jjg = yrank
!            end if
!            
!            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
!              iflgy = 1
!              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
!              if (iperiodicy .ne. 0) cycle
!              iig = xrank
!              jjg = yrank
!            end if
!            
!            ! Linear periodicity check in z-direction
!            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
!              iflgz = 1
!              kkg =modulo(kkg,ppiclf_n_bins(3))
!              if (iperiodicz .ne. 0) cycle
!              ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
!              ! rxdrng(3) = -1.0  if not periodic in Z
!              if (rxdrng(3) .gt. 0) then ! particle leaving from max z-face
!                if (iadd(3) .ge. ppiclf_n_bins(3)) then
!                  rxnew(3) = rxnew(3) - rxdrng(3)
!                elseif (iadd(3) .lt. 0) then
!                  rxnew(3) = rxnew(3) + rxdrng(3)
!                end if ! iadd
!              end if ! rxrdrng
!            else ! z-linear periodicity not applicable
!              kkg = zrank
!            end if ! kkg
!            
!            iflgsum = iflgx + iflgy + iflgz
!            ndumn  = iig + ppiclf_n_bins(1)*jjg + 
!     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
!             nrank = ndumn
!
!            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle
!
!            ! 08/27/24 - Thierry - modification for angular periodicty ends here
!!-----------------------------------------------------------------------
!
!            do i=1,isave
!               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 222
!            enddo
!            isave = isave + 1
!            gpsave(isave) = nrank
!
!            ibctype = iflgx+iflgy+iflgz
!
!            rxnew(1) = xrot(1)
!            rxnew(2) = xrot(2)
!            ppiclf_cp_map(4,ip) = vrot(1)
!            ppiclf_cp_map(5,ip) = vrot(2)
!                 
!            ppiclf_npart_gp = ppiclf_npart_gp + 1
!
!            ppiclf_iprop_gp(1,ppiclf_npart_gp) = ppiclf_iprop(1,ip)
!            ppiclf_iprop_gp(2,ppiclf_npart_gp) = ppiclf_iprop(2,ip)
!            ppiclf_iprop_gp(3,ppiclf_npart_gp) = nrank
!            ppiclf_iprop_gp(4,ppiclf_npart_gp) = iig
!            ppiclf_iprop_gp(5,ppiclf_npart_gp) = jjg
!            ppiclf_iprop_gp(6,ppiclf_npart_gp) = kkg
!            ppiclf_iprop_gp(7,ppiclf_npart_gp) = nrank
!
!            ! Thierry - we don't need ppiclf_comm_CheckPeriodicBC anymore for the angular periodic ghost algorithm
!            !           as this is now taken care of when anticipating where the particle might be when calling
!            !           ppiclf_comm_InvokeAngularPeriodic
!            !           we only need to assign xr and vr to ppiclf_rprop_gp
!
!            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1) ! ppiclf_y(PPICLF_JX, ip) for the periodic ghost particle
!            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2) ! JY
!            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3) ! JZ
!            
!            do k=4,PPICLF_LRP_GP
!               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
!            enddo
!  222 continue
!         enddo
!
!         ! corners
!         do ifc=1,ncornergp
!            ist = (ifc-1)*3
!            ii1 = iip + el_corner_num(ist+1) 
!            jj1 = jjp + el_corner_num(ist+2)
!            kk1 = kkp + el_corner_num(ist+3)
!
!            iig = ii1
!            jjg = jj1
!            kkg = kk1
!
!            distchk = 0.0d0
!            dist = 0.0d0
!            if (ii1-iip .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
!               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
!            endif
!            if (jj1-jjp .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
!               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
!            endif
!            if (ppiclf_ndim .gt. 2) then
!            if (kk1-kkp .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
!               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
!            endif
!            endif
!            distchk = sqrt(distchk)
!            dist = sqrt(dist)
!
!            if (ang_case==1) then  ! for wedge geometry
!            
!               call ppiclf_solve_InitAngularPlane(ip,
!     >                                 ang_per_rin  , ang_per_rout  ,
!     >                                 ang_per_angle, ang_per_xangle,
!     >                                 dist1, dist2)
!               if ((dist .gt. distchk).and.(dist1.gt.distchk)
!     >           .and.(dist2.gt.distchk)) cycle
!            else
!               if (dist .gt. distchk) cycle
!            endif
!
!            iflgx = 0
!            iflgy = 0
!            iflgz = 0
!
!!-----------------------------------------------------------------------
!            ! 08/27/24 - Thierry - modification for angular periodicty starts here
!
!               ! angle between particle and x-axis
!                alpha = atan2(ppiclf_y(PPICLF_JY,ip), 
!     >                        ppiclf_y(PPICLF_JX,ip))
!                
!
!                call ppiclf_solve_InvokeAngularPeriodic(ip, 
!     >                                                  ang_per_flag,
!     >                                                  alpha,         
!     >                                                  ang_per_angle,  
!     >                                                  ang_per_xangle, 
!     >                                                  0)
!
!              ! Thierry - this is how FindParticle implements it
!              ! need to find a way to make the code deal with negative xrot values
!
!            xrank = iig ; yrank=jjg; zrank = kkg
!            ! Thierry - previously placed before the CheckPeriodicBC call, had to move them for the periodic check
!            iadd(1) = ii1
!            iadd(2) = jj1
!            iadd(3) = kk1
!            rxnew(1) = rxval
!            rxnew(2) = ryval
!            rxnew(3) = rzval ! z-coordinate does not change when angular periodicity is invoked
!            
!            ! Angular periodicity check in x- and y-directions
!            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
!              iflgx = 1
!              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
!              if (iperiodicx .ne. 0) cycle
!              iig = xrank
!              jjg = yrank
!            end if
!            
!            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
!              iflgy = 1
!              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
!              if (iperiodicy .ne. 0) cycle
!              iig = xrank
!              jjg = yrank
!            end if
!            
!            ! Linear periodicity check in z-direction
!            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
!              iflgz = 1
!              kkg =modulo(kkg,ppiclf_n_bins(3))
!              if (iperiodicz .ne. 0) cycle
!              ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
!              ! rxdrng(3) = -1.0  if not periodic in Z
!              if (rxdrng(3) .gt. 0) then ! particle leaving from max z-face
!                if (iadd(3) .ge. ppiclf_n_bins(3)) then
!                  rxnew(3) = rxnew(3) - rxdrng(3)
!                elseif (iadd(3) .lt. 0) then
!                  rxnew(3) = rxnew(3) + rxdrng(3)
!                end if ! iadd
!              end if ! rxrdrng
!            else ! z-linear periodicity not applicable
!              kkg = zrank
!            end if ! kkg
!            
!            iflgsum = iflgx + iflgy + iflgz
!            ndumn  = iig + ppiclf_n_bins(1)*jjg + 
!     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
!             nrank = ndumn
!
!            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle
!
!            ! 08/27/24 - Thierry - modification for angular periodicty ends here
!!-----------------------------------------------------------------------
!            do i=1,isave
!               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 333
!            enddo
!            isave = isave + 1
!            gpsave(isave) = nrank
!
!            ibctype = iflgx+iflgy+iflgz
!
!            rxnew(1) = xrot(1)
!            rxnew(2) = xrot(2)
!            ppiclf_cp_map(4,ip) = vrot(1)
!            ppiclf_cp_map(5,ip) = vrot(2)
!
!            ppiclf_npart_gp = ppiclf_npart_gp + 1
!
!            ppiclf_iprop_gp(1,ppiclf_npart_gp) = ppiclf_iprop(1,ip)
!            ppiclf_iprop_gp(2,ppiclf_npart_gp) = ppiclf_iprop(2,ip)
!            ppiclf_iprop_gp(3,ppiclf_npart_gp) = nrank
!            ppiclf_iprop_gp(4,ppiclf_npart_gp) = iig
!            ppiclf_iprop_gp(5,ppiclf_npart_gp) = jjg
!            ppiclf_iprop_gp(6,ppiclf_npart_gp) = kkg
!            ppiclf_iprop_gp(7,ppiclf_npart_gp) = nrank
!
!            ! Thierry - we don't need ppiclf_comm_CheckPeriodicBC anymore for the angular periodic ghost algorithm
!            !           as this is now taken care of when anticipating where the particle might be when calling
!            !           ppiclf_comm_InvokeAngularPeriodic
!            !           we only need to assign xr and vr to ppiclf_rprop_gp
!
!            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1) ! ppiclf_y(PPICLF_JX, ip) for the periodic ghost particle
!            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2) ! JY
!            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3) ! JZ
!
!            do k=4,PPICLF_LRP_GP
!               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
!            enddo
!  333 continue
!         enddo
!
!      enddo ! ip 
!
!      return
!      end
!!----------------------------------------------------------------------
!      subroutine ppiclf_comm_CheckAngularBC(xrank, yrank, zrank)
!!
!      implicit none
!!
!      include "PPICLF"
!!
!! Local:
!!
!      integer*4 xrank, yrank, zrank
!!
!! Output:
!!
!
!      SELECT CASE (ang_case)
!        CASE(1) ! general wedge ; 0 <= angle < 90
!!          print*, "Wedge CheckAngularBC"
!          xrank  = FLOOR((xrot(1)-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
!          yrank  = FLOOR((xrot(2)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
!          zrank  = FLOOR((xrot(3)-ppiclf_binb(5))/ppiclf_bins_dx(3))
!
!        CASE(2) ! quarter cylinder ; angle = 90
!!          print*, "Quarter Cylinder CheckAngularBC"
!          xrank  = FLOOR((abs(xrot(1))-ppiclf_binb(1))
!     >                    /ppiclf_bins_dx(1)) 
!          yrank  = FLOOR((abs(xrot(2))-ppiclf_binb(3))
!     >                   /ppiclf_bins_dx(2)) 
!          zrank  = FLOOR((xrot(3)-ppiclf_binb(5))/ppiclf_bins_dx(3))
!
!        CASE(3) ! half cylinder ; angle = 180
!          print*, "Half Cylinder CheckAngularBC"
!
!        CASE DEFAULT
!            call ppiclf_exittr('Invalid Ghost Rotational Case!$',
!     >       0.0d0 ,ppiclf_nid)
!          END SELECT
!
!      return
!      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! David's old binning method left below for now
!      finished(1) = 0
!      finished(2) = 0
!      finished(3) = 0
!      total_bin = 1 
!
!      BinMinLen(1) = MAX(BinMinLen(1),BinMinLen(2),BinMinLen(3))
!      do i=1,ppiclf_ndim
!         finished(i) = 0
!         exit_1_array(i) = ppiclf_bins_set(i)
!         exit_2_array(i) = 0
!         if (ppiclf_bins_set(i) .ne. 1) ppiclf_n_bins(i) = 1
!         ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                        ppiclf_binb(2*(i-1)+1)  ) / 
!     >                       ppiclf_n_bins(i)
!         ! Make sure exit_2 is not violated by user input
!         if (ppiclf_bins_dx(i) .lt. BinMinLen(1)) then
!            do while (ppiclf_bins_dx(i) .lt. BinMinLen(1))
!               ppiclf_n_bins(i) = max(1, ppiclf_n_bins(i)-1)
!               ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                              ppiclf_binb(2*(i-1)+1)  ) / 
!     >                             ppiclf_n_bins(i)
!         WRITE(*,*) "Inf. loop in CreateBin", i, 
!     >              ppiclf_bins_dx(i), BinMinLen(1)
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
!                if (ppiclf_bins_dx(i) .lt. BinMinLen(1)) then
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


c-----------------------------------------------------------------------
      subroutine ppiclf_gop( x, w, op, n)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      real*8 x(n), w(n)
      character*3 op
      integer*4 n
!
! Internal:
!
      integer*4 i, ie
!
      if (op.eq.'+  ') then
      call mpi_allreduce
     >        (x,w,n,MPI_DOUBLE_PRECISION,mpi_sum,ppiclf_comm,ie)
      elseif (op.EQ.'M  ') then
      call mpi_allreduce
     >        (x,w,n,MPI_DOUBLE_PRECISION,mpi_max,ppiclf_comm,ie)
      elseif (op.EQ.'m  ') then
      call mpi_allreduce
     >        (x,w,n,MPI_DOUBLE_PRECISION,mpi_min,ppiclf_comm,ie)
      elseif (op.EQ.'*  ') then
      call mpi_allreduce
     >        (x,w,n,MPI_DOUBLE_PRECISION,mpi_prod,ppiclf_comm,ie)
      endif

      do i=1,n
         x(i) = w(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_igop( x, w, op, n)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      integer*4 x(n), w(n)
      character*3 op
      integer*4 n
!
! Internal:
!
      integer*4 i, ierr
!
      if     (op.eq.'+  ') then
        call MPI_Allreduce (x,w,n,mpi_integer,mpi_sum ,ppiclf_comm,ierr)
      elseif (op.EQ.'M  ') then
        call MPI_Allreduce (x,w,n,mpi_integer,mpi_max ,ppiclf_comm,ierr)
      elseif (op.EQ.'m  ') then
        call MPI_Allreduce (x,w,n,mpi_integer,mpi_min ,ppiclf_comm,ierr)
      elseif (op.EQ.'*  ') then
        call MPI_Allreduce (x,w,n,mpi_integer,mpi_prod,ppiclf_comm,ierr)
      endif

      do i=1,n
         x(i) = w(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      integer*4 function ppiclf_iglsum(a,n)
! 
      implicit none
! 
! Input:
! 
      integer*4 a(1)
      integer*4 n
! 
! Internal:
! 
      integer*4 tsum
      integer*4 tmp(1),work(1)
      integer*4 i
!
      tsum= 0
      do i=1,n
         tsum=tsum+a(i)
      enddo
      tmp(1)=tsum
      call ppiclf_igop(tmp,work,'+  ',1)
      ppiclf_iglsum=tmp(1)
      return
      end
C-----------------------------------------------------------------------
      real*8 function ppiclf_glsum (x,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of A
      real*8 x(n),tsum
      real*8 tmp(1),work(1)
      integer*4 i,n
!
      TSUM = 0.0d0
      DO 100 I=1,N
         TSUM = TSUM+X(I)
 100  CONTINUE
      TMP(1)=TSUM
      CALL ppiclf_GOP(TMP,WORK,'+  ',1)
      ppiclf_GLSUM = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      real*8 function ppiclf_glmax(a,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of A
      REAL*8 A(n),tmax
      real*8 TMP(1),WORK(1)
      integer*4 i,n
!
      TMAX=-99.0e20
      DO 100 I=1,N
         TMAX=MAX(TMAX,A(I))
  100 CONTINUE
      TMP(1)=TMAX
      CALL ppiclf_GOP(TMP,WORK,'M  ',1)
      ppiclf_GLMAX=TMP(1)
      return
      END
c-----------------------------------------------------------------------
      integer*4 function ppiclf_iglmax(a,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of A
      integer*4 a(n),tmax
      integer*4 tmp(1),work(1)
      integer*4 i,n
!
      tmax= -999999999
      do i=1,n
         tmax=max(tmax,a(i))
      enddo
      tmp(1)=tmax
      call ppiclf_igop(tmp,work,'M  ',1)
      ppiclf_iglmax=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      real*8 function ppiclf_glmin(a,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of A
      REAL*8 A(n),tmin
      real*8 TMP(1),WORK(1)
      integer*4 i,n
!
      TMIN=99.0e20
      DO 100 I=1,N
         TMIN=MIN(TMIN,A(I))
  100 CONTINUE
      TMP(1)=TMIN
      CALL ppiclf_GOP(TMP,WORK,'m  ',1)
      ppiclf_GLMIN = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      integer*4 function ppiclf_iglmin(a,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of a
      integer*4 a(n),tmin
      integer*4 tmp(1),work(1)
      integer*4 i, n
!
      tmin=  999999999
      do i=1,n
         tmin=min(tmin,a(i))
      enddo
      tmp(1)=tmin
      call ppiclf_igop(tmp,work,'m  ',1)
      ppiclf_iglmin=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      real*8 function ppiclf_vlmin(vec,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of VEC
      REAL*8 VEC(n),tmin
      integer*4 i, n
!
      TMIN = 99.0E20
      DO 100 I=1,N
         TMIN = MIN(TMIN,VEC(I))
 100  CONTINUE
      ppiclf_VLMIN = TMIN
      return
      END
c-----------------------------------------------------------------------
      real*8 function ppiclf_vlmax(vec,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of VEC
      REAL*8 VEC(n),tmax
      integer*4 i, n
!
      TMAX =-99.0E20
      do i=1,n
         TMAX = MAX(TMAX,VEC(I))
      enddo
      ppiclf_VLMAX = TMAX
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_copy(a,b,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of A and B
      real*8 a(n),b(n)
      integer*4 i,n
!

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_icopy(a,b,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of A and B
      INTEGER*4 A(n), B(n)
      integer*4 i,n
!
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_chcopy(a,b,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed A and B dimenions
      CHARACTER*1 A(n), B(n)
      integer*4 i,n
!
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_exittr(stringi,rdata,idata)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
! 
! Vars:
! 
      character*1 stringi(132)
      character*1 stringo(132)
      character*25 s25
      integer*4 ilen, ierr, k, idata
      integer*4 ppiclf_indx1
      real*8 rdata
      external ppiclf_indx1
!
      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')
      write(s25,25) rdata,idata
   25 format(1x,1p1e14.6,i10)
      call ppiclf_chcopy(stringo(ilen),s25,25)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen+24)
    1 format('PPICLF: ERROR ',132a1)

c     call mpi_finalize (ierr)
      call mpi_abort(ppiclf_comm, 1, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_printsri(stringi,rdata,idata)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!

      character*1 stringi(132)
      character*1 stringo(132)
      character*25 s25
      integer*4 ilen, idata, k, ierr
      integer*4 ppiclf_indx1
      real*8 rdata
      external ppiclf_indx1
#ifdef TEST
      RETURN
#endif!
      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')
      write(s25,25) rdata,idata
   25 format(1x,1p1e14.6,i10)
      call ppiclf_chcopy(stringo(ilen),s25,25)

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen+24)
    1 format('PPICLF: ',132a1)

      call mpi_barrier(ppiclf_comm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_printsi(stringi,idata)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character*1 stringi(132)
      character*1 stringo(132)
      character*10 s10
      integer*4 ilen, idata, k, ierr
      integer*4 ppiclf_indx1
      external ppiclf_indx1
!
#ifdef TEST
      RETURN
#endif      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')
      write(s10,10) idata
   10 format(1x,i9)
      call ppiclf_chcopy(stringo(ilen),s10,10)

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen+9)
    1 format('PPICLF: ',132a1)

      call mpi_barrier(ppiclf_comm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_printsr(stringi,rdata)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character*1 stringi(132)
      character*1 stringo(132)
      character*15 s15
      integer*4 ilen, k, ierr
      integer*4 ppiclf_indx1
      real*8 rdata
      external ppiclf_indx1
!
#ifdef TEST
      RETURN
#endif      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')
      write(s15,15) rdata
   15 format(1x,1p1e14.6)
      call ppiclf_chcopy(stringo(ilen),s15,15)

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen+14)
    1 format('PPICLF: ',132a1)

      call mpi_barrier(ppiclf_comm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_prints(stringi)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character*1 stringi(132)
      character*1 stringo(132)
      integer*4 ilen, k, ierr
      integer*4 ppiclf_indx1
      external ppiclf_indx1
!
#ifdef TEST
      RETURN
#endif      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen-1)
    1 format('PPICLF: ',132a1)

      call mpi_barrier(ppiclf_comm,ierr)

      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE PPICLF_BLANK(A,N)
! 
      implicit none
! 
! Vars:
!
      ! TLJ changed dimension of A
      CHARACTER*1 A(N)
      CHARACTER*1 BLNK
      SAVE        BLNK
      DATA        BLNK /' '/
      integer*4 i,n
!
C
      DO 10 I=1,N
         A(I)=BLNK
   10 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      INTEGER*4 FUNCTION PPICLF_INDX1(S1,S2)
! 
      implicit none
! 
! Vars:
!
      CHARACTER*1 S1(132),S2
      integer*4 n1, i
!
      N1=132
      PPICLF_INDX1=0
      IF (N1.LT.1) return
C
      DO 100 I=1,N1
         IF (S1(I).EQ.S2) THEN
            PPICLF_INDX1=I
            return
         ENDIF
  100 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      character*132 FUNCTION PPICLF_CHSTR(S1,indx1)
! 
      implicit none
! 
! Vars:
!
      ! TLJ modified, but not sure why I had to
      CHARACTER S1
      INTEGER indx1
!
      PPICLF_CHSTR = S1(1:indx1)
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_byte_open_mpi(fnamei,mpi_fh,ifro,ierr)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character fnamei*(*)
      logical ifro
      CHARACTER*1 BLNK
      DATA BLNK/' '/
      character*132 fname
      character*1   fname1(132)
      equivalence  (fname1,fname)
      integer*4 imode, ierr, mpi_fh
!
      imode = MPI_MODE_WRONLY+MPI_MODE_CREATE
      if(ifro) then
        imode = MPI_MODE_RDONLY 
      endif

      call MPI_file_open(ppiclf_comm,fnamei,imode,
     &                   MPI_INFO_NULL,mpi_fh,ierr)

      return
      end
C--------------------------------------------------------------------------
      subroutine ppiclf_byte_read_mpi(buf,icount,mpi_fh,ierr)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      real*4 buf(1)          ! buffer
      integer*4 iout, icount, mpi_fh, ierr
!

      iout = icount ! icount is in 4-byte words
      call MPI_file_read_all(mpi_fh,buf,iout,MPI_REAL,
     &                       MPI_STATUS_IGNORE,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine ppiclf_byte_write_mpi(buf,icount,iorank,mpi_fh,ierr)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      real*4 buf(1)          ! buffer
      integer*4 icount, iorank, mpi_fh, ierr, iout
!

      iout = icount ! icount is in 4-byte words
      if(iorank.ge.0 .and. ppiclf_nid.ne.iorank) iout = 0
      call MPI_file_write_all(mpi_fh,buf,iout,MPI_REAL,
     &                        MPI_STATUS_IGNORE,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine ppiclf_byte_close_mpi(mpi_fh,ierr)
! 
      implicit none
! 
      include 'mpif.h'
!
! Vars:
!
      integer*4 mpi_fh, ierr
!

      call MPI_file_close(mpi_fh,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine ppiclf_byte_set_view(ioff_in,mpi_fh)
! 
      implicit none
! 
      include 'mpif.h'
!
! Vars:
!
      integer*8 ioff_in
      integer*4 mpi_fh, ierr
!
      call MPI_file_set_view(mpi_fh,ioff_in,MPI_BYTE,MPI_BYTE,
     &                       'native',MPI_INFO_NULL,ierr)

      return
      end
C--------------------------------------------------------------------------
      subroutine ppiclf_bcast(buf,len)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      real*4 buf(1)
      integer*4 len, ierr
!

      call mpi_bcast (buf,len,mpi_byte,0,ppiclf_comm,ierr)

      return
      end
C--------------------------------------------------------------------------
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_io_ReadParticleVTU(filein1,istartoutin,
     > npart,dp_max)
     > bind(C, name="ppiclc_io_ReadParticleVTU")
#else
      subroutine ppiclf_io_ReadParticleVTU(filein1,istartoutin,
     > npart,dp_max)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      character*1 filein1(132)
!
! Internal:
!
      real*4  rout_pos(3      *PPICLF_LPART) 
     >       ,rout_sln(PPICLF_LRS*PPICLF_LPART)
     >       ,rout_lrp(PPICLF_LRP*PPICLF_LPART)
     >       ,rout_lip(3      *PPICLF_LPART)
      real*8 dp_max
      character*1 dum_read
      character*132 filein2
      integer*8 idisp_pos,idisp_sln,idisp_lrp,idisp_lip,stride_len
      integer*4 vtu, isize, jx, jy, jz, ivtu_size, ifound, i, npt_total,
     >          npart_min, npmax, npart, ndiff, iorank, icount_pos,
     >          icount_sln, icount_lrp, icount_lip, j, ic_pos, ic_sln,
     >          ic_lrp, ic_lip, pth, ierr, indx1
      integer*4 ppiclf_indx1
      external ppiclf_indx1
      character*132 PPICLF_CHSTR
      EXTERNAL PPICLF_CHSTR
      ! Sam - this modifies the interface for ppiclC. Will now need to
      ! include NULL for optional arguments
      integer*4, optional :: istartoutin
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!
      call ppiclf_prints(' *Begin ReadParticleVTU$')
      
      call ppiclf_solve_InitZero
      PPICLF_RESTART = .true.

      if (present(istartoutin)) then
        istartout = istartoutin
      else
        istartout = 0
      end if


      indx1 = ppiclf_indx1(filein1,'.')
      indx1 = indx1 + 3 ! v (1) t (2) u (3)
      filein2 = ppiclf_chstr(filein1(1:indx1),indx1)
      if (ppiclf_nid == 0) then
         print*,'   * ParticleVTU filein2 ',trim(filein2)
      endif

      isize = 4
      jx    = 1
      jy    = 2
      jz    = 3

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=trim(filein2)
     >    ,access='stream',form="unformatted")

      ivtu_size = -1
      ifound = 0
      do i=1,1000000
      read(vtu) dum_read
      if (dum_read == '_') ifound = ifound + 1
      if (ifound .eq. 2) then
         ivtu_size = i
         exit
      endif
      enddo
      read(vtu) npt_total
      close(vtu)
      npt_total = npt_total/isize/3
      endif

      call ppiclf_bcast(npt_total,isize)
      call ppiclf_bcast(ivtu_size,isize)


      npart_min = npt_total/ppiclf_np+1
      if (npart_min*ppiclf_np .gt. npt_total) npart_min = npart_min-1

      if (npt_total .gt. PPICLF_LPART*ppiclf_np) 
     >   call ppiclf_exittr('Increase LPART to at least$',0.0d0
     >    ,npart_min)


      npmax = min(npt_total/PPICLF_LPART+1,ppiclf_np)
      stride_len = 0
      if (ppiclf_nid .le. npmax-1 .and. ppiclf_nid. ne. 0) 
     >stride_len = ppiclf_nid*PPICLF_LPART

      npart = PPICLF_LPART
      if (ppiclf_nid .gt. npmax-1) npart = 0

      ndiff = npt_total - (npmax-1)*PPICLF_LPART
      if (ppiclf_nid .eq. npmax-1) npart = ndiff

      iorank = -1

      call ppiclf_byte_open_mpi(trim(filein2),pth,.true.,ierr)

      idisp_pos = ivtu_size + isize*(3*stride_len + 1)
      icount_pos = npart*3   
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_read_mpi(rout_pos,icount_pos,pth,ierr)

      do i=1,PPICLF_LRS
         idisp_sln = ivtu_size + isize*(3*npt_total 
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + i)
         j = 1 + (i-1)*npart
         icount_sln = npart

         call ppiclf_byte_set_view(idisp_sln,pth)
         call ppiclf_byte_read_mpi(rout_sln(j),icount_sln
     >                            ,pth,ierr)
      enddo

      do i=1,PPICLF_LRP
         idisp_lrp = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + i)

         j = 1 + (i-1)*npart
         icount_lrp = npart

         call ppiclf_byte_set_view(idisp_lrp,pth)
         call ppiclf_byte_read_mpi(rout_lrp(j),icount_lrp
     >                            ,pth,ierr)
      enddo

      do i=1,3
         idisp_lip = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + PPICLF_LRP*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + PPICLF_LRP + i)

         j = 1 + (i-1)*npart
         icount_lip = npart

         call ppiclf_byte_set_view(idisp_lip,pth)
         call ppiclf_byte_read_mpi(rout_lip(j),icount_lip
     >                            ,pth,ierr)
      enddo

      call ppiclf_byte_close_mpi(pth,ierr)

      ic_pos = 0
      ic_sln = 0
      ic_lrp = 0
      ic_lip = 0
      do i=1,npart
         ic_pos = ic_pos + 1
         ppiclf_y(jx,i) = rout_pos(ic_pos)
         ic_pos = ic_pos + 1
         ppiclf_y(jy,i) = rout_pos(ic_pos)
         ic_pos = ic_pos + 1
         if (ppiclf_ndim .eq. 3) then
         ppiclf_y(jz,i) = rout_pos(ic_pos)
         endif
      enddo
      do j=1,PPICLF_LRS
      do i=1,npart
         ic_sln = ic_sln + 1
         ppiclf_y(j,i) = rout_sln(ic_sln)
      enddo
      enddo
      do j=1,PPICLF_LRP
      do i=1,npart
         ic_lrp = ic_lrp + 1
         ppiclf_rprop(j,i) = rout_lrp(ic_lrp)
      enddo
!*** need to add ppiclf_rprop2 ppiclf_rprop3, rprop4, & prop5
      enddo
      ! This reads the particle tag infomation.
      do j=1,3
      do i=1,npart
         ic_lip = ic_lip + 1
         ppiclf_iprop(j,i) = int(rout_lip(ic_lip))
      enddo
      enddo

      ppiclf_npart = npart

      dp_max = MAXVAL(ppiclf_rprop(PPICLF_R_JDP,:))
      call ppiclf_printsi('  End ReadParticleVTU$',npt_total)

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_io_ReadWallVTK(filein1)
     > bind(C, name="ppiclc_io_ReadWallVTK")
#else
      subroutine ppiclf_io_ReadWallVTK(filein1)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      character*1 filein1(132)
!
! Internal:
!
      real*8 points(3,4*PPICLF_LWALL)
      integer*4 fid, nmax, i, j, i1, i2, i3, isize, irsize, ierr,
     >          npoints, nwalls, indx1
      character*1000 text
      character*132 filein2
      integer*4 ppiclf_indx1
      external ppiclf_indx1
      character*132 PPICLF_CHSTR
      external PPICLF_CHSTR
!
      !! THROW ERRORS HERE IN FUTURE
      call ppiclf_prints(' *Begin ReadWallVTK$')
      
      indx1 = ppiclf_indx1(filein1,'.')
      indx1 = indx1 + 3 ! v (1) t (2) k (3)
      !filein2 = ppiclf_chstr(filein1(1:indx1))
      filein2 = 'filein.vtk' !ppiclf_chstr(filein1(1:indx1))

      if (ppiclf_nid .eq. 0) then

      fid = 432

      open (unit=fid,file=trim(filein2),action="read")
      
      nmax = 10000
      do i=1,nmax
         read (fid,*,iostat=ierr) text,npoints

         do j=1,npoints
            read(fid,*) points(1,j),points(2,j),points(3,j)
         enddo

         read (fid,*,iostat=ierr) text,nwalls

         do j=1,nwalls
            if (ppiclf_ndim .eq. 2) then
               read(fid,*) i1,i2
    
               i1 = i1 + 1
               i2 = i2 + 1

               call ppiclf_solve_InitWall( 
     >                 (/points(1,i1),points(2,i1)/),
     >                 (/points(1,i2),points(2,i2)/),
     >                 (/points(1,i1),points(2,i1)/))  ! dummy 2d
    
            elseif (ppiclf_ndim .eq. 3) then
               read(fid,*) i1,i2,i3

               i1 = i1 + 1
               i2 = i2 + 1
               i3 = i3 + 1
    
               call ppiclf_solve_InitWall( 
     >                 (/points(1,i1),points(2,i1),points(3,i1)/),
     >                 (/points(1,i2),points(2,i2),points(3,i2)/),
     >                 (/points(1,i3),points(2,i3),points(3,i3)/))

            endif
         enddo
            
         exit
      enddo

      close(fid)

      endif

      isize  = 4
      irsize = 8
      call ppiclf_bcast(ppiclf_nwall, isize)
      call ppiclf_bcast(ppiclf_wall_c,9*PPICLF_LWALL*irsize)
      call ppiclf_bcast(ppiclf_wall_n,4*PPICLF_LWALL*irsize)

      call ppiclf_printsi('  End ReadWallVTK$',nwalls)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteBinVTU(filein1)
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) filein1
!
! Internal:
!
      character*3 filein
      character*12 vtufile
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth, nvtx_total, ncll_total
      integer*8 idisp_pos,idisp_cll,stride_lenv(8),stride_lenc
      integer*4 iint, nnp, nxx, ndxgpp1, ndygpp1, ndxygpp1, if_sz, ibin,
     >          jbin, kbin, il, ir, jl, jr, kl, kr, nbinpa, nbinpb,
     >          nbinpc, nbinpd, nbinpe, nbinpf, nbinpg, nbinph, i, j, k,
     >          ii, npa, npb, npc, npd, npe, npf, npg, nph, 
     >          ioff_dum, itype, iorank, if_cll, if_pos, icount_pos, 
     >          icount_cll, ierr, isize, ivtu_size
      real*4 rpoint(3)
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!


      if (icalld1 .eq. 0) icalld1 = istartout

      icalld1 = icalld1+1

      !ppiclf_printbinvtu set true in ppiclf_comm_CreateBin
      IF(ppiclf_printbinvtu .OR. ppiclf_time .EQ. 0.0) THEN
        call ppiclf_printsi(' *Begin WriteBinVTU$',ppiclf_cycle)
        ! Set false for next iteration
        ppiclf_printbinvtu = .FALSE. 
      ELSE
        RETURN
      END IF

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      nvtx_total = (ppiclf_n_bins(1)+1)*(ppiclf_n_bins(2)+1)
      if (ppiclf_ndim .gt. 2) 
     >    nvtx_total = nvtx_total*(ppiclf_n_bins(3)+1)
      ncll_total = ppiclf_n_bins(1)*ppiclf_n_bins(2)
      if (ppiclf_ndim .gt. 2) ncll_total = ncll_total*ppiclf_n_bins(3)

      ndxgpp1 = ppiclf_n_bins(1) + 1
      ndygpp1 = ppiclf_n_bins(2) + 1
      ndxygpp1 = ndxgpp1*ndygpp1

      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'bin'
      else 
         write(filein,'(A3)') filein1
      endif

      isize  = 4

      ! get which bin this processor holds
      ibin = modulo(ppiclf_nid,ppiclf_n_bins(1))
      jbin = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
      kbin = 0
      if (ppiclf_ndim .eq. 3)
     >kbin = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))

      il = ibin
      ir = ibin+1
      jl = jbin
      jr = jbin+1
      kl = kbin
      kr = kbin
      if (ppiclf_ndim .eq. 3) then
         kl = kbin
         kr = kbin+1
      endif

      nbinpa = il + ndxgpp1*jl + ndxygpp1*kl
      nbinpb = ir + ndxgpp1*jl + ndxygpp1*kl
      nbinpc = il + ndxgpp1*jr + ndxygpp1*kl
      nbinpd = ir + ndxgpp1*jr + ndxygpp1*kl
      if (ppiclf_ndim .eq. 3) then
         nbinpe = il + ndxgpp1*jl + ndxygpp1*kr
         nbinpf = ir + ndxgpp1*jl + ndxygpp1*kr
         nbinpg = il + ndxgpp1*jr + ndxygpp1*kr
         nbinph = ir + ndxgpp1*jr + ndxygpp1*kr
      endif

      stride_lenv(1) = 0
      stride_lenv(2) = 0
      stride_lenv(3) = 0
      stride_lenv(4) = 0
      stride_lenv(5) = 0
      stride_lenv(6) = 0
      stride_lenv(7) = 0
      stride_lenv(8) = 0
 
      stride_lenc = 0
      if (ppiclf_nid .le. 
     >      ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         stride_lenv(1) = nbinpa
         stride_lenv(2) = nbinpb
         stride_lenv(3) = nbinpc
         stride_lenv(4) = nbinpd
         if (ppiclf_ndim .eq. 3) then
            stride_lenv(5) = nbinpe
            stride_lenv(6) = nbinpf
            stride_lenv(7) = nbinpg
            stride_lenv(8) = nbinph
         endif
         
         stride_lenc    = ppiclf_nid
      endif

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (ppiclf_iendian .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (ppiclf_iendian .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') ppiclf_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') ppiclf_cycle
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='yes') '  </FieldData>'
      write(vtu,'(A)',advance='no') '  <Piece '
      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
      write(vtu,'(I0)',advance='no') nvtx_total
      write(vtu,'(A)',advance='no') '" NumberOfCells="'
      write(vtu,'(I0)',advance='no') ncll_total
      write(vtu,'(A)',advance='yes') '"> '

! -----------
! COORDINATES 
! -----------
      iint = 0
      write(vtu,'(A)',advance='yes') '   <Points>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
      iint = iint + 3*isize*nvtx_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'
      write(vtu,'(A)',advance='yes') '   </PointData> '
      write(vtu,'(A)',advance='yes') '   <CellData>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"PPR",1,iint)
      iint = iint + 1   *isize*ncll_total + isize
      write(vtu,'(A)',advance='yes') '   </CellData> '

! ----------
! END MATTER
! ----------
      write(vtu,'(A)',advance='yes') '   <Cells> '
      write(vtu,'(A)',advance='no')  '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="connectivity" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write connectivity here
      do ii=0,ncll_total-1
         i = modulo(ii,ppiclf_n_bins(1))
         j = modulo(ii/ppiclf_n_bins(1),ppiclf_n_bins(2))
         k = 0
         if (ppiclf_ndim .eq. 3)
     >   k = ii/(ppiclf_n_bins(1)*ppiclf_n_bins(2))
          
c     do K=0,ppiclf_n_bins(3)-1
         kl = K
         kr = K
         if (ppiclf_ndim .eq. 3) then
            kl = K
            kr = K+1
         endif
c     do J=0,ppiclf_n_bins(2)-1
         jl = J
         jr = J+1
c     do I=0,ppiclf_n_bins(1)-1
         il = I
         ir = I+1

         if (ppiclf_ndim .eq. 3) then
            npa = il + ndxgpp1*jl + ndxygpp1*kl
            npb = ir + ndxgpp1*jl + ndxygpp1*kl
            npc = il + ndxgpp1*jr + ndxygpp1*kl
            npd = ir + ndxgpp1*jr + ndxygpp1*kl
            npe = il + ndxgpp1*jl + ndxygpp1*kr
            npf = ir + ndxgpp1*jl + ndxygpp1*kr
            npg = il + ndxgpp1*jr + ndxygpp1*kr
            nph = ir + ndxgpp1*jr + ndxygpp1*kr
            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0,A)',advance='no')  npd, ' '
            write(vtu,'(I0,A)',advance='no')  npe, ' '
            write(vtu,'(I0,A)',advance='no')  npf, ' '
            write(vtu,'(I0,A)',advance='no')  npg, ' '
            write(vtu,'(I0)'  ,advance='yes') nph
         else
            npa = il + ndxgpp1*jl + ndxygpp1*kl
            npb = ir + ndxgpp1*jl + ndxygpp1*kl
            npc = il + ndxgpp1*jr + ndxygpp1*kl
            npd = ir + ndxgpp1*jr + ndxygpp1*kl
            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0)'  ,advance='yes') npd
         endif
c     enddo
c     enddo
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write offsetts here
      ioff_dum = 4
      if (ppiclf_ndim .eq. 3) ioff_dum = 8
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') ioff_dum*i
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="UInt8" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      itype = 8
      if (ppiclf_ndim .eq. 3) itype = 11
      ! write types here
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') itype
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='yes') '   </Cells> '
      write(vtu,'(A)',advance='yes') '  </Piece> '
      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
      close(vtu)

c1511 continue

      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >    ,position='append')
      write(vtu) '_'
      close(vtu)

      inquire(file=vtufile,size=ivtu_size)
      endif

      call ppiclf_bcast(ivtu_size, isize)

      iorank = -1

      if_pos = 3*isize*nvtx_total


      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif


      call mpi_barrier(ppiclf_comm,ierr)

      ! write points first
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      ! point A
      icount_pos = 0
      if (ppiclf_nid .le. 
     >    ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         icount_pos = 3
      endif
      idisp_pos  = ivtu_size + isize*(3*stride_lenv(1) + 1)
      rpoint(1)  = sngl(ppiclf_bin_pos(1,1))
      rpoint(2)  = sngl(ppiclf_bin_pos(1,2))
      rpoint(3)  = 0.0
      if (ppiclf_ndim .eq. 3)
     >rpoint(3)  = sngl(ppiclf_bin_pos(1,3))
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)

      ! 3d
      if (ppiclf_ndim .gt. 2) then

         ! point B
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(2) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(2,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(1,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(1,3))
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point C
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(3) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(1,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(2,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(1,3))
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point E
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(5) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(1,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(1,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(2,3))
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point D
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(4) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(2,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(2,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(1,3))
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point F
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(6) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(2,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(1,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(2,3))
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point G
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(7) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(1,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(2,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(2,3))
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point H
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(8) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(2,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(2,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(2,3))
         endif
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)


      ! 2d
      else

         ! point B
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3*stride_lenv(2) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(2,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(1,2))
            rpoint(3)  = 0.0
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point C
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3*stride_lenv(3) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(1,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(2,2))
            rpoint(3)  = 0.0
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point D
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3*stride_lenv(4) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(2,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(2,2))
            rpoint(3)  = 0.0
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)

      endif

      call ppiclf_byte_close_mpi(pth,ierr)

      call mpi_barrier(ppiclf_comm,ierr)

      if_cll = 1*isize*ncll_total

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_cll
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write points first
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      !
      ! cell values
      idisp_cll = ivtu_size + isize*(3*(nvtx_total) 
     >     + 1*stride_lenc + 2)
      icount_cll = 0
      if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         icount_cll = 1
      endif
      rpoint(1)  = real(nxx)
      call ppiclf_byte_set_view(idisp_cll,pth)
      call ppiclf_byte_write_mpi(rpoint,icount_cll,iorank,pth,ierr)

      call ppiclf_byte_close_mpi(pth,ierr)

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi('  End WriteBinVTU$',ppiclf_cycle)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteParticleVTU(filein1)
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) filein1
!
! Internal:
!
      real*4  rout_pos(3      *PPICLF_LPART) 
     >       ,rout_sln(PPICLF_LRS*PPICLF_LPART)
     >       ,rout_lrp(PPICLF_LRP*PPICLF_LPART)
     >       ,rout_lip(3      *PPICLF_LPART)
      character*3 filein
      character*12 vtufile
      character*6  prostr
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth,prevs(2,ppiclf_np)
      integer*8 idisp_pos,idisp_sln,idisp_lrp,idisp_lip,stride_len
      integer*4 iint, nnp, nxx, npt_total, jx, jy, jz, if_sz, isize,
     >          iadd, if_pos, if_sln, if_lrp, if_lip, ic_pos, ic_sln,
     >          ic_lrp, ic_lip, i, j, ie, nps, nglob, nkey, ndum,
     >          icount_pos, icount_sln, icount_lrp, icount_lip, iorank,
     >          ierr, ivtu_size
      integer*4 ppiclf_iglsum
      external ppiclf_iglsum
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!

      call ppiclf_printsi(' *Begin WriteParticleVTU$',ppiclf_cycle)

      if (icalld1 .eq. 0) icalld1 = istartout

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      npt_total = ppiclf_iglsum(nxx,1)

      jx    = 1
      jy    = 2
      jz    = 1
      if (ppiclf_ndim .eq. 3)
     >jz    = 3

      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'par'
      else 
         write(filein,'(A3)') filein1
      endif

! --------------------------------------------------
! COPY PARTICLES TO OUTPUT ARRAY
! --------------------------------------------------

      isize = 4

      iadd = 0
      if_pos = 3*isize*npt_total
      if_sln = 1*isize*npt_total
      if_lrp = 1*isize*npt_total
      if_lip = 1*isize*npt_total

      ic_pos = iadd
      ic_sln = iadd
      ic_lrp = iadd
      ic_lip = iadd
      do i=1,nxx

         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(ppiclf_y(jx,i))
         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(ppiclf_y(jy,i))
         ic_pos = ic_pos + 1
         if (ppiclf_ndim .eq. 3) then
            rout_pos(ic_pos) = sngl(ppiclf_y(jz,i))
         else
            rout_pos(ic_pos) = 0.0
         endif
      enddo
      do j=1,PPICLF_LRS
      do i=1,nxx
         ic_sln = ic_sln + 1
         rout_sln(ic_sln) = sngl(ppiclf_y(j,i))
      enddo
      enddo
      do j=1,PPICLF_LRP
      do i=1,nxx
         ic_lrp = ic_lrp + 1
         rout_lrp(ic_lrp) = sngl(ppiclf_rprop(j,i))
      enddo
      enddo
      do j=1,3
      do i=1,nxx
         ! This prints out the particle tag info
         ic_lip = ic_lip + 1
         rout_lip(ic_lip) = real(ppiclf_iprop(j,i))
      enddo
      enddo

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------
      do i=1,nnp
         prevs(1,i) = i-1
         prevs(2,i) = nxx
      enddo

      nps   = 1 ! index of new proc for doing stuff
      nglob = 1 ! unique key to sort by
      nkey  = 1 ! number of keys (just 1 here)
      ndum = 2
      call pfgslib_crystal_ituple_transfer(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nnp,nps)
      call pfgslib_crystal_ituple_sort(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nglob,nkey)

      stride_len = 0
      if (ppiclf_nid .ne. 0) then
      do i=1,ppiclf_nid
         stride_len = stride_len + prevs(2,i)
      enddo
      endif

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (ppiclf_iendian .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (ppiclf_iendian .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') ppiclf_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') ppiclf_cycle
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='yes') '  </FieldData>'
      write(vtu,'(A)',advance='no') '  <Piece '
      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
      write(vtu,'(I0)',advance='no') npt_total
      write(vtu,'(A)',advance='yes') '" NumberOfCells="0"> '

! -----------
! COORDINATES 
! -----------
      iint = 0
      write(vtu,'(A)',advance='yes') '   <Points>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
      iint = iint + 3*isize*npt_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'


      do ie=1,PPICLF_LRS
         write(prostr,'(A1,I2.2)') "y",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      do ie=1,PPICLF_LRP
         write(prostr,'(A4,I2.2)') "rprop",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      do ie=1,3
         write(prostr,'(A3,I2.2)') "tag",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      write(vtu,'(A)',advance='yes') '   </PointData> '

! ----------
! END MATTER
! ----------
      write(vtu,'(A)',advance='yes') '   <Cells> '
      write(vtu,'(A)',advance='no')  '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="connectivity" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='yes') '   </Cells> '
      write(vtu,'(A)',advance='yes') '  </Piece> '
      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
      close(vtu)

      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >    ,position='append')
      write(vtu) '_'
      close(vtu)

      inquire(file=vtufile,size=ivtu_size)
      endif

      call ppiclf_bcast(ivtu_size, isize)

      ! byte-displacements
      idisp_pos = ivtu_size + isize*(3*stride_len + 1)

      ! how much to write
      icount_pos = 3*nxx
      icount_sln = 1*nxx
      icount_lrp = 1*nxx
      icount_lip = 1*nxx

      iorank = -1

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_write_mpi(rout_pos,icount_pos,iorank,pth,ierr)
      call ppiclf_byte_close_mpi(pth,ierr)

      call mpi_barrier(ppiclf_comm,ierr)

      do i=1,PPICLF_LRS
         idisp_sln = ivtu_size + isize*(3*npt_total 
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + i)

         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_sln
           close(vtu)
         endif
   
         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*ppiclf_npart + 1
   
         ! write
         call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
         call ppiclf_byte_set_view(idisp_sln,pth)
         call ppiclf_byte_write_mpi(rout_sln(j),icount_sln,iorank,pth
     >                             ,ierr)
         call ppiclf_byte_close_mpi(pth,ierr)
      enddo

      do i=1,PPICLF_LRP
         idisp_lrp = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + i)

         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_lrp
           close(vtu)
         endif
   
         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*ppiclf_npart + 1
   
         ! write
         call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
         call ppiclf_byte_set_view(idisp_lrp,pth)
         call ppiclf_byte_write_mpi(rout_lrp(j),icount_lrp,iorank,pth
     >                             ,ierr)
         call ppiclf_byte_close_mpi(pth,ierr)
      enddo

      do i=1,3
         idisp_lip = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + PPICLF_LRP*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + PPICLF_LRP + i)
         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_lip
           close(vtu)
         endif

         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*ppiclf_npart + 1
   
         ! write
         call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
         call ppiclf_byte_set_view(idisp_lip,pth)
         call ppiclf_byte_write_mpi(rout_lip(j),icount_lip,iorank,pth
     >                             ,ierr)
         call ppiclf_byte_close_mpi(pth,ierr)
      enddo

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi(' *End WriteParticleVTU$',ppiclf_cycle)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteDataArrayVTU(vtu,dataname,ncomp,idist)
!
      implicit none
!
! Input:
!
      integer*4 vtu,ncomp
      integer*4 idist
      character (len = *) dataname
!
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="'
      write(vtu,'(A)',advance='no') dataname
      write(vtu,'(A)',advance='no') '" NumberOfComponents="'
      write(vtu,'(I0)',advance='no') ncomp
      write(vtu,'(A)',advance='no') '" format="append" '
      write(vtu,'(A)',advance='no') 'offset="'
      write(vtu,'(I0)',advance='no') idist
      write(vtu,'(A)',advance='yes') '"/>'

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagAll
!
      implicit none
!
      include "PPICLF"
!
      call ppiclf_prints('*********** PPICLF OUTPUT *****************$')
      call ppiclf_io_OutputDiagGen
      call ppiclf_io_OutputDiagGhost
      if (ppiclf_overlap) call ppiclf_io_OutputDiagGrid

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGen
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 npart_max, npart_min, npart_tot, npart_ide, nbin_total
      integer*4 ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
      external ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
!
      call ppiclf_prints(' *Begin General Info$')
         npart_max = ppiclf_iglmax(ppiclf_npart,1)
         npart_min = ppiclf_iglmin(ppiclf_npart,1)
         npart_tot = ppiclf_iglsum(ppiclf_npart,1)
         npart_ide = npart_tot/ppiclf_np

         nbin_total = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)

      call ppiclf_printsi('  -Cycle                  :$',ppiclf_cycle)
      call ppiclf_printsi('  -Output Freq.           :$',ppiclf_iostep)
      call ppiclf_printsr('  -Time                   :$',ppiclf_time)
      call ppiclf_printsr('  -dt                     :$',ppiclf_dt)
      call ppiclf_printsi('  -Global particles       :$',npart_tot)
      call ppiclf_printsi('  -Local particles (Max)  :$',npart_max)
      call ppiclf_printsi('  -Local particles (Min)  :$',npart_min)
      call ppiclf_printsi('  -Local particles (Ideal):$',npart_ide)
      call ppiclf_printsi('  -Total ranks            :$',ppiclf_np)
      call ppiclf_printsi('  -Problem dimensions     :$',ppiclf_ndim)
      call ppiclf_printsi('  -Integration method     :$',ppiclf_imethod)
      call ppiclf_printsi('  -Number of bins total   :$',nbin_total)
      call ppiclf_printsi('  -Number of bins (x)     :$',
     >                    ppiclf_n_bins(1))
      call ppiclf_printsi('  -Number of bins (y)     :$'
     >                    ,ppiclf_n_bins(2))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsi('  -Number of bins (z)     :$'
     >                    ,ppiclf_n_bins(3))
      call ppiclf_printsr('  -Bin xl coordinate      :$',ppiclf_binb(1))
      call ppiclf_printsr('  -Bin xr coordinate      :$',ppiclf_binb(2))
      call ppiclf_printsr('  -Bin yl coordinate      :$',ppiclf_binb(3))
      call ppiclf_printsr('  -Bin yr coordinate      :$',ppiclf_binb(4))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin zl coordinate      :$',ppiclf_binb(5))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin zr coordinate      :$',ppiclf_binb(6))

      call ppiclf_prints('  End General Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGrid
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 nel_max_orig, nel_min_orig, nel_total_orig, nel_max_map,
     >          nel_min_map, nel_total_map
      integer*4 ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
      external ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
!
      call ppiclf_prints(' *Begin Grid Info$')

         nel_max_orig   = ppiclf_iglmax(ppiclf_nFVCells,1)
         nel_min_orig   = ppiclf_iglmin(ppiclf_nFVCells,1)
         nel_total_orig = ppiclf_iglsum(ppiclf_nFVCells,1)

         nel_max_map   = ppiclf_iglmax(ppiclf_nCells_FV2PICL,1)
         nel_min_map   = ppiclf_iglmin(ppiclf_nCells_FV2PICL,1)
         nel_total_map = ppiclf_iglsum(ppiclf_nCells_FV2PICL,1)

      call ppiclf_printsi('  -Orig. Global cells     :$',nel_total_orig)
      call ppiclf_printsi('  -Orig. Local cells (Max):$',nel_max_orig)
      call ppiclf_printsi('  -Orig. Local cells (Min):$',nel_min_orig)
      call ppiclf_printsi('  -Map Global cells       :$',nel_total_map)
      call ppiclf_printsi('  -Map Local cells (Max)  :$',nel_max_map)
      call ppiclf_printsi('  -Map Local cells (Min)  :$',nel_min_map)

      call ppiclf_prints('  End Grid Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 npart_max, npart_min, npart_tot
      integer*4 ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
      external ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
!
      call ppiclf_prints(' *Begin Ghost Info$')

         npart_max = ppiclf_iglmax(ppiclf_npart_gp,1)
         npart_min = ppiclf_iglmin(ppiclf_npart_gp,1)
         npart_tot = ppiclf_iglsum(ppiclf_npart_gp,1)

      call ppiclf_printsi('  -Global ghosts          :$',npart_tot)
      call ppiclf_printsi('  -Local ghosts (Max)     :$',npart_max)
      call ppiclf_printsi('  -Local ghosts (Min)     :$',npart_min)

      call ppiclf_prints('  End Ghost Info$')

      return
      end
!-----------------------------------------------------------------------
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
      REAL*8    dSQl, dSQi, dSQ(28), xp(3), dl, 
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
      LOGICAL   remove
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
        IF(nnearest .LT. 1) remove = .TRUE.

        ! Remove particle outside of any cell
        IF(nnearest .LT. 27 .AND. nnearest .GE. 1) THEN
         ie = CellID_nearest(1)
         DO l = 1,3
          ! similar distance check as above, but only for the cell that is
          ! closest to the particle (1st index in nnearest)
          ! used ABS since distance isn't squared.
          ! ppiclf_interp_dchk is set to be 1.5xmax cell length per
          ! dimension (in SUBROUTINE ppiclf_solve_InterpTupleTransfer)
          IF(ppiclf_linperiodic(l) .AND.
     >                            ppiclf_EqualDomain(l)) THEN
            dl = ABS(MIN((ppiclf_picl_grid(l,ie) - xp(l)), 
     >             (binblength(l)-ABS(ppiclf_picl_grid(l,ie)
     >              - xp(l)))))
          ELSE
            dl = ABS(ppiclf_picl_grid(l,ie) - xp(l))
          END IF
          ! Ensure particle is within 1/2 cell distance of one cell.
          IF(dl .GT. ppiclf_interp_dchk(l)/1.5D0*0.5D0) remove = .TRUE.
         END DO
        END IF

        IF (remove) THEN
          ! Particle is outside of fluid domain.
          ! iprop(8,ip) set to -1 means it will be removed
          ! from ppiclf_y & ppiclf_rprop, rprop2, rprop3, rprop4, rprop5
          ppiclf_iprop(9,ip) = -1
          ppiclf_remove_particle = .TRUE.
          PRINT*, 'part # on proc # removed.', ppiclf_iprop(1,ip),
     >            ppiclf_nid
          remove = .FALSE.
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
c-----------------------------------------------------------------------
      subroutine pfgslib_userExitHandler()
!
      implicit none
!
! This has been added so when linked with another code that links
! to gslib, there are no naming conflicts. See USREXIT=1 flag in
! install script for gslib.
!
      return
      end
c-----------------------------------------------------------------------
      subroutine pfgslib_mxm(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3
      real*8 a(n1,n2),b(n2,n3),c(n1,n3)
!
      call ppiclf_mxmf2(a,n1,b,n2,c,n3)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxmf2(A,N1,B,N2,C,N3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3
      real*8 a(n1,n2),b(n2,n3),c(n1,n3)
!
      if (n2.le.8) then
         if (n2.eq.1) then
            call ppiclf_mxf1(a,n1,b,n2,c,n3)
         elseif (n2.eq.2) then
            call ppiclf_mxf2(a,n1,b,n2,c,n3)
         elseif (n2.eq.3) then
            call ppiclf_mxf3(a,n1,b,n2,c,n3)
         elseif (n2.eq.4) then
            call ppiclf_mxf4(a,n1,b,n2,c,n3)
         elseif (n2.eq.5) then
            call ppiclf_mxf5(a,n1,b,n2,c,n3)
         elseif (n2.eq.6) then
            call ppiclf_mxf6(a,n1,b,n2,c,n3)
         elseif (n2.eq.7) then
            call ppiclf_mxf7(a,n1,b,n2,c,n3)
         else
            call ppiclf_mxf8(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.16) then
         if (n2.eq.9) then
            call ppiclf_mxf9(a,n1,b,n2,c,n3)
         elseif (n2.eq.10) then
            call ppiclf_mxf10(a,n1,b,n2,c,n3)
         elseif (n2.eq.11) then
            call ppiclf_mxf11(a,n1,b,n2,c,n3)
         elseif (n2.eq.12) then
            call ppiclf_mxf12(a,n1,b,n2,c,n3)
         elseif (n2.eq.13) then
            call ppiclf_mxf13(a,n1,b,n2,c,n3)
         elseif (n2.eq.14) then
            call ppiclf_mxf14(a,n1,b,n2,c,n3)
         elseif (n2.eq.15) then
            call ppiclf_mxf15(a,n1,b,n2,c,n3)
         else
            call ppiclf_mxf16(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.24) then
         if (n2.eq.17) then
            call ppiclf_mxf17(a,n1,b,n2,c,n3)
         elseif (n2.eq.18) then
            call ppiclf_mxf18(a,n1,b,n2,c,n3)
         elseif (n2.eq.19) then
            call ppiclf_mxf19(a,n1,b,n2,c,n3)
         elseif (n2.eq.20) then
            call ppiclf_mxf20(a,n1,b,n2,c,n3)
         elseif (n2.eq.21) then
            call ppiclf_mxf21(a,n1,b,n2,c,n3)
         elseif (n2.eq.22) then
            call ppiclf_mxf22(a,n1,b,n2,c,n3)
         elseif (n2.eq.23) then
            call ppiclf_mxf23(a,n1,b,n2,c,n3)
         elseif (n2.eq.24) then
            call ppiclf_mxf24(a,n1,b,n2,c,n3)
         endif
      else
         call ppiclf_mxm44_0(a,n1,b,n2,c,n3)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf1(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,1),b(1,n3),c(n1,n3)
!
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf2(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,2),b(2,n3),c(n1,n3)
!
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf3(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,3),b(3,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf4(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,4),b(4,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf5(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,5),b(5,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf6(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,6),b(6,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf7(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,7),b(7,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf8(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,8),b(8,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf9(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,9),b(9,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf10(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,10),b(10,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf11(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,11),b(11,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf12(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,12),b(12,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf13(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,13),b(13,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf14(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,14),b(14,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf15(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,15),b(15,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf16(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,16),b(16,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf17(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,17),b(17,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf18(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,18),b(18,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf19(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,19),b(19,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf20(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,20),b(20,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf21(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,21),b(21,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf22(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,22),b(22,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf23(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,23),b(23,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
     $             + a(i,23)*b(23,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf24(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,24),b(24,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
     $             + a(i,23)*b(23,j)
     $             + a(i,24)*b(24,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxm44_0(a, m, b, k, c, n)
!
      implicit none
!
! Internal:
!
      integer*4 m, k, n, i, j, l, m1, n1, nresid, mresid
      real*8 a(m,k), b(k,n), c(m,n)
      real*8 s11, s12, s13, s14, s21, s22, s23, s24
      real*8 s31, s32, s33, s34, s41, s42, s43, s44
!
      mresid = iand(m,3) 
      nresid = iand(n,3) 
      m1 = m - mresid + 1
      n1 = n - nresid + 1

      do i=1,m-mresid,4
        do j=1,n-nresid,4
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s41 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          s42 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0
          s43 = 0.0d0
          s14 = 0.0d0
          s24 = 0.0d0
          s34 = 0.0d0
          s44 = 0.0d0
          do l=1,k
            s11 = s11 + a(i,l)*b(l,j)
            s12 = s12 + a(i,l)*b(l,j+1)
            s13 = s13 + a(i,l)*b(l,j+2)
            s14 = s14 + a(i,l)*b(l,j+3)

            s21 = s21 + a(i+1,l)*b(l,j)
            s22 = s22 + a(i+1,l)*b(l,j+1)
            s23 = s23 + a(i+1,l)*b(l,j+2)
            s24 = s24 + a(i+1,l)*b(l,j+3)

            s31 = s31 + a(i+2,l)*b(l,j)
            s32 = s32 + a(i+2,l)*b(l,j+1)
            s33 = s33 + a(i+2,l)*b(l,j+2)
            s34 = s34 + a(i+2,l)*b(l,j+3)

            s41 = s41 + a(i+3,l)*b(l,j)
            s42 = s42 + a(i+3,l)*b(l,j+1)
            s43 = s43 + a(i+3,l)*b(l,j+2)
            s44 = s44 + a(i+3,l)*b(l,j+3)
          enddo
          c(i,j)     = s11 
          c(i,j+1)   = s12 
          c(i,j+2)   = s13
          c(i,j+3)   = s14

          c(i+1,j)   = s21 
          c(i+2,j)   = s31 
          c(i+3,j)   = s41 

          c(i+1,j+1) = s22
          c(i+2,j+1) = s32
          c(i+3,j+1) = s42

          c(i+1,j+2) = s23
          c(i+2,j+2) = s33
          c(i+3,j+2) = s43

          c(i+1,j+3) = s24
          c(i+2,j+3) = s34
          c(i+3,j+3) = s44
        enddo
* Residual when n is not multiple of 4
        if (nresid .ne. 0) then
          if (nresid .eq. 1) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,n)
              s21 = s21 + a(i+1,l)*b(l,n)
              s31 = s31 + a(i+2,l)*b(l,n)
              s41 = s41 + a(i+3,l)*b(l,n)
            enddo
            c(i,n)     = s11 
            c(i+1,n)   = s21 
            c(i+2,n)   = s31 
            c(i+3,n)   = s41 
          elseif (nresid .eq. 2) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,j)
              s12 = s12 + a(i,l)*b(l,j+1)

              s21 = s21 + a(i+1,l)*b(l,j)
              s22 = s22 + a(i+1,l)*b(l,j+1)

              s31 = s31 + a(i+2,l)*b(l,j)
              s32 = s32 + a(i+2,l)*b(l,j+1)

              s41 = s41 + a(i+3,l)*b(l,j)
              s42 = s42 + a(i+3,l)*b(l,j+1)
            enddo
            c(i,j)     = s11 
            c(i,j+1)   = s12

            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 

            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42
          else
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            s13 = 0.0d0
            s23 = 0.0d0
            s33 = 0.0d0
            s43 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,j)
              s12 = s12 + a(i,l)*b(l,j+1)
              s13 = s13 + a(i,l)*b(l,j+2)

              s21 = s21 + a(i+1,l)*b(l,j)
              s22 = s22 + a(i+1,l)*b(l,j+1)
              s23 = s23 + a(i+1,l)*b(l,j+2)

              s31 = s31 + a(i+2,l)*b(l,j)
              s32 = s32 + a(i+2,l)*b(l,j+1)
              s33 = s33 + a(i+2,l)*b(l,j+2)

              s41 = s41 + a(i+3,l)*b(l,j)
              s42 = s42 + a(i+3,l)*b(l,j+1)
              s43 = s43 + a(i+3,l)*b(l,j+2)
            enddo
            c(i,j)     = s11 
            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 
            c(i,j+1)   = s12 
            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42
            c(i,j+2)   = s13
            c(i+1,j+2) = s23
            c(i+2,j+2) = s33
            c(i+3,j+2) = s43
          endif
        endif
      enddo

* Residual when m is not multiple of 4
      if (mresid .eq. 0) then
        return
      elseif (mresid .eq. 1) then
        do j=1,n-nresid,4
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          s14 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,j)
            s12 = s12 + a(m,l)*b(l,j+1)
            s13 = s13 + a(m,l)*b(l,j+2)
            s14 = s14 + a(m,l)*b(l,j+3)
          enddo
          c(m,j)     = s11 
          c(m,j+1)   = s12 
          c(m,j+2)   = s13
          c(m,j+3)   = s14
        enddo
* mresid is 1, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n)
          enddo
          c(m,n) = s11
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s12 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n-1)
            s12 = s12 + a(m,l)*b(l,n)
          enddo
          c(m,n-1) = s11
          c(m,n) = s12
          return
        else
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n-2)
            s12 = s12 + a(m,l)*b(l,n-1)
            s13 = s13 + a(m,l)*b(l,n)
          enddo
          c(m,n-2) = s11
          c(m,n-1) = s12
          c(m,n) = s13
          return
        endif          
      elseif (mresid .eq. 2) then
        do j=1,n-nresid,4
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          s14 = 0.0d0
          s21 = 0.0d0
          s22 = 0.0d0
          s23 = 0.0d0
          s24 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,j)
            s12 = s12 + a(m-1,l)*b(l,j+1)
            s13 = s13 + a(m-1,l)*b(l,j+2)
            s14 = s14 + a(m-1,l)*b(l,j+3)

            s21 = s21 + a(m,l)*b(l,j)
            s22 = s22 + a(m,l)*b(l,j+1)
            s23 = s23 + a(m,l)*b(l,j+2)
            s24 = s24 + a(m,l)*b(l,j+3)
          enddo
          c(m-1,j)   = s11 
          c(m-1,j+1) = s12 
          c(m-1,j+2) = s13
          c(m-1,j+3) = s14
          c(m,j)     = s21
          c(m,j+1)   = s22 
          c(m,j+2)   = s23
          c(m,j+3)   = s24
        enddo
* mresid is 2, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          s21 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n)
          enddo
          c(m-1,n) = s11
          c(m,n) = s21
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s21 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n-1)
            s12 = s12 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n-1)
            s22 = s22 + a(m,l)*b(l,n)
          enddo
          c(m-1,n-1) = s11
          c(m-1,n)   = s12
          c(m,n-1)   = s21
          c(m,n)     = s22
          return
        else
          s11 = 0.0d0
          s21 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n-2)
            s12 = s12 + a(m-1,l)*b(l,n-1)
            s13 = s13 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n-2)
            s22 = s22 + a(m,l)*b(l,n-1)
            s23 = s23 + a(m,l)*b(l,n)
          enddo
          c(m-1,n-2) = s11
          c(m-1,n-1) = s12
          c(m-1,n)   = s13
          c(m,n-2)   = s21
          c(m,n-1)   = s22
          c(m,n)     = s23
          return
        endif
      else
* mresid is 3
        do j=1,n-nresid,4
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0

          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0

          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0

          s14 = 0.0d0
          s24 = 0.0d0
          s34 = 0.0d0

          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,j)
            s12 = s12 + a(m-2,l)*b(l,j+1)
            s13 = s13 + a(m-2,l)*b(l,j+2)
            s14 = s14 + a(m-2,l)*b(l,j+3)

            s21 = s21 + a(m-1,l)*b(l,j)
            s22 = s22 + a(m-1,l)*b(l,j+1)
            s23 = s23 + a(m-1,l)*b(l,j+2)
            s24 = s24 + a(m-1,l)*b(l,j+3)

            s31 = s31 + a(m,l)*b(l,j)
            s32 = s32 + a(m,l)*b(l,j+1)
            s33 = s33 + a(m,l)*b(l,j+2)
            s34 = s34 + a(m,l)*b(l,j+3)
          enddo
          c(m-2,j)   = s11 
          c(m-2,j+1) = s12 
          c(m-2,j+2) = s13
          c(m-2,j+3) = s14

          c(m-1,j)   = s21 
          c(m-1,j+1) = s22
          c(m-1,j+2) = s23
          c(m-1,j+3) = s24

          c(m,j)     = s31 
          c(m,j+1)   = s32
          c(m,j+2)   = s33
          c(m,j+3)   = s34
        enddo
* mresid is 3, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n)
          enddo
          c(m-2,n) = s11
          c(m-1,n) = s21
          c(m,n) = s31
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n-1)
            s12 = s12 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n-1)
            s22 = s22 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n-1)
            s32 = s32 + a(m,l)*b(l,n)
          enddo
          c(m-2,n-1) = s11
          c(m-2,n)   = s12
          c(m-1,n-1) = s21
          c(m-1,n)   = s22
          c(m,n-1)   = s31
          c(m,n)     = s32
          return
        else
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n-2)
            s12 = s12 + a(m-2,l)*b(l,n-1)
            s13 = s13 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n-2)
            s22 = s22 + a(m-1,l)*b(l,n-1)
            s23 = s23 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n-2)
            s32 = s32 + a(m,l)*b(l,n-1)
            s33 = s33 + a(m,l)*b(l,n)
          enddo
          c(m-2,n-2) = s11
          c(m-2,n-1) = s12
          c(m-2,n)   = s13
          c(m-1,n-2) = s21
          c(m-1,n-1) = s22
          c(m-1,n)   = s23
          c(m,n-2)   = s31
          c(m,n-1)   = s32
          c(m,n)     = s33
          return
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
