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
      integer :: ppiclf_moveparticle
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

      real*8 factor, rmass_add

      real*8 gkern

! Needed for Added mass calculation
      integer*4 j, l
      real*8 SDrho

      real*8 vgradrhog
      integer*4 i, n, ic, k
      integer*4 store_forces

! Needed for heat transfer
      real*8 qq, rmass_therm, temp, qq_du

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
      integer*4 nstage, iStage
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
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >       + 1) - 1), (ppiclf_npart+ppiclf_npart_gp))
      INTEGER*4  SBin_counter( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >       + 1) - 1))
      INTEGER*4 i_Bin(3), n_SBin(3), tot_SBin

      INTEGER*4 nsubbin_size
      INTEGER*4 nbin_total
! 
!-----------------------------------------------------------------------
!
!
! Code:
!

      icallb = icallb + 1
      nstage = 3
      iStage = mod(icallb,nstage)
      if (iStage .eq. 0) iStage = 3  

      ! Count every iStage=1 for debug output
      if (iStage .eq. 1) idebug = idebug + 1

      ! Print dt and time every time step
      if (ppiclf_nid==0) then
      !if (iStage .eq. 1) then
        write(6,'(a,2x,i4,2x,2(1pe14.6),2x,i3)') 
     >      '*** PPICLF dt, time = ',
     >      iStage,ppiclf_dt,ppiclf_time
      !endif
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

      fac = ppiclf_rk3dtrk(iStage)*ppiclf_dt
      if (1==2) then
         if (ppiclf_nid==0) print*,'dt,fac=',
     >      iStage,ppiclf_dt,fac,
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
!     do i=1,ppiclf_npart
!        ppiclf_rprop(PPICLF_R_JDPe,i) = 
!     > (0.00005/ppiclf_rprop(PPICLF_R_JSPT,i))
!     > * ppiclf_rprop(PPICLF_R_JDP,i)  
        !if (ppiclf_npart .gt. 0) then
        !if ((i .eq. 1) .or. (i .eq. ppiclf_npart)) then
        !  write(*,*) "i,JSPT",i,ppiclf_rprop(PPICLF_R_JSPT,i)       
        !endif
        !endif
!      end do 
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

!
!-----------------------------------------------------------------------
! Avery added 10/10/2024 - Map particles to subbins if collisional force, 
! Briney Added Mass force, or QS fluctation force is flagged
!
      !nearest neighbor search is used for col, am_flag 2, qs_fluct
      if (sbNearest_flag == 1) then

         if ((am_flag==2).or.(collisional_flag>=1)
     >          .or.(qs_fluct_flag>=1)) then

            call ppiclf_user_subbinMap(i_Bin, n_SBin, tot_SBin 
     >                               ,SBin_counter ,SBin_map)

         endif ! Collisions, QS Fluct, or Briney AM flags on
      
         ! Print out relevant information about subbin
         if (ppiclf_nid==0) then
         if (iStage==1) then

         nbin_total = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)
         nsubbin_size =
     >     (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >        + 1) *
     >     (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >        + 1) *
     >     (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >       + 1) 

!         if (ppiclf_time .EQ. 0.0) then
!         write(*,*) 'Subbin Method used!'
!         write(6,*) 'SUBBIN ', 
!     >     ppiclf_time,
!     >     ppiclf_bins_dx(1:3),
!     >     nsubbin_size,
!     >     tot_SBin,n_SBin(1:3),
!     >     ppiclf_npart,ppiclf_npart_gp,
!     >     nsubbin_size*(ppiclf_npart+ppiclf_npart_gp),
!     >     ' GB: ',nsubbin_size*
!     >             (ppiclf_npart+ppiclf_npart_gp)*4/1e9
!         write(6,*) 'Viscous Unsteady',
!     >     ppiclf_nUnsteadyData,ppiclf_nTimeBH,
!     >     ppiclf_dt
!
!         endif ! end ppiclf_time = 0

         endif ! end iStage = 1
         endif ! end ppiclf_nid = 0

      endif ! end sbNearest_flag = 1

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
         rem    = (1.0d0-rphip)*rep   ! mean slip Reynolds Number
         asndf  = ppiclf_rprop(PPICLF_R_JCS,i)
         rmachp = vmag/asndf
         rhop   = ppiclf_rprop(PPICLF_R_JRHOP,i)

         ! TLJ - 04/03/2025; Do not calculate forces if vmag = 0
         !       Otherwise the particles might move before the 
         !       shock arrives
!         if (vmag <= 1.d-8) cycle
         
         ! 08/08/2025 - Thierry  - 1.d-8 is very small
         if (vmag <= 1.d-3) cycle

         ! Thierry - at initial times when rmachp and vmag are very
         ! small, we would get very large CD (>100)
         ! This leads to NaN projected force when PseudoTurbulence is 
         ! enabled. One fix I'm trying is to cycle when rmachp and vmag
         ! are small

         ! Thierry - seeing if that fixes the issue of very large CD
         ! initially 
         !if(vmag .lt. 1.0 .or. rmachp .lt. 1.d-3) cycle

         ! TLJ - redefined rprop(PPICLF_R_JSPT,i) to be the particle
         !   velocity magnitude for plotting purposes - 01/03/2025
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

         ! TLJ: Needed for viscous unsteady force
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
         qq=0.0d0; qq_du=0.0d0
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
               gkern = sqrt(rpi*ppiclf_filter**2/
     >                (4.0d0*log(2.0d0)))**(-ppiclf_ndim)
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
         end if ! qs_fluct_flag .and. vmag

         ! add neighbors
         IF ( sbNearest_flag .EQ. 1) THEN
            CALL ppiclf_solve_NearestNeighborSB(
     >           i,tot_SBin,SBin_counter,SBin_map,n_SBin,i_Bin)
         ELSE
             CALL ppiclf_solve_NearestNeighbor(i)
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
     >                                      xi_par,xi_perp,xi_T,
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
         ppiclf_rprop(PPICLF_R_XIT,i) = xi_T

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
! Step 7: Viscous-unsteady force with history kernel
!         Diffusive-unsteady heat transfer is also calculated in this call
         ! 1 == original Trapezoidal code for uniform dt
         ! 2 == original Trapezoidal code for variable dt
         ! 3 == modified Trapezoidal method for variable dt
         ! 4 == Hinsberg method for variable dt
         if (ViscousUnsteady_flag==1) then
            call ppiclf_user_history_uniformtrap
     >           (i,iStage,fvux,fvuy,fvuz,qq_du)
         elseif (ViscousUnsteady_flag==2) then
            call ppiclf_user_history_variabletrap
     >           (i,iStage,fvux,fvuy,fvuz,qq_du)
         elseif (ViscousUnsteady_flag==3) then
            call ppiclf_user_history_modifiedtrap
     >           (i,iStage,fvux,fvuy,fvuz,qq_du)
         elseif (ViscousUnsteady_flag==4) then
            call ppiclf_user_history_hinsberg
     >           (i,iStage,fvux,fvuy,fvuz,qq_du)
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
! Step 8c : Diffusive Unsteady Heat Transfer model
!           Calculated in same subroutine as Viscous Unsteady
         if (HTUnsteady_flag==0) qq_du = 0.0d0
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

         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JZ, i) = ppiclf_y(PPICLF_JVZ,i)
         ppiclf_ydot(PPICLF_JVX,i) = (fqsx+famx+fdpdx+fvux+liftx+fcx)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JVY,i) = (fqsy+famy+fdpdy+fvuy+lifty+fcy)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JVZ,i) = (fqsz+famz+fdpdz+fvuz+liftz+fcz)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JT,i)  = (qq+qq_du)/rmass_therm
         ppiclf_ydot(PPICLF_JOX,i) = taux/rmass_omega
         ppiclf_ydot(PPICLF_JOY,i) = tauy/rmass_omega
         ppiclf_ydot(PPICLF_JOZ,i) = tauz/rmass_omega
         ppiclf_ydot(PPICLF_JMETAL,i)  = mdot_me
         ppiclf_ydot(PPICLF_JOXIDE,i)  = mdot_ox

!
! Update data for viscous unsteady case at every RK3 stage
!
         if (ViscousUnsteady_flag>=1) then
            call ppiclf_user_UpdatePlag(i,iStage)
         endif

!
! Step 11: Feed Back force to the gas phase
!
         ! Comment: ydotc represented the collisional force in the
         !    particle eqautions above. Here, we over-write the
         !    ydotc vectors for the feedback force used in Rocflu.
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

         if (feedback_flag==0) then
            ppiclf_ydotc(PPICLF_JVX,i) = 0.0d0 
            ppiclf_ydotc(PPICLF_JVY,i) = 0.0d0 
            ppiclf_ydotc(PPICLF_JVZ,i) = 0.0d0 
            ppiclf_ydotc(PPICLF_JT,i)  = 0.0d0
         endif

         if (feedback_flag==1) then
            ! Momentum equations feedback terms
            ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_rprop(PPICLF_R_JSPL,i) *
     >         (ppiclf_ydot(PPICLF_JVX,i)*rmass - fcx)
            ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_rprop(PPICLF_R_JSPL,i) *
     >         (ppiclf_ydot(PPICLF_JVY,i)*rmass - fcy)
            ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_rprop(PPICLF_R_JSPL,i) *
     >         (ppiclf_ydot(PPICLF_JVZ,i)*rmass - fcz)

            ! Energy equation feedback term
            if(pseudoTurb_flag==1) then
            ! 09/02/2025 -  Addition of PTKE to Rocflu's Energy Equation
              ppiclf_ydotc(PPICLF_JT,i) = ppiclf_rprop(PPICLF_R_JSPL,i)*
     >         ( (fqsx+fvux+famx+liftx)*ppiclf_y(PPICLF_JVX,i) + 
     >           (fqsy+fvuy+famy+lifty)*ppiclf_y(PPICLF_JVY,i) + 
     >           (fqsz+fvuz+famz+liftz)*ppiclf_y(PPICLF_JVZ,i) +
     >           qq )

            else ! Pseudo Turbulence is OFF
            ! 09/19/2025 - Thierry - Added Lift force
!            ppiclf_ydotc(PPICLF_JT,i) = ppiclf_rprop(PPICLF_R_JSPL,i) *
!     >         ( (fqsx+fvux+liftx)*ppiclf_y(PPICLF_JVX,i) + 
!     >           (fqsy+fvuy+lifty)*ppiclf_y(PPICLF_JVY,i) + 
!     >           (fqsz+fvuz+liftz)*ppiclf_y(PPICLF_JVZ,i) +
!     >                  famx*ppiclf_rprop(PPICLF_R_JUX,i) +
!     >                  famy*ppiclf_rprop(PPICLF_R_JUY,i) +
!     >                  famz*ppiclf_rprop(PPICLF_R_JUZ,i) +
!     >                  taux_hydro*ppiclf_y(PPICLF_JOX,i) +
!     >                  tauy_hydro*ppiclf_y(PPICLF_JOY,i) +
!     >                  tauz_hydro*ppiclf_y(PPICLF_JOZ,i) +
!     >           qq )

            ! Testing Equation 19 of the projection notes
            ! F'_j . u_g + F'_qs . (u_p - u_g)
!            ppiclf_ydotc(PPICLF_JT,i) = ppiclf_rprop(PPICLF_R_JSPL,i) *
!     >         ( (fqsx+famx+fvux+liftx)*ppiclf_rprop(PPICLF_R_JUX,i) + 
!     >           (fqsy+famy+fvuy+lifty)*ppiclf_rprop(PPICLF_R_JUY,i) + 
!     >           (fqsz+famz+fvuz+liftz)*ppiclf_rprop(PPICLF_R_JUZ,i) +
!     >           (fqsx*(-vx) + fqsy*(-vy) + fqsz*(-vz))              +
!     >           qq )

            ! Work projection without dissipative term for Case-12
            ppiclf_ydotc(PPICLF_JT,i) = ppiclf_rprop(PPICLF_R_JSPL,i) *
     >         ( (fqsx+famx+fvux+liftx)*ppiclf_rprop(PPICLF_R_JUX,i) + 
     >           (fqsy+famy+fvuy+lifty)*ppiclf_rprop(PPICLF_R_JUY,i) + 
     >           (fqsz+famz+fvuz+liftz)*ppiclf_rprop(PPICLF_R_JUZ,i) +
     >           qq )

            !ppiclf_ydotc(PPICLF_JT,i) = -1.0d0*ppiclf_ydotc(PPICLF_JT,i)
            endif ! pseudoTurb_flag

      ! 07/21/2025 - Thierry - Added Reynolds Subgrid Stress Tensor
      ! We need to store the tensor calculations here in a certain array and then
      ! access that in the mapping subroutine for when projection is needed

            ppiclf_rprop4(PPICLF_R_JRSG11,i) = Rsg(1,1)
            ppiclf_rprop4(PPICLF_R_JRSG12,i) = Rsg(1,2)
            ppiclf_rprop4(PPICLF_R_JRSG13,i) = Rsg(1,3)
            ppiclf_rprop4(PPICLF_R_JRSG21,i) = Rsg(2,1)
            ppiclf_rprop4(PPICLF_R_JRSG22,i) = Rsg(2,2)
            ppiclf_rprop4(PPICLF_R_JRSG23,i) = Rsg(2,3)
            ppiclf_rprop4(PPICLF_R_JRSG31,i) = Rsg(3,1)
            ppiclf_rprop4(PPICLF_R_JRSG32,i) = Rsg(3,2)
            ppiclf_rprop4(PPICLF_R_JRSG33,i) = Rsg(3,3) 
            
            ppiclf_rprop4(PPICLF_R_JTSG1,i) = T_par(1)
            ppiclf_rprop4(PPICLF_R_JTSG2,i) = T_par(2)
            ppiclf_rprop4(PPICLF_R_JTSG3,i) = T_par(3)
            
         endif ! feedback_flag

!
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
!
         ! Thierry - just testing this for now to see if it makes a difference in R_perp
         ! commenting out for now since it's not being outputted anywhere
         ! I need to figure out how to output ppiclf_force to vtu files
!         ppiclf_force(PPICLF_R_FQSX,i)  = fqsx
!         ppiclf_force(PPICLF_R_FQSY,i)  = fqsy
!         ppiclf_force(PPICLF_R_FQSZ,i)  = fqsz
!         ppiclf_force(PPICLF_R_FAMX,i)  = famx-rmass_add*ppiclf_ydot(PPICLF_JVX,i)
!         ppiclf_force(PPICLF_R_FAMY,i)  = famy-rmass_add*ppiclf_ydot(PPICLF_JVY,i)
!         ppiclf_force(PPICLF_R_FAMZ,i)  = famz-rmass_add*ppiclf_ydot(PPICLF_JVZ,i)
!         ppiclf_force(PPICLF_R_FAMBX,i) = FamBinary(1)
!         ppiclf_force(PPICLF_R_FAMBY,i) = FamBinary(2)
!         ppiclf_force(PPICLF_R_FAMBZ,i) = FamBinary(3)
!         ppiclf_force(PPICLF_R_FCX,i)   = fcx
!         ppiclf_force(PPICLF_R_FCY,i)   = fcy
!         ppiclf_force(PPICLF_R_FCZ,i)   = fcz
!         ppiclf_force(PPICLF_R_FVUX,i)  = fvux
!         ppiclf_force(PPICLF_R_FVUY,i)  = fvuy
!         ppiclf_force(PPICLF_R_FVUZ,i)  = fvuz
!         ppiclf_force(PPICLF_R_QQ,i)    = qq
!         ppiclf_force(PPICLF_R_FPGX,i)  = fdpdx
!         ppiclf_force(PPICLF_R_FPGY,i)  = fdpdy
!         ppiclf_force(PPICLF_R_FPGZ,i)  = fdpdz

!
! Step 14: If debug mode is ON, calculate and print the max values.
!          The user should not have this ON for production runs.
!
         if (ppiclf_debug .ge. 1) then
            if (sbNearest_flag.eq.1 .and. ppiclf_debug.eq.2) then
               write(7001,*) ppiclf_time, ppiclf_bins_dx(1:3),
     >            nsubbin_size, tot_SBin,n_SBin(1:3),
     >            ppiclf_npart, ppiclf_npart_gp,
     >            nsubbin_size*(ppiclf_npart+ppiclf_npart_gp),
     >            nsubbin_size*(ppiclf_npart+ppiclf_npart_gp)*4/1e9 
                  ! last entry in GB; assuming 4 bytes for integer*4
            endif
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
     >                rep,rphip,rphif,rmachp,rhop,rmu,rnu,rkappa
                  endif
                  if (i==ppiclf_npart) then
                     write(7011,*) i,ppiclf_time,rmass,vmag,rhof,dp,
     >                rep,rphip,rphif,rmachp,rhop,rmu,rnu,rkappa
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
      if (mod(idebug,1) .eq. 0) then
      !if (mod(idebug,10) .eq. 0) then
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
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >       + 1) - 1), (ppiclf_npart+ppiclf_npart_gp))
      INTEGER*4  SBin_counter( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
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
            n_SBin(l) = FLOOR((ppiclf_bins_dx(l)+2*ppiclf_d2chk(3))
     >                       /ppiclf_d2chk(3))
          ELSE
            n_SBin(l) = 0
          END IF
        END DO

        ! Determine total number of subbins
        tot_SBin = (n_SBin(1)+1)*(n_SBin(2)+1)
        IF (ppiclf_ndim .EQ. 3) THEN
          tot_SBin = tot_SBin*(n_SBin(3)+1)
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
     >            - ppiclf_d2chk(3)))/ppiclf_d2chk(3)) 
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
          IF (xp(1) .GT. (bin_xMin(1)-ppiclf_d2chk(3))
     >  .AND. xp(2) .GT. (bin_xMin(2)-ppiclf_d2chk(3))
     >  .AND. xp(3) .GT. (bin_xMin(3)-ppiclf_d2chk(3))
     >  .AND. xp(1) .LT. (bin_xMin(1)+ppiclf_bins_dx(1)+ppiclf_d2chk(3))
     >  .AND. xp(2) .LT. (bin_xMin(2)+ppiclf_bins_dx(2)+ppiclf_d2chk(3))
     >  .AND. xp(3) .LT. (bin_xMin(3)+ppiclf_bins_dx(3)+ppiclf_d2chk(3))
     >        ) THEN
            ! Determine subbin
            DO l = 1,3
              IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
                i_SBin(l) = FLOOR((xp(l) - (bin_xMin(l) 
     >          - ppiclf_d2chk(3)))/ppiclf_d2chk(3)) 
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
     >         rphif,asndf,rmachp,rhop,rnu,fac,              ! 11-18
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
     >         rphif,asndf,rmachp,rhop,rnu,fac,              ! 11-18
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
      
      ! ppiclf_d2chk(3) is neighbor width - user defined
      dr_max = ppiclf_d2chk(3)

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
      elseif (heattransfer_flag == 5) then
         call HTModel_Sun(i,Nuss)
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
! Created October 16, 2025
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: define Nusselt number Nu = Nu(Pr,Re,M)
!
!    B. Sun, S. Tenneti, S. Subramaniam (2015)
!        "Modeling average gas-solid heat transfer using
!        particle-resolved direct numerical simulation"
!     International Journal of Heat and Mass Transfer
!        
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_Sun(i,Nuss)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 Nuss, theta, q
!
! Code:
!     
      ! Model is provided for 
      ! mean slip Reynolds number (1-100)
      ! particle volume fraction (0.1-0.5)

      ! Need to validate that the model works and doesn't blow up 
      ! in our kind of simulations

      ! define Nusselt number Nu = Nu(Pr,Re,M)
      Nuss = (-0.46d0 + 1.77d0*rphif + 0.69d0*rphif**2)*(rphif**(-3))
     > +(1.37d0 - 2.4d0*rphif + 1.2*rphif**2)*(rep**0.7)*(rpr**OneThird)

      theta = 1.0d0 - 1.6d0*rphip*(1.0d0 - rphip) 
     >        - (3.0d0*rphip*(1.0d0 - rphip)**4)*exp(-(rep**0.4)*rphip)

      ! Sun's formulation Eq. (29) of the paprt is per volume of grid cell
      ! This is transformed to per particle by multiplying by 
      ! (pi dp^3)/(6 phip) 
      ! This is done below in the Nusslet number with simplifications
      Nuss = rpi/4.0d0 * Nuss/theta


      return
      end
!-----------------------------------------------------------------------
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

      factor = 0.125d0*rpi*dp*dp*dp*rhof

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
! Created Feb.  1, 2024
! Modified Dec. 1, 2025     
!
! Subroutine for viscous unsteady force with history kernel
! Updated to includes diffusive-unsteady Heat Transfer
!
! Mei-Adrian history kernel for viscous-unsteady
!
! Original code copied from either files in rocintereact/
!   INRT_CalcDragUnsteady_AMImplicit.F90
!   INRT_CalcDragUnsteady_AMExplicit.F90
!
! The number of time steps kept for the history
!   kernel is set in libpicl/ppiclF/source/PPICLF_STD.h
!
!
! Below we include several schemes for the integration
!   of this viscous and diffusive history integrals
!   1. Original trapezoidal method for uniform time steps
!   2. Original trapezoidal method for variable time steps
!   3. Modified trapezoidal method for variable time steps
!   4. Hinsberg method for variable time steps
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Trapezoidal Method with uniform time steps only
!     Caution - Only use this if you manually set 
!               global%dtMin = constant inside 
!               rocflu/RFLU_TimeStepping.F90
!               see line 507
!
      subroutine ppiclf_user_history_uniformtrap
     >           (i,iStage,fvux,fvuy,fvuz,qq_du)
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
      real*8 factor_du,b_du,Pe,kernel_du,qq_du,alpha_fluid
      real*8 dti, theta_sangani
      integer*4 i1,i2

      i1 = 2
      i2 = 78

!
! Code:
!
      fvux = 0.0d0
      fvuy = 0.0d0
      fvuz = 0.0d0
      iT   = 1
      time = 0.0d0

      ! Thermal diffusivity and Peclet number
      alpha_fluid = rkappa/(rhof * rcp_fluid)
      Pe = rep * rpr 

      ! Viscous-unsteady
      !   Sangani's volume fraction correction for dilute random arrays
      !   Capping volume fraction at 0.5 
      theta_sangani = 1.0d0 + 2.28d0*min(rphip,0.5d0)
      factor = 3.0d0*rpi*rnu*dp*theta_sangani
      fH     = 0.75d0 + .1055d0*rep

      ! Diffusive-unsteady Heat Transfer
      factor_du = 2.0d0*(dp**2)*sqrt(rpi*rkappa*rhof*rcp_fluid)
      b_du = (1.63d0 - 0.92d0*erf(0.017*(Pe - 80.0d0)))
     >       *(1.0d0 - 0.4d0*exp(-Pe/16.0d0)-0.6d0*exp(-Pe**2/30.0d0))

      ! For uniform time step
      dti = ppiclf_dt

      if (ppiclf_nTimeBH > 1) then
         do iT = 2,ppiclf_nTimeBH-1
            time = ppiclf_timeBH(iT)

            ! Viscous-unsteady
            A  = (4.0d0*rpi*time*rnu/(dp**2))**(.25d0)
            B  = (rpi*(vmag**3)*(time**2)/ 
     >                 (dp*rnu*(fH**3)))**(.5d0)

            kernelVU = factor*(A+B)**(-2)

            fvux = fvux + dti * kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JX,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JX,iT,i) )
            fvuy = fvuy + dti * kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JY,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JY,iT,i) )
            fvuz = fvuz + dti * kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JZ,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JZ,iT,i) )
         
            ! Diffusive-unsteady
            kernel_du = exp(-b_du*vmag*time/dp)/sqrt(time)
            qq_du = qq_du + dti * factor_du*kernel_du*
     >                (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i))


!           if(ppiclf_iprop(5,i) .eq. 7 .and.
!    >         ppiclf_iprop(7,i) .eq. 1724) then
!              !open(unit=66,file='fort.66',position='append')   

!              write(66,*) ppiclf_time, iT, time,
!    >            factor_du, kernel_du,
!    >            ppiclf_dTdtMixt(iT,i),
!    >            ppiclf_dTdtPlag(iT,i),
!    >            (ppiclf_dTdtMixt(iT,i) - ppiclf_dTdtPlag(iT,i)),
!    >            factor_du*kernel_du*
!    >            (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i)),
!    >            qq_du,fvux,fvuy,fvuz,
!    >            kernelVU, A, B, factor 

!              flush(66)
!           endif
         enddo

         iT = ppiclf_nTimeBH
         time = ppiclf_timeBH(iT)

         ! Viscous-unsteady
         A  = (4.0d0*rpi*time*rnu/(dp**2))**(.25d0)
         B  = (rpi*(vmag**3)*(time**2)/ 
     >                 (dp*rnu*(fH**3)))**(.5d0)
         kernelVU = 0.5d0*factor*(A+B)**(-2)
         fvux = fvux + dti * kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JX,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JX,iT,i) )
         fvuy = fvuy + dti * kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JY,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JY,iT,i) )
         fvuz = fvuz + dti * kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JZ,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JZ,iT,i) )

         ! Diffusive-unsteady
         kernel_du = 0.5d0*exp(-b_du*vmag*time/dp)/sqrt(time)
         qq_du = qq_du + dti * factor_du*kernel_du*
     >             (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i))


!        if(ppiclf_iprop(5,i) .eq. 7 .and.
!    >      ppiclf_iprop(7,i) .eq. 1724) then
!           !open(unit=66,file='fort.66',position='append')   

!           write(66,*) ppiclf_time, iT, time,
!    >            factor_du, kernel_du,
!    >            ppiclf_dTdtMixt(iT,i),
!    >            ppiclf_dTdtPlag(iT,i),
!    >            (ppiclf_dTdtMixt(iT,i) - ppiclf_dTdtPlag(iT,i)),
!    >            factor_du*kernel_du*
!    >            (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i)),
!    >            qq_du,fvux,fvuy,fvuz,
!    >            kernelVU, A, B, factor 

!           flush(66)
!        endif
      endif

      if (ppiclf_iprop(5,i) .eq. i1 .and.
     >    ppiclf_iprop(7,i) .eq. i2) then
          if (istage == 3) then
             write(66,*) ppiclf_time,qq_du,fvux,fvuy,fvuz
          endif
      endif

      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Trapezoidal Method with non-uniform time steps
!
      subroutine ppiclf_user_history_variabletrap
     >           (i,iStage,fvux,fvuy,fvuz,qq_du)
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
      real*8 factor_du,b_du,Pe,kernel_du,qq_du,alpha_fluid
      real*8 dti, theta_sangani
      integer*4 i1,i2

      i1 = 2
      i2 = 78
!
! Code:
!
      fvux = 0.0d0
      fvuy = 0.0d0
      fvuz = 0.0d0
      iT   = 1
      time = 0.0d0

      qq_du = 0.0d0

      ! Thermal diffusivity and Peclet number
      alpha_fluid = rkappa/(rhof * rcp_fluid)
      Pe = rep * rpr

      ! Viscous-unsteady
      !   Sangani's volume fraction correction for dilute random arrays
      !   Capping volume fraction at 0.5 
      theta_sangani = 1.0d0 + 2.28d0*min(rphip,0.5d0)
      factor = 3.0d0*rpi*rnu*dp*theta_sangani
      fH     = 0.75d0 + .1055d0*rep

      ! Diffusive-unsteady Heat Transfer
      factor_du = 2.0d0*(dp**2)*sqrt(rpi*rkappa*rhof*rcp_fluid)
      b_du = (1.63d0 - 0.92d0*erf(0.017*(Pe - 80.0d0)))
     >    *(1.0d0 - 0.4d0*exp(-Pe/16.0d0)-0.6d0*exp(-Pe**2/30.0d0))

      if (ppiclf_nTimeBH > 1) then
         do iT = 3,ppiclf_nTimeBH
            dti = ppiclf_timeBH(iT)-ppiclf_timeBH(iT-1)

            time = ppiclf_timeBH(iT)
            ! Viscous-unsteady
            A  = (4.0d0*rpi*time*rnu/(dp**2))**(.25d0)
            B  = (rpi*(vmag**3)*(time**2)/ 
     >                 (dp*rnu*(fH**3)))**(.5d0)
            kernelVU = 0.5d0 * factor*(A+B)**(-2)
            fvux = fvux + dti * kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JX,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JX,iT,i) )
            fvuy = fvuy + dti * kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JY,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JY,iT,i) )
            fvuz = fvuz + dti * kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JZ,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JZ,iT,i) )
            ! Diffusive-unsteady
            kernel_du = 0.5d0 * exp(-b_du*vmag*time/dp)/sqrt(time)
            qq_du = qq_du + dti * factor_du*kernel_du*
     >                (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i))

            time = ppiclf_timeBH(iT-1)
            ! Viscous-unsteady
            A  = (4.0d0*rpi*time*rnu/(dp**2))**(.25d0)
            B  = (rpi*(vmag**3)*(time**2)/ 
     >                 (dp*rnu*(fH**3)))**(.5d0)
            kernelVU = 0.5d0 * factor*(A+B)**(-2)
            fvux = fvux + dti * kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JX,iT-1,i) -
     >                     ppiclf_drudtPlag(PPICLF_JX,iT-1,i) )
            fvuy = fvuy + dti * kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JY,iT-1,i) -
     >                     ppiclf_drudtPlag(PPICLF_JY,iT-1,i) )
            fvuz = fvuz + dti * kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JZ,iT-1,i) -
     >                     ppiclf_drudtPlag(PPICLF_JZ,iT-1,i) )
            ! Diffusive-unsteady
            kernel_du = 0.5d0 * exp(-b_du*vmag*time/dp)/sqrt(time)
            qq_du = qq_du + dti * factor_du*kernel_du*
     >                (ppiclf_dTdtMixt(iT-1,i) -ppiclf_dTdtPlag(iT-1,i))


!           if(ppiclf_iprop(5,i) .eq. i1 .and.
!    >         ppiclf_iprop(7,i) .eq. i2) then
!           if (iStage==3) then
!             !stop
!              !open(unit=66,file='fort.66',position='append')   
!
!              time = ppiclf_timeBH(iT)
!              write(66,*) ppiclf_time, iT, time,
!    >            factor_du, kernel_du,
!    >            ppiclf_dTdtMixt(iT,i),
!    >            ppiclf_dTdtPlag(iT,i),
!    >            (ppiclf_dTdtMixt(iT,i) - ppiclf_dTdtPlag(iT,i)),
!    >            factor_du*kernel_du*
!    >            (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i)),
!    >            qq_du,fvux,fvuy,fvuz,
!    >            kernelVU, A, B, factor 
!
!              flush(66)
!           endif
!           endif
         enddo

         iT = 2
         time = ppiclf_timeBH(2)
         dti = ppiclf_timeBH(2)-ppiclf_timeBH(1)
         ! Viscous-unsteady
         A  = (4.0d0*rpi*time*rnu/(dp**2))**(.25d0)
         B  = (rpi*(vmag**3)*(time**2)/ 
     >                 (dp*rnu*(fH**3)))**(.5d0)
         kernelVU = 0.5d0*factor*(A+B)**(-2)
         fvux = fvux + dti * kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JX,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JX,iT,i) )
         fvuy = fvuy + dti * kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JY,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JY,iT,i) )
         fvuz = fvuz + dti * kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JZ,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JZ,iT,i) )
         ! Diffusive-unsteady
         kernel_du = 0.5d0*exp(-b_du*vmag*time/dp)/sqrt(time)
         qq_du = qq_du + dti * factor_du*kernel_du*
     >             (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i))


!        if(ppiclf_iprop(5,i) .eq. i1 .and.
!    >      ppiclf_iprop(7,i) .eq. i2) then
!        if (istage == 3) then
!           !open(unit=66,file='fort.66',position='append')   

!           iT = 2
!           time = ppiclf_timeBH(2)
!           write(66,*) ppiclf_time, iT, time,
!    >            factor_du, kernel_du,
!    >            ppiclf_dTdtMixt(iT,i),
!    >            ppiclf_dTdtPlag(iT,i),
!    >            (ppiclf_dTdtMixt(iT,i) - ppiclf_dTdtPlag(iT,i)),
!    >            factor_du*kernel_du*
!    >            (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i)),
!    >            qq_du,fvux,fvuy,fvuz,
!    >            kernelVU, A, B, factor 

!           flush(66)
!        endif
!        endif
      endif

      if (ppiclf_iprop(5,i) .eq. i1 .and.
     >    ppiclf_iprop(7,i) .eq. i2) then
          if (istage == 3) then
             write(66,*) ppiclf_time,qq_du,fvux,fvuy,fvuz
          endif
      endif


      return
      end
!
!-----------------------------------------------------------------------
!
!     Modified Trapezodial Method with non-uniform time steps
!       This subroutines factors out 1/sqrt(s) integrable singularity
!
      subroutine ppiclf_user_history_modifiedtrap
     >           (i,iStage,fvux,fvuy,fvuz,qq_du)
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
      real*8 factor_du,b_du,Pe,kernel_du,qq_du,alpha_fluid
      real*8 s1,s2,dsi
      real*8 sx1,sx2,sy1,sy2,sz1,sz2,g1,g2
      real*8 theta_sangani
      integer*4 i1,i2

      i1 = 2
      i2 = 78

!
! Code:
!
      fvux = 0.0d0
      fvuy = 0.0d0
      fvuz = 0.0d0
      iT   = 1
      time = 0.0d0

      qq_du = 0.0d0

      ! Thermal diffusivity and Peclet number
      alpha_fluid = rkappa/(rhof * rcp_fluid)
      Pe = rep * rpr

      ! Viscous-unsteady
      !   Sangani's volume fraction correction for dilute random arrays
      !   Capping volume fraction at 0.5 
      theta_sangani = 1.0d0 + 2.28d0*min(rphip,0.5d0)
      factor = 3.0d0*rpi*rnu*dp*theta_sangani
      fH     = 0.75d0 + .1055d0*rep

      ! Diffusive-unsteady Heat Transfer
      factor_du = 2.0d0*(dp**2)*sqrt(rpi*rkappa*rhof*rcp_fluid)
      b_du = (1.63d0 - 0.92d0*erf(0.017*(Pe - 80.0d0)))
     >       *(1.0d0 - 0.4d0*exp(-Pe/16.0d0)-0.6d0*exp(-Pe**2/30.0d0))

      if (ppiclf_nTimeBH > 1) then
         ! Do interval [t_i,t_i-1],  i=2,...,N-1
         do iT = 2,ppiclf_nTimeBH
            s1 = ppiclf_timeBH(iT-1)
            s2 = ppiclf_timeBH(iT)
            dsi = (s2 - s1)/(sqrt(s1)+sqrt(s2))

            ! Viscous-unsteady
            ! at s1 = s(i-1)
            A  = (4.0d0*rpi*rnu/(dp**2))**(.25d0)
            B  = (rpi*(vmag**3)/(dp*rnu*(fH**3)))**(.5d0)
     >           * (s1**0.75d0)
            kernelVU = factor*(A+B)**(-2)
            sx1 = kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JX,iT-1,i) -
     >                     ppiclf_drudtPlag(PPICLF_JX,iT-1,i) )
            sy1 = kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JY,iT-1,i) -
     >                     ppiclf_drudtPlag(PPICLF_JY,iT-1,i) )
            sz1 = kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JZ,iT-1,i) -
     >                     ppiclf_drudtPlag(PPICLF_JZ,iT-1,i) )

            ! at s2 = s(i)
            A  = (4.0d0*rpi*rnu/(dp**2))**(.25d0)
            B  = (rpi*(vmag**3)/(dp*rnu*(fH**3)))**(.5d0)
     >           * (s2**0.75d0)
            kernelVU = factor*(A+B)**(-2)
            sx2 = kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JX,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JX,iT,i) )
            sy2 = kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JY,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JY,iT,i) )
            sz2 = kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JZ,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JZ,iT,i) )

            fvux = fvux + (sx1 + sx2)*dsi
            fvuy = fvuy + (sy1 + sy2)*dsi
            fvuz = fvuz + (sz1 + sz2)*dsi

            ! Diffusive-unsteady Heat Transfer
            ! at s1 = s(i-1)
            kernel_du = exp(-b_du*vmag*s1/dp)
            g1 = factor_du*kernel_du*
     >                (ppiclf_dTdtMixt(iT-1,i) -ppiclf_dTdtPlag(iT-1,i))
            ! at s2 = s(i)
            kernel_du = exp(-b_du*vmag*s2/dp)
            g2 = factor_du*kernel_du*
     >                (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i))

            qq_du = qq_du + (g1+g2)*dsi


!           if(ppiclf_iprop(5,i) .eq. 7 .and.
!    >           ppiclf_iprop(7,i) .eq. 1724) then
!              !open(unit=66,file='fort.66',position='append')   

!              time = ppiclf_timeBH(iT)
!              write(66,*) ppiclf_time, iT, time,
!    >            factor_du, kernel_du,
!    >            ppiclf_dTdtMixt(iT,i),
!    >            ppiclf_dTdtPlag(iT,i),
!    >            (ppiclf_dTdtMixt(iT,i) - ppiclf_dTdtPlag(iT,i)),
!    >            factor_du*kernel_du*
!    >            (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i)),
!    >            qq_du,fvux,fvuy,fvuz,
!    >            kernelVU, A, B, factor 

!              flush(66)
!           endif
         enddo

      endif

      if (ppiclf_iprop(5,i) .eq. i1 .and.
     >    ppiclf_iprop(7,i) .eq. i2) then
          if (istage == 3) then
             write(66,*) ppiclf_time,qq_du,fvux,fvuy,fvuz
          endif
      endif

      return
      end
!
!-----------------------------------------------------------------------
!
!     Hinsberg Method with non-uniform time steps
!
      subroutine ppiclf_user_history_hinsberg
     >           (i,iStage,fvux,fvuy,fvuz,qq_du)
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
      real*8 factor_du,b_du,Pe,kernel_du,qq_du,alpha_fluid
      real*8 TwoThirds
      real*8 s1,s2,dsinv,dsi1,dsi2,dsi3,dsi4
      real*8 sx1,sx2,sy1,sy2,sz1,sz2
      real*8 gx,gy,gz,g1,g2,g
      real*8 theta_sangani
      integer*4 i1,i2

      i1 = 2
      i2 = 78

!
! Code:
!
      fvux = 0.0d0
      fvuy = 0.0d0
      fvuz = 0.0d0
      iT   = 1
      time = 0.0d0
      TwoThirds = 2.0d0/3.0d0

      qq_du = 0.0d0

      ! Thermal diffusivity and Peclet number
      alpha_fluid = rkappa/(rhof * rcp_fluid)
      Pe = rep * rpr

      ! Viscous-unsteady
      !   Sangani's volume fraction correction for dilute random arrays
      !   Capping volume fraction at 0.5 
      theta_sangani = 1.0d0 + 2.28d0*min(rphip,0.5d0)
      factor = 3.0d0*rpi*rnu*dp*theta_sangani
      fH     = 0.75d0 + .1055d0*rep

      ! Diffusive-unsteady Heat Transfer
      factor_du = 2.0d0*(dp**2)*sqrt(rpi*rkappa*rhof*rcp_fluid)
      b_du = (1.63d0 - 0.92d0*erf(0.017*(Pe - 80.0d0)))
     >       *(1.0d0 - 0.4d0*exp(-Pe/16.0d0)-0.6d0*exp(-Pe**2/30.0d0))

      if (ppiclf_nTimeBH > 1) then
         do iT = 3,ppiclf_nTimeBH
            ! Do interval [t_i,t_i-1],  i=2,...,N-1
            s1 = ppiclf_timeBH(iT-1)
            s2 = ppiclf_timeBH(iT)
            dsinv = 1.0d0/(sqrt(s1)+sqrt(s2))
            dsi1 = s1*dsinv
            dsi2 = s2*dsinv
            dsi3 = (s2+sqrt(s1*s2)+s1)*dsinv

            ! Viscous-unsteady
            A  = (4.0d0*rpi*rnu/(dp**2))**(.25d0)
            B  = (rpi*(vmag**3)/(dp*rnu*(fH**3)))**(.5d0)
     >           * (s1**0.75d0)
            kernelVU = factor*(A+B)**(-2)
            sx1 = kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JX,iT-1,i) -
     >                     ppiclf_drudtPlag(PPICLF_JX,iT-1,i) )
            sy1 = kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JY,iT-1,i) -
     >                     ppiclf_drudtPlag(PPICLF_JY,iT-1,i) )
            sz1 = kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JZ,iT-1,i) -
     >                     ppiclf_drudtPlag(PPICLF_JZ,iT-1,i) )

            A  = (4.0d0*rpi*rnu/(dp**2))**(.25d0)
            B  = (rpi*(vmag**3)/(dp*rnu*(fH**3)))**(.5d0)
     >           * (s2**0.75d0)
            kernelVU = factor*(A+B)**(-2)
            sx2 = kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JX,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JX,iT,i) )
            sy2 = kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JY,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JY,iT,i) )
            sz2 = kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JZ,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JZ,iT,i) )

            gx = 2.0d0*(sx1*dsi2-sx2*dsi1)+TwoThirds*(sx2-sx1)*dsi3
            gy = 2.0d0*(sy1*dsi2-sy2*dsi1)+TwoThirds*(sx2-sx1)*dsi3
            gz = 2.0d0*(sz1*dsi2-sz2*dsi1)+TwoThirds*(sx2-sx1)*dsi3

            fvux = fvux + gx
            fvuy = fvuy + gy
            fvuz = fvuy + gz
         
            ! Diffusive-unsteady
            kernel_du = exp(-b_du*vmag*s1/dp)
            g1 = factor_du*kernel_du*
     >                (ppiclf_dTdtMixt(iT-1,i) -ppiclf_dTdtPlag(iT-1,i))
            kernel_du = exp(-b_du*vmag*s2/dp)
            g2 = factor_du*kernel_du*
     >                (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i))

            g = 2.0d0*(g1*dsi2-g2*dsi1)+TwoThirds*(g2-g1)*dsi3
            qq_du = qq_du + g

!           if(ppiclf_iprop(5,i) .eq. 7 .and.
!    >         ppiclf_iprop(7,i) .eq. 1724) then
!              !open(unit=66,file='fort.66',position='append')   

!              write(66,*) ppiclf_time, iT, time,
!    >            factor_du, kernel_du,
!    >            ppiclf_dTdtMixt(iT,i),
!    >            ppiclf_dTdtPlag(iT,i),
!    >            (ppiclf_dTdtMixt(iT,i) - ppiclf_dTdtPlag(iT,i)),
!    >            factor_du*kernel_du*
!    >            (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i)),
!    >            qq_du,fvux,fvuy,fvuz,
!    >            kernelVU, A, B, factor 

!              flush(66)
!           endif
         enddo

         ! Do interval [t_{N-1},t_{N-1}+dnk(iStage)]  
         !s1 = ppiclf_timeBH(1)
         s1 = ppiclf_timeBH(1)*(1.0d0+ppiclf_rk3dtrk(iStage))
         s2 = ppiclf_timeBH(2)
         dsinv = 1.0d0/(sqrt(s1)+sqrt(s2))
         dsi1 = s1*dsinv
         dsi2 = s2*dsinv
         dsi3 = (s2+sqrt(s1*s2)+s1)*dsinv
         dsi4 = (s2-s1)*dsinv

         ! Viscous-unsteady
         iT = 2
         A  = (4.0d0*rpi*rnu/(dp**2))**(.25d0)
         B  = (rpi*(vmag**3)/(dp*rnu*(fH**3)))**(.5d0)
     >        * (s2**0.75d0)
         kernelVU = factor*(A+B)**(-2)
         sx2 = kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JX,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JX,iT,i) )
         sy2 = kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JY,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JY,iT,i) )
         sz2 = kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JZ,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JZ,iT,i) )
         
         A  = (4.0d0*rpi*rnu/(dp**2))**(.25d0)
         B  = (rpi*(vmag**3)/(dp*rnu*(fH**3)))**(.5d0)
     >        * (s1**0.75d0)
         kernelVU = factor*(A+B)**(-2)
         sx1 = kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JX,iT-1,i) -
     >                  ppiclf_drudtPlag(PPICLF_JX,iT-1,i) )
         sy1 = kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JY,iT-1,i) -
     >                  ppiclf_drudtPlag(PPICLF_JY,iT-1,i) )
         sz1 = kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JZ,iT-1,i) -
     >                  ppiclf_drudtPlag(PPICLF_JZ,iT-1,i) )

         gx = (sx2+sx1)*dsi4
         gy = (sy2+sy1)*dsi4
         gz = (sz2+sz1)*dsi4

         fvux = fvux + gx
         fvuy = fvuy + gy
         fvuz = fvuy + gz
         
         ! Diffusive-unsteady Heat Transfer
         kernel_du = exp(-b_du*vmag*s2/dp)
         g2 = factor_du*kernel_du*
     >           (ppiclf_dTdtMixt(iT,i)-ppiclf_dTdtPlag(iT,i))

         kernel_du = exp(-b_du*vmag*s1/dp)
         g1 = factor_du*kernel_du*
     >           (ppiclf_dTdtMixt(iT-1,i)-ppiclf_dTdtPlag(iT-1,i))

         g = (g2+g1)*dsi4
         qq_du = qq_du + g


!        if(ppiclf_iprop(5,i) .eq. 7 .and.
!    >      ppiclf_iprop(7,i) .eq. 1724) then
!           !open(unit=66,file='fort.66',position='append')   

!           iT = 2
!           time = ppiclf_timeBH(iT)
!           write(66,*) ppiclf_time, iT, time,
!    >            factor_du, kernel_du,
!    >            ppiclf_dTdtMixt(iT,i),
!    >            ppiclf_dTdtPlag(iT,i),
!    >            (ppiclf_dTdtMixt(iT,i) - ppiclf_dTdtPlag(iT,i)),
!    >            factor_du*kernel_du*
!    >            (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i)),
!    >            qq_du,fvux,fvuy,fvuz,
!    >            kernelVU, A, B, factor 

!           flush(66)
!        endif
      endif

      if (ppiclf_iprop(5,i) .eq. i1 .and.
     >    ppiclf_iprop(7,i) .eq. i2) then
          if (istage == 3) then
             write(66,*) ppiclf_time,qq_du,fvux,fvuy,fvuz
          endif
      endif

      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
! Modified Dec. 1, 2025
!
! Shift arrays for Viscous Unsteady Force and Diffusive Unsteady
!    Heat Transfer
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

            ppiclf_dTdtMixt(iT,i) = ppiclf_dTdtMixt(iT-1,i)
            ppiclf_dTdtPlag(iT,i) = ppiclf_dTdtPlag(iT-1,i)
            ppiclf_TMixt(iT,i)    = ppiclf_TMixt(iT-1,i)
            ppiclf_TPlag(iT,i)    = ppiclf_TPlag(iT-1,i)
          
         enddo
      enddo


      if (ppiclf_nTimeBH < ppiclf_nUnsteadyData) then
            ppiclf_nTimeBH = ppiclf_nTimeBH + 1
      endif

      ! Note that ppiclf_timeBH stores s=t-tau, not tau
      ! i.e., ppiclf_timeBH(i) == s(i) == t(N-i+2), i=2,...,M
      do iT = ppiclf_nTimeBH,2,-1
            ppiclf_timeBH(it) = ppiclf_timeBH(iT-1) + ppiclf_dt
      enddo
      ppiclf_timeBH(1) = 0.0d0


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
! Note: index 1 is for current time
!
! See libpicl/user_files/ppiclf_user_AddedMass.f
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_UpdatePlag(i,iStage)
!
      implicit none
!
      include "PPICLF"
!
      integer*4 i, iStage, jj
      real*8 SDrho
      real*8 ug,vg,wg
      real*8 up,vp,wp
      real*8 vgradrho
      real*8 dt, time_plot

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

      ppiclf_TMixt(1,i) = ppiclf_rprop(PPICLF_R_JT,i)
      ppiclf_TPlag(1,i) = ppiclf_y(PPICLF_JT,i)

      ! Set the old step equal to the first one initially
      if(ppiclf_TMixt(2,i) .eq. 0.0d0) then
        ppiclf_TMixt(2,i) = ppiclf_TMixt(1,i)
      endif

      if(ppiclf_TPlag(2,i) .eq. 0.0d0) then
        ppiclf_TPlag(2,i) = ppiclf_TPlag(1,i)
      endif

      if(iStage==1) then
        dt = 5.0/15.0*ppiclf_dt
        time_plot = ppiclf_time
      elseif(iStage==2) then
        dt = (5.0/15.0 + 8.0/15.0)*ppiclf_dt
        time_plot = 5.0/15.0*dt + ppiclf_time
      elseif(iStage==3) then
        dt = ppiclf_dt
        time_plot = 13.0/15.0*dt + ppiclf_time
      else
        print*, "Unsupported beyond RK3 in Viscous-Unsteady"
        call ppiclf_exittr('Unsupported iStage', 0.0, 0)
      endif

      ! dT/dt update for Diffusive Unsteady HT
      ! first-order backward difference (explicit backward Euler form)
        ppiclf_dTdtMixt(1,i) = (ppiclf_TMixt(1,i)-ppiclf_TMixt(2,i))/dt
        ppiclf_dTdtPlag(1,i) = (ppiclf_TPlag(1,i)-ppiclf_TPlag(2,i))/dt

        if(ppiclf_iprop(5,i) .eq. 7 .and.
     >     ppiclf_iprop(7,i) .eq. 1724) then
          !open(unit=67,file='fort.67',position='append')   
          write(67,*) ppiclf_time, time_plot, dt, ppiclf_dt,
     >                ppiclf_dTdtMixt(1,i),     
     >                ppiclf_TMixt(1,i),
     >                ppiclf_TMixt(2,i),
     >                ppiclf_dTdtPlag(1,i),
     >                ppiclf_ydot(PPICLF_JT,i),
     >                ppiclf_TPlag(1,i),
     >                ppiclf_TPlag(2,i)
          flush(67)
      endif
!----------------------------------------------------------      

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
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_dTdtMixt(iT,i) = ppiclf_rprop3(k,i)
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_dTdtPlag(iT,i) = ppiclf_rprop3(k,i)
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_TMixt(iT,i) = ppiclf_rprop3(k,i)
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_TPlag(iT,i) = ppiclf_rprop3(k,i)
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
! Sets rprop3 from drudtMixt, drudtPlag, dTdtMixt, and dTdtPlag
! Needed for proper particle tracking
! Load particle data into communication buffers rprop3
! See rocpart/PLAG_RFLU_ModComm.F90:
!     SUBROUTINE PLAG_RFLU_UnloadBuffersRecv(pRegion)
!
! Modified Oct. 24, 2025 - Added dTdtMixt and dTdtPlag for diffusive
!                          unsteady heat transfer model
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
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_dTdtMixt(iT,i)
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_dTdtPlag(iT,i)
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_TMixt(iT,i)
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_TPlag(iT,i)
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
      ! Thierry - for PseudoTurbulence
      real*8 vmagj, asndfj, rhofj, dpj, phij, rej, mpj

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

         ! Filter widths are input values from *.inp
         ! ppiclf_d2chk(1 & 2) are 1/2 filter width - user defined
         !   = FILTERWIDTH / 2
         ! ppiclf_d2chk(3) is neighbor width - user defined
         !   = max(NEIGHBORWIDTH,4*Dp)
        
         if (am_flag == 2 .and. rdiff <= ppiclf_d2chk(3)) then
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
         end if ! am_flag==2 .and. rdiff <= ppiclf_d2chk(3)

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

         dist2 = ppiclf_d2chk(2)

         ! Box filter half-width dist2
         if (qs_fluct_filter_adapt_flag==0) then
            dist2 = ppiclf_d2chk(2)
         else if (qs_fluct_filter_adapt_flag>=1) then
            ! Adaptive filter defined wrt particle i
            ! Used for adaptive box or gaussian
            dpl = rpropi(PPICLF_R_JDP)
            phip = rpropi(PPICLF_R_JPHIP)
            adptfilter = ( 10.*(dpl**3)/max(1.e-4,phip) )**(1./3.)
            adptfilter = adptfilter/2.0
            dist2 = max(ppiclf_d2chk(2),adptfilter)
         endif

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
            gkern = sqrt(pi*ppiclf_filter**2/
     >              (4.0d0*log(2.0d0)))**(-ppiclf_ndim) * 
     >              exp(-dist**2/(ppiclf_filter**2/(4.0d0*log(2.0d0))))

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

         
         ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_ydotc(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_ydotc(PPICLF_JVY,i)
     >                              + rnmag*rn_12y
         ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_ydotc(PPICLF_JVZ,i)
     >                              + rnmag*rn_12z
        
         !write(*,*) "Wall NEAR",i,ppiclf_ydotc(PPICLF_JVY,i)  
      endif


      return
      end
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for projection map
!
!
! This subroutine maps certain ppiclf terms onto the fluid mesh
!
!    This routine is needed for
!       (i)  feedback_flag = 1 (rocpicl/PICL_TEMP_Runge.F90)
!       (ii) # PROBE is in the *.inp file (libfloflu/WriteProbe.F90)

!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop,rprop4)
!
      implicit none
      include "PPICLF"
!
! Input:
!
      real*8 y    (PPICLF_LRS)
      real*8 ydot (PPICLF_LRS)
      real*8 ydotc(PPICLF_LRS)
      real*8 rprop(PPICLF_LRP)
      real*8 rprop4(PPICLF_LRP4)
!
! Output:
!
      real*8 map  (PPICLF_LRP_PRO)
!
! Code:
!
!
      !Add here sploading factor         
      map(PPICLF_P_JPHIP) = rprop(PPICLF_R_JVOLP)*rprop(PPICLF_R_JSPL)
      map(PPICLF_P_JFX)   = ydotc(PPICLF_JVX) * rprop(PPICLF_R_JSPL)
      map(PPICLF_P_JFY)   = ydotc(PPICLF_JVY) * rprop(PPICLF_R_JSPL)
      map(PPICLF_P_JFZ)   = ydotc(PPICLF_JVZ) * rprop(PPICLF_R_JSPL)
      map(PPICLF_P_JE)    = ydotc(PPICLF_JT)  * rprop(PPICLF_R_JSPL)

      ! TLJ - modified 12/21/2024
      map(PPICLF_P_JPHIPD) = rprop(PPICLF_R_JVOLP)*rprop(PPICLF_R_JRHOP)
      map(PPICLF_P_JPHIPU) = rprop(PPICLF_R_JVOLP)*y(PPICLF_JVX)
      map(PPICLF_P_JPHIPV) = rprop(PPICLF_R_JVOLP)*y(PPICLF_JVY)
      map(PPICLF_P_JPHIPW) = rprop(PPICLF_R_JVOLP)*y(PPICLF_JVZ)
      map(PPICLF_P_JPHIPT) = rprop(PPICLF_R_JVOLP)*y(PPICLF_JT)

      ! 07/21/2025 - Thierry - Added Reynolds Subgrid Stress Tensor
      ! Still thinking whether I should expand the map array to include the projected Reynolds
      ! Stress Tensor values or just take them from rprop2 

      map(PPICLF_P_JRSG11) = rprop4(PPICLF_R_JRSG11)
      map(PPICLF_P_JRSG12) = rprop4(PPICLF_R_JRSG12)
      map(PPICLF_P_JRSG13) = rprop4(PPICLF_R_JRSG13)
      map(PPICLF_P_JRSG21) = rprop4(PPICLF_R_JRSG21)
      map(PPICLF_P_JRSG22) = rprop4(PPICLF_R_JRSG22)
      map(PPICLF_P_JRSG23) = rprop4(PPICLF_R_JRSG23)
      map(PPICLF_P_JRSG31) = rprop4(PPICLF_R_JRSG31)
      map(PPICLF_P_JRSG32) = rprop4(PPICLF_R_JRSG32)
      map(PPICLF_P_JRSG33) = rprop4(PPICLF_R_JRSG33)

      map(PPICLF_P_JTSG1)  = rprop4(PPICLF_R_JTSG1)
      map(PPICLF_P_JTSG2)  = rprop4(PPICLF_R_JTSG2)
      map(PPICLF_P_JTSG3)  = rprop4(PPICLF_R_JTSG3)

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
      integer*4 i,j,k

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
            fH     = 0.75d0 + .1055d0*rep
            factor = 3.0d0*rpi*rmu*dp*fac
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
      print*, "stringo", ppiclf_nid, ilen, stringo
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
!

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
      call ppiclf_blank(stringo,132)
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
      call ppiclf_blank(stringo,132)
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
      call ppiclf_blank(stringo,132)
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
      jz    = 1
      if (ppiclf_ndim .eq. 3)
     >jz    = 3

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
      enddo
      do j=5,7
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
      subroutine ppiclf_io_WriteSubBinVTU(filein1)
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
      character*6  prostr
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth, nvtx_total, ncll_total
      integer*8 idisp_pos
      integer*8 stride_lenv
      integer*4 icount_pos(PPICLF_BX1, PPICLF_BY1, PPICLF_BZ1)
      real*4 rpoint(3)
      integer*4 iint, nnp, nxx, ndxgpp1, ndygpp1, ndxygpp1, if_sz,
     >          isize, ibin, jbin, kbin, ndumx, ndumy, i, j, k,
     >          ie, ndum, kmax, ii, jj, kk, itmp, jtmp, ktmp, il, ir,
     >          jl, jr, kl, kr, npa, npb, npc, npd, npe, npf, npg, nph,
     >          ioff_dum, itype, ivtu_size, if_pos, ierr, icount_dum,
     >          iorank
!

      call ppiclf_printsi(' *Begin WriteSubBinVTU$',ppiclf_cycle)

      call ppiclf_prints(' *Begin InitSolve$')
         call ppiclf_solve_InitSolve
      call ppiclf_prints('  End InitSolve$')

      call ppiclf_prints(' *Begin CreateSubBin$')
         call ppiclf_comm_CreateSubBin
      call ppiclf_prints('  End CreateSubBin$')

      call ppiclf_prints(' *Begin ProjectParticleSubBin$')
         call ppiclf_solve_ProjectParticleSubBin
      call ppiclf_prints('  End ProjectParticleSubBin$')

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      ncll_total = ppiclf_n_bins(1)*(ppiclf_bx-1)
     >            *ppiclf_n_bins(2)*(ppiclf_by-1)
      if (ppiclf_ndim.eq.3) ncll_total = ncll_total
     >            *ppiclf_n_bins(3)*(ppiclf_bz-1)
      nvtx_total = (ppiclf_n_bins(1)*(ppiclf_bx-1)+1)
     >            *(ppiclf_n_bins(2)*(ppiclf_by-1)+1)
      if (ppiclf_ndim.eq.3) nvtx_total = nvtx_total
     >            *(ppiclf_n_bins(3)*(ppiclf_bz-1)+1)

      ndxgpp1 = ppiclf_n_bins(1) + 1
      ndygpp1 = ppiclf_n_bins(2) + 1
      ndxygpp1 = ndxgpp1*ndygpp1

      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'grd'
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

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

! test skip
c     goto 1511

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
      write(vtu,'(A)',advance='no') 'type="ProjectParticleSubBinInt32" '
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
      do ie=1,PPICLF_LRP_PRO
         write(prostr,'(A4,I2.2)') "PRO-",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*nvtx_total + isize
      enddo
c     call ppiclf_io_WriteDataArrayVTU(vtu,"PPR",1,iint)
c     iint = iint + 1*isize*nvtx_total + isize
c     call ppiclf_io_WriteDataArrayVTU(vtu,"NID",1,iint)
c     iint = iint + 1*isize*nvtx_total + isize
      write(vtu,'(A)',advance='yes') '   </PointData> '
      write(vtu,'(A)',advance='yes') '   <CellData>'
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
      ndumx = ppiclf_n_bins(1)*(ppiclf_bx-1) + 1
      ndumy = ppiclf_n_bins(2)*(ppiclf_by-1) + 1
      ndum = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)
      do ie=0,ndum-1
         i = modulo(ie,ppiclf_n_bins(1))
         j = modulo(ie/ppiclf_n_bins(1),ppiclf_n_bins(2))
         k = ie/(ppiclf_n_bins(1)*ppiclf_n_bins(2))

      kmax = 1
      if (ppiclf_ndim .eq. 3) kmax = ppiclf_bz-1

      do kk=1,kmax
      do jj=1,ppiclf_by-1
      do ii=1,ppiclf_bx-1

         itmp = i*(ppiclf_bx-1) + (ii-1)
         jtmp = j*(ppiclf_by-1) + (jj-1)
         ktmp = 0
         if (ppiclf_ndim .eq. 3)
     >   ktmp = k*(ppiclf_bz-1) + (kk-1)

         kl = ktmp
         kr = ktmp+1
         jl = jtmp
         jr = jtmp+1
         il = itmp
         ir = itmp+1

c        ndum = itmp + ndumx*jtmp + ndumx*ndumy*ktmp

         if (ppiclf_ndim .eq. 3) then
            npa = il + ndumx*jl + ndumx*ndumy*kl
            npb = ir + ndumx*jl + ndumx*ndumy*kl
            npc = il + ndumx*jr + ndumx*ndumy*kl
            npd = ir + ndumx*jr + ndumx*ndumy*kl
            npe = il + ndumx*jl + ndumx*ndumy*kr
            npf = ir + ndumx*jl + ndumx*ndumy*kr
            npg = il + ndumx*jr + ndumx*ndumy*kr
            nph = ir + ndumx*jr + ndumx*ndumy*kr

            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0,A)',advance='no')  npd, ' '
            write(vtu,'(I0,A)',advance='no')  npe, ' '
            write(vtu,'(I0,A)',advance='no')  npf, ' '
            write(vtu,'(I0,A)',advance='no')  npg, ' '
            write(vtu,'(I0)'  ,advance='yes') nph
         else
            npa = il + ndumx*jl + ndumx*ndumy*kl
            npb = ir + ndumx*jl + ndumx*ndumy*kl
            npc = il + ndumx*jr + ndumx*ndumy*kl
            npd = ir + ndumx*jr + ndumx*ndumy*kl

            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0)'  ,advance='yes') npd
         endif
      enddo
      enddo
      enddo
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ioff_dum = 4
      if (ppiclf_ndim .eq. 3) ioff_dum = 8
      ! write offsetts here
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

      do k=1,ppiclf_bz
      do j=1,ppiclf_by
      do i=1,ppiclf_bx
         icount_pos(i,j,k) = 0
      enddo
      enddo
      enddo
      if (ppiclf_nid .le. 
     >      ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ppiclf_ndim .eq. 3) then
         do k=1,ppiclf_bz
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            if (i .ne. ppiclf_bx .and.
     >          j .ne. ppiclf_by .and.
     >          k .ne. ppiclf_bz) then
                  icount_pos(i,j,k) = 3
            endif

            if (i .eq. ppiclf_bx) then
            if (j .ne. ppiclf_by) then
            if (k .ne. ppiclf_bz) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (j .eq. ppiclf_by) then
            if (i .ne. ppiclf_bx) then
            if (k .ne. ppiclf_bz) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (k .eq. ppiclf_bz) then
            if (i .ne. ppiclf_bx) then
            if (j .ne. ppiclf_by) then
            if (kbin .eq. ppiclf_n_bins(3)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (j .eq. ppiclf_by) then
            if (k .ne. ppiclf_bz) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (k .eq. ppiclf_bz) then
            if (j .ne. ppiclf_by) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
            if (kbin .eq. ppiclf_n_bins(3)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (j .eq. ppiclf_by) then
            if (k .eq. ppiclf_bz) then
            if (i .ne. ppiclf_bx) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
            if (kbin .eq. ppiclf_n_bins(3)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (j .eq. ppiclf_by) then
            if (k .eq. ppiclf_bz) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
            if (kbin .eq. ppiclf_n_bins(3)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif
            endif

         enddo
         enddo
         enddo
         ! 2D
         else
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            if (i .ne. ppiclf_bx .and.
     >          j .ne. ppiclf_by) then
                  icount_pos(i,j,1) = 3
            endif

            if (i .eq. ppiclf_bx) then
            if (j .ne. ppiclf_by) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
               icount_pos(i,j,1) = 3
            endif
            endif
            endif

            if (j .eq. ppiclf_by) then
            if (i .ne. ppiclf_bx) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
               icount_pos(i,j,1) = 3
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (j .eq. ppiclf_by) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
               icount_pos(i,j,1) = 3
            endif
            endif
            endif
            endif
         enddo
         enddo

         endif

      endif

      do k=1,ppiclf_bz
      do j=1,ppiclf_by
      do i=1,ppiclf_bx
         stride_lenv = ppiclf_grid_i(i,j,k)
         idisp_pos   = ivtu_size + isize*(3*stride_lenv + 1)
         icount_dum  = icount_pos(i,j,k)
         rpoint(1)   = ppiclf_grid_x(i,j,k)
         rpoint(2)   = ppiclf_grid_y(i,j,k)
         rpoint(3)   = 0.0d0
         if (ppiclf_ndim .eq. 3)
     >   rpoint(3)   = ppiclf_grid_z(i,j,k)
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_dum,iorank,pth,ierr)

      enddo
      enddo
      enddo

      call ppiclf_byte_close_mpi(pth,ierr)


      ! projected fields
      do ie=1,PPICLF_LRP_PRO

      if_pos = 1*isize*nvtx_total

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      do k=1,ppiclf_bz
      do j=1,ppiclf_by
      do i=1,ppiclf_bx
           stride_lenv = ppiclf_grid_i(i,j,k)
           idisp_pos   = ivtu_size + isize*(3*nvtx_total ! position fld
     >                   + (ie-1)*nvtx_total ! prev fields
     >                   + 1*stride_lenv    ! this fld
     >                   + 1 + ie)          ! ints
           icount_dum  = icount_pos(i,j,k)/3 ! either zero or 1
           rpoint(1)   = ppiclf_grid_fld(i,j,k,ie)
           call ppiclf_byte_set_view(idisp_pos,pth)
           call ppiclf_byte_write_mpi(rpoint,icount_dum,iorank,pth,ierr)
      enddo
      enddo
      enddo

      call ppiclf_byte_close_mpi(pth,ierr)

      enddo

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi('  End WriteSubBinVTU$',ppiclf_cycle)

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

      call ppiclf_printsi(' *Begin WriteBinVTU$',ppiclf_cycle)

      if (icalld1 .eq. 0) icalld1 = istartout

      icalld1 = icalld1+1

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
      rpoint(1)  = sngl(ppiclf_binx(1,1))
      rpoint(2)  = sngl(ppiclf_biny(1,1))
      rpoint(3)  = 0.0
      if (ppiclf_ndim .eq. 3)
     >rpoint(3)  = sngl(ppiclf_binz(1,1))
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
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(1,1))
            rpoint(3)  = sngl(ppiclf_binz(1,1))
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
            rpoint(1)  = sngl(ppiclf_binx(1,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = sngl(ppiclf_binz(1,1))
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
            rpoint(1)  = sngl(ppiclf_binx(1,1))
            rpoint(2)  = sngl(ppiclf_biny(1,1))
            rpoint(3)  = sngl(ppiclf_binz(2,1))
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
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = sngl(ppiclf_binz(1,1))
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
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(1,1))
            rpoint(3)  = sngl(ppiclf_binz(2,1))
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
            rpoint(1)  = sngl(ppiclf_binx(1,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = sngl(ppiclf_binz(2,1))
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
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = sngl(ppiclf_binz(2,1))
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
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(1,1))
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
            rpoint(1)  = sngl(ppiclf_binx(1,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
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
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
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
     >       ,rout_lrp5(PPICLF_LRP5*PPICLF_LPART)
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
     >          ierr, ivtu_size, ic_lrp5
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
      ic_lrp5 = iadd
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
      do j=5,7
      do i=1,nxx
         ic_lip = ic_lip + 1
         rout_lip(ic_lip) = real(ppiclf_iprop(j,i))
      enddo
      enddo

!      ! testing if writing of forces works
!      do j=1,PPICLF_LRP5
!      do i=1,nxx
!         ic_lrp5 = ic_lrp5 + 1
!         rout_lrp5(ic_lrp5) = sngl(ppiclf_force(j,i))
!      enddo
!      enddo

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

!      do ie=1,PPICLF_LRP5
!         write(prostr,'(A4,I2.2)') "force",ie
!         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
!         iint = iint + 1*isize*npt_total + isize
!      enddo

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

! Thierry - added here for testing of force writing
! ----------------------------------------------------------------------
!      do i=1,PPICLF_LRP5
!         idisp_lip = ivtu_size + isize*(3*npt_total  
!     >                         + PPICLF_LRS*npt_total
!     >                         + PPICLF_LRP*npt_total
!     >                         + 3*npt_total    
!     >                         + (i-1)*npt_total
!     >                         + (1)*stride_len
!     >                         + 1 + PPICLF_LRS + PPICLF_LRP + 3 + i)
!
!         ! integer write
!         if (ppiclf_nid .eq. 0) then
!           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
!     >         ,position='append')
!           write(vtu) if_lrp
!           close(vtu)
!         endif
!   
!         call mpi_barrier(ppiclf_comm,ierr)
!
!         j = (i-1)*ppiclf_npart + 1
!   
!         ! write
!         call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
!         call ppiclf_byte_set_view(idisp_lrp,pth)
!         call ppiclf_byte_write_mpi(rout_lrp(j),icount_lrp,iorank,pth
!     >                             ,ierr)
!         call ppiclf_byte_close_mpi(pth,ierr)
!      enddo
! ----------------------------------------------------------------------

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
      if (ppiclf_lsubbin) call ppiclf_io_OutputDiagSubBin
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

         nel_max_orig   = ppiclf_iglmax(ppiclf_nee,1)
         nel_min_orig   = ppiclf_iglmin(ppiclf_nee,1)
         nel_total_orig = ppiclf_iglsum(ppiclf_nee,1)

         nel_max_map   = ppiclf_iglmax(ppiclf_neltb,1)
         nel_min_map   = ppiclf_iglmin(ppiclf_neltb,1)
         nel_total_map = ppiclf_iglsum(ppiclf_neltb,1)

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
      subroutine ppiclf_io_OutputDiagSubBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 nbin_total
!

      call ppiclf_prints(' *Begin SubBin Info$')

      nbin_total = ppiclf_bx*ppiclf_by*ppiclf_bz

      call ppiclf_printsi('  -Number of local bin    :$',nbin_total)
      call ppiclf_printsi('  -Number of local bin (x):$',PPICLF_BX)
      call ppiclf_printsi('  -Number of local bin (y):$',PPICLF_BY)
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsi('  -Number of local bin (z):$',PPICLF_BZ)
      call ppiclf_printsr('  -Bin width (x)          :$',ppiclf_rdx)
      call ppiclf_printsr('  -Bin width (y)          :$',ppiclf_rdy)
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin width (z)          :$',ppiclf_rdz)
      call ppiclf_printsr('  -Filter width           :$',ppiclf_filter)
      call ppiclf_printsr('  -Filter cut-off         :$'
     >                                                 ,ppiclf_d2chk(2))
      call ppiclf_printsi('  -SubBins per filter res.:$',ppiclf_ngrids)

      call ppiclf_prints('  End SubBin Info$')

      return
      end
!-----------------------------------------------------------------------
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
      do ie=1,ppiclf_neltb
         call ppiclf_copy(ppiclf_xm1bi(1,1,1,ie,1)
     >                   ,ppiclf_xm1b(1,1,1,1,ie),n)
         call ppiclf_copy(ppiclf_xm1bi(1,1,1,ie,2)
     >                   ,ppiclf_xm1b(1,1,1,2,ie),n)
         call ppiclf_copy(ppiclf_xm1bi(1,1,1,ie,3)
     >                   ,ppiclf_xm1b(1,1,1,3,ie),n)
      ENDdo

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
      do ie=1,ppiclf_neltbbb
         call ppiclf_icopy(ppiclf_er_mapc(1,ie),ppiclf_er_maps(1,ie)
     >             ,PPICLF_LRMAX)
      ENDdo
      !PRINT*, 'Processor ID, ppiclf_neltbbb', ppiclf_nid, ppiclf_neltbbb
      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InterpField(j)
!
      implicit none
!
      include "PPICLF"
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
      do ie=1,ppiclf_neltbbb
         iee = ppiclf_er_mapc(1,ie)
         call ppiclf_copy(ppiclf_int_fld (1,1,1,j  ,ie)
     >                   ,ppiclf_int_fldu(1,1,1,iee,j ),n)
      ENDdo

      RETURN
      END
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_FinalizeInterp
!
      implicit none
!
      include "PPICLF"
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

         do ie=1,ppiclf_neltbbb
            call ppiclf_copy(fld(1,1,1,ie)
     >                      ,ppiclf_int_fld(1,1,1,i,ie),nxyz)
         ENDdo

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
