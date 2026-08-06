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
#ifdef PERF
      real*8 ppiclf_pt0
      real*8 MPI_WTIME
      external MPI_WTIME
#endif
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
         !rphip = min(rphip,0.62d0)
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
         taux_hydro = 0.0d0; tauy_hydro = 0.0d0;
         tauz_hydro = 0.0d0;
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
#ifdef PERF
          ppiclf_pt0 = MPI_WTIME()
#endif
          CALL ppiclf_solve_NearestNeighborSB(i)
#ifdef PERF
          PPICLF_TPPNNSearch = PPICLF_TPPNNSearch
     >        + (MPI_WTIME() - ppiclf_pt0)
#endif

         end if ! end Step 1b; nearestneighbor


!
! Step 2: Force component quasi-steady
#ifdef PERF
      ppiclf_pt0 = MPI_WTIME()
#endif
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
#ifdef PERF
      PPICLF_TQuasiSteady = PPICLF_TQuasiSteady
     >     + (MPI_WTIME() - ppiclf_pt0)
#endif
! Step 4: Force component added mass
#ifdef PERF
      ppiclf_pt0 = MPI_WTIME()
#endif
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
#ifdef PERF
      PPICLF_TAddedMass = PPICLF_TAddedMass
     >     + (MPI_WTIME() - ppiclf_pt0)
#endif
! Step 5: Force component pressure gradient
#ifdef PERF
      ppiclf_pt0 = MPI_WTIME()
#endif
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
#ifdef PERF
      PPICLF_TPresGrad = PPICLF_TPresGrad
     >     + (MPI_WTIME() - ppiclf_pt0)
#endif
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
#ifdef PERF
      ppiclf_pt0 = MPI_WTIME()
#endif
!
         if (heattransfer_flag >= 1) then
            call ppiclf_user_HT_driver(i,qq)
         endif ! heattransfer_flag >= 1


!
#ifdef PERF
      PPICLF_THeatTransfer = PPICLF_THeatTransfer
     >     + (MPI_WTIME() - ppiclf_pt0)
#endif
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
          ! With added mass treated implicitly on the particle LHS,
          ! (rmass+rmass_add)*ydot = fqs+fam+fdpd+fvu+flift+fc exactly,
          ! so subtracting fc returns the full hydrodynamic force per
          ! Daoud et al. (2026) Eq. (23). Using rmass alone underweights
          ! the feedback by m/(m+m_am) and leaks collision force into
          ! the gas whenever am_flag > 0.
          ppiclf_feedbk(PPICLF_P_JFX,i) = ppiclf_rprop(PPICLF_R_JSPL,i)*
     >      (ppiclf_ydot(PPICLF_JVX,i)*(rmass+rmass_add) - fcx)
          ppiclf_feedbk(PPICLF_P_JFY,i) = ppiclf_rprop(PPICLF_R_JSPL,i)*
     >      (ppiclf_ydot(PPICLF_JVY,i)*(rmass+rmass_add) - fcy)
          ppiclf_feedbk(PPICLF_P_JFZ,i) = ppiclf_rprop(PPICLF_R_JSPL,i)*
     >      (ppiclf_ydot(PPICLF_JVZ,i)*(rmass+rmass_add) - fcz)

          ! Energy equation feedback term
          ! Work coupling per Daoud, Jackson & Balachandar, IJMF (2026)
          ! Eq. (18)/(21): dissipative forces (QS incl. fluctuations,
          ! VU, lift) do work at the PARTICLE velocity; non-dissipative
          ! forces (added mass, pressure gradient) do work at the GAS
          ! velocity, so they cancel in the dissipation and conserve
          ! total energy. 
          ! 09/19/2025 - Thierry - Added Lift force
          ppiclf_feedbk(PPICLF_P_JE,i) = ppiclf_rprop(PPICLF_R_JSPL,i)
     >     * ( (fqsx+fvux+liftx)*ppiclf_y(PPICLF_JVX,i) + 
     >         (fqsy+fvuy+lifty)*ppiclf_y(PPICLF_JVY,i) + 
     >         (fqsz+fvuz+liftz)*ppiclf_y(PPICLF_JVZ,i) +
     >                famx*ppiclf_rprop(PPICLF_R_JUX,i) +
     >                famy*ppiclf_rprop(PPICLF_R_JUY,i) +
     >                famz*ppiclf_rprop(PPICLF_R_JUZ,i) +
              ! PG work deliberately EXCLUDED: Rocflu's MJWL energy
              ! flux is vFracF*rho*u*H = phig*(rho u E + u P), so the
              ! particle share of pressure work never enters the gas
              ! energy equation and must not be removed by feedback.
              ! (Daoud et al. Eq. 21 puts PG work at u_g because THEIR
              ! energy equation carries the unscaled div(P u); the
              ! feedback must match the discretization, not the paper.)
              ! The fdpd momentum feedback stays: the momentum flux
              ! carries the UNSCALED grad(P), so that share IS
              ! over-applied to the gas and must be returned.
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
          ! WARNING: this branch routes the added-mass slip work
          ! Fam.(ug-up) to the PTKE (Rsg) projections, which are
          ! currently commented out below -- that energy is lost
          ! until the Rsg/Tsg feedback path is restored.
          ELSE
            Rsg   = 0.0D0
            T_par = 0.0D0
          END IF ! pseudoTurb_flag


          !  ****************************************************
          ! Avery added for troubleshooting negative internal energy.
          !ppiclf_feedbk(PPICLF_P_JE,i) = 0.0!qq
          !  ****************************************************

          ! 07/21/2025 - Thierry - Added Reynolds Subgrid Stress Feedback
!          ppiclf_feedbk(PPICLF_P_JRSG11,i) = Rsg(1,1) 
!     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
!          ppiclf_feedbk(PPICLF_P_JRSG12,i) = Rsg(1,2) 
!     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
!          ppiclf_feedbk(PPICLF_P_JRSG13,i) = Rsg(1,3) 
!     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
!          ppiclf_feedbk(PPICLF_P_JRSG21,i) = Rsg(2,1) 
!     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
!          ppiclf_feedbk(PPICLF_P_JRSG22,i) = Rsg(2,2) 
!     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
!          ppiclf_feedbk(PPICLF_P_JRSG23,i) = Rsg(2,3) 
!     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
!          ppiclf_feedbk(PPICLF_P_JRSG31,i) = Rsg(3,1) 
!     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
!          ppiclf_feedbk(PPICLF_P_JRSG32,i) = Rsg(3,2) 
!     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
!          ppiclf_feedbk(PPICLF_P_JRSG33,i) = Rsg(3,3)  
!     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
!
!          ppiclf_feedbk(PPICLF_P_JTSG1,i) = T_par(1) 
!     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
!          ppiclf_feedbk(PPICLF_P_JTSG2,i) = T_par(2) 
!     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)
!          ppiclf_feedbk(PPICLF_P_JTSG3,i) = T_par(3) 
!     >                                   * ppiclf_rprop(PPICLF_R_JSPL,i)

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
