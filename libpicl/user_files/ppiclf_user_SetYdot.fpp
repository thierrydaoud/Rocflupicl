#:include "PPICLF_PARTMACROS.fypp"
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
submodule (ppiclf_user) ppiclf_user_SetYdot_imp
    ! particle data
    use ppiclf_data, only: ppiclf_npart
    use ppiclf_m_particledata, only: @{USEMODVAR(PPICLF_t_particle, ppiclf_parts)}@, @{USEMODVAR(PPICLF_t_ghostParticle, ppiclf_gparts)}@
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

    use mpi, only: MPI_WTIME
    use ppiclf_op, only: ppiclf_exittr
    use ppiclf_solve, only: ppiclf_solve_NearestNeighborSB

    use ppiclf_m_user_ForceModels, only: ppiclf_user_BR_driver, ppiclf_user_Lift_driver, ppiclf_user_Torque_driver, ppiclf_user_HT_driver
    use ppiclf_m_user_ForceModels, only: ppiclf_user_QS_Parmar, ppiclf_user_QS_ModifiedParmar, ppiclf_user_QS_Osnes, ppiclf_user_QS_Gidaspow
    use ppiclf_m_user_ForceModels, only: ppiclf_user_QS_fluct_Lattanzi, ppiclf_user_QS_fluct_Osnes
    use ppiclf_m_user_ForceModels, only: ppiclf_user_VU_Rocflu, ppiclf_user_UpdatePlag, ppiclf_user_ShiftUnsteadyData, ppiclf_user_plag2prop, ppiclf_user_prop2plag
    use ppiclf_m_user_ForceModels, only: ppiclf_user_AM_Briney_Unary, ppiclf_user_AM_Briney_Binary, ppiclf_user_AM_Parmar
    use ppiclf_m_user_SubbinMap, only: ppiclf_user_subbinMap
    implicit none
    contains

module procedure ppiclf_user_SetYdot
        
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
    real*8 taux, tauy, tauz, rmass_omega, taux_hydro, tauy_hydro, tauz_hydro 
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
    
    INTEGER*4 SBin_map( 0 : ( (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_nndist)/ppiclf_nndist)         + 1) * (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_nndist)/ppiclf_nndist)        + 1) * (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_nndist)/ppiclf_nndist)        + 1) - 1), (ppiclf_npart+ppiclf_npart_gp))
    INTEGER*4  SBin_counter( 0 : ( (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_nndist)/ppiclf_nndist)         + 1) * (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_nndist)/ppiclf_nndist)        + 1) * (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_nndist)/ppiclf_nndist)        + 1) - 1))
    INTEGER*4 i_Bin(3), n_SBin(3), tot_SBin

#ifdef PERF
    REAL*8 tstart, tfinal
#endif

    ! Unit Test only Code:
    !-----------------------------------------------------------------------
    !
#ifdef TEST
    ! unit_test only tests the SB nearest neighbor search in this
    ! subroutine.  The full subroutine is called to ensure that
    ! the array initialization is correct.
    CALL ppiclf_user_subbinMap(i_Bin, n_SBin, tot_SBin,SBin_counter ,SBin_map)
    DO i = 1,ppiclf_npart
        CALL ppiclf_solve_NearestNeighborSB(i,tot_SBin,SBin_counter,SBin_map,n_SBin,i_Bin)
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
            write(6,'(a,2x,2(1pe14.6),2x,i3)') '*** PPICLF dt, time = ', ppiclf_dt,ppiclf_time
        endif
    endif

    burnrate_model = 0
    if (burnrate_flag .gt. 0) then
        if ( TRIM(ppiclf_matname)=='AL'.or. TRIM(ppiclf_matname)=='Al' ) then
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
        if (ppiclf_nid==0) print*,'dt,fac=',istage,ppiclf_dt,fac,stationary, qs_flag, am_flag, pg_flag,collisional_flag, heattransfer_flag, feedback_flag,qs_fluct_flag, ppiclf_debug, ppiclf_nTimeBH,ppiclf_nUnsteadyData
    endif

    

    !
    !-----------------------------------------------------------------------
    !
    ! Reapply axi-sym collision correction
    ! Right now hard coding smallest radius  
    do i=1,ppiclf_npart
       @{USEPARTICLE(ppiclf_parts(i)%rprop%JDPe)}@ = (0.00005/(@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPT)}@)) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%DP)}@) 
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
    if(PPICLF_PPInteractions) then
#ifdef PERF
        tstart = MPI_WTIME()
#endif
        call ppiclf_user_subbinMap(i_Bin, n_SBin, tot_SBin,SBin_counter ,SBin_map)
#ifdef PERF
        tfinal = MPI_WTIME()
        PPICLF_TParticleParticleModels = tfinal - tstart
#endif
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
            temp    = @{USEPARTICLE(ppiclf_parts(i)%rprop%JT)}@
            rmu     = rmu_ref*sqrt(temp/tref)*(1.0d0+suth/tref)/(1.0d0+suth/temp)
        else
            call ppiclf_exittr('Unknown viscosity law$', 0.0d0, 0)
        endif
        rkappa = rcp_fluid*rmu/rpr


        ! Useful values
        rmass  = @{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@ * @{USEPARTICLE(ppiclf_parts(i)%rprop%RHOP)}@
        vx     = @{USEPARTICLE(ppiclf_parts(i)%rprop%U%X)}@ - @{USEPARTICLE(ppiclf_parts(i)%y%vel%x)}@
        vy     = @{USEPARTICLE(ppiclf_parts(i)%rprop%U%Y)}@ - @{USEPARTICLE(ppiclf_parts(i)%y%vel%y)}@
        vz     = @{USEPARTICLE(ppiclf_parts(i)%rprop%U%Z)}@ - @{USEPARTICLE(ppiclf_parts(i)%y%vel%z)}@
        vmag   = sqrt(vx*vx + vy*vy + vz*vz)
        rhof   = @{USEPARTICLE(ppiclf_parts(i)%rprop%RHOF)}@  
        dp     = @{USEPARTICLE(ppiclf_parts(i)%rprop%DP)}@
        rep    = vmag*dp*rhof/rmu
        rphip  = @{USEPARTICLE(ppiclf_parts(i)%rprop%PHIP)}@
        rphif  = 1.0d0-(@{USEPARTICLE(ppiclf_parts(i)%rprop%PHIP)}@)
        asndf  = @{USEPARTICLE(ppiclf_parts(i)%rprop%CS)}@
        rmachp = vmag/asndf
        rhop   = @{USEPARTICLE(ppiclf_parts(i)%rprop%RHOP)}@

        ! TLJ - 04/03/2025; Do not calculate forces if vmag = 0
        !       Otherwise the particles might move before the 
        !       shock arrives
        if (vmag <= 1.d-3) cycle

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

        ! ***This is bad practice and leads to difficult code to debug -Avery***
        @{USEPARTICLE(ppiclf_parts(i)%rprop%JSPT)}@ = sqrt((@{USEPARTICLE(ppiclf_parts(i)%y%vel%x)}@)**2 + (@{USEPARTICLE(ppiclf_parts(i)%y%vel%y)}@)**2 + (@{USEPARTICLE(ppiclf_parts(i)%y%vel%z)}@)**2)

        rep = max(0.1d0,rep)

        ! Redefine volume fractions
        ! Need to make sure phi_p + phi_f = 1
        rphip = @{USEPARTICLE(ppiclf_parts(i)%rprop%PHIP)}@
        rphif = 1.0d0-rphip

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
#ifdef PERF
        tstart = MPI_WTIME()
#endif
        if (am_flag == 2) then 
            ! 07/14/24 - Thierry - If Briney Algorithm flag and fluct_flag
            !   are ON -> evaluate added-mass unary term before evaluating
            !   neighbor-induced acceleration in EvalNearestNeighbor
            call ppiclf_user_AM_Briney_Unary(i,iStage,famx,famy,famz,rmass_add)
        endif ! end am_flag = 2


        !
        ! Step 1b: Call NearestNeighbor if particles i and j interact
        !
        if(PPICLF_PPInteractions) THEN

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
                phipmean = @{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@
                upmean   = @{USEPARTICLE(ppiclf_parts(i)%y%vel%x)}@
                vpmean   = @{USEPARTICLE(ppiclf_parts(i)%y%vel%y)}@
                wpmean   = @{USEPARTICLE(ppiclf_parts(i)%y%vel%z)}@
                u2pmean  = upmean**2
                v2pmean  = vpmean**2
                w2pmean  = wpmean**2
                icpmean  = 1
                else if (qs_fluct_filter_flag==1) then
                ! gaussian kernel
                ! r = 0
                !***AVERY - verify max filter dimension should be used***
                maxFilter = MAX(ppiclf_filter(1),ppiclf_filter(2), ppiclf_filter(3))
                gkern = sqrt(rpi*maxFilter**2/ (4.0d0*log(2.0d0)))**(-ppiclf_ndim)
                !***AVERY - I think phipmean is using wrong index ***
                ! Also, none of the below are means...
                phipmean = gkern * @{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@
                upmean   = gkern * (@{USEPARTICLE(ppiclf_parts(i)%y%vel%x)}@) * @{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@
                vpmean   = gkern * (@{USEPARTICLE(ppiclf_parts(i)%y%vel%y)}@) * @{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@
                wpmean   = gkern * (@{USEPARTICLE(ppiclf_parts(i)%y%vel%z)}@) * @{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@
                u2pmean  = gkern * ((@{USEPARTICLE(ppiclf_parts(i)%y%vel%x)}@)**2)* (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@)
                v2pmean  = gkern * ((@{USEPARTICLE(ppiclf_parts(i)%y%vel%y)}@)**2)* (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@)
                w2pmean  = gkern * ((@{USEPARTICLE(ppiclf_parts(i)%y%vel%z)}@)**2)* (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@)
                icpmean = 1
                end if
            end if

            ! add neighbors
            ! AVERY - always use subbins
            IF ( 1 .EQ. 1 ) THEN !sbNearest_flag .EQ. 1) THEN
                CALL ppiclf_solve_NearestNeighborSB(i,tot_SBin,SBin_counter,SBin_map,n_SBin,i_Bin)
                !         ELSE
                !             CALL ppiclf_solve_NearestNeighbor(i)
            END IF

        end if ! end Step 1b; nearestneighbor
#ifdef PERF
        tfinal = MPI_WTIME()
        PPICLF_TParticleParticleModels = PPICLF_TParticleParticleModels + (tfinal - tstart)
        ! TODO: Should this be tstart
        tfinal = MPI_WTIME()
#endif 

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
            call ppiclf_exittr('Wrong QS Model Choice$', 0.0d0, qs_flag)
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
            call ppiclf_user_QS_fluct_Osnes(i,iStage,fqs_fluct, xi_par,xi_perp,xi_T,fqsx, fqsy, fqsz)
        endif

        ! Add fluctuation part to quasi-steady force
        fqsx = fqsx + fqs_fluct(1)
        fqsy = fqsy + fqs_fluct(2)
        fqsz = fqsz + fqs_fluct(3)

        ! Store quasi-steady fluctuating force
        @{USEPARTICLE(ppiclf_parts(i)%rprop%FLUCTF%X)}@ = fqs_fluct(1)
        @{USEPARTICLE(ppiclf_parts(i)%rprop%FLUCTF%Y)}@ = fqs_fluct(2)
        @{USEPARTICLE(ppiclf_parts(i)%rprop%FLUCTF%Z)}@ = fqs_fluct(3)
        
        ! Store normally distributed random variables xi for PseudoTurbulence
        @{USEPARTICLE(ppiclf_parts(i)%rprop%XIPAR)}@  = xi_par
        @{USEPARTICLE(ppiclf_parts(i)%rprop%XIPERP)}@ = xi_perp
        @{USEPARTICLE(ppiclf_parts(i)%rprop%XIT)}@    = xi_T


        !
        ! Step 4: Force component added mass
        !
        if (am_flag == 1) then 
            call ppiclf_user_AM_Parmar(i,iStage,famx,famy,famz,rmass_add)

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
                call ppiclf_user_AM_Briney_Binary(i,iStage,famx,famy,famz,rmass_add)
                FamBinary(1) = famx - FamUnary(1)
                FamBinary(2) = famy - FamUnary(2)
                FamBinary(3) = famz - FamUnary(3)
            else
                ! if particle has no neighbors, need to multiply added mass forces
                ! by volume, as this is taken care of in Binary subroutine
                famx = famx*@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@
                famy = famy*@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@
                famz = famz*@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@
            endif
        endif


        !-----------------------------------------------------------------------

        !
        ! Step 5: Force component pressure gradient
        !
        if (pg_flag == 1) then
            fdpdx = -(@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%DPDX%X)}@)
            fdpdy = -(@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%DPDX%Y)}@)
            fdpdz = -(@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%DPDX%Z)}@)

            if (flow_model == 1) then ! Navier-Stokes Flow Model
                fdpvdx = (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%DPVDX%X)}@)
                fdpvdy = (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%DPVDX%Y)}@)
                fdpvdz = (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%DPVDX%Z)}@)
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
            
            fcx  = @{USEPARTICLE(ppiclf_parts(i)%ydotc%vel%x)}@
            fcy  = @{USEPARTICLE(ppiclf_parts(i)%ydotc%vel%y)}@
            fcz  = @{USEPARTICLE(ppiclf_parts(i)%ydotc%vel%z)}@

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
            call ppiclf_user_BR_driver(i,iStage,burnrate_model,qq,mdot_me,mdot_ox)
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
            taux  = @{USEPARTICLE(ppiclf_parts(i)%ydotc%ang_vel%x)}@
            tauy  = @{USEPARTICLE(ppiclf_parts(i)%ydotc%ang_vel%y)}@
            tauz  = @{USEPARTICLE(ppiclf_parts(i)%ydotc%ang_vel%z)}@
            call ppiclf_user_Torque_driver(i,iStage,taux,tauy,tauz,taux_hydro,tauy_hydro,tauz_hydro)
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
        @{USEPARTICLE(ppiclf_parts(i)%ydot%pos%x)}@ = @{USEPARTICLE(ppiclf_parts(i)%y%vel%x)}@
        @{USEPARTICLE(ppiclf_parts(i)%ydot%pos%y)}@ = @{USEPARTICLE(ppiclf_parts(i)%y%vel%y)}@
        @{USEPARTICLE(ppiclf_parts(i)%ydot%pos%z)}@ = @{USEPARTICLE(ppiclf_parts(i)%y%vel%z)}@
        @{USEPARTICLE(ppiclf_parts(i)%ydot%vel%x)}@ = (fqsx+famx+fdpdx+fvux+liftx+fcx)/(rmass+rmass_add)
        @{USEPARTICLE(ppiclf_parts(i)%ydot%vel%x)}@ = (fqsy+famy+fdpdy+fvuy+lifty+fcy)/(rmass+rmass_add)
        @{USEPARTICLE(ppiclf_parts(i)%ydot%vel%z)}@ = (fqsz+famz+fdpdz+fvuz+liftz+fcz)/(rmass+rmass_add)
        @{USEPARTICLE(ppiclf_parts(i)%ydot%t)}@  = qq/rmass_therm
        @{USEPARTICLE(ppiclf_parts(i)%ydot%ang_vel%x)}@ = taux/rmass_omega
        @{USEPARTICLE(ppiclf_parts(i)%ydot%ang_vel%y)}@ = tauy/rmass_omega
        @{USEPARTICLE(ppiclf_parts(i)%ydot%ang_vel%z)}@ = tauz/rmass_omega
        @{USEPARTICLE(ppiclf_parts(i)%ydot%metal)}@  = mdot_me
        @{USEPARTICLE(ppiclf_parts(i)%ydot%oxide)}@  = mdot_ox

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
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JFX)}@ = 0.0d0 
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JFY)}@ = 0.0d0 
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JFZ)}@ = 0.0d0 
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JE)}@  = 0.0d0
        END IF

        IF(feedback_flag==1) THEN
            ! Momentum equations feedback terms
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JFX)}@ = (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@) * ((@{USEPARTICLE(ppiclf_parts(i)%ydot%vel%x)}@)*rmass - fcx)
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JFY)}@ = (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@) * ((@{USEPARTICLE(ppiclf_parts(i)%ydot%vel%y)}@)*rmass - fcy)
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JFZ)}@ = (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@) * ((@{USEPARTICLE(ppiclf_parts(i)%ydot%vel%z)}@)*rmass - fcz)

            ! Energy equation feedback term
            ! 09/19/2025 - Thierry - Added Lift force
            ! Still need to add Torue \cdot angular velocity
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JE)}@ = (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)           &
                * ( (fqsx+fvux+liftx)*(@{USEPARTICLE(ppiclf_parts(i)%y%vel%X)}@)    +               &
                    (fqsy+fvuy+lifty)*(@{USEPARTICLE(ppiclf_parts(i)%y%vel%Y)}@)    +               &
                    (fqsz+fvuz+liftz)*(@{USEPARTICLE(ppiclf_parts(i)%y%vel%Z)}@)    +               &
                    famx*(@{USEPARTICLE(ppiclf_parts(i)%rprop%U%X)}@)               +               &
                    famy*(@{USEPARTICLE(ppiclf_parts(i)%rprop%U%Y)}@)               +               &
                    famz*(@{USEPARTICLE(ppiclf_parts(i)%rprop%U%Z)}@)               +               &
                    taux_hydro*(@{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%X)}@)       +               &
                    tauy_hydro*(@{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%Y)}@)       +               &
                    tauz_hydro*(@{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%Z)}@)       +               &
                    qq )

            IF(pseudoTurb_flag==1) THEN
                ! 09/02/2025 -  Addition of PTKE to Rocflu's Energy Equation
                @{USEPARTICLE(ppiclf_parts(i)%feedback%JE)}@ = (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)       &
                    * (                                                                             &
                    (fqsx+fvux+famx+liftx) * (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%X)}@)     +       &
                    (fqsy+fvuy+famy+lifty) * (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%Y)}@)     +       &
                    (fqsz+fvuz+famz+liftz) * (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%Z)}@)     +       &
                    qq )
            ELSE
                Rsg   = 0.0D0
                T_par = 0.0D0
            END IF ! pseudoTurb_flag
            ! 07/21/2025 - Thierry - Added Reynolds Subgrid Stress Feedback
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JRSG11)}@ = Rsg(1,1) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JRSG12)}@ = Rsg(1,2) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JRSG13)}@ = Rsg(1,3) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JRSG21)}@ = Rsg(2,1) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JRSG22)}@ = Rsg(2,2) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JRSG23)}@ = Rsg(2,3) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JRSG31)}@ = Rsg(3,1) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JRSG32)}@ = Rsg(3,2) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JRSG33)}@ = Rsg(3,3) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)

            @{USEPARTICLE(ppiclf_parts(i)%feedback%JTSG1)}@ = T_par(1) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JTSG2)}@ = T_par(2) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)
            @{USEPARTICLE(ppiclf_parts(i)%feedback%JTSG3)}@ = T_par(3) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)

        END IF ! Feedback flag

        ! Update volume fraction feedback quantities with feedback on or off
        @{USEPARTICLE(ppiclf_parts(i)%feedback%P_JPHIP)}@  = (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%JSPL)}@)
        @{USEPARTICLE(ppiclf_parts(i)%feedback%JPHIPD)}@   = (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%RHOP)}@)
        @{USEPARTICLE(ppiclf_parts(i)%feedback%JPHIPU)}@   = (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%X)}@)
        @{USEPARTICLE(ppiclf_parts(i)%feedback%JPHIPV)}@   = (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%Y)}@)
        @{USEPARTICLE(ppiclf_parts(i)%feedback%JPHIPW)}@   = (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%Z)}@)
        @{USEPARTICLE(ppiclf_parts(i)%feedback%JPHIPT)}@   = (@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@) * (@{USEPARTICLE(ppiclf_parts(i)%y%T)}@)


        ! Step 12: If stationary, don't move particles. Feedback can still be on
        ! though.
        !
        if (stationary .gt. 0) then
            if (stationary==1) then
                @{USEPARTICLE(ppiclf_parts(i)%ydot%pos%X)}@   = 0.0d0
                @{USEPARTICLE(ppiclf_parts(i)%ydot%pos%Y)}@   = 0.0d0
                @{USEPARTICLE(ppiclf_parts(i)%ydot%pos%Z)}@   = 0.0d0
                @{USEPARTICLE(ppiclf_parts(i)%ydot%vel%X)}@  = 0.0d0
                @{USEPARTICLE(ppiclf_parts(i)%ydot%vel%Y)}@  = 0.0d0
                @{USEPARTICLE(ppiclf_parts(i)%ydot%vel%Z)}@  = 0.0d0
                @{USEPARTICLE(ppiclf_parts(i)%ydot%T)}@   = 0.0d0
                @{USEPARTICLE(ppiclf_parts(i)%ydot%ang_vel%X)}@  = 0.0d0
                @{USEPARTICLE(ppiclf_parts(i)%ydot%ang_vel%Y)}@  = 0.0d0
                @{USEPARTICLE(ppiclf_parts(i)%ydot%ang_vel%Z)}@  = 0.0d0
            else
                call ppiclf_exittr('Unknown stationary flag$', 0.0d0, 0)
            endif
        elseif(stationary .lt. 0) then
            call ppiclf_user_unit_tests(i,iStage,famx,famy,famz)
            @{USEPARTICLE(ppiclf_parts(i)%ydot%pos%X)}@     = @{USEPARTICLE(ppiclf_parts(i)%y%vel%X)}@
            @{USEPARTICLE(ppiclf_parts(i)%ydot%pos%Y)}@     = @{USEPARTICLE(ppiclf_parts(i)%y%vel%Y)}@
            @{USEPARTICLE(ppiclf_parts(i)%ydot%pos%Z)}@     = @{USEPARTICLE(ppiclf_parts(i)%y%vel%Z)}@
            @{USEPARTICLE(ppiclf_parts(i)%ydot%vel%X)}@     = (@{USEPARTICLE(ppiclf_parts(i)%ydot%vel%X)}@)+famx
            @{USEPARTICLE(ppiclf_parts(i)%ydot%vel%Y)}@     = (@{USEPARTICLE(ppiclf_parts(i)%ydot%vel%Y)}@)+famy
            @{USEPARTICLE(ppiclf_parts(i)%ydot%vel%Z)}@     = (@{USEPARTICLE(ppiclf_parts(i)%ydot%vel%Z)}@)+famz
            @{USEPARTICLE(ppiclf_parts(i)%ydot%T)}@         = 0.0d0
            @{USEPARTICLE(ppiclf_parts(i)%ydot%ang_vel%X)}@ = 0.0d0
            @{USEPARTICLE(ppiclf_parts(i)%ydot%ang_vel%Y)}@ = 0.0d0
            @{USEPARTICLE(ppiclf_parts(i)%ydot%ang_vel%Z)}@ = 0.0d0
        endif



        !
        ! Step 13: Store forces

        ! TODO: update this conditional
        ! IF(PPICLF_LRP5 .EQ. 19) THEN
        !     {USEPARTICLE(i, rprop, FQSX)}@  = fqsx
        !     {USEPARTICLE(i, rprop, FQSY)}@  = fqsy
        !     {USEPARTICLE(i, rprop, FQSZ)}@  = fqsz
        !     {USEPARTICLE(i, rprop, FAMX)}@  = famx - rmass_add * ({USEPARTICLE(i, ydot, VX)}@)
        !     {USEPARTICLE(i, rprop, FAMY)}@  = famy - rmass_add * ({USEPARTICLE(i, ydot, VY)}@)
        !     {USEPARTICLE(i, rprop, FAMZ)}@  = famz - rmass_add * ({USEPARTICLE(i, ydot, VZ)}@)
        !     {USEPARTICLE(i, rprop, FAMBX)}@ = FamBinary(1)
        !     {USEPARTICLE(i, rprop, FAMBY)}@ = FamBinary(2)
        !     {USEPARTICLE(i, rprop, FAMBZ)}@ = FamBinary(3)
        !     {USEPARTICLE(i, rprop, FCX)}@   = fcx
        !     {USEPARTICLE(i, rprop, FCY)}@   = fcy
        !     {USEPARTICLE(i, rprop, FCZ)}@   = fcz
        !     {USEPARTICLE(i, rprop, FVUX)}@  = fvux
        !     {USEPARTICLE(i, rprop, FVUY)}@  = fvuy
        !     {USEPARTICLE(i, rprop, FVUZ)}@  = fvuz
        !     {USEPARTICLE(i, rprop, QQ)}@    = qq
        !     {USEPARTICLE(i, rprop, FPGX)}@  = fdpdx
        !     {USEPARTICLE(i, rprop, FPGY)}@  = fdpdy
        !     {USEPARTICLE(i, rprop, FPGZ)}@  = fdpdz
        ! END IF
        !
        ! Step 14: If debug mode is ON, calculate and print the max values.
        !          The user should not have this ON for production runs.
        !
        if (ppiclf_debug .ge. 1) then
            phimax = max(phimax,abs(rphip))

            fqsx_max = max(fqsx_max,abs(fqsx))
            fqsy_max = max(fqsy_max,abs(fqsy))
            fqsz_max = max(fqsz_max,abs(fqsz))
            fqs_mag  = max(fqs_mag,sqrt(fqsx*fqsx+fqsy*fqsy+fqsz*fqsz))

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
            fam_mag  = max(fam_mag,sqrt(famx*famx+famy*famy+famz*famz))

            fdpdx_max = max(fdpdx_max,abs(fdpdx))
            fdpdy_max = max(fdpdy_max,abs(fdpdy))
            fdpdz_max = max(fdpdz_max,abs(fdpdz))
            fdp_mag   = max(fdp_mag,sqrt(fdpdx*fdpdx+fdpdy*fdpdy+fdpdz*fdpdz))

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
                        write(7010,*) i,ppiclf_time,rmass,vmag,rhof,dp,rep,rphip,rphif,rmachp,rhop,rhoMixt,reyL,rmu,rnu,rkappa
                    endif
                    if (i==ppiclf_npart) then
                        write(7011,*) i,ppiclf_time,rmass,vmag,rhof,dp,rep,rphip,rphif,rmachp,rhop,rhoMixt,reyL,rmu,rnu,rkappa
                    endif
                endif
            endif

        endif ! ppiclf_debug .ge. 1
            
        ! write out for debug
        if (ppiclf_debug==3) then
            if (ppiclf_nid==0 .and. iStage==1) then
                if (mod(idebug,1)==0) then
                    if (i<=5) then
                        write(7020+i,*) i, ppiclf_time, rhof,   &
                            (@{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%X)}@),        &
                            (@{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%Y)}@),        &
                            (@{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%Z)}@),        &
                            (@{USEPARTICLE(ppiclf_parts(i)%ydot%vel%X)}@),     &
                            (@{USEPARTICLE(ppiclf_parts(i)%ydot%vel%Y)}@),     &
                            (@{USEPARTICLE(ppiclf_parts(i)%ydot%vel%Z)}@),     &
                            (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%X)}@),        &
                            (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%Y)}@),        &
                            (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%Z)}@),        &
                            (@{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%X)}@),        &
                            (@{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%Y)}@),        &
                            (@{USEPARTICLE(ppiclf_parts(i)%y%ang_vel%Z)}@)

                        write(7040+i,*) i, ppiclf_time,                                         &
                            @{USEPARTICLE(ppiclf_parts(i)%rprop%SDR)}@,    &   ! Du/Dt
                            @{USEPARTICLE(ppiclf_parts(i)%rprop%SDO)}@        ! DOmega/Dt

                        write(7050+i,*) i, ppiclf_time, fqs_mag,fam_mag,fdp_mag,fc_mag,tau_max

                        write(7060+i,*) i, ppiclf_time, fcx,fcy,fcz,liftx,lifty,liftz,taux,tauy,tauz
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
            SDrho = (@{USEPARTICLE(ppiclf_parts(i)%rprop%RHSR)}@)                         +   & 
            (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%X)}@) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%PGC%X)}@)   +   & 
            (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%Y)}@) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%PGC%Y)}@)   +   & 
            (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%Z)}@) * (@{USEPARTICLE(ppiclf_parts(i)%rprop%PGC%Z)}@)
            
            ! material derivative is phi weighted in Rocflu
            ! drho/dt
            SDrho = SDrho / (rphif)  
            vgradrhog = vx * (@{USEPARTICLE(ppiclf_parts(i)%rprop%RHOG%X)}@) +  &
                        vy * (@{USEPARTICLE(ppiclf_parts(i)%rprop%RHOG%Y)}@) +  &
                        vz * (@{USEPARTICLE(ppiclf_parts(i)%rprop%RHOG%Z)}@)
      
            ! Fluid density
            rhof   = @{USEPARTICLE(ppiclf_parts(i)%rprop%RHOF)}@

            vx = (@{USEPARTICLE(ppiclf_parts(i)%rprop%U%X)}@) - (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%X)}@)
            vy = (@{USEPARTICLE(ppiclf_parts(i)%rprop%U%Y)}@) - (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%Y)}@)
            vz = (@{USEPARTICLE(ppiclf_parts(i)%rprop%U%Z)}@) - (@{USEPARTICLE(ppiclf_parts(i)%y%Vel%Z)}@)
            ug = (@{USEPARTICLE(ppiclf_parts(i)%rprop%U%X)}@)
            vg = (@{USEPARTICLE(ppiclf_parts(i)%rprop%U%Y)}@)
            wg = (@{USEPARTICLE(ppiclf_parts(i)%rprop%U%Z)}@)
            ! Unary added mass solves rho^g d(u^p)/dt implicitly
            ! Binary added mass solves it explicitly and not implicitly
            ! WDOTX = D(rho^g u^g)/Dt - d(rho^g u^p)/dt)
            ! X-acceleration
            @{USEPARTICLE(ppiclf_parts(i)%rprop%WDOT%X)}@ = vx*SDrho + rhof*(@{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%X)}@) + ug*vgradrhog - rhof*(@{USEPARTICLE(ppiclf_parts(i)%ydot%Vel%X)}@)
          
            ! Y-acceleration
            @{USEPARTICLE(ppiclf_parts(i)%rprop%WDOT%Y)}@ = vy*SDrho + rhof*(@{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%Y)}@) + vg*vgradrhog - rhof*(@{USEPARTICLE(ppiclf_parts(i)%ydot%Vel%Y)}@)

            ! Z-acceleration
            @{USEPARTICLE(ppiclf_parts(i)%rprop%WDOT%Z)}@ = vz*SDrho + rhof*(@{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%Z)}@) + wg*vgradrhog - rhof*(@{USEPARTICLE(ppiclf_parts(i)%ydot%Vel%Z)}@)
          
            ! write out for debug
            if (ppiclf_debug==2) then
                if (ppiclf_nid==0 .and. iStage==1) then
                    if (mod(idebug,10)==0) then
                        if (i<=3) then
                            write(7020+i,*) i, ppiclf_time, rhof,   &
                                @{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%X)}@,          &
                                @{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%Y)}@,          &
                                @{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%Z)}@,          &
                                @{USEPARTICLE(ppiclf_parts(i)%ydot%vel%X)}@,       &
                                @{USEPARTICLE(ppiclf_parts(i)%ydot%vel%Y)}@,       &
                                @{USEPARTICLE(ppiclf_parts(i)%ydot%vel%Z)}@,       &
                                @{USEPARTICLE(ppiclf_parts(i)%y%vel%X)}@,          &
                                @{USEPARTICLE(ppiclf_parts(i)%y%vel%Y)}@,          &
                                @{USEPARTICLE(ppiclf_parts(i)%y%vel%Z)}@

                            write(7030+i,*) i, ppiclf_time, @{USEPARTICLE(ppiclf_parts(i)%rprop%WDOT)}@ 

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
#ifdef PERF
    tfinal = MPI_WTIME()
    PPICLF_TParticleParticleModels = PPICLF_TParticleParticleModels + (tfinal - tstart)
#endif
    ! ----------------------------------------------------------------------

    return
end procedure ppiclf_user_SetYdot
    
subroutine ppiclf_user_unit_tests(i,iStage,famx,famy,famz)
    integer*4 i, iStage
    real*8 famx, famy, famz
end subroutine ppiclf_user_unit_tests

subroutine ppiclf_user_debug
end subroutine ppiclf_user_debug

end submodule ppiclf_user_SetYdot_imp