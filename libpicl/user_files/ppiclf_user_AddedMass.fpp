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
#:include "PPICLF_PARTMACROS.fypp"
#include "PPICLF_STD.h"
submodule (ppiclf_m_user_ForceModels) ppiclf_m_user_ForceModels_AddedMass
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
    use ppiclf_data, only: ppiclf_nid, ppiclf_np
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
    
    use ppiclf_user_AM_functions, only : IA_analytical, IA_numerical, II_analytical, II_numerical

    use ppiclf_op, only: ppiclf_exittr
    implicit none


    contains

    module procedure ppiclf_user_AM_Parmar
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
        rmass_add = rhof*(@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@)*rcd_am

        !NEW Added mass, using how rocflu does it
        !1st Derivative, substantial how rocflu does it
        SDrho = @{USEPARTICLE(ppiclf_parts(i)%rprop%RHSR)}@                                             &
            + @{USEPARTICLE(ppiclf_parts(i)%y%Vel%X)}@ * @{USEPARTICLE(ppiclf_parts(i)%rprop%PGC%X)}@   &
            + @{USEPARTICLE(ppiclf_parts(i)%y%Vel%Y)}@ * @{USEPARTICLE(ppiclf_parts(i)%rprop%PGC%Y)}@   &
            + @{USEPARTICLE(ppiclf_parts(i)%y%Vel%Z)}@ * @{USEPARTICLE(ppiclf_parts(i)%rprop%PGC%Z)}@


        ! 03/11/2025 - Thierry - substantial derivative from Rocflu is 
        !              weighted by \phi^g.
        ! d(rho^g phi^g)/dt = rho^g * d(phi^g)/dt + phi^g * d(rho^g)/dt
        !                   = phi^g * d(rho^g)/dt
        !  
        !     d(rho^g)/dt   = SDrho = d(rho phi^g)/dt / phi^g
        SDrho = SDrho / (rphif) 

        ! 03/23/2025 - TLJ - added extra term involving grad(rhog)
        vgradrho =  vx*(@{USEPARTICLE(ppiclf_parts(i)%rprop%RHOG%X)}@) + &
                    vy*(@{USEPARTICLE(ppiclf_parts(i)%rprop%RHOG%Y)}@) + &
                    vz*(@{USEPARTICLE(ppiclf_parts(i)%rprop%RHOG%Z)}@)

        ug = @{USEPARTICLE(ppiclf_parts(i)%rprop%U%X)}@
        vg = @{USEPARTICLE(ppiclf_parts(i)%rprop%U%Y)}@
        wg = @{USEPARTICLE(ppiclf_parts(i)%rprop%U%Z)}@

        famx = rcd_am*@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@ * (vx*SDrho + rhof*(@{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%X)}@) + ug*vgradrho)

        famy = rcd_am*@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@ * (vy*SDrho + rhof*(@{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%Y)}@) + vg*vgradrho)

        famz = rcd_am*@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@ * (vz*SDrho + rhof*(@{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%Z)}@) + wg*vgradrho)


        ! if (1==2) then
        !     ! This is Ling's 2012 formulation where he replaced
        !     !   D(rhog*ug) with -grad(pg)
        !     famx_Ling = rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i)*(-ppiclf_rprop(PPICLF_R_JDPDX,i) - ppiclf_y(PPICLF_JVX,i)*SDrho)

        !     ! Original version
        !     ! /home/tlj/Codes_Rocflu/Rocflu_picl_tlj/ppiclf/source/ppiclf_user.f
        !     ! In the original version JSDRX was assumed to be D(rhog*ug)/Dt,
        !     !    but this was before we realized that the conserved Rocflu
        !     !    variable was phig*rhog*ug
        !     famx_Brad = rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i) *
        !     >   (ppiclf_rprop(PPICLF_R_JSDRX,i) - ppiclf_y(PPICLF_JVX,i)*SDrho)

        !     ! writing data only for the median particle
        !     if((ppiclf_iprop(5,i).eq.29.0) .and. (ppiclf_iprop(6,i).eq.0.0).and. (ppiclf_iprop(7,i).eq.151.0)) then
            
        !         open(unit=20,file='fort.20',position='append') 
        !         write(20,*) ppiclf_time,famx-rmass_add*ppiclf_ydot(PPICLF_JVX,i),
        !         >            rcd_am, ppiclf_rprop(PPICLF_R_JVOLP,i),
        !         >            rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i),
        !         >            vx, SDrho, vx*SDrho,
        !         >            rhof, ppiclf_rprop(PPICLF_R_JSDRX,i),
        !         >            rhof*ppiclf_rprop(PPICLF_R_JSDRX,i),
        !         >            ug, vgradrho,
        !         >            ug*vgradrho,
        !         >            famx_Ling-rmass_add*ppiclf_ydot(PPICLF_JVX,i),
        !         >            famx_Brad-rmass_add*ppiclf_ydot(PPICLF_JVX,i)
        !         flush(20)
        !     endif
        ! endif


        return
    end procedure ppiclf_user_AM_Parmar
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

    module procedure ppiclf_user_AM_Briney_Unary

        integer j, k, l, n, jj
        real*8 rad
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

        rmass_add = rhof*(@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@)*rcd_am

        !NEW Added mass, using how rocflu does it
        !1st Derivative, substantial how rocflu does it
        SDrho = @{USEPARTICLE(ppiclf_parts(i)%rprop%RHSR)}@                                             &
            + @{USEPARTICLE(ppiclf_parts(i)%y%Vel%X)}@ * @{USEPARTICLE(ppiclf_parts(i)%rprop%PGC%X)}@   &
            + @{USEPARTICLE(ppiclf_parts(i)%y%Vel%Y)}@ * @{USEPARTICLE(ppiclf_parts(i)%rprop%PGC%Y)}@   &
            + @{USEPARTICLE(ppiclf_parts(i)%y%Vel%Z)}@ * @{USEPARTICLE(ppiclf_parts(i)%rprop%PGC%Z)}@
        ! material derivative is phi weighted in Rocflu
        ! drho/dt
        SDrho = SDrho / (rphif) 

        ! 03/23/2025 - TLJ - added extra term involving grad(rhog)
        vgradrho =  vx*(@{USEPARTICLE(ppiclf_parts(i)%rprop%RHOG%X)}@) +    &
                    vy*(@{USEPARTICLE(ppiclf_parts(i)%rprop%RHOG%Y)}@) +    &
                    vz*(@{USEPARTICLE(ppiclf_parts(i)%rprop%RHOG%Z)}@)

        ug = @{USEPARTICLE(ppiclf_parts(i)%rprop%U%X)}@
        vg = @{USEPARTICLE(ppiclf_parts(i)%rprop%U%Y)}@
        wg = @{USEPARTICLE(ppiclf_parts(i)%rprop%U%Z)}@

        ! Take care of volume in Binary subroutine
        famx = rcd_am*(vx*SDrho + rhof*(@{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%X)}@) + ug*vgradrho)

        famy = rcd_am*(vy*SDrho + rhof*(@{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%Y)}@) + vg*vgradrho)

        famz = rcd_am*(vz*SDrho + rhof*(@{USEPARTICLE(ppiclf_parts(i)%rprop%SDR%Z)}@) + wg*vgradrho)

        ! Multiply by neighbors here for storing
        FamUnary(1) = famx*(@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@)
        FamUnary(2) = famy*(@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@)
        FamUnary(3) = famz*(@{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@)

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
                    write(7050+i,*) i, ppiclf_nid, ppiclf_np, ppiclf_time,   & ! 0-3
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%RHSR)}@,          & ! 4
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%PGC%X)}@,         & ! 5
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%PGC%Y)}@,         & ! 6
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%PGC%Z)}@,         & ! 7
                        SDrho,                                                & ! 8
                        rhof, rphip, rmachp, SDrho,                           & ! 9-12
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%WDOT%X)}@,        & ! 13
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%WDOT%Y)}@,        & ! 14
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%WDOT%Z)}@,        & ! 15
                        Wdot_neighbor_mean(1:3), nneighbors,                  & ! 16-19
                        famx, famy, famz,                                     & ! 20-22
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@*Fam(1),   & ! 23
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@*Fam(2),   & ! 24
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@*Fam(3)      ! 25
                    flush(7051)
                    flush(7052)
                    flush(7053)
                    flush(7054)
                    flush(7055)
                end if  
            end if  
        end if  

        return
    end procedure ppiclf_user_AM_Briney_Unary
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
    module procedure ppiclf_user_AM_Briney_Binary
        integer j, k, l, n, jj
        real*8 rad
        real*8 dr_max
        real*8 IA, II
        real*8 alpha

        ! ! Declare functions
        ! real*8 IA_analytical, IA_numerical, II_analytical, II_numerical
        !
        ! Code:
        !
        ! particle radius
        rad = @{USEPARTICLE(ppiclf_parts(i)%rprop%DP)}@ * 0.5d0
        
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
                print*, "***WARNING*** - NUMERICAL FUNCTIONS USED IN ADDED MASS"
            endif
            IA = IA_numerical(dr_max, rad, alpha)
            II = II_numerical(dr_max, rad, alpha)
        end if

        do j=1,3
            Fam(j) = Fam(j) + IA*(@{USEPARTICLE(ppiclf_parts(i)%rprop%wdot, skipIndex)}@(j))  ! added mass
            Fam(j) = Fam(j) + II*Wdot_neighbor_mean(j) / nneighbors ! induced added mass
        end do

        ! multiply by volume before adding unary term
        ! doing so here implies that the particle volume is
        ! the same for all particles; i.e., monodisperse packs
        do j=1,3
            Fam(j) = @{USEPARTICLE(ppiclf_parts(i)%rprop%VOLP)}@*Fam(j) 
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
                    write(7060+i,*) i,iStage, ppiclf_time,                  & ! 0-2
                        rphip, rmachp,                                      & ! 3-4
                        IA, II,                                             & ! 5-6
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%WDOT%X)}@,      & ! 7
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%WDOT%Y)}@,      & ! 8
                        @{USEPARTICLE(ppiclf_parts(i)%rprop%WDOT%Z)}@,      & ! 9
                        Wdot_neighbor_mean(1:3), nneighbors,                & ! 10-13
                        famx, famy, famz,                                   & ! 14-16
                        Fam(1), Fam(2), Fam(3)                                ! 17-19
                    flush(7061)
                    flush(7062)
                    flush(7063)
                    flush(7064)
                    flush(7065)
                end if
            end if
        end if
        

        return
    end procedure ppiclf_user_AM_Briney_Binary

end submodule ppiclf_m_user_ForceModels_AddedMass