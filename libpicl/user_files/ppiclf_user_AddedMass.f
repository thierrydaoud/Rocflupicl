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
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag, flow_model
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag, flow_model
      integer*4 i, iStage
      real*8 famx, famy, famz, rmass_add
      real*8 rcd_am
      real*8 SDrho
      real*8 ug,vg,wg,vgradrho

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

      ! Sangani's volume fraciton correction for dilute random arrays
      ! Capping volume fraction at 0.5
      rcd_am = rcd_am*(1.0+3.32*MIN(rphip,0.5))

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
! Implementing Added Mass Algorithm from S.Briney (2024)
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
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag, flow_model
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag, flow_model
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
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag, flow_model
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag, flow_model
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
      
      ! ppiclf_nndist is neighbor width - user defined
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
