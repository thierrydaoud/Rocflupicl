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
