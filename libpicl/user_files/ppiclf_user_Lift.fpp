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
