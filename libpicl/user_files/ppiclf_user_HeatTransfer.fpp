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
