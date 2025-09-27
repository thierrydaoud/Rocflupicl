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
      subroutine ppiclf_user_QS_Parmar(i,beta,cd)
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

      integer*4 i
      real*8 gamma,mp,phi,re
      real*8 rcd1,rmacr,rcd_mcr,rcd_std,rmach_rat,rcd_M1,
     >   rcd_M2,C1,C2,C3,f1M,f2M,f3M,lrep,factor,cd,beta

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
      !phi_corr = (1.0+5.94*MIN(rphip,0.5))

      cd = (24.0/re)*rcd1*(1.+2.*phi)/((1.-phi)**3)

      beta = rcd1*3.0*rpi*rmu*dp

      beta = beta*(1.+2.*phi)/((1.-phi)**3)


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
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
      subroutine ppiclf_user_QS_Osnes(i,beta,cd)
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

      integer*4 i
      real*8 gamma,mp,phi,re,Knp,fKn,CD1,s,JM,CD2,
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
