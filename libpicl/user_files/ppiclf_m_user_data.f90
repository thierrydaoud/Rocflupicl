#include "PPICLF_STD.h"
module ppiclf_m_user_data

    implicit none
    save

    !
    ! General useage
    !
    integer*4 :: stationary, qs_flag, am_flag, pg_flag,             &
        collisional_flag, heattransfer_flag, feedback_flag,         &
        qs_fluct_flag, ppiclf_debug, rmu_flag,                      &
        rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,      &
        qs_fluct_filter_adapt_flag,                                 &
        ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,  &
        sbNearest_flag, burnrate_flag, flow_model, pseudoTurb_flag

    real*8 :: rmu_ref, tref, suth, ksp, erest


    real*8 rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,             &
        rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac,               &
        vx,vy,vz,                                                   &
        rcp_part,rpr,                                               &
        phi, mp, re, rem


    !
    ! For misc values
    !
    real*8, parameter ::  OneThird = 1.0d0/3.0d0

 

    !
    ! For ppiclf_user_Fluctuations.f
    !
    integer*4 icpmean
    real*8 upmean, vpmean, wpmean, phipmean
    real*8 u2pmean, v2pmean, w2pmean


    real*8 UnifRnd(6), Rsg(3,3), T_par(3)


    !
    ! For ppiclf_user_debug.f
    !
    real*8 phimax,                                      &
        fqsx_max,fqsy_max,fqsz_max,                     &
        famx_max,famy_max,famz_max,                     &
        fdpdx_max,fdpdy_max,fdpdz_max,                  &
        fcx_max,fcy_max,fcz_max,                        &
        umean_max,vmean_max,wmean_max,                  &
        fqs_mag,fam_mag,fdp_mag,fc_mag,                 &
        fqsx_fluct_max,fqsy_fluct_max,fqsz_fluct_max,   &
        fqsx_total_max,fqsy_total_max,fqsz_total_max,   &
        fvux_max,fvuy_max,fvuz_max,                     &
        qq_max,tau_max,lift_max


    !
    ! For ppiclf_user_AddedMass.f
    !
    integer*4 nneighbors
    real*8 Fam(3), FamUnary(3), FamBinary(3),Wdot_neighbor_mean(3), R_pair(6,6)


    !
    ! For ppiclf_solve_InitAngularPeriodic
    !
        
    ! MOVED TO ppiclf_data.f90

end module ppiclf_m_user_data