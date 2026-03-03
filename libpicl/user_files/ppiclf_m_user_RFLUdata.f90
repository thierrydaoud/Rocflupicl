#include "PPICLF_STD.h"
module ppiclf_m_user_RFLUdata

    implicit none
    save
    
    !
    ! General useage
    !
    
    ! COMMON Name: RFLU_ppiclf

    integer*4 :: stationary, qs_flag, am_flag, pg_flag,             &
        collisional_flag, heattransfer_flag, feedback_flag,         &
        qs_fluct_flag, ppiclf_debug, rmu_flag,                      &
        rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,      &
        qs_fluct_filter_adapt_flag,                                 &
        ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,  &
        sbNearest_flag, burnrate_flag, flow_model, pseudoTurb_flag
    real*8 :: rmu_ref, tref, suth, ksp, erest

end module ppiclf_m_user_RFLUdata