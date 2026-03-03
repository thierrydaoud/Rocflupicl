!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine to set user-defined values at time t=0
!
!-----------------------------------------------------------------------
!

submodule (ppiclf_user) ppiclf_user_InitZero_imp
    use ppiclf_data, only : PPICLF_TIMEBH, PPICLF_DRUDTPLAG, PPICLF_DRUDTMIXT
    implicit none
    contains

    module procedure ppiclf_user_InitZero
        !
        ! Internal:
        !
        !
        ! Code:
        !
        ppiclf_TimeBH = 0.0d0

        ppiclf_drudtMixt = 0.0d0
        ppiclf_drudtPlag = 0.0d0
        return
    end procedure ppiclf_user_InitZero
    
end submodule ppiclf_user_InitZero_imp

