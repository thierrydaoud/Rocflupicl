module ppiclf_m_user_ForceModels
    implicit none
    interface
        module subroutine ppiclf_user_BR_driver(i,iStage,burnrate_model,qq,mdot_me,mdot_ox)
            integer*4 i, iStage, burnrate_model
            real*8 qq
            real*8 mdot_me, mdot_ox
        end subroutine ppiclf_user_BR_driver
        module subroutine ppiclf_user_Lift_driver(i,iStage,liftx,lifty,liftz)
            integer*4 i,iStage
            real*8 liftx, lifty, liftz
        end subroutine ppiclf_user_Lift_driver
        module subroutine ppiclf_user_Torque_driver(i,iStage,taux,tauy,tauz,taux_hydro,tauy_hydro,tauz_hydro)
            integer*4 i,iStage
            real*8 taux, tauy, tauz
            real*8 taux_hydro, tauy_hydro, tauz_hydro
        end subroutine ppiclf_user_Torque_driver
        module subroutine ppiclf_user_HT_driver(i,qq)
            integer*4 i
            real*8 qq
        end subroutine ppiclf_user_HT_driver

        module subroutine ppiclf_user_QS_Parmar(i,beta, cd)
            integer*4 i
            real*8 beta, cd
        end subroutine ppiclf_user_QS_Parmar
        module subroutine ppiclf_user_QS_ModifiedParmar(i,beta, cd)
            integer*4 i
            real*8 beta, cd
        end subroutine ppiclf_user_QS_ModifiedParmar
        module subroutine ppiclf_user_QS_Osnes(i,beta, cd)
            integer*4 i
            real*8 beta, cd
        end subroutine ppiclf_user_QS_Osnes
        module subroutine ppiclf_user_QS_Gidaspow(i,beta, cd)
            integer*4 i
            real*8 beta, cd
        end subroutine ppiclf_user_QS_Gidaspow

        module subroutine ppiclf_user_QS_fluct_Lattanzi(i,iStage,fqs_fluct)
            integer*4 i, iStage
            real*8 fqs_fluct(3)
        end subroutine ppiclf_user_QS_fluct_Lattanzi

        module subroutine ppiclf_user_QS_fluct_Osnes(i,iStage,fqs_fluct,xi_par,xi_perp,xi_T,fqsx,fqsy,fqsz)
            !
            ! Input:
            integer*4 i, iStage
            real*8 fqsx, fqsy, fqsz
            !
            ! Output:
            real*8 xi_par, xi_perp, xi_T
            real*8 fqs_fluct(3)
        end subroutine ppiclf_user_QS_fluct_Osnes

        module subroutine ppiclf_user_VU_Rocflu(i,iStage,fvux,fvuy,fvuz)
            integer*4 i, iStage
            real*8 fvux,fvuy,fvuz
        end subroutine ppiclf_user_VU_Rocflu
        module subroutine ppiclf_user_ShiftUnsteadyData
        end subroutine ppiclf_user_ShiftUnsteadyData

        module subroutine ppiclf_user_UpdatePlag(i)
            integer*4 i
        end subroutine ppiclf_user_UpdatePlag
        module subroutine ppiclf_user_prop2plag
        end subroutine ppiclf_user_prop2plag

        module subroutine ppiclf_user_plag2prop
        end subroutine ppiclf_user_plag2prop

        ! Added mass functions
        module subroutine ppiclf_user_AM_Parmar(i,iStage,famx,famy,famz,rmass_add)
            integer*4 i, iStage
            real*8 famx, famy, famz, rmass_add
        end subroutine ppiclf_user_AM_Parmar
        module subroutine ppiclf_user_AM_Briney_Unary(i,iStage, famx,famy,famz,rmass_add)
            integer i
            integer*4 iStage
            real*8 famx, famy, famz, rmass_add
        end subroutine ppiclf_user_AM_Briney_Unary

        module subroutine ppiclf_user_AM_Briney_Binary(i,iStage,famx,famy,famz,rmass_add)
            integer i
            integer*4 iStage
            real*8 famx, famy, famz, rmass_add
        end subroutine ppiclf_user_AM_Briney_Binary
    end interface
end module ppiclf_m_user_ForceModels