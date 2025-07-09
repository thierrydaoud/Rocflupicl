!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for projection map
!
!
! This subroutine maps certain ppiclf terms onto the fluid mesh
!
!    This routine is needed for
!       (i)  feedback_flag = 1 (rocpicl/PICL_TEMP_Runge.F90)
!       (ii) # PROBE is in the *.inp file (libfloflu/WriteProbe.F90)

!-----------------------------------------------------------------------
!
      SUBROUTINE ppiclf_user_MapProjPart(i)
!
      IMPLICIT NONE
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4 i
!
! Code:
!
!
       ppiclf_feedbk(PPICLF_P_JPHIP,i) = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >    *ppiclf_rprop(PPICLF_R_JSPL,i)
       ppiclf_feedbk(PPICLF_P_JPHIPD,i) =
     >    ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_rprop(PPICLF_R_JRHOP,i)
       ppiclf_feedbk(PPICLF_P_JPHIPU,i) = 
     >    ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JVX,i)
       ppiclf_feedbk(PPICLF_P_JPHIPV,i) =
     >    ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JVY,i)
       ppiclf_feedbk(PPICLF_P_JPHIPW,i) = 
     >    ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JVZ,i)
       ppiclf_feedbk(PPICLF_P_JPHIPT,i) =
     >    ppiclf_rprop(PPICLF_R_JVOLP,i)*ppiclf_y(PPICLF_JT,i)

      RETURN
      END
