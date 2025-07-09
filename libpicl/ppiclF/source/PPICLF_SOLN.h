#include "PPICLF_USER.h"
#include "PPICLF_STD.h"
! Computational particles
      REAL*8 PPICLF_Y     (PPICLF_LRS ,PPICLF_LPART)  ! Solution
     >      ,PPICLF_YDOT  (PPICLF_LRS ,PPICLF_LPART)  ! Total solution RHS
     >      ,PPICLF_YDOTC (PPICLF_LRS ,PPICLF_LPART)  ! Coupled solution RHS
     >      ,PPICLF_RPROP (PPICLF_LRP ,PPICLF_LPART)  ! Real particle properties
     >      ,PPICLF_RPROP2(PPICLF_LRP2,PPICLF_LPART)  ! Secondary real particle properties
     >      ,PPICLF_RPROP3(PPICLF_LRP3,PPICLF_LPART)  ! Third real particle properties
     >      ,PPICLF_FEEDBK(PPICLF_LRP_PRO,PPICLF_LPART) !Feedback particle terms
      COMMON /PPICLF_SLN_CURRENT_R/ PPICLF_Y
     >                             ,PPICLF_YDOT
     >                             ,PPICLF_YDOTC
     >                             ,PPICLF_RPROP
     >                             ,PPICLF_RPROP2
     >                             ,PPICLF_RPROP3
     >                             ,PPICLF_FEEDBK

      INTEGER*4 PPICLF_IPROP(PPICLF_LIP,PPICLF_LPART) ! Integer particle properties
      COMMON /PPICLF_SLN_CURRENT_I/  PPICLF_IPROP

      COMMON /PPICLF_SLN_CURRENT_N/ PPICLF_NPART
      INTEGER*4 PPICLF_NPART

! Previous time step solutions, may grow later
      REAL*8 PPICLF_Y1(PPICLF_LRS, PPICLF_LPART)
      COMMON /PPICLF_SLN_PREVIOUS_R/  PPICLF_Y1

! Previous time step solutions, may grow later
      REAL*8 PPICLF_TIMEBH(PPICLF_VU)
      REAL*8 PPICLF_DRUDTPLAG(3,PPICLF_VU,PPICLF_LPART)
      REAL*8 PPICLF_DRUDTMIXT(3,PPICLF_VU,PPICLF_LPART)
      COMMON /PPICLF_SLN_UNSTEADY/  PPICLF_TIMEBH
     >                             ,PPICLF_DRUDTPLAG
     >                             ,PPICLF_DRUDTMIXT

