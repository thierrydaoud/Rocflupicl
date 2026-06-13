#include "PPICLF_USER.h"
#include "PPICLF_STD.h"
c Particle options
      LOGICAL PPICLF_RESTART, PPICLF_OVERLAP, PPICLF_LCOMM
     >       ,PPICLF_LINIT, PPICLF_LINTP, PPICLF_LPROJ
     >       ,PPICLF_LSUBSUBBIN,PPICLF_EQUALDOMAIN(3)
     >       ,PPICLF_LINPERIODIC(3), PPICLF_REMOVE_PARTICLE
     >       ,PPICLF_BINCHANGED, PPICLF_PRINTBINVTU
     >       ,PPICLF_READYTOSOLVE,PPICLF_BINORDERSET
      COMMON /PPICLF_OPT_PARAM_L/ PPICLF_RESTART, PPICLF_OVERLAP
     >                           ,PPICLF_LCOMM, PPICLF_LINIT
     >                           ,PPICLF_LINTP
     >                           ,PPICLF_LPROJ 
     >                           ,PPICLF_LSUBSUBBIN
     >                           ,PPICLF_EQUALDOMAIN
     >                           ,PPICLF_LINPERIODIC
     >                           ,PPICLF_REMOVE_PARTICLE
     >                           ,PPICLF_BINCHANGED
     >                           ,PPICLF_PRINTBINVTU
     >                           ,PPICLF_READYTOSOLVE
     >                           ,PPICLF_BINORDERSET
      DATA PPICLF_LCOMM /.false./
      DATA PPICLF_RESTART /.false./

      INTEGER*4 PPICLF_NDIM, PPICLF_IMETHOD 
     >         ,PPICLF_NGRIDS, PPICLF_CYCLE, PPICLF_IOSTEP
     >         ,PPICLF_IENDIAN, PPICLF_IWALLM
      COMMON /PPICLF_OPT_PARAM_I/ PPICLF_NDIM, PPICLF_IMETHOD
     >                           ,PPICLF_NGRIDS
     >                           ,PPICLF_CYCLE, PPICLF_IOSTEP
     >                           ,PPICLF_IENDIAN, PPICLF_IWALLM
      REAL*8 PPICLF_FILTER(3), PPICLF_RK3COEF(3,3),
     >       PPICLF_DT, PPICLF_TIME, PPICLF_NNDIST,
     >       PPICLF_INTERP_DCHK(3), PPICLF_TOTNNDIST(PPICLF_LPART)
      REAL*8 PPICLF_RK3ARK(3)
      COMMON /PPICLF_OPT_PARAM_R/ PPICLF_FILTER
     >                           ,PPICLF_RK3COEF, PPICLF_DT
     >                           ,PPICLF_TIME, PPICLF_NNDIST
     >                           ,PPICLF_RK3ARK, PPICLF_TOTNNDIST
     >                           ,PPICLF_INTERP_DCHK
