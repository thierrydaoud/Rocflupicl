#include "PPICLF_STD.h" 
! PPICLF_STD has been modified to include PPICLF_USER.h

module ppiclf_data
    implicit none

    save
    
    ! originally PPICLF_SOLN.h
    ! computational particl data
    ! REAL*8 PPICLF_Y     (PPICLF_LRS ,PPICLF_LPART)  &   ! Solution
    !     ,PPICLF_YDOT  (PPICLF_LRS ,PPICLF_LPART)    &   ! Total solution RHS
    !     ,PPICLF_YDOTC (PPICLF_LRS ,PPICLF_LPART)    &   ! Coupled solution RHS
    !     ,PPICLF_RPROP (PPICLF_LRP ,PPICLF_LPART)    &   ! Real particle properties
    !     ,PPICLF_RPROP2(PPICLF_LRP2,PPICLF_LPART)    &   ! Secondary real particle properties
    !     ,PPICLF_RPROP3(PPICLF_LRP3,PPICLF_LPART)    &   ! Third real particle properties
    !     ,PPICLF_RPROP4(PPICLF_LRP4,PPICLF_LPART)    &   ! Fourth real particle properties
    !     ,PPICLF_RPROP5(PPICLF_LRP5,PPICLF_LPART)    &   ! Fifth real particle properties
    !     ,PPICLF_FEEDBK(PPICLF_LRP_PRO,PPICLF_LPART)     !Feedback particle terms

    ! INTEGER*4 PPICLF_IPROP(PPICLF_LIP,PPICLF_LPART) ! Integer particle properties

    INTEGER*4 PPICLF_NPART 

    ! Previous time step solutions, may grow later
    ! REAL*8 PPICLF_Y1(PPICLF_LRS, PPICLF_LPART)

    ! Previous time step solutions, may grow later
    REAL*8 PPICLF_TIMEBH(PPICLF_VU)
    REAL*8 PPICLF_DRUDTPLAG(3,PPICLF_VU,PPICLF_LPART)
    REAL*8 PPICLF_DRUDTMIXT(3,PPICLF_VU,PPICLF_LPART)

    

    ! Originally PPICLF_GRID.h

    ! Grid
    REAL*8  PPICLF_XDRANGE(2,3)

    REAL*8 PPICLF_PRO_FLD(PPICLF_LEE,PPICLF_LRP_PRO),       &
        PPICLF_PRO_FLD_PICL(PPICLF_LRP_PRO,PPICLF_LEE),     &
        PPICLF_INT_FLD_INPUT(PPICLF_LEE,PPICLF_LRP_INT),    &
        PPICLF_INT_FLD(PPICLF_LRP_INT,PPICLF_LEE),          &
        PPICLF_FLUID_GRID (7,PPICLF_LEE),                   &
        PPICLF_PICL_GRID  (7,PPICLF_LEE),                   &
        PPICLF_PART2CELL_DIST(PPICLF_LPART,27)


                            
    INTEGER*4 PPICLF_CELL_MAP(PPICLF_LRMAX,PPICLF_LEE),     &
        PPICLF_CELL_MAP_ORIG(PPICLF_LRMAX,PPICLF_LEE),      &
        PPICLF_CELL_MAP_INTERP(PPICLF_LRMAX,PPICLF_LEE),    &
        PPICLF_CELL_MAP_PROJ(PPICLF_LRMAX,PPICLF_LEE),      &
        PPICLF_NCELLS_FV2PICL, PPICLF_NCELLS_FV2PICL_ORIG,  &
        PPICLF_NCELLS_INTERP, PPICLF_NCELLS_PROJ,           &
        PPICLF_PART2CELL_MAP(PPICLF_LPART,27),              &
        PPICLF_NPART2CELL(PPICLF_LPART)

    INTEGER*4 PPICLF_NFVCELLS

    INTEGER*4 PPICLF_INT_ICNT, PPICLF_INT_MAP(PPICLF_LRP_INT)


    ! Originally PPICLF_OPT.h
    ! Particle options
    LOGICAL PPICLF_RESTART, PPICLF_OVERLAP, PPICLF_LCOMM, PPICLF_LINIT, PPICLF_LINTP, PPICLF_LPROJ, PPICLF_LSUBSUBBIN,PPICLF_EQUALDOMAIN(3), PPICLF_LINPERIODIC(3), &
        PPICLF_REMOVE_PARTICLE, PPICLF_BINCHANGED, PPICLF_PRINTBINVTU, PPICLF_READYTOSOLVE

    DATA PPICLF_LCOMM /.false./
    DATA PPICLF_RESTART /.false./

    INTEGER*4 PPICLF_NDIM, PPICLF_IMETHOD, PPICLF_NGRIDS, PPICLF_CYCLE, PPICLF_IOSTEP, PPICLF_IENDIAN, PPICLF_IWALLM

    REAL*8 PPICLF_FILTER(3), PPICLF_RK3COEF(3,3), PPICLF_DT, PPICLF_TIME, PPICLF_NNDIST, PPICLF_INTERP_DCHK(3), PPICLF_TOTNNDIST(PPICLF_LPART)
    REAL*8 PPICLF_RK3ARK(3)


    ! Originally PPICLF_PARALLEL.h
    ! Communication

    INTEGER*4 PPICLF_COMM, PPICLF_NP, PPICLF_NID, PPICLF_CR_HNDL,PPICLF_FP_HNDL, PPICLF_COMM_NID
    DATA PPICLF_NID /0/

    ! Bins
    INTEGER*4 PPICLF_N_BINS(3) 

    REAL*8 PPICLF_BINS_DX(3), PPICLF_BINB(6), PPICLF_BIN_POS(2,3),PPICLF_PREVIOUSBINB(6)
      

    ! Ghost particles
    ! REAL*8 PPICLF_RPROP_GP(PPICLF_LRP_GP,PPICLF_LPART_GP),PPICLF_CP_MAP(PPICLF_LRP_GP,PPICLF_LPART)


    ! INTEGER*4 PPICLF_IPROP_GP(PPICLF_LIP_GP,PPICLF_LPART_GP),
    INTEGER*4 PPICLF_PARTICLEMOVED

    

    INTEGER*4  PPICLF_NPART_GP

    LOGICAL PPICLF_GPREQUIRED

    INTEGER*4 PARTICLE_NN(PPICLF_LPART)


    ! Originally PPICLF_GEOM.h

    ! Wall support
    REAL*8 PPICLF_WALL_C(9,PPICLF_LWALL),PPICLF_WALL_N(4,PPICLF_LWALL)

    INTEGER*4  PPICLF_NWALL


    ! Originally PPICLF_USER_COMMON.h
    ! this seems like really poor design, I dont think the user code should be defining variables used by 
    ! PPICLF core (configuring things like array sizes via defines are not the same)

    ! For ppiclf_solve_InitAngularPeriodic
    integer*4 x_per_flag, y_per_flag, z_per_flag, ang_per_flag,ang_case 
    real*8 ang_per_angle, ang_per_xangle,ang_per_rin, ang_per_rout,xrot(3) , vrot(3)
    real*8 x_per_min, x_per_max,y_per_min, y_per_max, z_per_min, z_per_max

end module ppiclf_data