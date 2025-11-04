!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Quasi-steady force fluctuations
! Osnes, Vartdal, Khalloufi, Capecelatro (2023)
!   Comprehensive quasi-steady force correlations
!   for compressible flow through random particle
!   suspensions.
!   International Journal of Multiphase Flows,
!   Vo. 165, 104485.
! Lattanzi, Tavanashad, Subramaniam, Capecelatro (2022)
!   Stochastic model for the hydrodynamic force in
!   Euler-Lagrange silumations of particle-laden flows.
!   Physical Review Fluids, Vol. 7, 014301.
! Note: To compute the granular temperature, we assume
!   the velocity fluctuations are uncorrelated.
! Note: The means are computed using a box filter with an
!   adaptive filter width.
! Compute mean using box filter for langevin model - not for fedback
!
! The mean is calcuated according to Lattanzi etal,
!   Physical Review Fluids, 2022.
!
! Sam - TODO: either couple the projection filter to the
! fluctuation filter or completely decouple them.
! Right now they use the same filter width and are assumed to
! be the same. 
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_fluct_Lattanzi(i,iStage,fqs_fluct)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 fqs_fluct(3)

      real*8 aSDE,bq,bSDE,chi,denum,dW1,dW2,dW3,fq,Fs,gkern,
     >   sigF,tF_inv,theta,upflct,vpflct,wpflct,Z1,Z2,Z3
      real*8 TwoPi

!
! Code:
!
      TwoPi = 2.0d0*acos(-1.0d0)


      if (qs_fluct_filter_flag==0) then
         denum = max(dfloat(icpmean),1.d0)  ! for arithmetic mean
      else if (qs_fluct_filter_flag==1) then
         denum = max(phipmean,1.0e-6)  ! for volume mean
      endif
      

      upmean = upmean / denum
      vpmean = vpmean / denum
      wpmean = wpmean / denum
      u2pmean = u2pmean / denum
      v2pmean = v2pmean / denum
      w2pmean = w2pmean / denum

      if (ppiclf_debug==2) then
         if (ppiclf_nid==0 .and. iStage==1) then
            write(6,"(2x,E16.8,i5,16(1x,F13.8))") ppiclf_time,
     >                denum,upmean,vpmean,wpmean
         endif
      endif

      ! Lattenzi is valid only for incompressible flows,
      !   so here we use the compressible correction of Osnes
      !   Equations (9-12) of Osnes paper, with eqn (9) corrected
      !   Note that sigD has units of force; N = Pa-m^2
      fq = 6.52*rphip - 22.56*(rphip**2) + 49.90*(rphip**3)
      Fs = 3.0*rpi*rmu*dp*(1.0+0.15*((rep*(1.0-rphip))**0.687))
     >          *(1.0-rphip)*vmag
      bq = min(sqrt(20.0*rmachp),1.0)*0.55*(rphip**0.7) 
     >          *(1.0+tanh((rmachp-0.5)/0.2))
      sigF = (fq + bq)*Fs

      chi = (1.0+2.50*rphip+4.51*(rphip**2)+4.52*(rphip**3))
     >         /((1.0-(rphip/0.64)**3)**0.68)

      fqs_fluct = 0.0d0

! Particle velocity fluctuation and Granular Temperature
! Need particle velocity mean
! Though the theory assumes granular temperature to be an 
!    average over neighboring particles, here it is approximated 
!    as that of the chosen particle - Comment 3/6/24
!    This is now fixed - Comment 4/12/24
!
      ! Particle velocity fluctuation
      upflct = ppiclf_y(PPICLF_JVX,i) - upmean
      vpflct = ppiclf_y(PPICLF_JVY,i) - vpmean
      wpflct = ppiclf_y(PPICLF_JVZ,i) - wpmean

      ! Granular temperature
      ! theta = (upflct*upflct + vpflct*vpflct + wpflct*wpflct)/3.0
      ! This is averaged over neighboring particles
      theta  = ((u2pmean + v2pmean + w2pmean) - 
     >          (upmean**2 + vpmean**2 + wpmean**2))/3.0d0

      ! 11/21/24 - Thierry - prevent NaN variables
      if(theta.le.1.d-12) then
        theta = 0.0d0
      endif

      tF_inv = (24.0*rphip*chi)/dp * sqrt(theta/rpi)

      aSDE = tF_inv
      bSDE = sigF*sqrt(2.0*tF_inv)

      call RANDOM_NUMBER(UnifRnd)

      Z1 = sqrt(-2.0d0*log(UnifRnd(1)))*cos(TwoPi*UnifRnd(2))
      Z2 = sqrt(-2.0d0*log(UnifRnd(3)))*cos(TwoPi*UnifRnd(4))
      Z3 = sqrt(-2.0d0*log(UnifRnd(5)))*cos(TwoPi*UnifRnd(6))

      dW1 = sqrt(fac)*Z1
      dW2 = sqrt(fac)*Z2
      dW3 = sqrt(fac)*Z3

      fqs_fluct(1) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFX,i)
     >           + bSDE*dW1
      fqs_fluct(2) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFY,i)
     >           + bSDE*dW2
      fqs_fluct(3) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFZ,i)
     >           + bSDE*dW3


      if (ppiclf_debug==2 .and. (iStage==1 .and. ppiclf_nid==0)) then
         if (ppiclf_time.gt.2.d-8) then
         if (i<=10) then
            write(7350+(i-1)*1,*) i,ppiclf_time,             ! 0-1
     >         rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,  ! 2-10
     >         rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac, ! 11-18
     >         vx,vy,vz,ppiclf_dt,                           ! 19-22
     >         ppiclf_npart,ppiclf_n_bins(1:3),              ! 23-26
     >         ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3), ! 27
     >         ppiclf_binb(1:6),                             ! 28-33
     >         upmean,vpmean,wpmean,phipmean,                ! 34-37
     >         ppiclf_y(PPICLF_JVX:PPICLF_JVZ,i),            ! 38-40
     >         upflct,vpflct,wpflct,icpmean,                 ! 41-44
     >         fq,Fs,bq,theta,chi,tF_inv,                    ! 45-50
     >         aSDE,bSDE,sigF,                               ! 51-53
     >         fqs_fluct(1:3),                               ! 54-56
     >         Z1,Z2,Z3,dW1,dW2,dW3,                         ! 57-62
     >         ppiclf_np
         endif
         endif
      endif


      return
      end
!
!
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024 - T.L. Jackson
! Modified 3/6/24 - Balachandar
!
! Quasi-steady force fluctuations
! Osnes, Vartdal, Khalloufi, Capecelatro (2023)
!   Comprehensive quasi-steady force correlations
!   for compressible flow through random particle
!   suspensions.
!   International Journal of Multiphase Flows,
!   Vo. 165, 104485.
! Lattanzi, Tavanashad, Subramaniam, Capecelatro (2022)
!   Stochastic model for the hydrodynamic force in
!   Euler-Lagrange silumations of particle-laden flows.
!   Physical Review Fluids, Vol. 7, 014301.
! Note: To compute the granular temperature, we assume
!   the velocity fluctuations are uncorrelated.
! Note: The means are computed using a box filter with an
!   adaptive filter width.
! Compute mean using box filter for langevin model - not for fedback
!
! The mean is calcuated according to Lattanzi etal,
!   Physical Review Fluids, 2022.
!
! Sam - TODO: either couple the projection filter to the
! fluctuation filter or completely decouple them.
! Right now they use the same filter width and are assumed to
! be the same. 
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_fluct_Osnes(i,iStage,fqs_fluct,
     >                                      xi_par,xi_perp,xi_T,
     >                                      fqsx,fqsy,fqsz)
!                                                    
      implicit none
!
      include "PPICLF"
!
! Input:
      integer*4 i, iStage
      real*8 fqsx, fqsy, fqsz
!
! Output:
      real*8 xi_par, xi_perp, xi_T
      real*8 fqs_fluct(3)
!
! Internal:
!
      real*8 aSDE,bq,chi,denum,dW1,dW2,dW3,fq,Fs,gkern,
     >   sigD,tF_inv,theta,upflct,vpflct,wpflct,Z1,Z2,Z3
      real*8 TwoPi
      real*8 bSDE_CD, bSDE_CL, bSDE_CT, CD_frac, CD_prime
      real*8 sigT,sigCT
      real*8 sigmoid_cf, f_CF
      real*8 avec(3)
      real*8 bvec(3)
      real*8 cvec(3)
      real*8 dvec(3)
      real*8 cosrand,sinrand
      real*8 eunit(3)

      integer*4 m
      real*8 s_par, s_perp, s_T, Rmean_par, Rmean_perp, R_par, R_perp
      real*8 R(3,3), Q(3,3), Qt(3,3), Tmean_par(3)
      real*8 CD_average
      real*8 k_tilde, k_Mach, b_tilde, b_Mach, b_par, b_perp,
     >       k_Osnes, b_Osnes
      real*8 C1, C2, C3, C4, C5,
     >       D1, D2, D3, D4, D5, D6, D7, D8,
     >       E1, E2, E3, E4
      real*8 F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11,
     >       G1, G2, G3, G4, G5, G6, G7, G8,
     >       H1, H2, H3, H4, H5, H6, H7, H8,
     >       A1, A2, A3, A4
      real*8 D9, D10, D11
      real*8 fit_func
!
! Code:
!
      TwoPi = 2.0d0*acos(-1.0d0)

      if (qs_fluct_filter_flag==0) then
         denum = max(dfloat(icpmean),1.d0)  ! for arithmetic mean
      else if (qs_fluct_filter_flag==1) then
         denum = max(phipmean,1.0e-6)  ! for volume mean
      endif
      
      upmean = upmean / denum
      vpmean = vpmean / denum
      wpmean = wpmean / denum
      u2pmean = u2pmean / denum
      v2pmean = v2pmean / denum
      w2pmean = w2pmean / denum

      if (ppiclf_debug==2) then
         if (ppiclf_nid==0 .and. iStage==1) then
            write(6,*) 'FLUC1 ',
     >          ppiclf_time,denum,upmean,vpmean,wpmean,
     >          abs(upmean),abs(vpmean),abs(wpmean),
     >          u2pmean,v2pmean,w2pmean
         endif
      endif

      ! Equations (9-12) of Osnes paper, with eqn (9) corrected
      ! Note that sigD has units of force; N = Pa-m^2
      fq = 6.52*rphip - 22.56*(rphip**2) + 49.90*(rphip**3)
      Fs = 3.0*rpi*rmu*dp*(1.0+0.15*((rep*(1.0-rphip))**0.687))
     >          *(1.0-rphip)*vmag
      bq = min(sqrt(20.0*rmachp),1.0)*0.55*(rphip**0.7) 
     >          *(1.0+tanh((rmachp-0.5)/0.2))
      sigD = (fq + bq)*Fs

! Particle velocity fluctuation and Granular Temperature
! Need particle velocity mean
! Though the theory assumes granular temperature to be an 
!    average over neighboring particles, here it is approximated 
!    as that of the chosen particle - Comment 3/6/24
!    This is now fixed - Comment 4/12/24
!
      upflct = ppiclf_y(PPICLF_JVX,i) - upmean
      vpflct = ppiclf_y(PPICLF_JVY,i) - vpmean
      wpflct = ppiclf_y(PPICLF_JVZ,i) - wpmean

      ! Granular temperature
      ! theta = (upflct*upflct + vpflct*vpflct + wpflct*wpflct)/3.0
      ! This is averaged over neighboring particles
      theta  = ((u2pmean + v2pmean + w2pmean) - 
     >          (upmean**2 + vpmean**2 + wpmean**2))/3.0d0

      ! 11/21/24 - Thierry - prevent NaN variables
      if(theta.le.1.d-12) then
        theta = 0.0d0
      endif

      chi = (1.0 + 2.50*rphip + 4.51*(rphip**2) + 4.52*(rphip**3))
     >         /((1.0-(rphip/0.64)**3)**0.68)

      tF_inv = (24.0*rphip*chi/dp) * sqrt(theta/rpi)

      aSDE = tF_inv
      bSDE_CD = sigD*sqrt(2.0*tF_inv)  ! Modified 3/6/24

! Fluctuating perpendicular force
! Compare CD' (units of force) against sigma_CD (units of force) 
!    to determine which of 5 bins to use for sigCT
!
! Added 3/6/24 
! Modified 3/14/24 
!
      ! 03/13/2025 - Thierry - if velocity is very small, don't impose fluctuations
      if(vmag > 1.d-8) then
        avec = [vx,vy,vz]/vmag

        CD_prime = ppiclf_rprop(PPICLF_R_FLUCTFX,i)*avec(1) +
     >             ppiclf_rprop(PPICLF_R_FLUCTFY,i)*avec(2) +
     >             ppiclf_rprop(PPICLF_R_FLUCTFZ,i)*avec(3)
        CD_frac  = CD_prime/sigD

      else
        avec     = [1.0, 0.0, 0.0]
        CD_prime = 0.0
        sigD     = 0.0
        CD_frac  = 0.0
      endif

      ! Thierry Daoud - Updated June 2, 2024
      sigmoid_cf = 1.0 / (1.0 + exp(-CD_frac))
      f_CF = 0.39356905*sigmoid_cf + 0.43758848
      sigT  = f_CF*sigD
      bSDE_CL = sigT*sqrt(2.0*tF_inv)

      if (ppiclf_debug==2) then
      if (i<=4) then
         if (ppiclf_nid==0 .and. iStage==1) then
            write(6,*) 'FLUC1 ',i,
     >        CD_prime,CD_frac,sigD,theta,bSDE_CD,bSDE_CL
         endif
      endif
      endif


! Calculate the three orthogonal unit vectors
! The first vector (avec) is vx/vmag, vy/vmag, and vz/vmag
! The second (bvec) is constructued by taking cross-product with eunit
!   Note: if avec is in dir. of e_x=(1,0,0), use e_y=(0,1,0) to get e_z
!       : if avec is in dir. of e_y=(0,1,0), use e_z=(0,0,1) to get e_x
!       : if avec is in dir. of e_z=(0,0,1), use e_x=(1,0,0) to get e_y
! The third (cvec) is cross product of the first two
! written 3/6/24
!
! avec : unit vector in main direction
! bvec, cvec: two orthogonal vectors to avec
      eunit = [1,0,0]
      if (abs(avec(2))+abs(avec(3)) <= 1.d-8) then
         eunit = [0,1,0]
      elseif (abs(avec(1))+abs(avec(3)) <= 1.d-8) then
         eunit = [0,0,1]
      endif

      bvec(1) = avec(2)*eunit(3) - avec(3)*eunit(2)
      bvec(2) = avec(3)*eunit(1) - avec(1)*eunit(3)
      bvec(3) = avec(1)*eunit(2) - avec(2)*eunit(1)
      denum   = max(1.d-8,sqrt(bvec(1)**2 + bvec(2)**2 + bvec(3)**2))
      bvec    = bvec / denum

      cvec(1) = avec(2)*bvec(3) - avec(3)*bvec(2)
      cvec(2) = avec(3)*bvec(1) - avec(1)*bvec(3)
      cvec(3) = avec(1)*bvec(2) - avec(2)*bvec(1)
      denum   = max(1.d-8,sqrt(cvec(1)**2 + cvec(2)**2 + cvec(3)**2))
      cvec    = cvec / denum

      ! Generate  Gaussian Random Values
      call RANDOM_NUMBER(UnifRnd)

! Box-Muller transform for generating two independent standard normal
! (Gaussian) random variables
! Z1 & Z2 are standard normal random variables
      Z1 = sqrt(-2.0d0*log(UnifRnd(1)))*cos(TwoPi*UnifRnd(2))
      Z2 = sqrt(-2.0d0*log(UnifRnd(3)))*cos(TwoPi*UnifRnd(4))

! dW1 & dW2 are scaled stochastic amplitudes       
      dW1 = sqrt(fac)*Z1
      dW2 = sqrt(fac)*Z2

! Random mixture of bvec and cvec - make sure the new one is a unit vector
! Added 3/6/24
! dvec : Random unit direction in perpendicular plane 
      cosrand = cos(TwoPi*UnifRnd(5))
      sinrand = sin(TwoPi*UnifRnd(5)) 
      dvec(1) = bvec(1)*cosrand + cvec(1)*sinrand
      dvec(2) = bvec(2)*cosrand + cvec(2)*sinrand
      dvec(3) = bvec(3)*cosrand + cvec(3)*sinrand
      denum   = max(1.d-8,sqrt(dvec(1)**2 + dvec(2)**2 + dvec(3)**2))
      dvec    = dvec/denum

      ! Store and project back fqs_fluct values only if QSFLUCT = 2
      ! Otherwise, need to call subroutine for PseudoTurbulence
      if(qs_fluct_flag .eq. 0) then
        fqs_fluct = 0.0d0
      else
      fqs_fluct(1) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFX,i)
     >           + bSDE_CD*dW1*avec(1) + bSDE_CL*dW2*dvec(1)
      fqs_fluct(2) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFY,i)
     >           + bSDE_CD*dW1*avec(2) + bSDE_CL*dW2*dvec(2)
      fqs_fluct(3) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFZ,i)
     >           + bSDE_CD*dW1*avec(3) + bSDE_CL*dW2*dvec(3)
      endif

      if (ppiclf_debug==2 .and. (iStage==1 .and. ppiclf_nid==0)) then
         if (ppiclf_time.gt.2.d-8) then
         if (i<=10) then
            write(7350+(i-1)*1,*) i,ppiclf_time,             ! 0-1
     >         rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,  ! 2-10
     >         rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac, ! 11-18
     >         vx,vy,vz,ppiclf_dt,                           ! 19-22
     >         ppiclf_npart,ppiclf_n_bins(1:3),              ! 23-26
     >         ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3), ! 27
     >         ppiclf_binb(1:6),                             ! 28-33
     >         upmean,vpmean,wpmean,phipmean,                ! 34-37
     >         ppiclf_y(PPICLF_JVX:PPICLF_JVZ,i),            ! 38-40
     >         upflct,vpflct,wpflct,icpmean,                 ! 41-44
     >         fq,Fs,bq,theta,chi,tF_inv,                    ! 45-50
     >         aSDE,bSDE_CD,bSDE_CL,                         ! 51-53
     >         sigD,sigT,                                    ! 54-55
     >         CD_prime,CD_frac,sigmoid_cf,f_CF,             ! 56-59
     >         eunit,avec,bvec,                              ! 60-68
     >         cvec,dvec,rpi,                                ! 69-75
     >         fqs_fluct(1:3),                               ! 76-78
     >         Z1,Z2,dW1,dW2,                                ! 79-82
     >         ppiclf_np,                                    ! 83
     >         ppiclf_y(PPICLF_JX:PPICLF_JZ,i)               ! 84-86
         endif
         endif
      endif

!---------------------------------------------------------------------------
! Pseudo-Turbulence Calculations starts here 
      if(pseudoTurb_flag==1) then
        !phi = max(0.05d0, min(0.3d0, rphip))
        !mp = max(0.0d0, min(0.87d0, rmachp))
        !re = max(0.0d0, min(266.0d0, rep))
        phi = rphip
        mp = rmachp
        re = rep
        rem = (1.0-phi)*re

      ! Constants taken from Osnes PseudoTurbulent paper, Table 1
        C1 = -1.2152; D1 = -0.0462; E1 = -0.2906;
        C2 = -7.6314; D2 = -0.1068; E2 =  1.1899; 
        C3 =  0.2889; D3 =  0.6793; E3 =  0.5218; 
        C4 =  0.6143; D4 =  1.1461; E4 =  0.0699;
        C5 =  0.3082; D5 = -2.6886; 
                      D6 = -2.1376;
                      D7 =  0.4873;
                      D8 =  0.2395;
                    ! Dr. Bala's terms
                      D9  =  0.716
                      D10 = -2.14
                      D11 =  1.6

      ! Constants taken from Osnes PseudoTurbulent paper, Table 2
        F1  = -0.0022; G1 = -0.2867; H1 =  0.4992
        F2  = -0.0219; G2 =  0.2176; H2 = -1.3528
        F3  =  0.0932; G3 =  0.2826; H3 = -0.1358
        F4  = -0.0135; G4 = -0.0644; H4 = -0.1463
        F5  =  0.0361; G5 =  0.0466; H5 =  0.2583
        F6  =  0.0403; G6 =  0.0973; H6 = -0.3339
        F7  = -0.0761; G7 = -0.0081; H7 = -0.0407
        F8  =  0.0599; G8 = -0.0235; H8 = -0.0806
        F9  =  0.0164;
        F10 =  0.0453;
        F11 = -0.0265;
     
        ! zero out variables  at first
        Rmean_par = 0.0d0 ; Rmean_perp = 0.0d0
        R = 0.0d0; Rsg = 0.0d0
        Tmean_par = 0.0d0; T_par = 0.0d0
        
        ! CD_average is zero at early time steps

        avec = [vx,vy,vz]/vmag

        CD_average = fqsx*avec(1) +
     >               fqsy*avec(2) +
     >               fqsz*avec(3)

        ! avoiding singularity
        if(CD_average .lt. 1.d-8) return

        ! Reynolds Subgrid Stress Tensor - Eulerian Mean Model 
                                                                       
        ! Reynolds number and vol fraction dependent k^tilde and b_par 
        ! Mehrabadi's terms
        k_tilde = 2.0*phi + 2.5*phi*((1.0-phi)**3) * 
     >         exp(-phi*(rem**0.5))                                 
                                                                    
!        b_par = 0.523/(1.0+0.305*exp(-0.114*rem)) *
!     >          exp(-3.511*phi/(1.0+1.801*exp(-0.005*rem)))

        ! 08/25/2025 - Thierry - Fitted phip function to better match
        ! formulation with Osnes's low Mach number data
        fit_func = -10.18530152*phi**3 + 10.94163073*phi**2
     >              -7.07374862*phi +  0.38424203
        b_par = 0.523/(1.0+0.305*exp(-0.114*rem)) *
     >          exp(fit_func/(1.0+1.801*exp(-0.005*rem)))
                                                                       
        ! Mach number correction provided by Osnes                            
!        k_Mach = phi*(C1 + C2*phi + re**C3) * 
!     >        (tanh(C4/C5) + tanh((mp - C4)/C5))

        ! Mach number correction at Re=100, coeff taken from Osnes
        ! cap vol fraction here at 0.3
        k_Mach = min(phi,0.3)*(-6.918*min(phi,0.3) + 2.238) *
     >        (tanh(C4/C5) + tanh((mp - C4)/C5))

!        b_Mach = (D1 + (re/300.0)*(D2 + D3*re/300.0) +
!     >         phi*(D4 + D5*(re**2/300.0**2) + D6*phi)) *
!     >         (tanh(-D7/D8) - tanh((mp-D7)/D8))

        ! First term was corrected by Dr. Bala (D9, D10, D11)
        b_Mach = phi*(D9 + D10*phi + D11*phi**2) *
     >         (tanh(-D7/D8) - tanh((mp-D7)/D8))

        ! Corrected k^tilde and b_par components                       
        k_Osnes =  k_tilde*(1.0d0 + k_Mach)
        b_Osnes =  b_par  *(1.0d0 + b_Mach)
        b_perp  = -b_Osnes/2.0d0
                                                                       
        ! Mean Eulerian Reynolds Subgrid Stress - Parallel Component   
        Rmean_par  = 2.0d0*k_Osnes*(b_Osnes  + 1.0d0/3.0d0)
                                                                       
        ! Mean Eulerian Reynolds Subgrid Stress - Perpendicular Component
        Rmean_perp = 2.0d0*k_Osnes*(b_perp + 1.0d0/3.0d0)

        ! Mean Pseudo Turbulent Kinetic Energy Model
        Tmean_par = E1 + E2*phi/(E3 + re/300.0) + E4*mp
        
c--  Multiply by the mean relative flow kinetic energy to dimentionalize      
        Rmean_par  = Rmean_par  * 0.5d0 * vmag**2
        Rmean_perp = Rmean_perp * 0.5d0 * vmag**2

c--  Multiply by the mean relative velocity & flow kinetic energy to dimentionalize      
        Tmean_par(1)  = Tmean_par(1) * vx * k_Osnes * 0.5d0 * vmag**2
        Tmean_par(2)  = Tmean_par(2) * vy * k_Osnes * 0.5d0 * vmag**2
        Tmean_par(3)  = Tmean_par(3) * vz * k_Osnes * 0.5d0 * vmag**2

c------ Lagrangian Model
  
        ! 08/15/2025 - Ditch A1 per Dr. Bala, and set A2 to this constant per Osnes
        A1 = 0.0d0
        A2 = 0.064

        A3 = G1 + G2/(min(phi,0.3) + G3) + G4 * mp
        
        ! 09/02/2025 - Cap according to Osnes model range
        A4 = H1 + H2*max(0.0d0, min(0.3d0, phi))
     >          + H3*max(0.0d0, min(0.87d0, mp)) 
     >          + H4*max(30.0d0, min(266.0d0, re))/300.0d0
  
        s_par = F8 + F9/(phi + F10) + F11 * mp
        !s_perp = G5/(phi + G6) + (G7*re)/(300.0*phi) + G8
c---     We ditch Osnes's expression for s_perp and assume it as big as s_par
        s_perp = s_par
        
        ! 09/02/2025 - Cap according to Osnes model range
        s_T = H5 + H6*max(0.0d0, min(0.3d0, phi))
     >           + H7*max(0.0d0, min(0.87d0, mp)) 
     >           + H8*max(30.0d0, min(266.0d0, re))/300.0d0

        tF_inv = (24.0*phi*chi/dp) * sqrt(theta/rpi)
        aSDE = tF_inv
        bSDE_CD = s_par *sqrt(2.0*tF_inv)
        bSDE_CL = s_perp*sqrt(2.0*tF_inv)
        bSDE_CT = s_T *sqrt(2.0*tF_inv)
  
        call RANDOM_NUMBER(UnifRnd)
        
        ! Box-Muller transform for generating two independent standard normal
        ! (Gaussian) random variables
        ! Z1 & Z2 are standard normal random variables
        Z1 = sqrt(-2.0d0*log(UnifRnd(1))) * cos(TwoPi*UnifRnd(2))
        Z2 = sqrt(-2.0d0*log(UnifRnd(3))) * sin(TwoPi*UnifRnd(4))
        Z3 = sqrt(-2.0d0*log(UnifRnd(5))) * cos(TwoPi*UnifRnd(6))

        ! dW1 & dW2 are scaled stochastic amplitudes       
        dW1 = sqrt(fac)*Z1
        dW2 = sqrt(fac)*Z2
        dW3 = sqrt(fac)*Z3
  
        ! Langevin Model implemented for xi_par, xi_perp, xi_T
        xi_par = (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_XIPAR,i)
     >            + bSDE_CD*dW1
        xi_perp = (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_XIPERP,i)
     >            + bSDE_CL*dW2
        xi_T = (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_XIT,i)
     >            + bSDE_CT*dW3

        ! CD_prime has unit of Force
        ! CD_average has unit of Force

        ! Lagrangian Reynolds Subgrid Stress - Parallel Component
        R_par = 1.0 + A1 + A2 * CD_prime / CD_average + xi_par
  
        ! Lagrangian Reynolds Subgrid Stress - Perpendicular Component
        R_perp = 1.0 + A3 * CD_prime / CD_average + xi_perp

c--  Multiply Lagrangian Model by the Eulerian Mean Model
        R_par  = R_par  * Rmean_par 
        R_perp = R_perp * Rmean_perp

c---  Q = [avec | bvec | cvec], 3x3 matrix
c---  avec : unit vector in main direction
c---  bvec, cvec: two orthogonal vectors to avec

        ! 08/15/2025 - Thierry - still need to finalize how to do the
        ! rotation of the R tensor 
        do m=1,3
          Q(m,1) = avec(m)
          Q(m,2) = bvec(m)
          Q(m,3) = cvec(m)
        enddo
  
        Qt = transpose(Q)
  
c--- R = |R_par,   0   ,   0   |
c---     | 0   , R_perp,   0   |
c---     | 0       0   , R_perp|
  
c---  R matrix only has diagonal components
        R(1,1) = R_par
        R(2,2) = R_perp
        R(3,3) = R_perp
  
c--- Now Rotate the matrix, Rsg = Q . R . Q^T
  
       Rsg = matmul(Q, matmul(R,Qt))

c--- Osnes Formulation for PTKE

      T_par = A4 * CD_prime/CD_average + xi_T

c--  Multiply by the mean relative velocity & flow kinetic energy to dimentionalize      
c--  then add mean PTKE
       T_par(1) = T_par(1) * vx * k_Osnes * 0.5d0 * vmag**2 
     >            + Tmean_par(1)

       T_par(2) = T_par(2) * vy * k_Osnes * 0.5d0 * vmag**2 
     >            + Tmean_par(2)

       T_par(3) = T_par(3) * vz * k_Osnes * 0.5d0 * vmag**2 
     >            + Tmean_par(3)

      endif ! pseudoTurb_flag

      return
      end
