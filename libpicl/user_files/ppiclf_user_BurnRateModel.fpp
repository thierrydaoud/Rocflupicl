!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for burn rate models for reactive particles
!
! if heattransfer_flag = 0  ignore heat transfer
!                      = 1  Stokes
!                      = 2  Ranz-Marshall (1952)
!                      = 3  Gunn (1977)
!                      = 4  Fox (1978)
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_BR_driver(i,iStage,burnrate_model,
     >   qq,mdot_me,mdot_ox)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage, burnrate_model
      real*8 qq
      real*8 mdot_me, mdot_ox

!
! Code:
!
      if (burnrate_model == 1) then
         call AL_CombModel(i,iStage,qq,mdot_me,mdot_ox)
      elseif (burnrate_model == 2) then
         !call Carbon_CombModel(i,iStage,qq,mdot_me,mdot_ox)
      elseif (burnrate_model == 3) then
         !call Mg_CombModel(i,iStage,qq,mdot_me,mdot_ox)
      else
         call ppiclf_exittr('Unknown combustion model$', 0.0d0, 0)
      endif


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
!-----------------------------------------------------------------------
!
      subroutine AL_CombModel(i,iStage,qq,mdot_me,mdot_ox)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i,iStage
      real*8 qq
      real*8 mdot_me, mdot_ox

      REAL*8 emi, A, E_a, hr, Ks, r_gas
      REAL*8 Q_HSR, Q_comb, T_part, T_ign, T_evap
      REAL*8 Dp_me, m_me, h_comb, h_evap, rho_me

      !
      real*8 rho_ox, m_ox
      real*8 MW_me, MW_o, MW_me2o3
      real*8 oxide_t, V_me, V_ox
      real*8 vol_avg, rho_p, Dia, D0
      real*8 Cp_me, Cp_ox, t
      real*8 rmass_therm

      REAL*8 T_flame
      REAL*8 C, Cs
      REAL*8 xi_co2, xi_h2, xi_h2o, xi_o2, xi_eff
      REAL*8 psi_me, V_p
      REAL*8 Pres
      INTEGER*4 IVALUE, OxideFilm


!===============================================================
!     PARTICLE PROPERTIES
!===============================================================

      T_part = PPICLF_y(PPICLF_JT,i)
      Pres   = PPICLF_RPROP(PPICLF_R_JP,i)

      m_me = PPICLF_y(PPICLF_JMETAL,i)
      m_ox = PPICLF_y(PPICLF_JOXIDE,i)

      Dia    = PPICLF_RPROP(PPICLF_R_JDP,i)
      V_p    = PPICLF_RPROP(PPICLF_R_JVOLP,i)

      D0 = PPICLF_RPROP(PPICLF_R_JIDP,i)

      V_me = m_me / rho_me
      psi_me = V_me / V_p
      Dp_me = ( 6.0d0 * m_me / rpi / rho_me )**(1.0d0/3.0d0)


!-----Al2O3 ignition as a function of diameter, Sundaram et al.--

      if (D0 .le. 100.0d-6) then
         T_ign = 368.0d0 * (D0*1.0d6)**0.268d0 + 780.0d0
      elseif (D0 .le. 1000.0d-6 .and. D0 .gt. 100.0d-6) then
         T_ign = 0.1617d0 * (D0*1.0d6) + 2040.0d0
      else
         T_ign = 3000.0d0
      endif

!-----Updating density--------------------------------------------

      !metal oxide
      if (T_part .le. 3250.0d0 .and. T_part .gt. 2345.0d0) then
         rho_ox = 5632.0d0 - 1.127d0*T_part
      elseif (T_part .le. 2345.0d0) then
         rho_ox = 3970.0d0 * (1.0d0 - 8.0d-6 * (T_part - 300.0d0))
      else
         rho_ox = 5632.0d0 - 1.127d0*3250.0d0
      endif

      !metal
      if (T_part .le. 2743.0d0 .and. T_part .gt. 933.0d0) then
         rho_me = 2385.0d0 - 0.280d0 * (T_part - 933.0d0)
      elseif (T_part .le. 933.0d0) then
         rho_me = 2700.0d0 * (1.0d0 - 23.1d-6 * (T_part - 300.0d0))
      else
         rho_me = 2385.0d0 - 0.280d0 * (2743.0d0 - 933.0d0)
      endif

      !----Updating particle properties--------------------------------

      V_me = m_me/rho_me
      V_ox = m_ox/rho_ox
      vol_avg  = V_me + V_ox

!      phi_me = m_me / (m_me + m_ox)
!      phi_ox = 1.0d0 - phi_me
!      psi_me = ( m_me/rho_me ) / vol_avg
!      psi_ox = 1.0d0 - psi_me

      rho_p = ( m_me + m_ox ) / vol_avg
      Dia = ( 6.0d0 * vol_avg / rpi )**( 1.0d0 / 3.0d0 )

      !-----Molar Weight------------------------------------------------

      MW_me = 0.0269815384d0       ! kg/mol
      MW_o = 0.015999d0            ! kg/mol
      MW_me2o3 = MW_me * 2.0d0 + MW_o * 3.0d0 ! kg/mol

!----------------------------------------------------------------
!     Heat capacity for metal and oxide
!----------------------------------------------------------------

      t = T_part / 1000.0d0

      !Metal heat capacity
      if (T_part .le. 933.0d0) then
         !solid phase heat capacity (Shomate Equation)
         Cp_me = 28.08920d0 - 5.414849d0 * t + 8.560423d0 * t**2.0
     >      + 3.427370d0 * t**3.0d0 - 0.277375d0 / t**2.0

         Cp_me = Cp_me / MW_me ! J/mol*K to J/kg*K
      else
         !liquid phase heat capacity (Shomate Equation)
         Cp_me = 31.75104d0 + 3.935826d-8 * t - 1.786515d-8 * t**2.0
     >       + 2.694171d-9 * t**3.0d0 + 5.480037d-9 / t**2.0

         Cp_me = Cp_me / MW_me ! J/mol*K to J/kg*K
      endif


      !Metal oxide heat capacity
      if (T_part .le. 2327.0d0) then
         !solid phase heat capacity (Shomate Equation)
         Cp_ox = 102.4290d0 + 38.7498d0 * t - 15.91090d0 * t**2.0
     >      + 2.628181d0 * t**3.0d0 - 3.007551d0 / t**2.0

         Cp_ox = Cp_ox / MW_me2o3 ! J/mol*K to J/kg*K
      else
         !liquid phase heat capacity (Shomate Equation)
         Cp_ox = 192.4640d0 + 9.519856d-8 * t - 2.858928d-8 * t**2.0
     >      + 2.929147d-9 * t**3.0d0 + 5.599405d-8 / t**2.0

         Cp_ox = Cp_ox / MW_me2o3 ! J/mol*K to J/kg*K
      endif

      rmass_therm = m_me * Cp_me + m_ox * Cp_ox

!----Updating PPICLF_RPOP values---------------------------------

      ! TEMP FIX
      if (Dia.gt.PPICLF_RPROP(PPICLF_R_JIDP,i)) then
!         print*,'Warning Dia too big',
!     >    i, PPICLF_RPROP(PPICLF_R_JIDP,i),Dia
         Dia = PPICLF_RPROP(PPICLF_R_JIDP,i)
      endif
      PPICLF_RPROP(PPICLF_R_JDP,i)   = Dia
      PPICLF_RPROP(PPICLF_R_JRHOP,i) = rho_p
      PPICLF_RPROP(PPICLF_R_JVOLP,i) = vol_avg

!===============================================================
!     COMBUSTION HEAT TRANSFER
!===============================================================

      !constants - Al
      emi=0.1d0
      A=200.0d0
      E_a=95395.0d0
      hr=3.1d7
      h_evap=1.183d7
      T_evap=2743.d0

      r_gas=8.314d0 ! Gas Constant

!----------------------------------------------------------------
!     Source terms for heat transfer
!----------------------------------------------------------------

      !initialize source terms to 0

      Q_HSR = 0.0d0
      Q_comb = 0.0d0

      !preheat - HSR
      if (T_part .ge. 933.0d0 .and. T_part .lt. T_ign) then
         Ks = A * exp(-E_a / (r_gas * T_part))
         Q_HSR = hr * Ks * Dp_me**2.0d0 * rpi
      endif

      !combustion - comb
      if (T_part .ge. T_ign) then
            if (T_part .le. 3300.d0) then
               h_comb = -4.3334d7
            else
               h_comb = 3.3d6
            endif
         Q_comb = -mdot_me * (h_comb + h_evap)
      endif

      qq = qq + Q_HSR + Q_comb


!===============================================================
!     COMBUSTION MODEL
!===============================================================

!----------------------------------------------------------------
!     Setup information
!----------------------------------------------------------------

      !particle
      IVALUE = 0                   ! 0 for B&W; 1 for Maggi

      !constants
      C = 5.9427624643167829d-7    ! (m**1.8,s**-1,K**-0.2,Pa**-0.1)
      Cs = 0.0d0                   ! kg/m^3

!----------------------------------------------------------------
!         WIDENER and Beckstead, AIAA-98-382
!         With Fadi Najjar's correction, J. Spacecraft & Rockets,
!           43(6), 2006, pp.1258--1270
!----------------------------------------------------------------

!-------------MOLAR FRACTIONS IN THE GAS PHASE-------------------

      if (IVALUE==0) then  ! Widener and Beckstead (1998)
          xi_o2  = 0.02d0
          xi_h2o = 0.2d0
          xi_co2 = 0.2d0
          xi_h2  = 0.2d0
      endif
      if (IVALUE==1) then  ! Filippo Maggi (CSAR, 2007)
          xi_o2  = 0.00002d0
          xi_h2o = 0.04809d0
          xi_co2 = 0.00394d0
          xi_h2  = 0.33582d0
      endif

!-----Relative oxidizer concentration----------------------------

      xi_eff = 0.0291d0
!      xi_eff = xi_o2 + 0.58d0*xi_h2o + 0.22d0*xi_co2



!----------------------------------------------------------------
!             Aluminum mass burning rate (UNITS: SI)
!             Weighted by aluminum volume fraction
!             this formula takes diameter in meters
!----------------------------------------------------------------

!-------------Check if particle is burning-----------------------

              if ((T_part .ge. T_ign .or. PPICLF_RPROP(PPICLF_R_JBRNT,i)
     >           .gt. 0.0) .and. Dp_me .gt. 5.0d-6) then

                 !Burn law
                 mdot_me = C*rho_me*(Dia**1.2d0) *
     >              xi_eff*(T_part**0.2d0) *
     >              (pres**0.1d0) * psi_me

                 mdot_ox = 0.25d0 * rpi * Dia**2.0d0 *
     >              abs(vmag) * 0.25d0 * Cs

                 !update total burn time
                 PPICLF_RPROP(PPICLF_R_JBRNT,i) =
     >              PPICLF_RPROP(PPICLF_R_JBRNT,i) + ppiclf_dt
              else
                 mdot_me = 0.0d0
                 mdot_ox = 0.0d0
              endif

      return
      end
