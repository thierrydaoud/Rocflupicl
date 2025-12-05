!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for viscous unsteady force with history kernel
!
! Mei-Adrian history kernel
!
! Copied from either files in rocintereact/
!   INRT_CalcDragUnsteady_AMImplicit.F90
!   INRT_CalcDragUnsteady_AMExplicit.F90
!
! The number of time steps kept for the history
!   kernel is set in libpicl/ppiclF/source/PPICLF_STD.h
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_VU_Rocflu(i,iStage,fvux,fvuy,fvuz,qq_du)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage, iT
      real*8 fvux,fvuy,fvuz
      real*8 time,fH,factor,A,B,kernelVU
      real*8 factor_du,b_du,Pe,kernel_du,qq_du,alpha_fluid

!
! Code:
!
      fvux = 0.0d0
      fvuy = 0.0d0
      fvuz = 0.0d0
      iT   = 1
      time = 0.0d0

      qq_du = 0.0d0

      Pe = rep * rpr ! Peclet number

      fH     = 0.75d0 + .1055d0*rep
      ! Sangani's volume fraction correction for dilute random arrays
      ! Capping volume fraction at 0.5 
      factor = 3.0d0*rpi*rnu*dp*ppiclf_dt*(1.0+2.28*min(rphip,0.5))

      ! Thermal diffusivity
      alpha_fluid = rkappa/(rhof * rcp_fluid)
      ! diffusive unsteady Heat Transfer
      factor_du = sqrt(rpi*(dp**3)*rkappa*rhof*rcp_fluid*vmag)*ppiclf_dt

      if (ppiclf_nTimeBH > 1) then
         do iT = 2,ppiclf_nTimeBH-1
            time = ppiclf_timeBH(iT)

            A  = (4.0d0*rpi*time*rnu/(dp**2))**(.25d0)
            B  = (rpi*(vmag**3)*(time**2)/ 
     >                 (dp*rnu*(fH**3)))**(.5d0)

            kernelVU = factor*(A+B)**(-2)

            fvux = fvux + kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JX,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JX,iT,i) )
            fvuy = fvuy + kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JY,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JY,iT,i) )
            fvuz = fvuz + kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JZ,iT,i) -
     >                     ppiclf_drudtPlag(PPICLF_JZ,iT,i) )
         
            ! Diffusive Unsteady Heat Transfer starts here
            b_du = (1.63d0 - 0.92d0*erf(0.017*(Pe - 80.0d0)))
     >       *(1.0d0 - 0.4d0*exp(-Pe/16.0d0)-0.6d0*exp(-Pe**2/30.0d0))
         
            ! Kernel for diffusive unsteady Heat Transfer
            kernel_du=2.0d0*sqrt(dp/(vmag*time))*exp(-b_du*vmag*time/dp)

            qq_du = qq_du + factor_du*kernel_du*
     >                (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i))

        if(ppiclf_iprop(5,i) .eq. 7 .and.
     >     ppiclf_iprop(7,i) .eq. 1724) then
      open(unit=66,file='fort.66',position='append')   

      write(66,*) ppiclf_time, iT, time,
     >            factor_du, kernel_du,
     >            ppiclf_dTdtMixt(iT,i),
     >            ppiclf_dTdtPlag(iT,i),
     >            (ppiclf_dTdtMixt(iT,i) - ppiclf_dTdtPlag(iT,i)),
     >            factor_du*kernel_du*
     >            (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i)),
     >            qq_du,
     >            kernelVU, A, B, factor 

      flush(66)

      endif
          enddo

         iT = ppiclf_nTimeBH
         time = ppiclf_timeBH(iT)

         A  = (4.0d0*rpi*time*rnu/(dp**2))**(.25d0)
         B  = (rpi*(vmag**3)*(time**2)/ 
     >                 (dp*rnu*(fH**3)))**(.5d0)

         kernelVU = 0.5d0*factor*(A+B)**(-2)

         fvux = fvux + kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JX,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JX,iT,i) )
         fvuy = fvuy + kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JY,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JY,iT,i) )
         fvuz = fvuz + kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JZ,iT,i) -
     >                  ppiclf_drudtPlag(PPICLF_JZ,iT,i) )

         ! Diffusive Unsteady Heat Transfer starts here
         b_du = (1.63d0 - 0.92d0*erf(0.017*(Pe - 80.0d0)))
     >    *(1.0d0 - 0.4d0*exp(-Pe/16.0d0)-0.6d0*exp(-Pe**2/30.0d0))
         
         ! Kernel for diffusive unsteady Heat Transfer
         kernel_du = 0.5d0*2.0d0*sqrt(dp/(vmag*time))
     >                      *exp(-b_du*vmag*time/dp)
         qq_du = qq_du + factor_du*kernel_du*
     >             (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i))

        if(ppiclf_iprop(5,i) .eq. 7 .and.
     >     ppiclf_iprop(7,i) .eq. 1724) then
      open(unit=66,file='fort.66',position='append')   

      write(66,*) ppiclf_time, iT, time,
     >            factor_du, kernel_du,
     >            ppiclf_dTdtMixt(iT,i),
     >            ppiclf_dTdtPlag(iT,i),
     >            (ppiclf_dTdtMixt(iT,i) - ppiclf_dTdtPlag(iT,i)),
     >            factor_du*kernel_du*
     >            (ppiclf_dTdtMixt(iT,i) -ppiclf_dTdtPlag(iT,i)),
     >            qq_du,
     >            kernelVU, A, B, factor 

      flush(66)

      endif

      endif


      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Shift arrays for Viscous Unsteady Force
!
! See rocpart/PLAG_RFLU_ShiftUnsteadyData.F90
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_ShiftUnsteadyData
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iT

!
! Code:
!
      do i=1,ppiclf_npart
         do iT = ppiclf_nUnsteadyData,2,-1
            ppiclf_drudtMixt(PPICLF_JX,iT,i) = 
     >                      ppiclf_drudtMixt(PPICLF_JX,iT-1,i)
            ppiclf_drudtMixt(PPICLF_JY,iT,i) = 
     >                      ppiclf_drudtMixt(PPICLF_JY,iT-1,i)
            ppiclf_drudtMixt(PPICLF_JZ,iT,i) = 
     >                      ppiclf_drudtMixt(PPICLF_JZ,iT-1,i)

            ppiclf_drudtPlag(PPICLF_JX,iT,i) = 
     >                      ppiclf_drudtPlag(PPICLF_JX,iT-1,i)
            ppiclf_drudtPlag(PPICLF_JY,iT,i) = 
     >                      ppiclf_drudtPlag(PPICLF_JY,iT-1,i)
            ppiclf_drudtPlag(PPICLF_JZ,iT,i) = 
     >                      ppiclf_drudtPlag(PPICLF_JZ,iT-1,i)

            ppiclf_dTdtMixt(iT,i) = ppiclf_dTdtMixt(iT-1,i)
            ppiclf_dTdtPlag(iT,i) = ppiclf_dTdtPlag(iT-1,i)
            ppiclf_TMixt(iT,i)    = ppiclf_TMixt(iT-1,i)
            ppiclf_TPlag(iT,i)    = ppiclf_TPlag(iT-1,i)
          
         enddo
      enddo


      if (ppiclf_nTimeBH < ppiclf_nUnsteadyData) then
            ppiclf_nTimeBH = ppiclf_nTimeBH + 1
      endif

      do iT = ppiclf_nTimeBH,2,-1
            ppiclf_timeBH(it) = ppiclf_timeBH(iT-1) + ppiclf_dt
      enddo


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Update arrays for Viscous Unsteady Force for JT=1 (current time step)
!
! See libpicl/user_files/ppiclf_user_AddedMass.f
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_UpdatePlag(i,iStage)
!
      implicit none
!
      include "PPICLF"
!
      integer*4 i, iStage, jj
      real*8 SDrho
      real*8 ug,vg,wg
      real*8 up,vp,wp
      real*8 vgradrho
      real*8 dt, time_plot

!
! Code:
!
      SDrho = ppiclf_rprop(PPICLF_R_JRHSR,i)
     >         + ppiclf_y(PPICLF_JVX,i) * ppiclf_rprop(PPICLF_R_JPGCX,i)
     >         + ppiclf_y(PPICLF_JVY,i) * ppiclf_rprop(PPICLF_R_JPGCY,i)
     >         + ppiclf_y(PPICLF_JVZ,i) * ppiclf_rprop(PPICLF_R_JPGCZ,i)

      ! 03/11/2025 - Thierry - substantial derivative from Rocflu is
      !              weighted by \phi^g.
      ! d(rho^g phi^g)/dt = rho^g * d(phi^g)/dt + phi^g * d(rho^g)/dt
      !                   = phi^g * d(rho^g)/dt
      !
      !     d(rho^g)/dt   = SDrho = d(rho phi^g)/dt / phi^g
      SDrho = SDrho / (rphif)

      ! 03/23/2025 - TLJ - added extra term involving grad(rhog)
      vgradrho = vx*ppiclf_rprop(PPICLF_R_JRHOGX,i) +
     >           vy*ppiclf_rprop(PPICLF_R_JRHOGY,i) +
     >           vz*ppiclf_rprop(PPICLF_R_JRHOGZ,i)

      ug = ppiclf_rprop(PPICLF_R_JUX,i)
      vg = ppiclf_rprop(PPICLF_R_JUY,i)
      wg = ppiclf_rprop(PPICLF_R_JUZ,i)
      up = ppiclf_y(PPICLF_JVX,i)
      vp = ppiclf_y(PPICLF_JVY,i)
      wp = ppiclf_y(PPICLF_JVZ,i)

      ! D(rhog*ug)/Dt
      ppiclf_drudtMixt(PPICLF_JX,1,i) =
     >   ug*(SDrho+vgradrho) + rhof*ppiclf_rprop(PPICLF_R_JSDRX,i)
      ppiclf_drudtMixt(PPICLF_JY,1,i) =
     >   vg*(SDrho+vgradrho) + rhof*ppiclf_rprop(PPICLF_R_JSDRY,i)
      ppiclf_drudtMixt(PPICLF_JZ,1,i) =
     >   wg*(SDrho+vgradrho) + rhof*ppiclf_rprop(PPICLF_R_JSDRZ,i)

      ! d(rhog*up)/dt
      ppiclf_drudtPlag(PPICLF_JX,1,i) =
     >   up*SDrho + rhof*ppiclf_ydot(PPICLF_JVX,i)
      ppiclf_drudtPlag(PPICLF_JY,1,i) =
     >   vp*SDrho + rhof*ppiclf_ydot(PPICLF_JVY,i)
      ppiclf_drudtPlag(PPICLF_JZ,1,i) =
     >   wp*SDrho + rhof*ppiclf_ydot(PPICLF_JVZ,i)

      ppiclf_TMixt(1,i) = ppiclf_rprop(PPICLF_R_JT,i)
      ppiclf_TPlag(1,i) = ppiclf_y(PPICLF_JT,i)

      ! Set the old step equal to the first one initially
      if(ppiclf_TMixt(2,i) .eq. 0.0d0) then
        ppiclf_TMixt(2,i) = ppiclf_TMixt(1,i)
      endif

      if(ppiclf_TPlag(2,i) .eq. 0.0d0) then
        ppiclf_TPlag(2,i) = ppiclf_TPlag(1,i)
      endif

      if(iStage==1) then
        dt = 5.0/15.0*ppiclf_dt
        time_plot = ppiclf_time
      elseif(iStage==2) then
        dt = (5.0/15.0 + 8.0/15.0)*ppiclf_dt
        time_plot = 5.0/15.0*dt + ppiclf_time
      elseif(iStage==3) then
        dt = ppiclf_dt
        time_plot = 13.0/15.0*dt + ppiclf_time
      else
        print*, "Unsupported beyond RK3 in Viscous-Unsteady"
        call ppiclf_exittr('Unsupported iStage', 0.0, 0)
      endif

      ! dT/dt update for Diffusive Unsteady HT
      ! first-order backward difference (explicit backward Euler form)
        ppiclf_dTdtMixt(1,i) = (ppiclf_TMixt(1,i)-ppiclf_TMixt(2,i))/dt
        ppiclf_dTdtPlag(1,i) = (ppiclf_TPlag(1,i)-ppiclf_TPlag(2,i))/dt

        if(ppiclf_iprop(5,i) .eq. 7 .and.
     >     ppiclf_iprop(7,i) .eq. 1724) then
          open(unit=67,file='fort.67',position='append')   
          write(67,*) ppiclf_time, time_plot, dt, ppiclf_dt,
     >                ppiclf_dTdtMixt(1,i),     
     >                ppiclf_TMixt(1,i),
     >                ppiclf_TMixt(2,i),
     >                ppiclf_dTdtPlag(1,i),
     >                ppiclf_ydot(PPICLF_JT,i),
     >                ppiclf_TPlag(1,i),
     >                ppiclf_TPlag(2,i)
          flush(67)
      endif
!----------------------------------------------------------      

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Sets drudtMixt and drudtPlag from rprop3
! Needed for proper particle tracking
! Load communication buffers rprop3 into particle data
! See rocpart/PLAG_RFLU_ModComm.F90:
!     SUBROUTINE PLAG_RFLU_LoadBuffersSend(pRegion)
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_prop2plag
!
      implicit none
!
      include "PPICLF"
!
      integer*4 i,k,ic,iT
!
! Code:
!
      do i=1,ppiclf_npart
         k = 0
         do ic = 1,3
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_drudtMixt(ic,iT,i) = ppiclf_rprop3(k,i)
         enddo
         enddo
         do ic = 1,3
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_drudtPlag(ic,iT,i) = ppiclf_rprop3(k,i)
         enddo
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_dTdtMixt(iT,i) = ppiclf_rprop3(k,i)
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_dTdtPlag(iT,i) = ppiclf_rprop3(k,i)
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_TMixt(iT,i) = ppiclf_rprop3(k,i)
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_TPlag(iT,i) = ppiclf_rprop3(k,i)
         enddo
      enddo

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Sets rprop3 from drudtMixt, drudtPlag, dTdtMixt, and dTdtPlag
! Needed for proper particle tracking
! Load particle data into communication buffers rprop3
! See rocpart/PLAG_RFLU_ModComm.F90:
!     SUBROUTINE PLAG_RFLU_UnloadBuffersRecv(pRegion)
!
! Modified Oct. 24, 2025 - Added dTdtMixt and dTdtPlag for diffusive
!                          unsteady heat transfer model
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_plag2prop
!
      implicit none
!
      include "PPICLF"
!
      integer*4 i,k,ic,iT
!
! Code:
!
      do i=1,ppiclf_npart
         k = 0
         do ic = 1,3
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_drudtMixt(ic,iT,i)
         enddo
         enddo
         do ic = 1,3
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_drudtPlag(ic,iT,i)
         enddo
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_dTdtMixt(iT,i)
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_dTdtPlag(iT,i)
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_TMixt(iT,i)
         enddo
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_TPlag(iT,i)
         enddo

      enddo


      return
      end
