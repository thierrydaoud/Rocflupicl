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
      subroutine ppiclf_user_VU_Rocflu(i,iStage,fvux,fvuy,fvuz)
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

!
! Code:
!
      fvux = 0.0d0
      fvuy = 0.0d0
      fvuz = 0.0d0
      iT   = 1
      time = 0.0d0

      fH     = 0.75d0 + .105d0*reyL
      ! Sangani's volume fraction correction for dilute random arrays
      ! Capping volume fraction at 0.5 
      factor = 3.0d0*rpi*rnu*dp*ppiclf_dt*(1.0+2.28*min(rphip,0.5))

      if (ppiclf_nTimeBH > 1) then
         do iT = 2,ppiclf_nTimeBH-1
            time = ppiclf_timeBH(iT)

            A  = (4.0d0*rpi*time*rnu/(dp**2))**(.25d0)
            B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >                 (0.5d0*dp*rnu*(fH**3)))**(.5d0)

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
         enddo

         iT = ppiclf_nTimeBH
         time = ppiclf_timeBH(iT)

         A  = (4.0d0*rpi*time*rnu/(dp**2))**(.25d0)
         B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >                 (0.5d0*dp*rnu*(fH**3)))**(.5d0)

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
      subroutine ppiclf_user_UpdatePlag(i)
!
      implicit none
!
      include "PPICLF"
!
      integer*4 i
      real*8 SDrho
      real*8 ug,vg,wg
      real*8 up,vp,wp
      real*8 vgradrho

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
! Sets rprop3 from drudtMixt and drudtPlag
! Needed for proper particle tracking
! Load particle data into communication buffers rprop3
! See rocpart/PLAG_RFLU_ModComm.F90:
!     SUBROUTINE PLAG_RFLU_UnloadBuffersRecv(pRegion)
!
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
      enddo


      return
      end
