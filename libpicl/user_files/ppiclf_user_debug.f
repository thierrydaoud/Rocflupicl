!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for output if ppiclf_debug=1
! fort.72## is reserved for debug
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_debug
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Internal:
!
      integer*4 i, n, ic, k, iStage

! Needed for allreduce
      integer*4 ngop
      parameter(ngop = 32)
      real*8 xin(ngop),wout(ngop)

! Needed for viscous unsteady
      integer*4 iT,ii
      real*8 time, factor, A, B, fH, kernelVU
      real*8 FVUoutput
      real*8 ppiclf_npart_sum
      integer*4 npart_tot, npart_max, npart_min
      integer*4 ppiclf_iglsum,ppiclf_iglmax,ppiclf_iglmin
      external  ppiclf_iglsum,ppiclf_iglmax,ppiclf_iglmin

      iStage = 1 ! for internal use only

!
! Code:
!
      ! Use ppiclf ALLREDUCE to compute values across processors
      ! Note that ALLREDUCE uses MPI_BARRIER, which is cpu expensive
      ! Print out every 10th iStage=1 counts

      !xin(1) = dfloat(ppiclf_npart)
      !call ppiclf_gop(xin, wout, '+  ', 1)
      !ppiclf_npart_sum = wout(1)
      npart_tot = ppiclf_iglsum(PPICLF_NPART,1)
      npart_max = ppiclf_iglmax(PPICLF_NPART,1)
      npart_min = ppiclf_iglmin(PPICLF_NPART,1)


      xin=(/phimax,
     >         fqsx_max,fqsy_max,fqsz_max,
     >         famx_max,famy_max,famz_max, 
     >         fdpdx_max,fdpdy_max,fdpdz_max, 
     >         fcx_max,fcy_max,fcz_max,
     >         umean_max,vmean_max,wmean_max,
     >         fqs_mag,fam_mag,fdp_mag,fc_mag,
     >         fqsx_fluct_max,fqsy_fluct_max,fqsz_fluct_max,
     >         fqsx_total_max,fqsy_total_max,fqsz_total_max,
     >         fvux_max,fvuy_max,fvuz_max,
     >         qq_max,tau_max,lift_max/)
      call ppiclf_gop(xin, wout, 'M  ', ngop)
      phimax     = wout(1)
      fqsx_max   = wout(2)
      fqsy_max   = wout(3)
      fqsz_max   = wout(4)
      famx_max   = wout(5)
      famy_max   = wout(6)
      famz_max   = wout(7)
      fdpdx_max  = wout(8)
      fdpdy_max  = wout(9)
      fdpdz_max  = wout(10)
      fcx_max    = wout(11)
      fcy_max    = wout(12)
      fcz_max    = wout(13)
      umean_max  = wout(14)
      vmean_max  = wout(15)
      wmean_max  = wout(16)
      fqs_mag    = wout(17)
      fam_mag    = wout(18)
      fdp_mag    = wout(19)
      fc_mag     = wout(20)
      fqsx_fluct_max = wout(21)
      fqsy_fluct_max = wout(22)
      fqsz_fluct_max = wout(23)
      fqsx_total_max = wout(24)
      fqsy_total_max = wout(25)
      fqsz_total_max = wout(26)
      fvux_max   = wout(27)
      fvuy_max   = wout(28)
      fvuz_max   = wout(29)
      qq_max     = wout(30)
      tau_max    = wout(31)
      lift_max   = wout(32)

      ! Sam - logging for debugging purposes
      ! TLJ - below is a mess I created, need to clean up
      if (ppiclf_nid.eq.0) then

         if (ViscousUnsteady_flag>=1) then
            fH     = 0.75d0 + .105d0*reyL
            factor = 3.0d0*rpi*rnu*dp*fac
            FVUoutput = 0.0
            if (ppiclf_nTimeBH>1) then
               do iT = 2,ppiclf_nTimeBH-1
                  time = ppiclf_timeBH(iT)
                  A  = (4.0d0*rpi*time*rnu/dp**2)**(.25d0)
                  B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >              (0.5d0*dp*rnu*fH**3))**(.5d0)
                  kernelVU = factor*(A+B)**(-2)
                  FVUoutput = FVUoutput + kernelVU*
     >               (ppiclf_drudtMixt(1,iT,1)-ppiclf_drudtPlag(1,iT,1))
                  if (abs(FVUoutput) < 1.d-20) FVUoutput = 0.0d0
               enddo
               iT = ppiclf_nTimeBH
               time = ppiclf_timeBH(iT)
               A  = (4.0d0*rpi*time*rnu/dp**2)**(.25d0)
               B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >              (0.5d0*dp*rnu*fH**3))**(.5d0)
               kernelVU = 0.5*factor*(A+B)**(-2)
               FVUoutput = FVUoutput + kernelVU*
     >             (ppiclf_drudtMixt(1,iT,1)-ppiclf_drudtPlag(1,iT,1))
            endif

            WRITE(7225,"(700(1x,E14.6))") ppiclf_time, 
     >        ((ppiclf_drudtMixt(1,i,1)-ppiclf_drudtPlag(1,i,1))
     >        ,i=1,6)
            WRITE(7227,"(2i5,2x,700(2x,E14.6))") iT,iStage,ppiclf_time,
     >        time,FVUoutput,A,B,kernelVU
            WRITE(7229,"(i5,2x,i5,2x,800(1x,E14.6))") ppiclf_nTimeBH,
     >          ppiclf_nUnsteadyData,ppiclf_dt,
     >          ppiclf_time,ppiclf_TimeBH(1:6)

         endif

         WRITE(7226,"(700(1x,E14.6))") ppiclf_time,
     >        ((ppiclf_drudtMixt(1,i,1)-ppiclf_drudtPlag(1,i,1))
     >        ,i=1,ppiclf_nUnsteadyData)
         WRITE(7228,"(70(1x,E14.6))") ppiclf_time,
     >        ((ppiclf_drudtMixt(3,i,1)-ppiclf_drudtPlag(3,i,1))
     >        ,i=1,ppiclf_nUnsteadyData)

         WRITE(7230,"(27(1x,E23.16))") ppiclf_time, ppiclf_y(1:12, 1)

         WRITE(7231,"(28(1x,E23.16))") ppiclf_time,phimax,
     >             fqsx_max, fqsy_max, fqsz_max,
     >             famx_max, famy_max, famz_max,
     >             fdpdx_max, fdpdy_max, fdpdz_max,
     >             fcx_max, fcy_max, fcz_max,
     >             qq_max,tau_max,lift_max,
     >             fqsx_total_max,fqsy_total_max,fqsz_total_max,
     >             fvux_max, fvuy_max, fvuz_max
         WRITE(7232,"(26(1x,F13.8))") ppiclf_time,
     >             umean_max,vmean_max,wmean_max
         WRITE(7233,"(i5,2x,28(1x,E23.16))")
     >             ppiclf_nid,ppiclf_dt,ppiclf_time,
     >             fac, phimax,
     >             fqsx_fluct_max, fqsy_fluct_max, fqsz_fluct_max
         WRITE(7234,*) ppiclf_nid,istage,PPICLF_LRS ,PPICLF_LPART,
     >             PPICLF_NPART,ppiclf_time,
     >             ppiclf_rprop(PPICLF_R_FLUCTFX:PPICLF_R_FLUCTFZ,1),
     >             ppiclf_ydotc(PPICLF_JVX:PPICLF_JT,1)
         WRITE(7235,"(26(1x,F13.8))") ppiclf_time,
     >             fqs_mag,fam_mag,fdp_mag,fc_mag
         WRITE(7236,"(26(1x,F13.8))") ppiclf_time,
     >             fcx_max, fcy_max, fcz_max
         WRITE(7237,"(26(1x,F13.8))") ppiclf_time,UnifRnd
         WRITE(7240,"(26(1x,F13.8))") ppiclf_time,
     >             fqsx_max, fqsy_max, fqsz_max,
     >             fqsx_fluct_max, fqsy_fluct_max, fqsz_fluct_max,
     >             fqsx_total_max,fqsy_total_max,fqsz_total_max
         WRITE(7241,"(5(1x,E23.16))") ppiclf_time,
     >             phipmean, upmean, vpmean, wpmean
         WRITE(7243,"(1x,E23.16,5(2x,I8))") ppiclf_time,
     >             npart_tot,npart_max,npart_min

         do i = 1,4
            write(7250+i,*) ppiclf_time,
     >              ppiclf_rprop(PPICLF_R_JSDRX:PPICLF_R_JSDRZ,i), ! Du/Dt
     >              ppiclf_rprop(PPICLF_R_JSDOX:PPICLF_R_JSDOZ,i)  ! DOmega/Dt
         enddo


      endif


      return
      end
