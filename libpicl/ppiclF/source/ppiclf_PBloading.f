!-----------------------------------------------------------------------
! The following subroutines are for the new particle-based load balance
! approach. The goal is to make ppiclf_npart approximately equal per
! processor.
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_CreateBinPartLB

      USE ppiclf_DynamicAllocation
      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
 
      INTEGER*4 i, j, ierr

      REAL*8    BinMinLen(3), local_extremes(6)
     >          ,BinBuffer(3), temp1, temp2, periodicDistCheck
     >          ,idum, jdum, kdum

      ! Right now, bin every RK stage. Infrequent binning probably
      ! a lot harder with particle-based load balancing
      ppiclf_binchanged = .TRUE.

      CALL MPI_ALLREDUCE(ppiclf_npart,ppiclf_glnpart,1
     >                   ,MPI_INTEGER4, MPI_SUM
     >                   ,ppiclf_comm, ierr)

      ! Bin must be larger than nearest neighbor search distance
      ! and the ppiclf_filter(1:3).  This makes a buffer around the bin
      ! domain. ppiclf_filter is 1.2x max cell length in each direction
      DO i = 1,3
        BinMinLen(i) = MAX(ppiclf_filter(i),ppiclf_nndist)
        ! For buffer: need ppiclf_filter to make sure you have 1 layer
        ! of outer fluid cells
        ! Need ppiclf_nndist/2 to ensure BinMinLen is never violated
        BinBuffer(i) = MAX(ppiclf_filter(i),ppiclf_nndist/2)
      END DO

      ! Looping through particles on this processor
      ! to find bin boundary locations
      DO i = 1,3
        local_extremes(2*i) = -1.0D10
        local_extremes(2*i-1) = 1.0D10
      END DO
      DO i=1,ppiclf_npart
        DO j = 1,3
          ! Finding min/max particle extremes.
          ! Add buffer so that layers of outer cells 
          ! are available for interpolation/projection.
          temp1 = ppiclf_y(j,i) - BinBuffer(j)
          temp2 = ppiclf_y(j,i) + BinBuffer(j)
          IF(temp1 .LT. local_extremes(2*j-1)) 
     >                               local_extremes(2*j-1) = temp1
          IF(temp2 .GT. local_extremes(2*j)) 
     >                               local_extremes(2*j) = temp2
        END DO
      END DO

      ! Finds global bin domain boundaries across MPI ranks
      DO i = 1,3
        CALL MPI_ALLREDUCE(local_extremes(2*i-1),ppiclf_binb(2*i-1), 1
     >                     ,MPI_DOUBLE_PRECISION, MPI_MIN
     >                     ,ppiclf_comm, ierr)

        CALL MPI_ALLREDUCE(local_extremes(2*i),ppiclf_binb(2*i), 1
     >                     ,MPI_DOUBLE_PRECISION, MPI_MAX
     >                     ,ppiclf_comm, ierr)
        ppiclf_BinDomLen(i) = ppiclf_binb(2*i) - ppiclf_binb(2*i-1)
      END DO

      ! Ensuring ppiclf_binb not greater than 
      ! cartesian fluid domain extremes.
      ! If dist within ppiclf_nndist, set ppiclf_binb
      ! equal to fluid domain for periodic ghost particles.
      ! Needed to know when to use linear periodic
      
      ppiclf_EqualDomain(1) = .FALSE.
      ppiclf_EqualDomain(2) = .FALSE.
      ppiclf_EqualDomain(3) = .FALSE.

      DO i = 1,3
        ! Check bin min domain
        periodicDistCheck = MAX(ppiclf_nndist,ppiclf_filter(i))
        IF(ppiclf_binb(i*2-1) - periodicDistCheck .LE. 
     >                          ppiclf_xdrange(1,i)) THEN
          ppiclf_binb(i*2-1) = ppiclf_xdrange(1,i)
          ppiclf_EqualDomain(i) = .TRUE.
        END IF
        ! Check bin max domain
        IF(ppiclf_binb(i*2)+periodicDistCheck .GE. 
     >                          ppiclf_xdrange(2,i)) THEN
          ppiclf_binb(i*2) = ppiclf_xdrange(2,i)
        ELSE
          ppiclf_EqualDomain(i) = .FALSE.
        END IF
      END DO

      ppiclf_totalBins = 1
      DO i = 1,3
        ppiclf_n_bins(i) = INT( (ppiclf_binb(2*i)-ppiclf_binb(2*i-1)) /
     >                          BinMinLen(i) )
        ppiclf_bins_dx(i) = (ppiclf_binb(2*i) - ppiclf_binb(2*i-1)) / 
     >                      ppiclf_n_bins(i)
        ppiclf_totalBins = ppiclf_totalBins*ppiclf_n_bins(i)
      END DO

      CALL ppiclf_dyn_alloc(ppiclf_totalBins,ppiclf_np)

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE ppiclf_comm_FindParticlePartLB

      USE ppiclf_DynamicAllocation
      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      INTEGER*4  i, ii, jj, kk, nbin, ierr, partcheck, NumBins
      EXTERNAL   ppiclf_iglmax
      INTEGER*4  ppiclf_iglmax

      DO i = 0,(ppiclf_totalBins - 1)
        ppiclf_ParticleCount(i) = 0
      END DO

      DO i = 1,ppiclf_npart
        ! Calculates particle's bin index
        ii  = FLOOR((ppiclf_y(1,i)-ppiclf_binb(1))/ppiclf_bins_dx(1))
        jj  = FLOOR((ppiclf_y(2,i)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
        kk  = FLOOR((ppiclf_y(3,i)-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
        
        ! Calculates particle's bin
        nbin = ii + ppiclf_n_bins(1)*jj + 
     >         ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk

        ! Maps particle to correct processor based on active bin number
        !ppiclf_iprop(4,i) = nrank ! Processor to send to
        ppiclf_iprop(5,i) = ii    ! x bin #
        ppiclf_iprop(6,i) = jj    ! y bin #
        ppiclf_iprop(7,i) = kk    ! z bin #
        ppiclf_iprop(8,i) = nbin ! total bin number
        ppiclf_ParticleCount(nbin) = ppiclf_ParticleCount(nbin) + 1
      END DO

      ! Now sum particles per bin across MPI Ranks
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_ParticleCount
     >                   ,ppiclf_totalBins ,MPI_INTEGER4, MPI_SUM
     >                   ,ppiclf_comm, ierr)

      CALL ppiclf_comm_partLoadBalance

      DO i =1,ppiclf_npart
        ! Now map particle to MPI Rank since we have a bin->rank map
        ii = ppiclf_iprop(8,i)
        ppiclf_iprop(4,i) = ppiclf_BinToRankMap(ii) ! Owning MPI Rank
      END DO


      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_partLoadBalance
   
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      INTEGER*4   ierr, i, j, k, targetParticleCnt, particleSum
     >           ,irank, bin, LbinCheck
     >           ,d2, d3, ii, jj, kk, nb1, nb2, nb3
     >           ,iloop, jloop, kloop
     >           ,iMinusFace, jMinusFace, kMinusFace
     >           ,iPlusFace, jPlusFace, kPlusFace
      REAL*8     rankBounds(6)

      CALL MPI_BARRIER(ppiclf_comm,ierr)
      ! Find the order to loop through dimensions. 
      ! Using largest dimension minimizes surface area between 
      ! processors, which likely minimizes ghost particle and
      ! overlap cell communication

      ! Sorting the domain lengths: ppiclf_dL = largest dimension,
      ! ppiclf_dM = medium dimension, ppiclf_dS = smallest dimension
      ppiclf_dL = MAXLOC(ppiclf_BinDomLen, DIM=1)
      ppiclf_dS = MINLOC(ppiclf_BinDomLen, DIM=1)
      ppiclf_dM = 6 - ppiclf_dL - ppiclf_dS

      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)

      targetParticleCnt = ppiclf_glnpart/ppiclf_np    

      particleSum   = 0
      irank = 0
      ! Save min rank boundary
      IF(irank .EQ. ppiclf_nid) THEN
        ppiclf_bin_pos(1,ppiclf_dL) = ppiclf_binb(2*ppiclf_dL-1) +
     >                                (0)*ppiclf_bins_dx(ppiclf_dL)
        ppiclf_bin_pos(1,ppiclf_dM) = ppiclf_binb(2*ppiclf_dM-1)  
        ppiclf_bin_pos(1,ppiclf_dS) = ppiclf_binb(2*ppiclf_dS-1)   
      END IF
      ! Iterate through loops of largest dimension.
      ! Increment loops by one and loop through all other dimensions.
      DO iloop = 0,(ppiclf_n_bins(ppiclf_dL) - 1)
        IF(ppiclf_dL .EQ. 1) THEN
          ii = iloop
        ELSE IF(ppiclf_dL .EQ. 2) THEN
          jj = iloop
        ELSE IF(ppiclf_dL .EQ. 3) THEN
          kk = iloop
        ELSE
          PRINT*, 'ERROR in ppiclf_comm_partLoadBalance'
          RETURN
        END IF
        DO jloop = 0,(ppiclf_n_bins(ppiclf_dM) - 1)
          IF(ppiclf_dM .EQ. 1) THEN
            ii = jloop
          ELSE IF(ppiclf_dM .EQ. 2) THEN
            jj = jloop
          ELSE IF(ppiclf_dM .EQ. 3) THEN
            kk = jloop
          ELSE
            PRINT*, 'ERROR in ppiclf_comm_partLoadBalance'
            RETURN
          END IF
          DO kloop = 0,(ppiclf_n_bins(ppiclf_dS) - 1)
            IF(ppiclf_dS .EQ. 1) THEN
              ii = kloop
            ELSE IF(ppiclf_dS .EQ. 2) THEN
              jj = kloop
            ELSE IF(ppiclf_dS .EQ. 3) THEN
              kk = kloop
            ELSE
              PRINT*, 'ERROR in ppiclf_comm_partLoadBalance'
              RETURN
            END IF
            bin  = ii + nb1*jj + nb1*nb2*kk
            particleSum = particleSum + ppiclf_ParticleCount(bin)
            ppiclf_LMapFluid(bin) = .FALSE.
            IF(ppiclf_ParticleCount(bin) .GT. 0) THEN
              ppiclf_LMapFluid(bin) = .TRUE.
            ELSE
        Loop: DO i = -1,1
                DO j = -1,1
                  DO k = -1,1
                    LbinCheck = (ii + i) + nb1*(jj + j) +
     >                          nb1*nb2*(kk + k)
                    IF(ppiclf_ParticleCount(LbinCheck) .GT. 0) THEN
                      ppiclf_LMapFluid(bin) = .TRUE.
                      EXIT Loop
                    END IF
                  END DO !k
                END DO !j
              END DO Loop !i
            END IF
            ! This maps the bin to the MPI-rank
            ppiclf_BinToRankMap(bin) = irank 
            IF(particleSum .GE. targetParticleCnt) THEN
              !ppiclf_RankBoundary(0:nbins_in-1,6)
              ! Set MPI Rank's maximum integers/reals
              IF(irank .EQ. ppiclf_nid) THEN
                ! Save max rank boundary
                ppiclf_bin_pos(2,ppiclf_dL) = ppiclf_binb(2*ppiclf_dL-1)
     >                             + (iloop+1)*ppiclf_bins_dx(ppiclf_dL)
                ppiclf_bin_pos(2,ppiclf_dM) = ppiclf_binb(2*ppiclf_dM)  
                ppiclf_bin_pos(2,ppiclf_dS) = ppiclf_binb(2*ppiclf_dS)
              END IF
              ! Save rank boundary indices
              ppiclf_IRankBoundary(irank,ppiclf_dL*2) = iloop
              ppiclf_IrankBoundary(irank,ppiclf_dM*2) = jloop
              ppiclf_IrankBoundary(irank,ppiclf_dS*2) = kloop

              ! Reset Particle counter and advance MPI rank iteration
              particleSum = 0
              irank = irank + 1
              IF(irank .GT. ppiclf_np - 1) THEN
                irank = ppiclf_np - 1
              ELSE
                ppiclf_IRankBoundary(irank,ppiclf_dL*2-1) = iloop
                ppiclf_IrankBoundary(irank,ppiclf_dM*2-1) = jloop
                ppiclf_IrankBoundary(irank,ppiclf_dS*2-1) = kloop

                IF(irank .EQ. ppiclf_nid) THEN
                  ppiclf_bin_pos(1,ppiclf_dL) = 
     >                           ppiclf_binb(2*ppiclf_dL-1) 
     >                           + (iloop+1)*ppiclf_bins_dx(ppiclf_dL)
                  ppiclf_bin_pos(1,ppiclf_dM) = 
     >                           ppiclf_binb(2*ppiclf_dM-1)
                  ppiclf_bin_pos(1,ppiclf_dS) = 
     >                           ppiclf_binb(2*ppiclf_dS-1)
                END IF
              END IF
            END IF
          END DO !kloop (shortest bin dimension)
        END DO !jloop (medium bin dimension)
      END DO !iloop (longest bin dimension)

      ppiclf_IRankBoundary(ppiclf_np-1,ppiclf_dL) =
     >                                       ppiclf_n_bins(ppiclf_dL) -1
      ppiclf_IRankBoundary(ppiclf_np-1,ppiclf_dM) =
     >                                       ppiclf_n_bins(ppiclf_dM) -1
      ppiclf_IRankBoundary(ppiclf_np-1,ppiclf_dS) =
     >                                       ppiclf_n_bins(ppiclf_dS) -1

      ! Store logical array which is True if bin face is on MPI Boundary
      ! Must loop through every bin
      iMinusFace = 0
      jMinusFace = 0
      kMinusFace = 0
      iPlusFace  = ppiclf_IRankBoundary(0,ppiclf_dL) 
      jPlusFace  = ppiclf_IRankBoundary(0,ppiclf_dM)
      jPlusFace  = ppiclf_IRankBoundary(0,ppiclf_dS)
      DO iloop = 0,(ppiclf_n_bins(ppiclf_dL) - 1)
        IF(ppiclf_dL .EQ. 1) THEN
          ii = iloop
        ELSE IF(ppiclf_dL .EQ. 2) THEN
          jj = iloop
        ELSE IF(ppiclf_dL .EQ. 3) THEN
          kk = iloop
        END IF
        DO jloop = 0,(ppiclf_n_bins(ppiclf_dM) - 1)
          IF(ppiclf_dM .EQ. 1) THEN
            ii = jloop
          ELSE IF(ppiclf_dM .EQ. 2) THEN
            jj = jloop
          ELSE IF(ppiclf_dM .EQ. 3) THEN
            kk = jloop
          END IF
          DO kloop = 0,(ppiclf_n_bins(ppiclf_dS) - 1)
            IF(ppiclf_dS .EQ. 1) THEN
              ii = kloop
            ELSE IF(ppiclf_dS .EQ. 2) THEN
              jj = kloop
            ELSE IF(ppiclf_dS .EQ. 3) THEN
              kk = kloop
            END IF
            bin  = ii + nb1*jj + nb1*nb2*kk
            irank = ppiclf_BinToRankMap(bin) 

            ! Shift Face index tracker if past MPI Boundary
            IF(iloop .GT. iPlusFace) THEN
              IF(ppiclf_IRankBoundary(irank,ppiclf_dL) 
     >                                 .NE. iPlusFace) THEN
                iMinusFace = iPlusFace + 1
              ELSE
                iMinusFace = iPlusFace
              END IF
              iPlusFace = ppiclf_IRankBoundary(irank,ppiclf_dL)
            END IF
            IF(jloop .GT. jPlusFace) THEN
              IF(ppiclf_IRankBoundary(irank,ppiclf_dM) 
     >                                 .NE. jPlusFace) THEN
                jMinusFace = jPlusFace + 1
              ELSE
                jMinusFace = jPlusFace
              END IF
              jMinusFace = jPlusFace + 1
              jPlusFace = ppiclf_IRankBoundary(irank,ppiclf_dM)
            END IF
            IF(kloop .GT. kPlusFace) THEN
              IF(ppiclf_IRankBoundary(irank,ppiclf_dS) 
     >                                 .NE. kPlusFace) THEN
                kMinusFace = kPlusFace + 1
              ELSE
                kMinusFace = kPlusFace
              END IF
              kMinusFace = kPlusFace + 1
              kPlusFace = ppiclf_IRankBoundary(irank,ppiclf_dS)
            END IF

            ! If + side of bin is a boundary
            IF(iloop .EQ. iMinusFace) THEN
              ppiclf_LRankBoundary(bin,2*ppiclf_dL-1) = .TRUE.
            ELSE
              ppiclf_LRankBoundary(bin,2*ppiclf_dL-1) = .FALSE.
            END IF
            IF(jloop .EQ. jMinusFace) THEN
              ppiclf_LRankBoundary(bin,2*ppiclf_dM-1) = .TRUE.
            ELSE
              ppiclf_LRankBoundary(bin,2*ppiclf_dM-1) = .FALSE.
            END IF
            IF(kloop .EQ. kMinusFace) THEN
              ppiclf_LRankBoundary(bin,2*ppiclf_dS-1) = .TRUE.
            ELSE
              ppiclf_LRankBoundary(bin,2*ppiclf_dS-1) = .FALSE.
            END IF

            ! If - side of bin is a boundary
            IF(iloop .EQ. iPlusFace) THEN
              ppiclf_LRankBoundary(bin,2*ppiclf_dL) = .TRUE.
            ELSE
              ppiclf_LRankBoundary(bin,2*ppiclf_dL) = .FALSE.
            END IF
            IF(jloop .EQ. jPlusFace) THEN
              ppiclf_LRankBoundary(bin,2*ppiclf_dM) = .TRUE.
            ELSE
              ppiclf_LRankBoundary(bin,2*ppiclf_dM) = .FALSE.
            END IF
            IF(kloop .EQ. kPlusFace) THEN
              ppiclf_LRankBoundary(bin,2*ppiclf_dS) = .TRUE.
            ELSE
              ppiclf_LRankBoundary(bin,2*ppiclf_dS) = .FALSE.
            END IF
          END DO !kloop (shortest bin dimension)
        END DO !jloop (medium bin dimension)
      END DO !iloop (longest bin dimension)

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_LBCheck

      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      INTEGER*4   ierr, i, rank
      INTEGER*4, ALLOCATABLE :: LB_Count(:)

      IF(ALLOCATED(LB_Count)) THEN
        DEALLOCATE(LB_Count)
      END IF

      ALLOCATE(LB_Count(0:ppiclf_np - 1))
 
      DO i = 0,ppiclf_np-1
        LB_Count(i) = 0
      END DO

      DO i = 1,ppiclf_npart
        rank = ppiclf_iprop(4,i)
        LB_Count(rank) = LB_Count(rank) + 1
      END DO 

      ! Now sum particles per bin across MPI Ranks
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, LB_Count
     >                   ,ppiclf_np ,MPI_INTEGER4, MPI_SUM
     >                   ,ppiclf_comm, ierr)

      IF(ppiclf_nid .EQ. 0) THEN
        DO i = 0,ppiclf_np-1
          PRINT*, 'Proc:',i,'Particles on Rank:',LB_Count(i)
        END DO
      END IF

      END SUBROUTINE

!----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_MoveParticlePartLB
!
! This subroutine is called from ppiclf_solve_InitSolve
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
!
! Internal:
!
      LOGICAL   partl ! dummy variable    
      INTEGER*4 rtempLim
      PARAMETER(rtempLim = PPICLF_LRS*4 + PPICLF_LRP + PPICLF_LRP2
     >       + PPICLF_LRP3 + PPICLF_LRP4 + PPICLF_LRP5 + PPICLF_LRP_PRO)
      REAL*8    rtemp(rtempLim,PPICLF_LPART)
      INTEGER*4 i, icount, j0
#ifdef PERF
      REAL*8    tstart, tfinal
#endif
!
      ! copy particle y, rprop, rprop2, rprop3 arrays into rtemp
      ! array for communication
      DO i=1,ppiclf_npart
        icount = 1
        CALL ppiclf_copy(rtemp(icount,i),ppiclf_y(1,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(rtemp(icount,i),ppiclf_y1(1,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(rtemp(icount,i),ppiclf_ydot(1,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(rtemp(icount,i),ppiclf_ydotc(1,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(rtemp(icount,i),ppiclf_rprop(1,i),PPICLF_LRP)
        icount = icount + PPICLF_LRP
        IF(PPICLF_LRP2 .GT. 1) THEN
          CALL ppiclf_copy(rtemp(icount,i),
     >                     ppiclf_rprop2(1,i),PPICLF_LRP2)
          icount = icount + PPICLF_LRP2
        END IF
        IF(PPICLF_LRP3 .GT. 1) THEN
          CALL ppiclf_copy(rtemp(icount,i),
     >                     ppiclf_rprop3(1,i),PPICLF_LRP3)
          icount = icount + PPICLF_LRP3
        END IF
        IF(PPICLF_LRP4 .GT. 1) THEN
          CALL ppiclf_copy(rtemp(icount,i),
     >                     ppiclf_rprop4(1,i),PPICLF_LRP4)
          icount = icount + PPICLF_LRP4
        END IF
        IF(PPICLF_LRP5 .GT. 1) THEN
          CALL ppiclf_copy(rtemp(icount,i),
     >                     ppiclf_rprop5(1,i),PPICLF_LRP5)
          icount = icount + PPICLF_LRP5
        END IF
        CALL ppiclf_copy(rtemp(icount,i),
     >                   ppiclf_feedbk(1,i),PPICLF_LRP_PRO)
      END DO
      
      j0 = 4 ! index of ppiclf_iprop that contains rank to send to

#ifdef PERF
      tstart = MPI_WTIME()
#endif

      CALL pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >             ,ppiclf_npart,PPICLF_LPART ! Setup
     >             ,ppiclf_iprop,PPICLF_LIP   ! Integer Comm
     >             ,partl,0                   ! Logical Comm
     >             ,rtemp,rtempLim            ! Real Comm
     >             ,j0)                       ! Receiver processor index

#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TDataTransfers = PPICLF_TDataTransfers + (tfinal - tstart)
#endif

      IF(ppiclf_npart .GT. PPICLF_LPART .OR. ppiclf_npart .LT. 0) THEN
        PRINT*,'Increase LPART. Processor:',ppiclf_nid,
     >   'LPART should be greater than:',ppiclf_npart
        CALL ppiclf_exittr('Increase LPART$',0.0d0,ppiclf_npart)
      END IF
 
      ! Update processor particle values with newly transfered rtemp
      ! array from communication
      DO i=1,ppiclf_npart
        icount = 1
        CALL ppiclf_copy(ppiclf_y(1,i),rtemp(icount,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(ppiclf_y1(1,i),rtemp(icount,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(ppiclf_ydot(1,i),rtemp(icount,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(ppiclf_ydotc(1,i),rtemp(icount,i),PPICLF_LRS)
        icount = icount + PPICLF_LRS
        CALL ppiclf_copy(ppiclf_rprop(1,i),rtemp(icount,i),PPICLF_LRP)
        icount = icount + PPICLF_LRP
        IF(PPICLF_LRP2 .GT. 1) THEN
        CALL ppiclf_copy(ppiclf_rprop2(1,i),rtemp(icount,i),
     >                   PPICLF_LRP2)
        icount = icount + PPICLF_LRP2
        END IF
        IF(PPICLF_LRP3 .GT. 1) THEN
          CALL ppiclf_copy(ppiclf_rprop3(1,i),rtemp(icount,i),
     >                     PPICLF_LRP3)
          icount = icount + PPICLF_LRP3
        END IF
        IF(PPICLF_LRP4 .GT. 1) THEN
          CALL ppiclf_copy(ppiclf_rprop4(1,i),rtemp(icount,i),
     >                     PPICLF_LRP4)
          icount = icount + PPICLF_LRP4
        END IF
        IF(PPICLF_LRP5 .GT. 1) THEN
          CALL ppiclf_copy(ppiclf_rprop5(1,i),rtemp(icount,i),
     >                     PPICLF_LRP5)
          icount = icount + PPICLF_LRP5
        END IF
        CALL ppiclf_copy(ppiclf_feedbk(1,i),rtemp(icount,i),
     >           PPICLF_LRP_PRO)
      END DO
        
      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_MapOverlapGridPartLB
!
! This subroutine is called from ppiclf_solve_InitSolve
!
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE 'mpif.h'
!
! Internal:
!
      INTEGER*4 icalld
      SAVE      icalld
      DATA      icalld /0/
      INTEGER*4 nkey(2), i, j, k, l, ie, iee, ii, jj, kk, irank,
     >          nl, nii, njj, nrr, iic, jjc, kkc, ierr 
      INTEGER*4 ix, iy, iz, ixLow, ixHigh, iyLow,
     >          iyHigh, izLow, izHigh, ibin, jbin, kbin, nbin,
     >          nRankMaps, RankMaps(27,4)
      REAL*8    rxval, ryval, rzval, 
     >          MinPoint(3), ppiclf_vlmin, ppiclf_vlmax,
     >          centeri(3), Max_CellLen(3) 
      LOGICAL   partl, ErrorFound, MapCell
      EXTERNAL  ppiclf_vlmin, ppiclf_vlmax
#ifdef PERF
      REAL*8    tstart, tfinal
#endif
!
! Code Start:
!
      ! Number of fluid finite volume cells that map to particle bins
      ppiclf_nCells_FV2PICL = 0 
      ! Loops through number of fluid FV cells on this processor
      DO ie=1,ppiclf_nFVCells  
        centeri(1) = ppiclf_fluid_grid(1,ie)
        centeri(2) = ppiclf_fluid_grid(2,ie)
        centeri(3) = ppiclf_fluid_grid(3,ie)
        ! Exits if fluid cell center is outside of any bin boundaries 
        IF (centeri(1) .GT. ppiclf_binb(2)) CYCLE
        IF (centeri(2) .GT. ppiclf_binb(4)) CYCLE
        IF (centeri(3) .GT. ppiclf_binb(6)) CYCLE
        IF (centeri(1) .LT. ppiclf_binb(1)) CYCLE
        IF (centeri(2) .LT. ppiclf_binb(3)) CYCLE
        IF (centeri(3) .LT. ppiclf_binb(5)) CYCLE
        ! Determines what bin the fluid cell is nominally mapped to
        ibin    = FLOOR((centeri(1)-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
        jbin    = FLOOR((centeri(2)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
        kbin    = FLOOR((centeri(3)-ppiclf_binb(5))/ppiclf_bins_dx(3))
        ! Calculates processor rank
        nbin  = ibin + ppiclf_n_bins(1)*jbin
     >          + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kbin
        ! Will loop through cell mapping once if all below stay as 0
        ixLow  = 0
        ixHigh = 0
        iyLow  = 0
        iyHigh = 0
        izLow  = 0
        izHigh = 0
        ! Change loop bounds if the bin is on a MPI Boundary Face
        IF(ppiclf_LRankBoundary(nbin,1)) ixLow  = -1
        IF(ppiclf_LRankBoundary(nbin,3)) iyLow  = -1
        IF(ppiclf_LRankBoundary(nbin,5)) izLow  = -1
        IF(ppiclf_LRankBoundary(nbin,2)) ixHigh =  1
        IF(ppiclf_LRankBoundary(nbin,4)) iyHigh =  1
        IF(ppiclf_LRankBoundary(nbin,5)) izHigh =  1
        nRankMaps = 0
        RankMaps  = -99
        DO ix=ixLow,ixHigh
          DO iy=iyLow,iyHigh
            DO iz=izLow,izHigh
              ! Adjust bin if on MPI Boundary
              ii = ibin + ix 
              jj = jbin + iy
              kk = kbin + iz
              ! This covers ghost exchanged cells for linear periodicity
              ! Maps cells greater than ppiclf bin domain to first bin
              ! Maps cells less than ppiclf bin domain to last bin
              IF(ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)
     >                               .AND. ppiclf_n_bins(1) .GT. 1) THEN
                IF(ii .EQ. ppiclf_n_bins(1)) ii = 0
                IF(ii .EQ. -1) ii = ppiclf_n_bins(1) - 1
              END IF
              IF(ppiclf_linperiodic(2) .AND. ppiclf_EqualDomain(2)
     >                               .AND. ppiclf_n_bins(2) .GT. 1) THEN
                IF(jj .EQ. ppiclf_n_bins(2)) jj = 0
                IF(jj .EQ. -1) jj = ppiclf_n_bins(2) - 1
              END IF
              IF(ppiclf_linperiodic(3) .AND. ppiclf_EqualDomain(3)
     >                               .AND. ppiclf_n_bins(3) .GT. 1) THEN
                IF(kk .EQ. ppiclf_n_bins(3)) kk = 0
                IF(kk .EQ. -1) kk = ppiclf_n_bins(3) - 1
              END IF
              
              ! Ensures duplicate cells don't get sent to same processor
              IF (ii .LT. 0 .OR. ii .GT. ppiclf_n_bins(1)-1) CYCLE
              IF (jj .LT. 0 .OR. jj .GT. ppiclf_n_bins(2)-1) CYCLE
              IF (kk .LT. 0 .OR. kk .GT. ppiclf_n_bins(3)-1) CYCLE
              nbin  = ii + ppiclf_n_bins(1)*jj
     >                + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
              irank = ppiclf_BinToRankMap(nbin)
              IF(ppiclf_LMapFluid(nbin)) THEN 
              ! Checks if bin is near particles (LMapFluid = True)
                MapCell = .TRUE.
                IF(nRankMaps .GT. 0) THEN
                ! Checks if bin is already mapped (irank not in list)
                  DO i = 1,nRankMaps
                    IF(RankMaps(i,1) .EQ. irank) THEN
                      MapCell = .FALSE.
                      EXIT
                    END IF 
                  END DO !i
                END IF !nRankMaps > 0
                IF(MapCell) THEN
                  nRankMaps = nRankMaps + 1
                  RankMaps(nRankMaps,1) = irank
                  RankMaps(nRankMaps,2) = ii
                  RankMaps(nRankMaps,3) = jj
                  RankMaps(nRankMaps,4) = kk
                END IF !MapCell
              END IF !ppiclf_LMapFluid
            END DO !iz
          END DO !iy
        END DO !ix
        IF(nRankMaps .GT. 0) THEN
          DO i = 1,nRankMaps
!            PRINT*, 'proc, i, ie, rank, ii, jj, kk, nmap',ppiclf_nid,
!     >              i,ie,
!     >              RankMaps(i,1),RankMaps(i,2),
!     >              RankMaps(i,3),RankMaps(i,4),
!     >              nRankMaps
            ppiclf_nCells_FV2PICL = ppiclf_nCells_FV2PICL + 1
            IF(ppiclf_nCells_FV2PICL .GT. PPICLF_LEE) THEN
              PRINT*, '***ERROR*** PPICLF_LEE',PPICLF_LEE, 'in', 
     >         'MapOverlapGrid must be greater than',
     >          ppiclf_nCells_FV2PICL 
              CALL ppiclf_exittr('Increase PPICLF_LEE$ (MapOverlap)',0.0D0
     >             ,ppiclf_nCells_FV2PICL)
            END IF
            ! Stores cells to rank mapping.
            ! Fluid solver cell ID
            ppiclf_cell_map(1,ppiclf_nCells_FV2PICL) = ie
            ! Fluid solver cell rank
            ppiclf_cell_map(2,ppiclf_nCells_FV2PICL) = ppiclf_nid
            ! Particle solver cell rank and bin indicies
            ppiclf_cell_map(3,ppiclf_nCells_FV2PICL) = RankMaps(i,1)
            ppiclf_cell_map(4,ppiclf_nCells_FV2PICL) = RankMaps(i,2)
            ppiclf_cell_map(5,ppiclf_nCells_FV2PICL) = RankMaps(i,3)
            ppiclf_cell_map(6,ppiclf_nCells_FV2PICL) = RankMaps(i,4)
          END DO !nRankMaps
        END IF
      END DO !ie

      DO ie=1,ppiclf_nCells_FV2PICL 
        ! These copy all indicies since Fortran is column-major
        iee = ppiclf_cell_map(1,ie)
        CALL ppiclf_copy(ppiclf_picl_grid(1,ie)
     >                 ,ppiclf_fluid_grid(1,iee),7)
 
        ! ppiclf_filter initially set in PICL_TEMP_InitSolver
        ! Want to only consider cells that reside in the particle domain
        ! Update ppiclf_filter for next binning cycle 1.2*dx since next
        ! layer of cells in a growing particle domain may be slightly larger
        ! and we want filter equal a minimum of 1 cell.
        DO l = 1,3
          IF(ppiclf_filter(l) .LT. ppiclf_fluid_grid(3+l,iee)) THEN
            ppiclf_filter(l) = 1.2D0*ppiclf_fluid_grid(3+l,iee)
          END IF
        END DO
      END DO

      ! Find max filter across processors
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_filter
     >                   ,3 ,MPI_DOUBLE_PRECISION
     >                   ,MPI_MAX ,ppiclf_comm, ierr)

      ! Copy mapping since it is need to send fluid properties in interp
      ppiclf_nCells_FV2PICL_Orig = ppiclf_nCells_FV2PICL
      DO ie=1,ppiclf_nCells_FV2PICL_Orig
         ! Copies cells to rank mapping (integer copy)
         CALL ppiclf_icopy(ppiclf_cell_map_Orig(1,ie)
     >            ,ppiclf_cell_map(1,ie),PPICLF_LRMAX)
      END DO

      ! GSLIB required info
      ! NumPiclCells - number of columns to transfer
      ! PPICLF_LEE - number of columns declared
      ! nl - partl row size (dummy logical variable)
      nl   = 0
      ! nii - ppiclf_cell_map row size declared
      nii  = PPICLF_LRMAX
      ! njj - Row index of ppiclf_cell_map with receiver processor/rank
      njj  = 3
      ! nrr - ppiclf_rocGrid row size declared
      nrr  = 7
      ! Defines sorting order
      nkey(1) = 2
      nkey(2) = 1

#ifdef PERF
      tstart = MPI_WTIME()
#endif

      CALL pfgslib_crystal_tuple_transfer(
     >        ppiclf_cr_hndl,ppiclf_nCells_FV2PICL,PPICLF_LEE !setup
     >        ,ppiclf_cell_map,nii ! Integer Comm
     >        ,partl,nl                 ! Logical Comm
     >        ,ppiclf_picl_grid,nrr      ! Real Comm
     >        ,njj)                      ! Receiver processor index
      CALL pfgslib_crystal_tuple_sort(
     >        ppiclf_cr_hndl,ppiclf_nCells_FV2PICL !setup
     >        ,ppiclf_cell_map,nii !Integer to sort
     >        ,partl,nl                 !Logical to sort
     >        ,ppiclf_picl_grid,nrr      !Real to sort
     >        ,nkey,2)                  !sorting method

#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TDataTransfers = PPICLF_TDataTransfers + (tfinal - tstart)
#endif

      ! Find distance check for interpolation.
      ! This is 1.5*MaxCellLength to ensure that at least
      ! 27 neighboring cells are mapped.
      Max_CellLen(1) = 0.0D0
      Max_CellLen(2) = 0.0D0
      Max_CellLen(3) = 0.0D0
      DO ie = 1,ppiclf_nCells_FV2PICL ! Loop through cells mapped to bin
        DO l = 1,3
          ! Find max cell lengths in all dimensions
          IF(ppiclf_picl_grid(3+l,ie) .GT. Max_CellLen(l))
     >      Max_CellLen(l) = ppiclf_picl_grid(3+l,ie)
        END DO !l
      END DO !ie
      DO l = 1,3
        ! Multiply by 1.5 so particle near face will
        ! find center one cell over in farthest direction
        ppiclf_interp_dchk(l) = Max_CellLen(l)*1.5D0
      END DO

      IF (icalld .EQ. 0) THEN 
         icalld = icalld + 1
         CALL ppiclf_prints('   *Begin mpi_comm_split$')
            CALL mpi_comm_split(ppiclf_comm
     >                         ,ppiclf_nid
     >                         ,0
     >                         ,ppiclf_comm_nid
     >                         ,ierr)
         CALL ppiclf_prints('    End mpi_comm_split$')
         CALL ppiclf_io_OutputDiagGrid
      END IF

      RETURN
      END SUBROUTINE 
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_CreateGhostPartLB
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
!
! Internal:
!
      REAL*8     GhostPos(3), PeriodicShift(3), 
     >           distSQ(3), distCheckSQ, buffer
      INTEGER*4  ip, idum, iip, jjp, kkp, iig, jjg, kkg, nrank, 
     >           j, k, l, GhostInc(3), ix, iy, iz
!
      ! Calculate the linear periodicity shift in each dimension
      DO l = 1,3
        IF(ppiclf_linperiodic(l)) THEN
          PeriodicShift(l) = ppiclf_xdrange(2,l) - ppiclf_xdrange(1,l)
        ELSE
          PeriodicShift(l) = 0.0D0
        END IF
      END DO

      ppiclf_npart_gp = 0
      
      DO ip=1,ppiclf_npart
        idum = 0
        ! Copy particle solution variables
        DO j=1,PPICLF_LRS
           idum = idum + 1
           ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip)
        END DO
        ! Copy particle property variables
        DO j=1,PPICLF_LRP
           idum = idum + 1
           ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip)
        END DO

        ! GP Bin Index
        iip    = ppiclf_iprop(5,ip)
        jjp    = ppiclf_iprop(6,ip)
        kkp    = ppiclf_iprop(7,ip)
   
        ! Found that buffer was needed in unit testing 
        ! due to round-off errors with periodicity
        buffer = 1.02 
        distCheckSQ = (ppiclf_nndist*(buffer**3))**2

        DO ix = 1,3
          distSQ = 0.0D0
          GhostPos(1) = ppiclf_cp_map(1,ip)
          IF(ix .LT. 3) THEN
            CALL ppiclf_comm_GhostDistCheck(ix,GhostPos(1),
     >                     ppiclf_nndist*buffer,GhostInc(1),1,distSQ(1))
            IF(GhostInc(1) .EQ. 0) CYCLE
          ELSE
            GhostInc(1) = 0 !For ghosts in other 2 dimensions only
          END IF
          iig = iip + GhostInc(1)

          ! Angular Periodicity Check
          ! *** Add here ***

          ! If ghost is outside of ppiclf domain:
          IF(iig .LT. 0 .OR. iig .GT. ppiclf_n_bins(1)-1) THEN
            IF(ppiclf_linperiodic(1).AND. ppiclf_EqualDomain(1)) THEN
              CALL ppiclf_comm_LinearPeriodicityGhost
     >                       (iig,1,GhostPos(1),PeriodicShift(1))
            ELSE
              CYCLE
            END IF
          END IF

          DO iy = 1,3
            GhostPos(2) = ppiclf_cp_map(2,ip)
            IF(iy .LT. 3) THEN
              CALL ppiclf_comm_GhostDistCheck(iy,GhostPos(2),
     >                    ppiclf_nndist*buffer,GhostInc(2),2,distSQ(2))
              IF(GhostInc(2) .EQ. 0.) CYCLE
              ! This corner/edge check caused issues when unit testing
              ! Removing this makes extra ghost particles
              !IF(distSQ(1)+distSQ(2) .GT. distCheckSQ) CYCLE !corner/edge check
            ELSE
              GhostInc(2) = 0 !For ghosts in other 2 dimensions only
            END IF
            jjg = jjp + GhostInc(2)

          ! Angular Periodicity Check
          ! *** Add here ***

            ! If ghost is outside of ppiclf domain:
            IF(jjg .LT. 0 .OR. jjg .GT. ppiclf_n_bins(2)-1) THEN
              IF(ppiclf_linperiodic(2) .AND.
     >                     ppiclf_EqualDomain(2)) THEN
                CALL ppiclf_comm_LinearPeriodicityGhost(jjg,2,
     >                           GhostPos(2),PeriodicShift(2))
              ELSE
                CYCLE
              END IF
            END IF

            DO iz = 1,3
              GhostPos(3) = ppiclf_cp_map(3,ip)
              IF(iz .LT. 3) THEN
                CALL ppiclf_comm_GhostDistCheck(iz,GhostPos(3),
     >                    ppiclf_nndist*buffer,GhostInc(3),3,distSQ(3))
                IF(GhostInc(3) .EQ. 0) CYCLE
                ! This corner/edge check caused issues when unit testing
                ! Removing this makes extra ghost particles
                !corner/edge check
                !IF(distSQ(1)+distSQ(2)+distSQ(3) .GT. distCheckSQ) CYCLE
              ELSE
                GhostInc(3) = 0
              END IF
              kkg = kkp + GhostInc(3)

          ! Angular Periodicity Check
          ! *** Add here ***              

              ! If ghost is outside of ppiclf domain:
              IF(kkg .LT. 0 .OR. kkg .GT. ppiclf_n_bins(3)-1) THEN
                IF(ppiclf_linperiodic(3) .AND.
     >                        ppiclf_EqualDomain(3)) THEN
                  CALL ppiclf_comm_LinearPeriodicityGhost
     >                           (kkg,3,GhostPos(3),PeriodicShift(3))
                ELSE
                  CYCLE
                END IF
              END IF

              ! This prevents ghost in same bin.
              IF(GhostInc(1) .EQ. 0 .AND. GhostInc(2) .EQ. 0 .AND.
     >           GhostInc(3) .EQ. 0 ) CYCLE
              
              ! Add ghost particle and map integer and real properties
              nrank = iig + ppiclf_n_bins(1)*jjg 
     >               + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
              !ghostsMade = ghostsMade + 1
              ppiclf_npart_gp = ppiclf_npart_gp + 1
              ! Copy particle ID info
              ppiclf_iprop_gp(1,ppiclf_npart_gp) = ppiclf_iprop(1,ip)
              ppiclf_iprop_gp(2,ppiclf_npart_gp) = ppiclf_iprop(2,ip)
              ppiclf_iprop_gp(3,ppiclf_npart_gp) = ppiclf_iprop(3,ip)
              ppiclf_iprop_gp(4,ppiclf_npart_gp) = nrank !*** change to processor
              ppiclf_iprop_gp(5,ppiclf_npart_gp) = iig
              ppiclf_iprop_gp(6,ppiclf_npart_gp) = jjg
              ppiclf_iprop_gp(7,ppiclf_npart_gp) = kkg
              ppiclf_iprop_gp(8,ppiclf_npart_gp) = nrank

              ppiclf_rprop_gp(1,ppiclf_npart_gp) = GhostPos(1)
              ppiclf_rprop_gp(2,ppiclf_npart_gp) = GhostPos(2)
              ppiclf_rprop_gp(3,ppiclf_npart_gp) = GhostPos(3)

              DO k=4,PPICLF_LRP_GP
                ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
              END DO
            END DO !iz = 1:3
          END DO !iy = 1:3
        END DO !ix = 1:3
      END DO !ip = 1:ppiclf_npart

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_GhostDistCheckPartLB(ix,Pos,distchk,
     >                                      GhostInc,l,dSQ)
      
      IMPLICIT NONE
      
      INCLUDE "PPICLF"

      ! ix: ghostcheck loop counter, GhostInc: bin +/-, l: dimenison
      INTEGER*4 ix, GhostInc, l
      ! Pos: Position of Ghost Particle, distchk: criteria to create
      ! ghost particle
      ! distSQ: used to evaluate distance ghost in edge & corner case
      REAL*8    Pos, distchk, dSQ

      ! ppiclf_bin_pos(1,1) is bin min position in x
      ! ppiclf_bin_pos(2,1) is bin max position in x
      IF(ABS(Pos - ppiclf_bin_pos(ix,l)) 
     >                          .LT. distchk) THEN
        dSQ = (Pos-ppiclf_bin_pos(ix,l))**2
        IF(ix .EQ. 1) GhostInc = -1 ! close to bin min
        IF(ix .EQ. 2) GhostInc =  1 ! clost to bin max
      ELSE
        GhostInc = 0
      END IF
      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_LinearPeriodicityGhostPartLB
     >           (iig,l,Pos,PerShift)
      
      IMPLICIT NONE
      
      INCLUDE "PPICLF"

      !iig: Ghost bin index, l: dimension Number (1:x,2:y,3:z)
      INTEGER*4 iig, l 
      ! Pos: GhostPos(l), PerShift: PeriodicShift(l)
      REAL*8    Pos, PerShift
      IF(iig .LT. 0) THEN
        iig = ppiclf_n_bins(l) - 1
        Pos = Pos + PerShift
      ELSE IF (iig .GT. ppiclf_n_bins(l) - 1) THEN
        iig = 0
        Pos = Pos - PerShift
      END IF

      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_MoveGhostPartLB
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
!
! Internal:
!
      INTEGER*4 iprop_proc_index
      LOGICAL   partl  ! Dummy variable       

#ifdef PERF
      REAL*8    tstart, tfinal
#endif

!
      iprop_proc_index = 4 ! since ppiclf_iprop(4,np) contains processor
                           ! that should receive ghost particle

#ifdef PERF
      tstart = MPI_WTIME()
#endif

      CALL pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >             ,ppiclf_npart_gp,PPICLF_LPART_GP ! Setup
     >             ,ppiclf_iprop_gp,PPICLF_LIP_GP   ! Integer Comm
     >             ,partl,0                         ! Logical Comm
     >             ,ppiclf_rprop_gp,PPICLF_LRP_GP   ! Real Comm
     >             ,iprop_proc_index)               ! Receiver processor index

#ifdef PERF
      tfinal = MPI_WTIME()
      PPICLF_TDataTransfers = PPICLF_TDataTransfers + (tfinal - tstart)
#endif

      RETURN
      END SUBROUTINE

