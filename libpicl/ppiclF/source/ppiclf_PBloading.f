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

      INTEGER*4 i, j, ierr, temp_dSize(3), itemp
      REAL*8    BinMinLen(3), local_extremes(6)
     >          ,global_extremes(6)
     >          ,BinBuffer(3), periodicDistCheck
     >          ,idum, jdum, kdum
#ifdef PERF
      REAL*8    tstart, tfinal
#endif

      DO i = 1,3
        BinMinLen(i) = MAX(ppiclf_filter(i), ppiclf_nndist)
        BinBuffer(i) = 2.0D0*MAX(ppiclf_filter(i), ppiclf_nndist)
      END DO
      IF(ppiclf_istage .NE. 2 .OR. ppiclf_istage .NE. 3) THEN
        local_extremes(1:3) =  1.0D10 ! Large "min" value
        local_extremes(4:6) = -1.0D10 ! Small "max" value

        DO i=1,ppiclf_npart
          local_extremes(1) = MIN(local_extremes(1), 
     >                            ppiclf_y(1,i) - BinBuffer(1))
          local_extremes(2) = MIN(local_extremes(2), 
     >                            ppiclf_y(2,i) - BinBuffer(2))
          local_extremes(3) = MIN(local_extremes(3), 
     >                            ppiclf_y(3,i) - BinBuffer(3))
          local_extremes(4) = MAX(local_extremes(4), 
     >                            ppiclf_y(1,i) + BinBuffer(1))
          local_extremes(5) = MAX(local_extremes(5), 
     >                            ppiclf_y(2,i) + BinBuffer(2))
          local_extremes(6) = MAX(local_extremes(6), 
     >                            ppiclf_y(3,i) + BinBuffer(3))
         END DO
        
        ! flip sign on min so that I only need to call ALLREDUCE
        ! for max values and not two expensive MPI calls.
        DO i = 1,3
          local_extremes(i) = - local_extremes(i)
        END DO
#ifdef PERF
      tstart = MPI_WTIME()
#endif
        ! Finds global bin domain boundaries across MPI ranks
        CALL MPI_ALLREDUCE(local_extremes, global_extremes, 6
     >                       ,MPI_DOUBLE_PRECISION, MPI_MAX
     >                       ,ppiclf_comm, ierr)
#ifdef PERF
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_allreduces = PPICLF_TMPI_allreduces + tfinal
      PPICLF_TCreateBin = PPICLF_TCreateBin - tfinal
#endif
        ! flip sign on min values back to positive
        DO i = 1,3
          global_extremes(i) = - global_extremes(i)
        END DO
        DO i = 1,3
          ppiclf_binb(2*i-1)  = global_extremes(i)
          ppiclf_binb(2*i)    = global_extremes(i+3)
          ppiclf_BinDomLen(i) = ppiclf_binb(2*i) - ppiclf_binb(2*i-1)
        END DO
      END IF

      IF(ppiclf_glnpart .LT. 1) THEN
        PRINT*, 'ERROR: PPICLF RAN WITH ZERO PARTICLES'
        CALL ppiclf_exittr('',0.0,0)
      END IF
      ppiclf_binchanged = .FALSE.

      ! If all particles are within last RK Stage binboundaries,
      ! do not calculate bins again
      DO i = 1,3
        IF((ppiclf_binb(2*i-1) + 0.4D0*BinMinLen(i)) .LT.
     >             ppiclf_previousbinb(2*i-1)) THEN
          ppiclf_binchanged = .TRUE.
          EXIT
        END IF
        IF((ppiclf_binb(2*i)   - 0.4D0*BinMinLen(i)) .GT.
     >             ppiclf_previousbinb(2*i))   THEN
          ppiclf_binchanged = .TRUE.
          EXIT
        END IF
      END DO

#ifdef TEST
      ppiclf_binchanged = .TRUE.
#endif

      IF(.NOT. ppiclf_binchanged) THEN
        ! Reset due to possible lack of full buffer
        DO i = 1,3
          ppiclf_binb(2*i-1)  = ppiclf_previousbinb(2*i-1)
          ppiclf_binb(2*i)    = ppiclf_previousbinb(2*i)
          ppiclf_BinDomLen(i) = ppiclf_binb(2*i) - ppiclf_binb(2*i-1)
        END DO
        ! Other data stays consistent - return from subroutine
        RETURN
      END IF

      ! Since bins changed, need make new BTRM and move particles
      ppiclf_particleMoved = .TRUE.
      ppiclf_rebalance     = .TRUE.

! Never print bins, since it isn't one per rank...
! the bin print out doesn't make much sense anymore
! Maybe re-attack later if we want to see the full binning picture?
! Will require printing new information in the printbin subroutine.
      ppiclf_printbinvtu = .FALSE.

      ppiclf_EqualDomain(1:3) = .FALSE.

      DO i = 1,3
        ! Check bin min domain
        periodicDistCheck = MAX(ppiclf_nndist, ppiclf_filter(i))
        IF( (ppiclf_binb(i*2-1) - periodicDistCheck .LE. 
     >                              ppiclf_xdrange(1,i))
     >      .AND.
     >      (ppiclf_binb(i*2) + periodicDistCheck .GE. 
     >                              ppiclf_xdrange(2,i)) ) THEN
          ppiclf_binb(i*2-1) = ppiclf_xdrange(1,i)
          ppiclf_binb(i*2)   = ppiclf_xdrange(2,i)
          ppiclf_EqualDomain(i) = .TRUE.
        END IF
        ! Ensure particle domain bin boundaries do not exceed fluid boundaries
        ppiclf_binb(i*2-1) = MAX(ppiclf_binb(i*2-1),ppiclf_xdrange(1,i))
        ppiclf_binb(i*2)   = MIN(ppiclf_binb(i*2),ppiclf_xdrange(2,i))
        ! Update Bin Domain Length
        ppiclf_BinDomLen(i) = ppiclf_binb(2*i) - ppiclf_binb(2*i-1) 
      END DO

      ! Set previous bin boundaries for next RK Stage check
      ppiclf_previousbinb(1:6) = ppiclf_binb(1:6)

      ppiclf_totalBins = 1
      DO i = 1,3
        ppiclf_n_bins(i)  = MAX( 1,
     >                      FLOOR(ppiclf_BinDomLen(i)/BinMinLen(i)))
        ppiclf_bins_dx(i) = ppiclf_BinDomLen(i)/DBLE(ppiclf_n_bins(i))
        ppiclf_totalBins  = ppiclf_totalBins * ppiclf_n_bins(i)
      END DO

      ! Perform on root processor only to ensure no rounding errors
      IF(ppiclf_binorderset) THEN
        ! Not the first time sorting    
        ! Sorting the domain lengths
        ! Don't want to constantly flip - do has to be 1.25x bigger
        temp_dSize(1) = ppiclf_dL
        temp_dSize(2) = ppiclf_dM
        temp_dSize(3) = ppiclf_dS
        DO i = 1,2
          DO j = i+1,3
            IF(ppiclf_BinDomLen(temp_dSize(j)) .GT.
     >         ppiclf_BinDomLen(temp_dSize(i))*1.25) THEN
              itemp = temp_dSize(i)
              temp_dSize(i) = temp_dSize(j)
              temp_dSize(j) = itemp
            END IF
          END DO !j
        END DO !i
       ELSE
        ! First time sorting    
        ! Sorting the domain lengths
        temp_dSize(1) = 1
        temp_dSize(2) = 2
        temp_dSize(3) = 3
        DO i = 1,2
          DO j = i+1,3
            IF(ppiclf_BinDomLen(temp_dSize(j)) .GT.
     >         ppiclf_BinDomLen(temp_dSize(i))) THEN
              itemp = temp_dSize(i)
              temp_dSize(i) = temp_dSize(j)
              temp_dSize(j) = itemp
            END IF
          END DO !j
        END DO !i
      END IF !binorderset

      ppiclf_dL = temp_dSize(1) 
      ppiclf_dM = temp_dSize(2) 
      ppiclf_dS = temp_dSize(3)
      ppiclf_binorderset = .TRUE.
      ! Allocate all arrays dependant on number of bins or processors
      CALL ppiclf_dyn_alloc(ppiclf_totalBins, ppiclf_np)

      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_FindParticlePartLB

      USE ppiclf_DynamicAllocation

      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      INTEGER*4  i, j, ii, jj, kk, bin, ierr
      INTEGER*4  nb1, nb2, nb3, nb1xnb2, prevBin
      INTEGER*4  iloop, jloop, kloop, nRank
      INTEGER*4  stride(3), stride_L, stride_M, stride_S
      INTEGER*4  bin_L, bin_M
      REAL*8     LB_iterationCount, LB_criteria, LB_target, inv_dx(3)
#ifdef PERF
      REAL*8     tstart, tfinal
#endif
     
      nb1  = ppiclf_n_bins(1)
      nb2  = ppiclf_n_bins(2)
      nb3  = ppiclf_n_bins(3)
      nb1xnb2 = nb1 * nb2

      ppiclf_ParticleCount = 0
      IF(.NOT. ppiclf_binchanged) THEN
        ppiclf_particleMoved = .FALSE.
      END IF

#ifdef TEST
      ppiclf_particleMoved = .TRUE.
#endif

      ! Pre-compute inverse of dx to avoid slow division in loop
      inv_dx(1) = 1.0D0 / ppiclf_bins_dx(1)
      inv_dx(2) = 1.0D0 / ppiclf_bins_dx(2)
      inv_dx(3) = 1.0D0 / ppiclf_bins_dx(3)

      ! Loop through particles
      DO i = 1,ppiclf_npart
        ! Use fast multiplication instead of division
        ii = FLOOR((ppiclf_y(1,i)-ppiclf_binb(1)) * inv_dx(1))
        jj = FLOOR((ppiclf_y(2,i)-ppiclf_binb(3)) * inv_dx(2)) 
        kk = FLOOR((ppiclf_y(3,i)-ppiclf_binb(5)) * inv_dx(3)) 

        ii = MAX(0, MIN(ii, nb1-1))
        jj = MAX(0, MIN(jj, nb2-1))
        kk = MAX(0, MIN(kk, nb3-1))

        bin = ii + nb1*jj + nb1xnb2*kk

        ! Maps particle to correct processor based on active bin number
        ppiclf_iprop(5,i) = ii    ! x bin #
        ppiclf_iprop(6,i) = jj    ! y bin #
        ppiclf_iprop(7,i) = kk    ! z bin #
        ppiclf_iprop(8,i) = bin   ! total bin number
        IF(.NOT. ppiclf_binchanged) THEN
          ! This is based on previous BTRM
          ! Skip if bins changed since BTRM will update
          nRank = ppiclf_BinToRankMap(bin) 
          IF(nRank .NE. ppiclf_iprop(4,i)) THEN
            ppiclf_particleMoved = .TRUE.
            ppiclf_iprop(4,i) = nRank
          END IF
        END IF
        ppiclf_ParticleCount(bin) = ppiclf_ParticleCount(bin) + 1
      END DO
      ! Sum particles per bin across MPI Ranks
#ifdef PERF
      tstart = MPI_WTIME()
#endif
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_ParticleCount
     >                   ,ppiclf_totalBins ,MPI_INTEGER4, MPI_SUM
     >                   ,ppiclf_comm, ierr)
#ifdef PERF
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_allreduces = PPICLF_TMPI_allreduces + tfinal
      PPICLF_TFindPart = PPICLF_TFindPart - tfinal
#endif

      ppiclf_glnpart = SUM(ppiclf_ParticleCount)

      IF(ppiclf_binchanged) THEN
        ppiclf_rebalance     = .TRUE.
        ppiclf_particleMoved = .TRUE.
        RETURN
      END IF

      ! Logical OR comparison across MPI Ranks
#ifdef PERF
      tstart = MPI_WTIME()
#endif
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_particleMoved
     >                   ,1, MPI_LOGICAL, MPI_LOR
     >                   ,ppiclf_comm, ierr)
#ifdef PERF
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_allreduces = PPICLF_TMPI_allreduces + tfinal
      PPICLF_TFindPart = PPICLF_TFindPart - tfinal
#endif

      ! Check if new BTRM required
      ! If NumPart > LB_criteria*TargetNumPart -> Reassign BTRM
      LB_criteria = 1.2D0
      LB_target   = CEILING(DBLE(ppiclf_glnpart)/DBLE(ppiclf_np))
      LB_criteria = LB_criteria*LB_target
      ppiclf_rebalance = .FALSE.

      stride(1) = 1
      stride(2) = nb1
      stride(3) = nb1xnb2

      stride_L = stride(ppiclf_dL)
      stride_M = stride(ppiclf_dM)
      stride_S = stride(ppiclf_dS)

      LB_iterationCount = 0.0D0 
      bin = 0 

      outer_loop: DO iloop = 0,(ppiclf_n_bins(ppiclf_dL) - 1)
        bin_L = iloop * stride_L
        
        DO jloop = 0,(ppiclf_n_bins(ppiclf_dM) - 1)
          bin_M = bin_L + jloop * stride_M
          
          DO kloop = 0,(ppiclf_n_bins(ppiclf_dS) - 1)
            prevBin = bin
            bin = bin_M + kloop * stride_S

            IF(ppiclf_BinToRankMap(prevBin) .NE.
     >         ppiclf_BinToRankMap(bin)       ) THEN
              LB_iterationCount = 0.0D0 ! new Rank
            END IF

            LB_iterationCount = LB_iterationCount 
     >                          + DBLE(ppiclf_ParticleCount(bin))

            IF(LB_iterationCount .GT. LB_criteria) THEN
              ppiclf_rebalance = .TRUE.
              EXIT outer_loop
            END IF

          END DO
        END DO
      END DO outer_loop

      ! Logical OR comparison across MPI Ranks
#ifdef PERF
      tstart = MPI_WTIME()
#endif
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_rebalance
     >                   ,1, MPI_LOGICAL, MPI_LOR
     >                   ,ppiclf_comm, ierr)
#ifdef PERF
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_allreduces = PPICLF_TMPI_allreduces + tfinal
      PPICLF_TFindPart = PPICLF_TFindPart - tfinal
#endif

      ! Ensure remap particles since BTRM will change
      IF(ppiclf_rebalance) ppiclf_particleMoved = .TRUE.

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_PartLoadBalance

      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
#ifdef PERF
      REAL*8 ppiclf_pt0, tstart, tfinal
#endif


      INTEGER*4   ierr, i, bin, irank, particleSum
     >           ,targetParticleCnt, prevParticleSum
     >           ,nb1, nb2, nb1xnb2, j, k
     >           ,iloop, jloop, kloop, remainingParticles
      INTEGER*4   stride(3), stride_L, stride_M, stride_S
      INTEGER*4   bin_L, bin_M
      REAL*8      max_dp

      ! Calculate BTRM on root processor and broadcast
      ! to all others. This ensures all processors have same
      ! key global mapping.
#ifdef PERF
      ppiclf_pt0 = MPI_WTIME()
#endif
      IF(ppiclf_nid. EQ. 0) THEN
        
        nb1 = ppiclf_n_bins(1)
        nb2 = ppiclf_n_bins(2)
        nb1xnb2 = nb1 * nb2

        ! Map dimensions to strides to avoid IF statements inside the loops.
        stride(1) = 1
        stride(2) = nb1
        stride(3) = nb1xnb2

        stride_L = stride(ppiclf_dL)
        stride_M = stride(ppiclf_dM)
        stride_S = stride(ppiclf_dS)

        targetParticleCnt =CEILING(DBLE(ppiclf_glnpart)/DBLE(ppiclf_np))
        remainingParticles = ppiclf_glnpart 
        particleSum = 0
        irank = 0
        !Iterate through all bins using the fast, branchless stride logic
        DO iloop = 0,(ppiclf_n_bins(ppiclf_dL) - 1)
          bin_L = iloop * stride_L
          DO jloop = 0,(ppiclf_n_bins(ppiclf_dM) - 1)
            bin_M = bin_L + jloop * stride_M
            DO kloop = 0,(ppiclf_n_bins(ppiclf_dS) - 1)
              bin = bin_M + kloop * stride_S
              prevParticleSum = particleSum
              particleSum = particleSum + ppiclf_ParticleCount(bin)
              IF(particleSum .GE. targetParticleCnt) THEN
                IF(irank .LT. ppiclf_np - 1) THEN
                  ! Check if the previous state (undershoot) was better
                  ! than the current state (overshoot)
                  IF((targetParticleCnt - prevParticleSum) .LE.
     >               (particleSum - targetParticleCnt)         ) THEN

                    remainingParticles = remainingParticles 
     >                                   - prevParticleSum
                    ! Assign the current bin to the NEXT rank.
                    irank = irank + 1
                    ppiclf_BinToRankMap(bin) = irank
                    ! Start the new rank's count with the current bin
                    particleSum = ppiclf_ParticleCount(bin) 
                  ELSE
                    ! Keep this bin on the current rank.
                    remainingParticles = remainingParticles
     >                                   - particleSum
                    ppiclf_BinToRankMap(bin) = irank
                    irank = irank + 1
                    ! Reset particle counter for the next rank
                    particleSum = 0 
                  END IF
                  targetParticleCnt = CEILING(DBLE(remainingParticles) /
     >                                        DBLE(ppiclf_np - irank))
                ELSE
                  ! This is the last rank, so can't increase ranks
                  ppiclf_BinToRankMap(bin) = irank
                END IF
              ELSE
                ! Ideal number of particles per rank not met yet
                ppiclf_BinToRankMap(bin) = irank
              END IF
            END DO !kloop
          END DO !jloop
        END DO !iloop

      END IF ! root Processor

      ! Share BTRM to all processors 
#ifdef PERF
      tstart = MPI_WTIME()
#endif
      CALL MPI_BCAST(ppiclf_BinToRankMap,ppiclf_totalBins,MPI_INTEGER4,
     >               0, ppiclf_comm, ierr)
#ifdef PERF
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_allreduces = PPICLF_TMPI_allreduces + tfinal
      PPICLF_TLoadBalance = PPICLF_TLoadBalance - tfinal
#endif

      ! Assign correct rank to each particle
      DO i = 1,ppiclf_npart
        ! Now map particle to MPI Rank since we have a bin->rank map
        bin = ppiclf_iprop(8,i)
        ppiclf_iprop(4,i) = ppiclf_BinToRankMap(bin) ! Owning MPI Rank
      END DO

#ifdef PERF
      PPICLF_TLoadBalance = PPICLF_TLoadBalance
     >     + (MPI_WTIME() - ppiclf_pt0)
#endif
#ifdef PERF
      ppiclf_pt0 = MPI_WTIME()
#endif
      CALL ppiclf_comm_setRankBoundaries
#ifdef PERF
      PPICLF_TRankBounds = PPICLF_TRankBounds
     >     + (MPI_WTIME() - ppiclf_pt0)
#endif
#ifdef PERF
      ppiclf_pt0 = MPI_WTIME()
#endif
      CALL ppiclf_comm_setEmptyIndicator
#ifdef PERF
      PPICLF_TEmptyInd = PPICLF_TEmptyInd
     >     + (MPI_WTIME() - ppiclf_pt0)
#endif
#ifdef PERF
      ppiclf_pt0 = MPI_WTIME()
#endif
      CALL ppiclf_comm_setInterfaceIndicator
#ifdef PERF
      PPICLF_TInterfaceInd = PPICLF_TInterfaceInd
     >     + (MPI_WTIME() - ppiclf_pt0)
#endif

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_setRankBoundaries
      
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"


      INTEGER*4 ii, jj, kk, i, iBin, nb1, nb2, nb3
     >         ,min_i, max_i, min_j, max_j, min_k, max_k
      LOGICAL owned_in_row, owned_in_layer

      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)

      ! Initialize with extreme impossible values
      min_i = nb1
      max_i = -1
      min_j = nb2
      max_j = -1
      min_k = nb3
      max_k = -1
      
      iBin = 0
      DO kk = 0, (nb3 - 1)
        owned_in_layer = .FALSE. ! Reset for this entire Z-layer
        DO jj = 0, (nb2 - 1)
          owned_in_row = .FALSE. ! Reset for this X-row
          DO ii = 0, (nb1 - 1)
            ! If THIS rank owns this bin, expand the i-bounds
            IF(ppiclf_BinToRankMap(iBin) .EQ. ppiclf_nid) THEN
              min_i = MIN(min_i, ii)
              max_i = MAX(max_i, ii)
              owned_in_row = .TRUE.
            END IF
            iBin = iBin + 1
          END DO
          ! Evaluate j-bounds ONLY once per row
          IF(owned_in_row) THEN
            min_j = MIN(min_j, jj)
            max_j = MAX(max_j, jj)
            owned_in_layer = .TRUE.
          END IF
        END DO
        ! Evaluate k-bounds ONLY once per layer
        IF(owned_in_layer) THEN
          min_k = MIN(min_k, kk)
          max_k = MAX(max_k, kk)
        END IF
      END DO

      ppiclf_binBIndex(1) = min_i
      ppiclf_binBIndex(2) = max_i
      ppiclf_binBIndex(3) = min_j
      ppiclf_binBIndex(4) = max_j
      ppiclf_binBIndex(5) = min_k
      ppiclf_binBIndex(6) = max_k

      ! Now set the physical continuous bounding box positions
      IF(max_i .GE. 0) THEN 
        ! This rank owns at least one bin
        ppiclf_bin_pos(1,1) = ppiclf_binb(1) + min_i*ppiclf_bins_dx(1)
        ppiclf_bin_pos(2,1) = ppiclf_binb(1) +
     >                        (max_i+1)*ppiclf_bins_dx(1)

        ppiclf_bin_pos(1,2) = ppiclf_binb(3) + min_j*ppiclf_bins_dx(2)
        ppiclf_bin_pos(2,2) = ppiclf_binb(3) +
     >                        (max_j+1)*ppiclf_bins_dx(2)

        ppiclf_bin_pos(1,3) = ppiclf_binb(5) + min_k*ppiclf_bins_dx(3)
        ppiclf_bin_pos(2,3) = ppiclf_binb(5) + 
     >                        (max_k+1)*ppiclf_bins_dx(3)
      ELSE
        ! This rank has zero particles/bins assigned (idle processor)
        DO i = 1,3
          ppiclf_bin_pos(1,i) = 0.0D0
          ppiclf_bin_pos(2,i) = 0.0D0
        
          ppiclf_binBIndex(2*i)   = 0
          ppiclf_binBIndex(2*i-1) = 0 
        END DO
      END IF

      ppiclf_total_SBin = 1
      DO i = 1,3
        ppiclf_binOffset(i) = MAX(0,ppiclf_binBIndex(2*i-1)-1)
        ppiclf_nSBin(i) = MIN(ppiclf_n_bins(i) - 1,
     >                        ppiclf_binBIndex(2*i) + 1)
     >                  - MAX(0,ppiclf_binBIndex(2*i-1)-1) + 1
        ppiclf_total_SBin = ppiclf_total_SBin*ppiclf_nSBin(i)
      END DO
      PPICLF_T_LocalBins = ppiclf_total_SBin 
      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_setEmptyIndicator
  
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"


      INTEGER*4  iBin, imin, imax, jmin, jmax, kmin, kmax
     >           ,ii, jj, kk, nb1, nb2, nb3, nb1xnb2 
     >           ,di, dj, dk, ni, nj, nk, neighborBin
      LOGICAL    wrap_i, wrap_j, wrap_k
      LOGICAL, ALLOCATABLE, SAVE :: LMapFluid_prev(:)

      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)
      nb1xnb2 = nb1*nb2
      
      wrap_i = ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)
      wrap_j = ppiclf_linperiodic(2) .AND. ppiclf_EqualDomain(2)
      wrap_k = ppiclf_linperiodic(3) .AND. ppiclf_EqualDomain(3)
      
      iBin = 0
      DO kk = 0, nb3-1
        DO jj = 0, nb2-1
          DO ii = 0, nb1-1
            ppiclf_LMapFluid(iBin) = .FALSE.
            IF(ppiclf_ParticleCount(iBin) .GT. 0) THEN
              ppiclf_LMapFluid(iBin) = .TRUE.
            ELSE
              ! See if neighboring bin has at least 1 particle
              search_loop: DO dk = -1,1
                nk = kk + dk          
                ! k out of bounds or periodicity check
                IF(nk .LT. 0 .OR. nk .GT. nb3-1) THEN
                  IF(wrap_k) THEN
                    IF(nk .GT. nb3-1) nk = 0
                    IF(nk .LT.  0)    nk = nb3-1
                  ELSE
                    CYCLE
                  END IF
                END IF          
                
                DO dj = -1,1
                  nj = jj + dj
                  ! j out of bounds or periodicity check
                  IF(nj .LT. 0 .OR. nj .GT. nb2-1) THEN
                    IF(wrap_j) THEN
                      IF(nj .GT. nb2-1) nj = 0
                      IF(nj .LT.  0)    nj = nb2-1
                    ELSE
                      CYCLE
                    END IF
                  END IF
                  DO di = -1,1
                    ni = ii + di
                    ! i out of bounds or periodicity check
                    IF(ni .LT. 0 .OR. ni .GT. nb1-1) THEN
                      IF(wrap_i) THEN
                        IF(ni .GT. nb1-1) ni = 0
                        IF(ni .LT.  0)    ni = nb1-1
                      ELSE
                        CYCLE
                      END IF
                    END IF

                    neighborBin = ni + nb1*nj + nb1xnb2*nk
                    IF(ppiclf_ParticleCount(neighborBin) .GT. 0) THEN
                      ppiclf_LMapFluid(iBin) = .TRUE.
                      EXIT search_loop  ! Instantly breaks dk, dj, di loops
                    END IF

                  END DO !di
                END DO !dj
              END DO search_loop !dk
            END IF
            iBin = iBin + 1 ! Fast linear tracking
          END DO
        END DO
      END DO

      ! --- Flag whether bin occupancy changed since the prev stage ---
      ! ppiclf_LMapFluid is globally identical on every rank (built from
      ! the already-allreduced ppiclf_ParticleCount), so this is a purely
      ! local comparison -- no extra collective. The flag lets
      ! ppiclf_solve_PostTimeStepPartLB skip the overlap-grid / cell
      ! remaps on stages where no bin flips occupied<->empty.
      IF(.NOT. ALLOCATED(LMapFluid_prev)) THEN
        ppiclf_emptyChanged = .TRUE.
      ELSE IF(SIZE(LMapFluid_prev) .NE. SIZE(ppiclf_LMapFluid)) THEN
        ppiclf_emptyChanged = .TRUE.
      ELSE
        ppiclf_emptyChanged =
     >        ANY(ppiclf_LMapFluid .NEQV. LMapFluid_prev)
      END IF

      ! Refresh the saved copy for next stage. Reallocate only when the
      ! bin count changed (a binchanged stage, where the remaps rerun
      ! regardless of this flag).
      IF(ALLOCATED(LMapFluid_prev)) THEN
        IF(SIZE(LMapFluid_prev) .NE. SIZE(ppiclf_LMapFluid))
     >    DEALLOCATE(LMapFluid_prev)
      END IF
      IF(.NOT. ALLOCATED(LMapFluid_prev))
     >  ALLOCATE(LMapFluid_prev(0:SIZE(ppiclf_LMapFluid)-1))
      LMapFluid_prev = ppiclf_LMapFluid

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_setInterfaceIndicator

      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"


      INTEGER*4 nb1, nb2, nb3, nb1xnb2 
      INTEGER*4 ii, jj, kk, iBin, iRank
      INTEGER*4 di, dj, dk, ni, nj, nk, neighborBin, nRank
      LOGICAL   wrap_i, wrap_j, wrap_k

      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)
      nb1xnb2 = nb1 * nb2

      wrap_i = ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)
      wrap_j = ppiclf_linperiodic(2) .AND. ppiclf_EqualDomain(2)
      wrap_k = ppiclf_linperiodic(3) .AND. ppiclf_EqualDomain(3)
 
      iBin = 0
      DO kk = 0, nb3-1
        DO jj = 0, nb2-1
          DO ii = 0, nb1-1
            iRank = ppiclf_BinToRankMap(iBin) 
            ppiclf_LRankBoundary(iBin, 1:6) = .FALSE.
            ! Single pass through all 26 neighbors
            DO dk = -1, 1
              nk = kk + dk
              ! Out of bounds and periodicity check
              IF(nk .LT. 0 .OR. nk .GT. nb3-1) THEN
                IF(wrap_k) THEN
                  IF(nk .GT. nb3-1) nk = 0
                  IF(nk .LT.  0)    nk = nb3-1
                ELSE
                  CYCLE
                END IF
              END IF

              DO dj = -1, 1
                nj = jj + dj
                ! Out of bounds and periodicity check
                IF(nj .LT. 0 .OR. nj .GE. nb2) THEN
                  IF(wrap_j) THEN
                    IF(nj .GT. nb2-1) nj = 0
                    IF(nj .LT.  0)    nj = nb2-1
                  ELSE
                    CYCLE
                  END IF
                END IF

                DO di = -1, 1
                  ! Skip the current bin
                  IF(di .EQ. 0 .AND. dj. EQ. 0 .AND. dk .EQ. 0) CYCLE
                  ni = ii + di
                  ! Out of bounds and periodicity check
                  IF(ni .LT. 0 .OR. ni .GE. nb1) THEN
                    IF(wrap_i) THEN
                      IF(ni .GT. nb1-1) ni = 0
                      IF(ni .LT.  0)    ni = nb1-1
                    ELSE
                      CYCLE
                    END IF
                  END IF

                  neighborBin = ni + nb1*nj + nb1xnb2*nk
                  nRank = ppiclf_BinToRankMap(neighborBin)
                  IF(nRank .NE. iRank) THEN
                    ! Identify which specific faces this neighbor touches
                    ! Diagonals will correctly trigger multiple faces
                    IF(di .EQ. -1) 
     >                ppiclf_LRankBoundary(iBin,1) = .TRUE.
                    IF(di .EQ.  1) 
     >                ppiclf_LRankBoundary(iBin,2) = .TRUE.
                    IF(dj .EQ. -1) 
     >                ppiclf_LRankBoundary(iBin,3) = .TRUE.
                    IF(dj .EQ.  1) 
     >                ppiclf_LRankBoundary(iBin,4) = .TRUE.
                    IF(dk .EQ. -1) 
     >                ppiclf_LRankBoundary(iBin,5) = .TRUE.
                    IF(dk .EQ.  1) 
     >              ppiclf_LRankBoundary(iBin,6) = .TRUE.
                  END IF

                END DO !di
              END DO !dj
            END DO !dk
            iBin = iBin + 1
          END DO !ii
        END DO !jj
      END DO !kk

      RETURN
      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_MoveParticlePartLB
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
      INTEGER*4 i, icount, j0, ierr
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

        IF(PPICLF_LRP2 .GT. 0) THEN
          CALL ppiclf_copy(rtemp(icount,i),
     >                     ppiclf_rprop2(1,i),PPICLF_LRP2)
          icount = icount + PPICLF_LRP2
        END IF
        IF(PPICLF_LRP3 .GT. 0) THEN
          CALL ppiclf_copy(rtemp(icount,i),
     >                     ppiclf_rprop3(1,i),PPICLF_LRP3)
          icount = icount + PPICLF_LRP3
        END IF
        IF(PPICLF_LRP4 .GT. 0) THEN
          CALL ppiclf_copy(rtemp(icount,i),
     >                     ppiclf_rprop4(1,i),PPICLF_LRP4)
          icount = icount + PPICLF_LRP4
        END IF
        IF(PPICLF_LRP5 .GT. 0) THEN
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
!
      CALL pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >             ,ppiclf_npart,PPICLF_LPART ! Setup
     >             ,ppiclf_iprop,PPICLF_LIP   ! Integer Comm
     >             ,partl,0                   ! Logical Comm
     >             ,rtemp,rtempLim            ! Real Comm
     >             ,j0)                       ! Receiver processor index
!
#ifdef PERF
      PPICLF_TMPI_moveRP = PPICLF_TMPI_moveRP
     >     + (MPI_WTIME() - tstart)
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
        IF(PPICLF_LRP2 .GT. 0) THEN
        CALL ppiclf_copy(ppiclf_rprop2(1,i),rtemp(icount,i),
     >                   PPICLF_LRP2)
        icount = icount + PPICLF_LRP2
        END IF
        IF(PPICLF_LRP3 .GT. 0) THEN
          CALL ppiclf_copy(ppiclf_rprop3(1,i),rtemp(icount,i),
     >                     PPICLF_LRP3)
          icount = icount + PPICLF_LRP3
        END IF
        IF(PPICLF_LRP4 .GT. 0) THEN
          CALL ppiclf_copy(ppiclf_rprop4(1,i),rtemp(icount,i),
     >                     PPICLF_LRP4)
          icount = icount + PPICLF_LRP4
        END IF
        IF(PPICLF_LRP5 .GT. 0) THEN
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
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE 'mpif.h'
!
!      INTEGER*4 icalld
!      SAVE      icalld
!      DATA      icalld /0/
      INTEGER*4 nkey(2), i, j, k, l, ie, iee, ii, jj, kk, irank,
     >          nl, nii, njj, nrr, iic, jjc, kkc, ierr, nb1, nb2,
     >          nb3, nb1xnb2 
      INTEGER*4 ix, iy, iz, ixLow, ixHigh, iyLow,
     >          iyHigh, izLow, izHigh, ibin, jbin, kbin, nbin, bin,
     >          nRankMaps, RankMaps(27,5)
      REAL*8    rxval, ryval, rzval, MinPoint(3),
     >          centeri(3), Max_CellLen(3), inv_dx(3)
      LOGICAL   partl, ErrorFound, MapCell 
      LOGICAL, ALLOCATABLE :: rank_is_mapped(:)
#ifdef PERF
      REAL*8    tstart, tfinal
#endif
!
! Code Start:
!
      ALLOCATE(rank_is_mapped(0:ppiclf_np-1))
      rank_is_mapped = .FALSE.
      ! Number of fluid finite volume cells that map to particle bins
      ppiclf_nCells_FV2PICL = 0 
      
      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)
      nb1xnb2 = nb1*nb2
      
      inv_dx(1) = 1.0D0 / ppiclf_bins_dx(1)
      inv_dx(2) = 1.0D0 / ppiclf_bins_dx(2)
      inv_dx(3) = 1.0D0 / ppiclf_bins_dx(3)
      
      ! Loops through number of fluid FV cells on this processor
      DO ie=1,ppiclf_nFVCells  
        centeri(1) = ppiclf_fluid_grid(1,ie)
        centeri(2) = ppiclf_fluid_grid(2,ie)
        centeri(3) = ppiclf_fluid_grid(3,ie)
        ! Cycles if fluid cell center is outside of any bin boundaries
        ! Bin buffer ensures cells just outside bin aren't needed 
        IF (centeri(1) .GT. ppiclf_binb(2)) CYCLE
        IF (centeri(2) .GT. ppiclf_binb(4)) CYCLE
        IF (centeri(3) .GT. ppiclf_binb(6)) CYCLE
        IF (centeri(1) .LT. ppiclf_binb(1)) CYCLE
        IF (centeri(2) .LT. ppiclf_binb(3)) CYCLE
        IF (centeri(3) .LT. ppiclf_binb(5)) CYCLE
        ! Determines what bin the fluid cell is nominally mapped to
        ibin = FLOOR((centeri(1)-ppiclf_binb(1))*inv_dx(1))
        jbin = FLOOR((centeri(2)-ppiclf_binb(3))*inv_dx(2)) 
        kbin = FLOOR((centeri(3)-ppiclf_binb(5))*inv_dx(3))
        ! Below accounts for round-off errors
        ibin = MAX(0, MIN(ibin, nb1-1))
        jbin = MAX(0, MIN(jbin, nb2-1))
        kbin = MAX(0, MIN(kbin, nb3-1))
        ! Calculates processor rank
        bin  = ibin + nb1*jbin + nb1xnb2*kbin
        ! Only map if fluid cell is needed for interpolation 
        ! or projection (when particle nearby).
        IF(.NOT. ppiclf_LMapFluid(bin)) CYCLE 
        ! Will loop through cell mapping once if all below stay as 0
        ixLow  = 0
        ixHigh = 0
        iyLow  = 0
        iyHigh = 0
        izLow  = 0
        izHigh = 0
        ! Change loop bounds if the bin is on a MPI Boundary Face
        IF(ppiclf_LRankBoundary(bin,1)) ixLow  = -1
        IF(ppiclf_LRankBoundary(bin,3)) iyLow  = -1
        IF(ppiclf_LRankBoundary(bin,5)) izLow  = -1
        IF(ppiclf_LRankBoundary(bin,2)) ixHigh =  1
        IF(ppiclf_LRankBoundary(bin,4)) iyHigh =  1
        IF(ppiclf_LRankBoundary(bin,6)) izHigh =  1
        nRankMaps = 0
        DO ix=ixLow,ixHigh
          ii = ibin + ix 
          ! Out of Bounds and Periodicity Check
          IF(ii .LT. 0 .OR. ii .GT. nb1-1) THEN
            IF(ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)
     >                               .AND. nb1 .GT. 1) THEN
              IF(ii .GT. nb1-1) ii = 0
              IF(ii .LT. 0)     ii = nb1 - 1
            ELSE
              CYCLE
            END IF
          END IF
          
          DO iy=iyLow,iyHigh
            jj = jbin + iy
            ! Out of Bounds and Periodicity Check
            IF (jj .LT. 0 .OR. jj .GT. nb2-1) THEN
              IF(ppiclf_linperiodic(2) .AND. ppiclf_EqualDomain(2)
     >                                 .AND. nb2 .GT. 1) THEN
                IF(jj .GT. nb2-1) jj = 0
                IF(jj .LT. 0)     jj = nb2 - 1
              ELSE
                CYCLE
              END IF
            END IF
            
            DO iz=izLow,izHigh
              kk = kbin + iz          
              ! Out of Bounds and Periodicity Check   
              IF (kk .LT. 0 .OR. kk .GT. nb3-1) THEN
                IF(ppiclf_linperiodic(3) .AND. ppiclf_EqualDomain(3)
     >                                   .AND. nb3 .GT. 1) THEN
                  IF(kk .GT. nb3-1) kk = 0
                  IF(kk .LT. 0)     kk = nb3 - 1
                ELSE
                  CYCLE
                END IF
              END IF

              nbin  = ii + nb1*jj + nb1xnb2*kk
              irank = ppiclf_BinToRankMap(nbin)
              IF(ppiclf_LMapFluid(nbin)) THEN 
                ! Only need cell once per rank.
                ! Min distance takes care of periodicity
                ! in particle to cell mapping 
                IF(.NOT. rank_is_mapped(irank)) THEN
                  rank_is_mapped(irank) = .TRUE.
                  nRankMaps             = nRankMaps + 1
                  RankMaps(nRankMaps,1) = irank ! Receiving Rank
                  RankMaps(nRankMaps,2) = bin   ! Original bin index
                  RankMaps(nRankMaps,3) = ibin  ! Original index
                  RankMaps(nRankMaps,4) = jbin  ! Original index
                  RankMaps(nRankMaps,5) = kbin  ! Original index
                END IF !rank_is_mapped
              END IF !ppiclf_LMapFluid(nbin)
            END DO !iz
          END DO !iy
        END DO !ix
        IF(nRankMaps .GT. 0) THEN
          DO i = 1,nRankMaps
            ! Counter for overlap cells created by this processor
            ppiclf_nCells_FV2PICL = ppiclf_nCells_FV2PICL + 1
            IF(ppiclf_nCells_FV2PICL .GT. PPICLF_LEE) THEN
              PRINT*, '***ERROR*** Issue when creating overlap',
     >                ' cell mapping. Either increase ppiclf cell',
     >                ' limit, decrease total cell in fluid grid, or',
     >                ' increase number of particles.'
              PRINT*, '***ERROR*** Due to rocflu domain cell partition',
     >                ' and not the ppiclf domain overlap cells.'
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
            ppiclf_cell_map(7,ppiclf_nCells_FV2PICL) = RankMaps(i,5)
          END DO !nRankMaps
        END IF
        DO i = 1,nRankMaps
          rank_is_mapped(RankMaps(i,1)) = .FALSE.
        END DO
      END DO !ie
      
      DEALLOCATE(rank_is_mapped)
      ! Save cell information in new array for those that
      ! are overlap cells
      DO ie=1,ppiclf_nCells_FV2PICL 
        ! These copy all indicies since Fortran is column-major
        iee = ppiclf_cell_map(1,ie)
        CALL ppiclf_copy(ppiclf_picl_grid(1,ie)
     >                 ,ppiclf_fluid_grid(1,iee),7)
      END DO

      ! Must copy because array will be updated after send/receive
      ! Must have later for sending cell flow data and projection
      ppiclf_nCells_FV2PICL_Orig = ppiclf_nCells_FV2PICL
      DO ie=1,ppiclf_nCells_FV2PICL_Orig
         ! Copies cells to rank mapping (integer copy)
         CALL ppiclf_icopy(ppiclf_cell_map_Orig(1,ie)
     >            ,ppiclf_cell_map(1,ie),PPICLF_LRMAX)
      END DO

      ! GSLIB required info
      ! nCells_FV2PICL - number of columns to transfer
      ! PPICLF_LEE - number of columns declared
      ! nl - partl row size (dummy logical variable)
      nl   = 0
      ! nii - ppiclf_cell_map row size declared
      nii  = PPICLF_LRMAX
      ! njj - Row index of ppiclf_cell_map with receiver processor/rank
      njj  = 3
      ! nrr - ppiclf_picl_grid row size declared
      nrr  = 7
      ! Defines sorting order
      nkey(1) = 2
      nkey(2) = 1

#ifdef PERF
      tstart = MPI_WTIME()
      PPICLF_T_FVCells       = ppiclf_nFVCells
      PPICLF_T_OverlapCells_sent  = ppiclf_nCells_FV2PICL
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
      PPICLF_T_OverlapCells_received  = ppiclf_nCells_FV2PICL
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_moveOvlp = PPICLF_TMPI_moveOvlp + tfinal
      PPICLF_TMapOverlap = PPICLF_TMapOverlap - tfinal
#endif

      ! Find distance check for interpolation.
      ! This is 1.5*MaxCellLength to ensure that at least
      ! 27 neighboring cells are mapped.
      Max_CellLen(1) = 0.0D0
      Max_CellLen(2) = 0.0D0
      Max_CellLen(3) = 0.0D0
      ! Loop through overlapcells mapped to bin
      DO ie = 1,ppiclf_nCells_FV2PICL 
        DO l = 1,3
          ! Find max cell lengths in all dimensions
          IF(ppiclf_picl_grid(3+l,ie) .GT. Max_CellLen(l))
     >      Max_CellLen(l) = ppiclf_picl_grid(3+l,ie)
        END DO !l
      END DO !ie
      ! Update ppiclf filter, which is important as bin 
      ! boundaries grow and captures larger cells
      DO l = 1,3
        ppiclf_filter(l) = Max_CellLen(l)*1.51D0
      END DO
      ! Find max filter across processors (all overlap cells considered)
#ifdef PERF
      tstart = MPI_WTIME()
#endif
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_filter
     >                   ,3 ,MPI_DOUBLE_PRECISION
     >                   ,MPI_MAX ,ppiclf_comm, ierr)
#ifdef PERF
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_allreduces = PPICLF_TMPI_allreduces + tfinal
      PPICLF_TMapOverlap = PPICLF_TMapOverlap - tfinal
#endif

      DO l = 1,3
        ! Multiply by 1.5 so particle near face will
        ! find center one cell over in farthest direction
        ! Only need cells on this processor - not all overlap cells
        ppiclf_interp_dchk(l) = Max_CellLen(l)*1.5D0
      END DO

      RETURN
      END SUBROUTINE 
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_subbinRealParticleMap
      
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"


      INTEGER*4  i, tempSBin, i_SBin(3) 
     
      ! Quadruple it for ghost particles in the future 
      ppiclf_maxParticlePerBin = MAX(ppiclf_maxParticlePerBin,
     >                               4*MAXVAL(ppiclf_ParticleCount))

      CALL ppiclf_allocate_BTP(ppiclf_total_SBin,
     >                         ppiclf_maxParticlePerBin)
       
      IF(ppiclf_npart .LT. 1) RETURN

      ! This creates a list of local particles contained within
      ! each bin that resides on this rank.
      ppiclf_binPartCount = 0
      DO i = 1,ppiclf_npart
        i_SBin(1) = ppiclf_iprop(5,i) - ppiclf_binOffset(1)
        i_SBin(2) = ppiclf_iprop(6,i) - ppiclf_binOffset(2)
        i_SBin(3) = ppiclf_iprop(7,i) - ppiclf_binOffset(3)
        tempSBin = i_SBin(1) + ppiclf_nSBin(1)*i_SBin(2) + 
     >             ppiclf_nSBin(1)*ppiclf_nSBin(2)*i_SBin(3)
        ppiclf_binPartCount(tempSBin) = 
     >                      ppiclf_binPartCount(tempSBin) + 1
        ppiclf_binPartList(tempSBin, 
     >                     ppiclf_binPartCount(tempSBin)) = i
      END DO ! i

      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------

      SUBROUTINE ppiclf_comm_CreateGhostPartLB

      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      REAL*8     PeriodicShift(3), distchk, xlo(3), xhi(3)
      REAL*8     gp_x, gp_y, gp_z
      INTEGER*4  ip, iip, jjp, kkp, iig, jjg, kkg, nrank, l, nbin
      INTEGER*4  nb1, nb2, nb3, nb1xnb2, ix, iy, iz, iBin
      INTEGER*4  sb_x, sb_y, sb_z, tempSBin, pIdx, p_id
      INTEGER*4  ierr, i, j
      INTEGER*4  sx, sy, sz, shiftKey
      LOGICAL    wrap_x, wrap_y, wrap_z, wrapped_x, wrapped_y, wrapped_z
      LOGICAL    bnd_x_neg, bnd_x_pos, bnd_y_neg, bnd_y_pos 
      LOGICAL    bnd_z_neg, bnd_z_pos, ghost_x_neg, ghost_x_pos
      LOGICAL    ghost_y_neg, ghost_y_pos, ghost_z_neg, ghost_z_pos
      INTEGER*4, ALLOCATABLE, SAVE :: sent_stamp(:,:)
      INTEGER*4, SAVE :: gpStamp = 0
  
      IF(ppiclf_npart .LT. 1) RETURN
      ! sent_stamp columns mark the SIGNED periodic shift applied to the
      ! ghost: shiftKey = (sx+1) + 3*(sy+1) + 9*(sz+1), sx/sy/sz in
      ! {-1,0,+1}. 13 == no wrap. Two ghosts to one rank are identical
      ! only if they share this shift vector, so this is the correct
      ! dedup key.
      IF(.NOT. ALLOCATED(sent_stamp)) THEN
        ALLOCATE(sent_stamp(0:ppiclf_np-1,0:26))
        sent_stamp = 0
      END IF 
      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)
      nb1xnb2 = nb1 * nb2
      ! This adds a buffer. Will make extra GPs, but won't miss any
      distchk = ppiclf_nndist * 1.03D0 
      ppiclf_npart_gp = 0

      wrap_x = ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)
      wrap_y = ppiclf_linperiodic(2) .AND. ppiclf_EqualDomain(2)
      wrap_z = ppiclf_linperiodic(3) .AND. ppiclf_EqualDomain(3)

      DO l = 1,3
        IF(ppiclf_linperiodic(l) .AND. ppiclf_EqualDomain(l)) THEN
          PeriodicShift(l) = ppiclf_binb(2*l) - ppiclf_binb(2*l-1)
        ELSE
          PeriodicShift(l) = 0.0D0
        END IF
      END DO

      DO sb_z = 0, ppiclf_nSBin(3)-1
        DO sb_y = 0, ppiclf_nSBin(2)-1
          DO sb_x = 0, ppiclf_nSBin(1)-1
            tempSBin = sb_x + ppiclf_nSBin(1)*sb_y +
     >                 ppiclf_nSBin(1)*ppiclf_nSBin(2)*sb_z
            iip  = sb_x + ppiclf_binOffset(1)
            jjp  = sb_y + ppiclf_binOffset(2)
            kkp  = sb_z + ppiclf_binOffset(3)
            iBin = iip + jjp*nb1 + nb1xnb2*kkp

            IF(iip .LT. 0 .OR. iip .GE. nb1) CYCLE
            IF(jjp .LT. 0 .OR. jjp .GE. nb2) CYCLE
            IF(kkp .LT. 0 .OR. kkp .GE. nb3) CYCLE

            ! Instant bailout for empty bins
            IF(ppiclf_ParticleCount(iBin) .EQ. 0) CYCLE
 
            bnd_x_neg = ppiclf_LRankBoundary(iBin, 1)
     >                  .OR. (wrap_x .AND. iip .EQ. 0)
            bnd_x_pos = ppiclf_LRankBoundary(iBin, 2)
     >                  .OR. (wrap_x .AND. iip .EQ. nb1-1)
            bnd_y_neg = ppiclf_LRankBoundary(iBin, 3)
     >                  .OR. (wrap_y .AND. jjp .EQ. 0)
            bnd_y_pos = ppiclf_LRankBoundary(iBin, 4)
     >                  .OR. (wrap_y .AND. jjp .EQ. nb2-1)
            bnd_z_neg = ppiclf_LRankBoundary(iBin, 5)
     >                  .OR. (wrap_z .AND. kkp .EQ. 0)
            bnd_z_pos = ppiclf_LRankBoundary(iBin, 6)
     >                  .OR. (wrap_z .AND. kkp .EQ. nb3-1)

            ! If sub-bin is interior or not on periodic boundary,
            ! skip ALL particles inside it instantly!
            IF (.NOT. (bnd_x_neg .OR. bnd_x_pos .OR. bnd_y_neg .OR. 
     >                 bnd_y_pos .OR. bnd_z_neg .OR. bnd_z_pos)) CYCLE

            ! Compute geometric walls for this sub-bin exactly ONCE
            xlo(1) = ppiclf_binb(1) + DBLE(iip)*ppiclf_bins_dx(1)
            xhi(1) = xlo(1) + ppiclf_bins_dx(1)
            xlo(2) = ppiclf_binb(3) + DBLE(jjp)*ppiclf_bins_dx(2)
            xhi(2) = xlo(2) + ppiclf_bins_dx(2)
            xlo(3) = ppiclf_binb(5) + DBLE(kkp)*ppiclf_bins_dx(3)
            xhi(3) = xlo(3) + ppiclf_bins_dx(3)

            ! Process ONLY the particles in this specific
            ! subbin on MPI or periodic boundary.
            DO pIdx = 1, ppiclf_binPartCount(tempSBin)
              ip = ppiclf_binPartList(tempSBin, pIdx)
              gpStamp = gpStamp + 1
              ghost_x_neg = bnd_x_neg .AND.
     >                      (ppiclf_y(1,ip) - xlo(1) .LT. distchk)
              ghost_x_pos = bnd_x_pos .AND. 
     >                      (xhi(1) - ppiclf_y(1,ip) .LT. distchk)
              ghost_y_neg = bnd_y_neg .AND. 
     >                      (ppiclf_y(2,ip) - xlo(2) .LT. distchk)
              ghost_y_pos = bnd_y_pos .AND. 
     >                      (xhi(2) - ppiclf_y(2,ip) .LT. distchk)
              ghost_z_neg = bnd_z_neg .AND. 
     >                      (ppiclf_y(3,ip) - xlo(3) .LT. distchk)
              ghost_z_pos = bnd_z_pos .AND. 
     >                      (xhi(3) - ppiclf_y(3,ip) .LT. distchk)

              ! Second Bailout: Particle is in a boundary bin,
              ! but physically too far from the wall
              IF (.NOT. (ghost_x_neg .OR. ghost_x_pos .OR. ghost_y_neg
     >            .OR. ghost_y_pos .OR. ghost_z_neg .OR. ghost_z_pos) )
     >            CYCLE

              ! Generate the required ghost particles
              DO iz = -1, 1
                ! See if particle is on a z MPI boundary
                IF (iz .EQ. -1 .AND. .NOT. ghost_z_neg) CYCLE
                IF (iz .EQ.  1 .AND. .NOT. ghost_z_pos) CYCLE
                
                DO iy = -1, 1
                  ! See if particle is on a y MPI boundary
                  IF (iy .EQ. -1 .AND. .NOT. ghost_y_neg) CYCLE
                  IF (iy .EQ.  1 .AND. .NOT. ghost_y_pos) CYCLE
                  
                  DO ix = -1, 1
                    ! See if particle is on a x MPI boundary 
                    IF (ix .EQ. -1 .AND. .NOT. ghost_x_neg) CYCLE
                    IF (ix .EQ.  1 .AND. .NOT. ghost_x_pos) CYCLE
                    ! This is the bin the real particle is in
                    IF (ix .EQ. 0 .AND. iy .EQ. 0 .AND. iz .EQ. 0) CYCLE

                    iig = iip + ix
                    jjg = jjp + iy
                    kkg = kkp + iz
                    
                    wrapped_x = .FALSE. 
                    wrapped_y = .FALSE. 
                    wrapped_z = .FALSE. 

                    gp_x = ppiclf_y(1,ip)
                    gp_y = ppiclf_y(2,ip)
                    gp_z = ppiclf_y(3,ip)

                    ! --- Periodic Wrapping Logic ---
                    IF(iig .LT. 0 .OR. iig .GE. nb1) THEN
                      IF(wrap_x) THEN
                        wrapped_x = .TRUE.
                        IF(iig .LT. 0) THEN
                          iig = nb1 - 1
                          gp_x = gp_x + PeriodicShift(1)
                        ELSE
                          iig = 0
                          gp_x = gp_x - PeriodicShift(1)
                        END IF
                      ELSE
                        CYCLE
                      END IF
                    END IF

                    IF(jjg .LT. 0 .OR. jjg .GE. nb2) THEN
                      IF(wrap_y) THEN
                        wrapped_y = .TRUE.
                        IF(jjg .LT. 0) THEN
                          jjg = nb2 - 1
                          gp_y = gp_y + PeriodicShift(2)
                        ELSE
                          jjg = 0
                          gp_y = gp_y - PeriodicShift(2)
                        END IF
                      ELSE
                        ! No GP needed
                        CYCLE
                      END IF
                    END IF

                    IF(kkg .LT. 0 .OR. kkg .GE. nb3) THEN
                      IF(wrap_z) THEN
                        wrapped_z = .TRUE.
                        IF(kkg .LT. 0) THEN
                          kkg = nb3 - 1
                          gp_z = gp_z + PeriodicShift(3)
                         ELSE
                          kkg = 0
                          gp_z = gp_z - PeriodicShift(3)
                         END IF
                      ELSE
                        ! No GP needed
                        CYCLE
                      END IF
                    END IF

                    ! This determines what rank(s) need the GP
                    nbin  = iig + nb1*jjg + nb1xnb2*kkg
                    nrank = ppiclf_BinToRankMap(nbin)

                    IF(.NOT. wrapped_x) THEN
                      iig = iip
                    END IF
                    IF(.NOT. wrapped_y) THEN
                      jjg = jjp
                    END IF
                    IF(.NOT. wrapped_z) THEN
                      kkg = kkp
                    END IF
                    ! This is the bin of the "real" particle,
                    ! which may be shifted when periodic.
                    nbin  = iig + nb1*jjg + nb1xnb2*kkg

                    sx = 0
                    sy = 0
                    sz = 0
                    IF(wrapped_x) sx = ix
                    IF(wrapped_y) sy = iy
                    IF(wrapped_z) sz = iz
                    shiftKey = (sx+1) + 3*(sy+1) + 9*(sz+1) ! 0..26

                    ! Skip if it belongs to this rank AND it didn't wrap
                    ! periodically (shiftKey==13 is the no-wrap center).
                    IF(nrank .EQ. ppiclf_nid .AND. shiftKey .EQ. 13) 
     >                CYCLE

                    IF(sent_stamp(nrank,shiftKey) .EQ. gpStamp) CYCLE

                    sent_stamp(nrank,shiftKey) = gpStamp
                    ppiclf_npart_gp = ppiclf_npart_gp + 1
                    ppiclf_iprop_gp(1, ppiclf_npart_gp) =
     >                                      ppiclf_iprop(1, ip)  
                    ppiclf_iprop_gp(2, ppiclf_npart_gp) =
     >                                      ppiclf_iprop(2, ip)  
                    ppiclf_iprop_gp(3, ppiclf_npart_gp) =
     >                                      ppiclf_iprop(3, ip)  
                    ppiclf_iprop_gp(4, ppiclf_npart_gp) = nrank
                    ppiclf_iprop_gp(5, ppiclf_npart_gp) = iig
                    ppiclf_iprop_gp(6, ppiclf_npart_gp) = jjg
                    ppiclf_iprop_gp(7, ppiclf_npart_gp) = kkg
                    ppiclf_iprop_gp(8, ppiclf_npart_gp) = nbin


                    ppiclf_rprop_gp(1, ppiclf_npart_gp) = gp_x
                    ppiclf_rprop_gp(2, ppiclf_npart_gp) = gp_y
                    ppiclf_rprop_gp(3, ppiclf_npart_gp) = gp_z
                    ppiclf_rprop_gp(4:PPICLF_LRS, ppiclf_npart_gp) =
     >                                         ppiclf_y(4:PPICLF_LRS,ip)
                    ppiclf_rprop_gp(1+PPICLF_LRS:PPICLF_LRP+PPICLF_LRS, 
     >                              ppiclf_npart_gp   ) = 
     >                                ppiclf_rprop(1:PPICLF_LRP, ip)
                  END DO !iz
                END DO !iy
              END DO !ix

            END DO ! pIdx
          END DO ! sb_x
        END DO ! sb_y
      END DO ! sb_z

      RETURN
      END SUBROUTINE
!      
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
      INTEGER*4 iprop_proc_index, ierr
      LOGICAL   partl  ! Dummy variable       
#ifdef PERF
      REAL*8    tstart, tfinal
#endif
!
      iprop_proc_index = 4 ! since ppiclf_iprop(4,np) contains processor
                           ! that should receive ghost particle
#ifdef PERF
      tstart = MPI_WTIME()
      PPICLF_T_GhostPartSent = ppiclf_npart_gp
#endif

      CALL pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >             ,ppiclf_npart_gp,PPICLF_LPART_GP ! Setup
     >             ,ppiclf_iprop_gp,PPICLF_LIP_GP   ! Integer Comm
     >             ,partl,0                         ! Logical Comm
     >             ,ppiclf_rprop_gp,PPICLF_LRP_GP   ! Real Comm
     >             ,iprop_proc_index)               ! Receiver processor index

#ifdef PERF
      PPICLF_T_GhostPartRec  = ppiclf_npart_gp
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_moveGP = PPICLF_TMPI_moveGP + tfinal
      PPICLF_TMoveGhost = PPICLF_TMoveGhost - tfinal
#endif
      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_subbinGhostParticleMap

      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"


      INTEGER*4  i, tempSBin, i_SBin(3), newMax

      IF(ppiclf_npart .LT. 1) RETURN
      IF(ppiclf_npart_gp .LT. 1) RETURN

      ! Single pass: Append ghosts
      DO i = 1, ppiclf_npart_gp
        i_SBin(1) = ppiclf_iprop_gp(5,i) - ppiclf_binOffset(1)
        i_SBin(2) = ppiclf_iprop_gp(6,i) - ppiclf_binOffset(2)
        i_SBin(3) = ppiclf_iprop_gp(7,i) - ppiclf_binOffset(3)

        ! These can be out-of-bounds since ppiclf_binOffsets on this
        ! rank is different from the values on the rank that sent the 
        ! ghost particle
        IF(i_SBin(1) .LT. -1 .OR. i_SBin(1) .GT. ppiclf_nSBin(1)) CYCLE
        IF(i_SBin(2) .LT. -1 .OR. i_SBin(2) .GT. ppiclf_nSBin(2)) CYCLE
        IF(i_SBin(3) .LT. -1 .OR. i_SBin(3) .GT. ppiclf_nSBin(3)) CYCLE

        i_SBin(1) = MAX(0, MIN(i_SBin(1), ppiclf_nSBin(1)-1))
        i_SBin(2) = MAX(0, MIN(i_SBin(2), ppiclf_nSBin(2)-1))
        i_SBin(3) = MAX(0, MIN(i_SBin(3), ppiclf_nSBin(3)-1))

        tempSBin = i_SBin(1) + ppiclf_nSBin(1)*i_SBin(2) + 
     >             ppiclf_nSBin(1)*ppiclf_nSBin(2)*i_SBin(3)

        IF(tempSBin .GT. SIZE(ppiclf_binPartCount)-1) THEN
          CYCLE
        END IF

        ppiclf_binPartCount(tempSBin) = 
     >                    ppiclf_binPartCount(tempSBin) + 1
        IF(ppiclf_binPartCount(tempSBin) 
     >     .GT. ppiclf_maxParticlePerBin) THEN
          ! Scale aggressively to prevent a cascade of slow reAllocs
          newMax = INT(DBLE(ppiclf_maxParticlePerBin) * 2.0D0) + 1
          CALL ppiclf_reallocate_BTP(ppiclf_total_SBin,
     >                               ppiclf_maxParticlePerBin,
     >                               newMax)
          ppiclf_maxParticlePerBin = newMax
        END IF

        ! Negative index indicates ghost particle
        ppiclf_binPartList(tempSBin, 
     >                     ppiclf_binPartCount(tempSBin)) = -i
      END DO 

      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_subbinFineParticleMap
!
!     Build a FINE sub-bin grid sized by the particle-particle search
!     cutoff (ppiclf_nndist) and map BOTH real and ghost particles into
!     it by POSITION. This is a separate, finer grid than the coarse
!     iprop sub-bins (sized by ppiclf_filter) used for communication and
!     ghost creation. Using it for the nearest-neighbor search shrinks
!     the 27-cell candidate set by roughly (ppiclf_filter/nndist)**3.
!
!     Storage convention (matches ppiclf_solve_NearestNeighborSB):
!       real  particle i -> stored as +i (position ppiclf_y(1:3,i))
!       ghost particle i -> stored as -i (position ppiclf_rprop_gp(1:3,i))
!
!     The fine/coarse DECISION (ppiclf_useFineGrid) is made by the
!     CALLER (the driver, e.g. ppiclf_solve_SetYdot) BEFORE this routine
!     is called; this routine only BUILDS the grid. The search then
!     locates a particle's own fine cell FROM ITS POSITION using
!       ppiclf_fineLo / ppiclf_fineInvLen / ppiclf_nFine
!     and sweeps the 27 neighbor cells of
!       ppiclf_finePartCount / ppiclf_finePartList.
!
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
!
! Internal:
!
      INTEGER*4  i, d, fi(3), tempBin, ng, icount, newMax
      REAL*8     xlo(3), xhi(3), boxLen(3), fineCut, fineLen(3), pad

      IF(ppiclf_npart .LT. 1) RETURN

      ng = ppiclf_npart_gp
      IF(ng .LT. 0) ng = 0

      ! Search cutoff = radius used in ppiclf_solve_NearestNeighborSB
      ! (driver guarantees ppiclf_nndist > 0 before calling).
      fineCut = ppiclf_nndist

      ! Built from actual positions (not the coarse halo) so every
      ! particle is guaranteed inside and no clamp can ever drop one.
      xlo(1) = ppiclf_y(1,1)
      xlo(2) = ppiclf_y(2,1)
      xlo(3) = ppiclf_y(3,1)
      xhi(1) = xlo(1)
      xhi(2) = xlo(2)
      xhi(3) = xlo(3)
      DO i = 1,ppiclf_npart
        DO d = 1,3
          xlo(d) = MIN(xlo(d), ppiclf_y(d,i))
          xhi(d) = MAX(xhi(d), ppiclf_y(d,i))
        END DO
      END DO
      DO i = 1,ng
        DO d = 1,3
          xlo(d) = MIN(xlo(d), ppiclf_rprop_gp(d,i))
          xhi(d) = MAX(xhi(d), ppiclf_rprop_gp(d,i))
        END DO
      END DO

      ! Build the fine grid. Cell length must be >= cutoff so the
      ! 27-cell stencil is complete (standard cell-list invariant).
      ! FLOOR(boxLen/cut) maximizes cells while keeping len >= cut.
      pad = 0.5D0*fineCut
      ppiclf_total_fineBin = 1
      DO d = 1,3
        xlo(d) = xlo(d) - pad
        xhi(d) = xhi(d) + pad
        boxLen(d) = xhi(d) - xlo(d)
        ppiclf_nFine(d) = MAX(1, FLOOR(boxLen(d)/fineCut))
        fineLen(d) = boxLen(d)/DBLE(ppiclf_nFine(d))
        ppiclf_fineLo(d) = xlo(d)
        ppiclf_fineInvLen(d) = 1.0D0/fineLen(d)
        ppiclf_total_fineBin = ppiclf_total_fineBin*ppiclf_nFine(d)
      END DO

      ! Per-fine-bin capacity. Keep the cheap data-driven estimate but
      ! never shrink it across calls; grow-on-demand below covers any
      ! residual overflow, so correctness no longer depends on it.
      ppiclf_maxPartPerFineBin =
     >     MAX(ppiclf_maxPartPerFineBin, 5*MAXVAL(ppiclf_ParticleCount))
      IF(ppiclf_maxPartPerFineBin .LT. 1) ppiclf_maxPartPerFineBin = 20
      CALL ppiclf_allocate_FineBTP(ppiclf_total_fineBin,
     >                            ppiclf_maxPartPerFineBin)
      ppiclf_finePartCount = 0

      DO i = 1,ppiclf_npart
        fi(1)=FLOOR((ppiclf_y(1,i)-ppiclf_fineLo(1))
     >              *ppiclf_fineInvLen(1))
        fi(2)=FLOOR((ppiclf_y(2,i)-ppiclf_fineLo(2))
     >              *ppiclf_fineInvLen(2))
        fi(3)=FLOOR((ppiclf_y(3,i)-ppiclf_fineLo(3))
     >              *ppiclf_fineInvLen(3))
        fi(1)=MAX(0,MIN(fi(1),ppiclf_nFine(1)-1))
        fi(2)=MAX(0,MIN(fi(2),ppiclf_nFine(2)-1))
        fi(3)=MAX(0,MIN(fi(3),ppiclf_nFine(3)-1))
        tempBin = fi(1) + ppiclf_nFine(1)*fi(2)
     >          + ppiclf_nFine(1)*ppiclf_nFine(2)*fi(3)
        icount = ppiclf_finePartCount(tempBin) + 1
        ppiclf_finePartCount(tempBin) = icount
        ! Grow on demand. icount rises by 1 per hit, so one doubling
        ! always restores capacity (newMax >= icount).
        IF(icount .GT. ppiclf_maxPartPerFineBin) THEN
          newMax = 2*ppiclf_maxPartPerFineBin + 1
          CALL ppiclf_reallocate_FineBTP(ppiclf_total_fineBin,
     >                                   ppiclf_maxPartPerFineBin,
     >                                   newMax)
          ppiclf_maxPartPerFineBin = newMax
        END IF
        ppiclf_finePartList(tempBin,icount) = i
      END DO

      ! Ghost particles -> -i (shifted position ppiclf_rprop_gp)
      DO i = 1,ng
        fi(1)=FLOOR((ppiclf_rprop_gp(1,i)-ppiclf_fineLo(1))
     >              *ppiclf_fineInvLen(1))
        fi(2)=FLOOR((ppiclf_rprop_gp(2,i)-ppiclf_fineLo(2))
     >              *ppiclf_fineInvLen(2))
        fi(3)=FLOOR((ppiclf_rprop_gp(3,i)-ppiclf_fineLo(3))
     >              *ppiclf_fineInvLen(3))
        fi(1)=MAX(0,MIN(fi(1),ppiclf_nFine(1)-1))
        fi(2)=MAX(0,MIN(fi(2),ppiclf_nFine(2)-1))
        fi(3)=MAX(0,MIN(fi(3),ppiclf_nFine(3)-1))
        tempBin = fi(1) + ppiclf_nFine(1)*fi(2)
     >          + ppiclf_nFine(1)*ppiclf_nFine(2)*fi(3)
        icount = ppiclf_finePartCount(tempBin) + 1
        ppiclf_finePartCount(tempBin) = icount
        ! Grow on demand -- no ghost is ever dropped.
        IF(icount .GT. ppiclf_maxPartPerFineBin) THEN
          newMax = 2*ppiclf_maxPartPerFineBin + 1
          CALL ppiclf_reallocate_FineBTP(ppiclf_total_fineBin,
     >                                   ppiclf_maxPartPerFineBin,
     >                                   newMax)
          ppiclf_maxPartPerFineBin = newMax
        END IF
        ppiclf_finePartList(tempBin,icount) = -i
      END DO

      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------
!
      SUBROUTINE ppiclf_comm_subbinCellMap

      USE ppiclf_DynamicAllocation

      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
!
! Input:
!
      INTEGER*4  ie, i, j, k, iTemp_SBin(3)
     >          ,tempSBin, i_SBin(3), icount, newMax
      
      IF(ppiclf_npart .LT. 1) RETURN

      IF(ppiclf_nCells_FV2PICL .LE. 0 .AND. ppiclf_npart .GT. 0) THEN
        PRINT*,'ERROR: ',ppiclf_npart, 'Particles mapped to bin:'
     >         ,ppiclf_nid
        PRINT*,'No cells mapped to bin for Interpolation/Projection.'
        CALL ppiclf_exittr('Failure in particle to cell mapping$',
     >                      0.D0,0)
      END IF

      ! 92 was chosen with the assumption that min_cell = 3*max_cell
      ! in all 3 dimension (3*3*3*1.5*1.5*1.5=91)
      IF(ppiclf_maxCellsPerBin .LT. 1) ppiclf_maxCellsPerBin = 92
      CALL ppiclf_allocate_BTC(ppiclf_total_SBin, ppiclf_maxCellsPerBin)
      ppiclf_binCellCount = 0

      ! NOTE: Overlap cells have not been shifted for periodicity
      !       Both the cell centroid and cell bin indices are of the
      !       original fluid cell.
      ! Creates a list of overlap cells for a given bin on this rank.
      ! In the i,j,k loops below, 0 takes care of non-periodic mapping
      ! and 1 takes care of periodic mapping.  If a cell is in corner,
      ! it can be mapped to a maximum of 2*2*2=8 bins.
      DO ie = 1,ppiclf_nCells_FV2PICL
        i_SBin(1) = ppiclf_cell_map(5,ie) - ppiclf_binOffset(1)
        i_SBin(2) = ppiclf_cell_map(6,ie) - ppiclf_binOffset(2)
        i_SBin(3) = ppiclf_cell_map(7,ie) - ppiclf_binOffset(3)
        DO i = 0,1
          IF(i .EQ. 0) THEN
            iTemp_SBin(1) = i_SBin(1)
            IF(iTemp_SBin(1) .LT. 0 .OR. 
     >         iTemp_SBin(1) .GT. ppiclf_nSBin(1) - 1) CYCLE 
          ELSE ! i .EQ. 1
            IF(ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)) THEN
              IF(i_SBin(1) .LE. 0) THEN
                iTemp_SBin(1) = ppiclf_nSBin(1) - 1
                IF(iTemp_SBin(1) .EQ. i_SBin(1)) CYCLE 
              ELSE IF(i_SBin(1) .GE. ppiclf_nSBin(1) - 1) THEN
                iTemp_SBin(1) = 0
                IF(iTemp_SBin(1) .EQ. i_SBin(1)) CYCLE
              ELSE
                CYCLE
              END IF
            ELSE 
              CYCLE
            END IF
          END IF
          DO j = 0,1
            IF(j .EQ. 0) THEN
              iTemp_SBin(2) = i_SBin(2)
              IF(iTemp_SBin(2) .LT. 0 .OR.  
     >           iTemp_SBin(2) .GT. ppiclf_nSBin(2) - 1) CYCLE
            ELSE ! j .EQ. 1
              ! This takes care of periodicity for single processor
              IF(ppiclf_linperiodic(2).AND.ppiclf_EqualDomain(2)) THEN
                IF(i_SBin(2) .LE. 0) THEN
                  iTemp_SBin(2) = ppiclf_nSBin(2) - 1
                  IF(iTemp_SBin(2) .EQ. i_SBin(2)) CYCLE
                ELSE IF(i_SBin(2) .GE. ppiclf_nSBin(2) - 1) THEN
                  iTemp_SBin(2) = 0
                  IF(iTemp_SBin(2) .EQ. i_SBin(2)) CYCLE
                ELSE
                  CYCLE
                END IF
              ELSE 
                CYCLE
              END IF
            END IF
            DO k = 0,1
              IF(k .EQ. 0) THEN
                iTemp_SBin(3) = i_SBin(3)
                IF(iTemp_SBin(3) .LT. 0 .OR. 
     >             iTemp_SBin(3) .GT. ppiclf_nSBin(3) - 1) CYCLE
              ELSE ! k .EQ. 1
                ! This takes care of periodicity for single processor
                IF(ppiclf_linperiodic(3) .AND. 
     >                                    ppiclf_EqualDomain(3)) THEN 
                  IF(i_SBin(3) .LE. 0) THEN
                    iTemp_SBin(3) = ppiclf_nSBin(3) - 1
                    IF(iTemp_SBin(3) .EQ. i_SBin(3)) CYCLE
                  ELSE IF(i_SBin(3) .GE. ppiclf_nSBin(3) - 1) THEN
                    iTemp_SBin(3) = 0
                    IF(iTemp_SBin(3) .EQ. i_SBin(3)) CYCLE
                  ELSE
                    CYCLE
                  END IF
                ELSE 
                  CYCLE
                END IF
              END IF
              ! Finally, add the cell to a subbin 
              tempSBin = iTemp_SBin(1) + ppiclf_nSBin(1)*iTemp_SBin(2)
     >                 + ppiclf_nSBin(1)*ppiclf_nSBin(2)*iTemp_SBin(3)
              IF(tempSBin .LT. 0 .OR. 
     >           tempSBin .GT. ppiclf_total_SBin-1) THEN
               PRINT*, 'ERROR:Bad Subbin Index in Overlap Cell Mapping'
     >                 , tempSBin
                CALL ppiclf_exittr('',0.0D0,0)
              END IF

              icount = ppiclf_binCellCount(tempSBin) + 1
              ppiclf_binCellCount(tempSBin) = icount

              ! Grow on demand. icount rises by 1 per hit, so a single
              ! doubling always restores capacity (newMax >= icount).
              IF(icount .GT. ppiclf_maxCellsPerBin) THEN
                newMax = 2*ppiclf_maxCellsPerBin + 1
                CALL ppiclf_reallocate_BTC(ppiclf_total_SBin,
     >                                     ppiclf_maxCellsPerBin,
     >                                     newMax)
                ppiclf_maxCellsPerBin = newMax
              END IF

              ppiclf_binCellList(tempSBin, icount) = ie
            END DO !k
          END DO !j 
        END DO !i
      END DO !ie

      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------
!
! This is for statistics gathering only
!
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_LBCheck(series,iteration)

      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      INTEGER*4   ierr, i, rank, Pmin, Pmax, iteration, series, gl_part
      INTEGER*4, ALLOCATABLE :: LB_Count(:)

      REAL*8      Pavg
 
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

      Pmax = 0
      gl_part = 0
      Pmin = 99000000
      Pavg = 0.0D0
      IF(ppiclf_nid .EQ. 0) THEN
        DO i = 0,ppiclf_np-1
          Pmax = MAX(Pmax,LB_Count(i))
          Pmin = MIN(Pmin,LB_Count(i))
          gl_part = gl_part + LB_Count(i)
        END DO
        IF(ppiclf_np .EQ. 1) THEN
          WRITE(series+iteration,*) 'nBins:',ppiclf_n_bins(1),
     >                              ppiclf_n_bins(2),ppiclf_n_bins(3)
          WRITE(series+iteration,*) 'nParticles:', gl_part 
      
          WRITE(series+iteration,*) '# of Processors | Max Overshoot',
     >                              ' | Max Overshoot Percent'
        END IF
        Pavg = gl_part/DBLE(ppiclf_np)
        WRITE(series+iteration,*) ppiclf_np, Pmax - Pavg, 
     >                            (Pmax - Pavg)/Pavg*100.0 
      END IF
    
      DEALLOCATE(LB_Count)
     
      RETURN

      END SUBROUTINE
!----------------------------------------------------------------------
