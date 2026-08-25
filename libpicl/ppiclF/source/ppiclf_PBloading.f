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
      REAL*8    tstart, tfinal

      DO i = 1,3
        BinMinLen(i) = MAX(ppiclf_filter(i), ppiclf_nndist)
        BinBuffer(i) = 2.0D0*MAX(ppiclf_filter(i), ppiclf_nndist)
      END DO
      IF(ppiclf_istage .NE. 2 .AND. ppiclf_istage .NE. 3) THEN
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
        
        ! Flip sign on min so that MPI_ALLREDUCE called once
        DO i = 1,3
          local_extremes(i) = - local_extremes(i)
        END DO
      tstart = MPI_WTIME()
        ! Finds global bin domain boundaries across MPI ranks
        CALL MPI_ALLREDUCE(local_extremes, global_extremes, 6
     >                       ,MPI_DOUBLE_PRECISION, MPI_MAX
     >                       ,ppiclf_comm, ierr)
#ifdef PERF
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_allreduces = PPICLF_TMPI_allreduces + tfinal
      PPICLF_TCreateBin = PPICLF_TCreateBin - tfinal
#endif
        ! Flip sign on min values back to positive
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
      ! Bin geometry changed: per-bin overlap cell counts used by the
      ! weighted load balance must be rebuilt (see FindCellPartLB).
      ppiclf_CellCountValid = .FALSE.

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
     >           ,aPP, bMAP, bPROJ, cCell, LB_curMax
      REAL*8     ppiclf_comm_BinWeight
      EXTERNAL   ppiclf_comm_BinWeight
      INTEGER*4  LB_pcount, capPart
      LOGICAL    docheck
      REAL*8     tstart, tfinal
     
      nb1  = ppiclf_n_bins(1)
      nb2  = ppiclf_n_bins(2)
      nb3  = ppiclf_n_bins(3)
      nb1xnb2 = nb1 * nb2

      ! ---- global-count / rebalance-criterion cadence gate ----
      ! The O(N_b) count zeroing, the N_b-word allreduce, and the
      ! O(N_b) weighted scan below are P-independent fixed per-stage
      ! costs. Run them every ppiclf_LB_checkfreq stages (default 3 =
      ! once per RK3 step); per-particle bin assignment and migration
      ! detection stay per-stage. docheck is identical on all ranks
      ! (module counter + globally consistent binchanged), so the
      ! gated collective remains collectively consistent. Staleness of
      ! the global counts on off-stages is benign: CreateGhostPartLB
      ! tests LOCAL per-bin counts (fresh every stage) and the EIB
      ! one-bin dilation covers any bin a particle can enter within a
      ! step. checkfreq = 1 restores the legacy behavior exactly.
      ppiclf_LB_stagectr = ppiclf_LB_stagectr + 1
      ppiclf_LB_countsfresh = .FALSE.
      docheck = ppiclf_binchanged .OR.
     >          ppiclf_LB_stagectr .LE. 1 .OR.
     >          MOD(ppiclf_LB_stagectr, ppiclf_LB_checkfreq) .EQ. 0

      IF(docheck) ppiclf_ParticleCount = 0
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
        IF(docheck)
     >    ppiclf_ParticleCount(bin) = ppiclf_ParticleCount(bin) + 1
      END DO
      ! Sum particles per bin across MPI Ranks (check stages only)
      IF(docheck) THEN
      tstart = MPI_WTIME()
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_ParticleCount
     >                   ,ppiclf_totalBins ,MPI_INTEGER4, MPI_SUM
     >                   ,ppiclf_comm, ierr)
#ifdef PERF
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_allreduces = PPICLF_TMPI_allreduces + tfinal
      PPICLF_TFindPart = PPICLF_TFindPart - tfinal
#endif

      ppiclf_glnpart = SUM(ppiclf_ParticleCount)
      ppiclf_LB_countsfresh = .TRUE.
      END IF

      IF(ppiclf_binchanged) THEN
        ppiclf_rebalance     = .TRUE.
        ppiclf_particleMoved = .TRUE.
        RETURN
      END IF

      ! Logical OR comparison across MPI Ranks
      tstart = MPI_WTIME()
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_particleMoved
     >                   ,1, MPI_LOGICAL, MPI_LOR
     >                   ,ppiclf_comm, ierr)
#ifdef PERF
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_allreduces = PPICLF_TMPI_allreduces + tfinal
      PPICLF_TFindPart = PPICLF_TFindPart - tfinal
#endif

      ! Off-stage: no fresh global counts, so skip the weighted scan
      ! and criterion; migration (above) has already been handled.
      IF(.NOT. docheck) THEN
        ppiclf_rebalance = .FALSE.
        RETURN
      END IF

      ! Check if new BTRM required using the same weighted cost model
      ! as the partitioner (see LBWeightCoef/BinWeight):
      !   W(bin) = Np + aPP*Np^2 + (bMAP*sMAP + bPROJ)*Np*Nc + cCell*Nc
      ! If RankWeight > LB_criteria*TargetWeight -> Reassign BTRM
      ! ppiclf_CellCount (centroid-in-bin overlap cell counts) is
      ! still valid here: the binchanged path returned above and the
      ! fluid grid partition is static, so cell->bin is unchanged.
      ! Analytic cost model coefficients (see LBWeightCoef)
      CALL ppiclf_comm_LBWeightCoef(aPP, bMAP, bPROJ, cCell)
      LB_target = 0.0D0
      DO i = 0,ppiclf_totalBins-1
        LB_target = LB_target
     >      + ppiclf_comm_BinWeight(i, aPP, bMAP, bPROJ, cCell)
      END DO
      LB_criteria = 1.2D0
      LB_target   = LB_target/DBLE(ppiclf_np)
      LB_criteria = LB_criteria*LB_target
      ppiclf_rebalance = .FALSE.

      stride(1) = 1
      stride(2) = nb1
      stride(3) = nb1xnb2

      stride_L = stride(ppiclf_dL)
      stride_M = stride(ppiclf_dM)
      stride_S = stride(ppiclf_dS)

      LB_iterationCount = 0.0D0 
      LB_curMax = 0.0D0
      bin = 0 
      ! Weight balance no longer bounds particle counts implicitly,
      ! so drift between rebalances could exceed the transfer/common
      ! buffer capacity without tripping the weight criterion.
      ! Trigger a rebalance when any rank's count reaches 95% of
      ! PPICLF_LPART (the partitioner cuts at 90%, leaving headroom).
      LB_pcount = 0
      capPart   = PPICLF_LPART - PPICLF_LPART/20

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
              LB_pcount         = 0
            END IF
            LB_pcount = LB_pcount + ppiclf_ParticleCount(bin)

            LB_iterationCount = LB_iterationCount
     >          + ppiclf_comm_BinWeight(bin, aPP, bMAP, bPROJ, cCell)
            LB_curMax = MAX(LB_curMax, LB_iterationCount)

            ! Count capacity is a hard safety limit: absolute trigger.
            IF(LB_pcount .GT. capPart) THEN
              ppiclf_rebalance = .TRUE.
              EXIT outer_loop
            END IF

          END DO
        END DO
      END DO outer_loop

      IF(.NOT. ppiclf_rebalance) THEN
        IF(ppiclf_lbJustBal) THEN
          ! First check after a rebalance: record the achieved max
          ! segment weight as the futility baseline instead of
          ! re-triggering. If the partition is granularity limited,
          ! this is the best the partitioner can do and the criterion
          ! must adapt to it rather than fire every step.
          ppiclf_lbBaseMax = LB_curMax
          ppiclf_lbJustBal = .FALSE.
        ELSE IF(LB_curMax .GT. LB_criteria) THEN
          ! Weight criterion exceeded: rebalance only if this is a
          ! genuine degradation beyond what the last rebalance
          ! achieved (5% margin), or no baseline exists yet.
          IF(ppiclf_lbBaseMax .LT. 0.0D0 .OR.
     >       LB_curMax .GT. 1.05D0*ppiclf_lbBaseMax) THEN
            ppiclf_rebalance = .TRUE.
          END IF
        END IF
      END IF

      ! Logical OR comparison across MPI Ranks
      tstart = MPI_WTIME()
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
      SUBROUTINE ppiclf_comm_FindCellPartLB
!
! Tallies overlap fluid cells per bin for the weighted load balance,
! using the centroid-in-bin approximation: each local fluid cell whose
! center lies inside the bin domain is counted once, in the bin that
! contains its center. The MPI-boundary-face duplication performed
! later in MapOverlapGridPartLB is intentionally NOT included here,
! because it depends on the bin-to-rank map that this count is used
! to build (circular dependency).
!
! The count depends only on bin geometry (binb, n_bins, bins_dx) and
! the static fluid grid partition, so it is recomputed only when
! ppiclf_CellCountValid is .FALSE. (set by CreateBinPartLB whenever
! the bins change). MUST be called by all ranks (collective).
!
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      INTEGER*4  ie, ii, jj, kk, bin, ierr
      INTEGER*4  nb1, nb2, nb3, nb1xnb2
      INTEGER*4  l
      REAL*8     xc, yc, zc, inv_dx(3)
      REAL*8     ppiclf_comm_BinMapStencil
      EXTERNAL   ppiclf_comm_BinMapStencil
      REAL*8     tstart, tfinal

      ! Cells contribute work only when ppiclF maps an overlap fluid
      ! grid. Without overlap, no cell is ever mapped to a bin, so all
      ! cells are outside the ppiclf domain of operations and must
      ! carry zero weight.
      IF(.NOT. ppiclf_overlap) THEN
        ppiclf_CellCount = 0
        ppiclf_CellCountValid = .FALSE.
        RETURN
      END IF

      ! Consistent across ranks: binchanged and dyn_alloc decisions
      ! are made from globally reduced data.
      IF(ppiclf_CellCountValid) RETURN

      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)
      nb1xnb2 = nb1*nb2

      inv_dx(1) = 1.0D0 / ppiclf_bins_dx(1)
      inv_dx(2) = 1.0D0 / ppiclf_bins_dx(2)
      inv_dx(3) = 1.0D0 / ppiclf_bins_dx(3)

      ppiclf_CellCount  = 0
      ppiclf_CellMaxLen = 0.0D0

      DO ie = 1,ppiclf_nFVCells
        xc = ppiclf_fluid_grid(1,ie)
        yc = ppiclf_fluid_grid(2,ie)
        zc = ppiclf_fluid_grid(3,ie)
        ! Skip cells whose center is outside the bin domain
        ! (same criterion as MapOverlapGridPartLB)
        IF (xc .GT. ppiclf_binb(2)) CYCLE
        IF (yc .GT. ppiclf_binb(4)) CYCLE
        IF (zc .GT. ppiclf_binb(6)) CYCLE
        IF (xc .LT. ppiclf_binb(1)) CYCLE
        IF (yc .LT. ppiclf_binb(3)) CYCLE
        IF (zc .LT. ppiclf_binb(5)) CYCLE

        ii = FLOOR((xc - ppiclf_binb(1))*inv_dx(1))
        jj = FLOOR((yc - ppiclf_binb(3))*inv_dx(2))
        kk = FLOOR((zc - ppiclf_binb(5))*inv_dx(3))
        ! Round-off protection
        ii = MAX(0, MIN(ii, nb1-1))
        jj = MAX(0, MIN(jj, nb2-1))
        kk = MAX(0, MIN(kk, nb3-1))

        bin = ii + nb1*jj + nb1xnb2*kk
        ppiclf_CellCount(bin) = ppiclf_CellCount(bin) + 1
        ! Track the largest cell edge lengths seen in this bin
        ! (fluid_grid components 4:6 are the cell dimensions)
        DO ii = 1,3
          IF(ppiclf_fluid_grid(3+ii,ie) .GT.
     >       ppiclf_CellMaxLen(ii,bin))
     >       ppiclf_CellMaxLen(ii,bin) = ppiclf_fluid_grid(3+ii,ie)
        END DO
      END DO

      ! Sum cell counts per bin across MPI Ranks
      tstart = MPI_WTIME()
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_CellCount
     >                   ,ppiclf_totalBins ,MPI_INTEGER4, MPI_SUM
     >                   ,ppiclf_comm, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_CellMaxLen
     >                   ,3*ppiclf_totalBins ,MPI_DOUBLE_PRECISION
     >                   ,MPI_MAX ,ppiclf_comm, ierr)

#ifdef PERF
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_allreduces = PPICLF_TMPI_allreduces + tfinal
      PPICLF_TLoadBalance = PPICLF_TLoadBalance - tfinal
#endif

      ! Do not latch an empty count as valid: at initialization the
      ! load balance can run before the fluid grid is registered, and
      ! a later call must then rebuild the counts.
      IF(SUM(INT(ppiclf_CellCount,8)) .GT. 0_8) THEN
        ppiclf_CellCountValid = .TRUE.
      END IF

      ! Cache the per-bin P2C stencil ratio sMAP: its inputs
      ! (CellMaxLen, bins_dx, filter) change only when this routine
      ! reruns, so the per-stage consumers (weighted scan in
      ! FindParticlePartLB, calibration basis in LBCalibAccum) can use
      ! a table lookup instead of re-deriving it for every bin.
      DO ie = 0,ppiclf_totalBins-1
        ppiclf_binSmap(ie) =
     >    ppiclf_comm_BinMapStencil(ie, DBLE(ppiclf_CellCount(ie)))
      END DO

      RETURN
      END SUBROUTINE


!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_LBWeightCoef(aPP, bMAP, bPROJ, cCell)
!
! Analytic per-bin cost model for the load balance, from operation
! counts of the per-step kernels. With uniform intra-bin density the
! time attributable to a bin is
!
!   T(bin) ~ C_L*Np + C_PP*Np*(rho_p*Vpp) + C_MAP*Np*(sig_c*Vmap)
!          + C_PROJ*Np*(sig_c*Vsup) + C_CELL*Nc
!
! rho_p = Np/Vbin, sig_c = Nc/Vbin, and the volumes follow from the
! search structures:
!   Vpp  : 3x3x3 fine-cell P2P stencil, fine cell ~ nndist
!   Vmap : 3x3x3 fluid sub-cell P2C stencil; the sub-cell size a
!          bin's owner builds is ~1.5x the largest cell edge in its
!          region, estimated here PER BIN from ppiclf_CellMaxLen
!          (mesh grading makes this ratio vary strongly across bins)
!   Vsup : Gaussian projection support, KAPPA*filter per dimension
! Normalizing by Vbin and C_L:
!
!   W(bin) = (1 + bPROJ)*Np + aPP*Np^2
!          + bMAP*sMAP(bin)*Np*Nc + cCell*Nc
!   sMAP(bin) = PROD_d 3/nsf_d(bin),
!   nsf_d(bin) = MAX(1, FLOOR(bins_dx_d/(1.5*CellMaxLen_d(bin))))
!
! computed by ppiclf_comm_BinWeight. Only coefficient RATIOS affect
! the partition; C_L is the normalizer. Flop-count estimates (tune
! if kernels change):
!   C_L    per particle-step: 27-cell interpolation of ~30 fields
!          + SetYdot force models O(1e3) + RK update
!   C_PP   per P2P candidate: dist^2 + compare
!   C_MAP  per P2C candidate: dist^2 + 27-nearest partial insertion
!   C_PROJ per projected pair: Gaussian eval + 2*LRP_PRO accumulates
!   C_CELL per overlap cell-step: map build/zero/reduce bookkeeping
! Flop counts ignore memory effects; those shift coefficients O(1),
! not the model form.
!
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      REAL*8, INTENT(OUT) :: aPP, bMAP, bPROJ, cCell
      INTEGER*4 d, nf
      REAL*8    sPP, sPROJ
      ! C_PROJ and C_CELL are EFFECTIVE costs, not flop counts:
      ! projection is a strided 22-field scatter-accumulate (memory
      ! bound, ~5-10x its flop count) and C_CELL carries the per-cell
      ! per-step communication (grid refresh in + projection fields
      ! out, ~30 doubles through the crystal router), which gates the
      ! step exactly like compute.
      ! The coefficients now live in ppiclf_LB_C(1:5) of the
      ! DynamicAllocation module, initialized to the analytic
      ! flop-count priors. In PERF builds they are periodically
      ! recalibrated online against the measured per-rank kernel
      ! timers (ppiclf_comm_LBCalibrate; the channel timers are
      ! unconditional source code, so this is active in every build);
      ! the offline script
      ! calibrate_lb_coefficients.py performs the same fit from the
      ! per-rank ppiclf_perf_<nid>.csv files for diagnostics.
      REAL*8    C_L, C_PP, C_MAP, C_PROJ, C_CELL, KAPPA
      PARAMETER(KAPPA  = 1.0D0)
      C_L    = ppiclf_LB_C(1)
      C_PP   = ppiclf_LB_C(2)
      C_MAP  = ppiclf_LB_C(3)
      C_PROJ = ppiclf_LB_C(4)
      C_CELL = ppiclf_LB_C(5)

      CALL ppiclf_comm_LBStencilRatios(sPP, sPROJ)

      aPP   = (C_PP/C_L)*sPP
      bMAP  = C_MAP/C_L
      ! C_PROJ is fitted against the per-PARTICLE basis X(4)=Np in
      ! LBCalibAccum (capped projection stencil). The weight must use
      ! the same basis: no sPROJ factor here, otherwise the fitted
      ! seconds-per-particle is multiplied by a stencil-volume ratio
      ! that was never in the regression.
      bPROJ = C_PROJ/C_L
      cCell = C_CELL/C_L

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_LBStencilRatios(sPP, sPROJ)
!
! Dimensionless stencil-to-bin volume ratios of the P2P search and
! the projection support, shared by the weight model (LBWeightCoef)
! and the online calibration accumulator (LBCalibAccum).
!
! P2P path check, mirroring ppiclf_solve (fine grid is built and
! used iff nndist>0 and MAX(filter) >= 2*nndist):
!   fine   -> 3x3x3 nndist-cell stencil, sPP = PROD min(1,3/nf)
!   coarse -> 3x3x3 BIN stencil over binPartList, sPP = 27
!   none   -> no NN search, sPP = 0
!
      IMPLICIT NONE
      INCLUDE "PPICLF"
      REAL*8, INTENT(OUT) :: sPP, sPROJ
      INTEGER*4 d, nf
      REAL*8    KAPPA
      PARAMETER(KAPPA = 1.0D0)

      sPP = 0.0D0
      IF(ppiclf_nndist .GT. 0.0D0) THEN
        IF(MAX(ppiclf_filter(1),ppiclf_filter(2),ppiclf_filter(3))
     >     .GE. 2.0D0*ppiclf_nndist) THEN
          sPP = 1.0D0
          DO d = 1,3
            nf  = MAX(1, FLOOR(ppiclf_bins_dx(d)/ppiclf_nndist))
            sPP = sPP * MIN(1.0D0, 3.0D0/DBLE(nf))
          END DO
        ELSE
          sPP = 2.7D1
        END IF
      END IF

      sPROJ = 1.0D0
      DO d = 1,3
        sPROJ = sPROJ * (KAPPA*ppiclf_filter(d)/ppiclf_bins_dx(d))
      END DO

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_LBCalibAccum
!
! Online cost-model calibration, part 1 of 2 (sampling). Called once
! per FindParticlePartLB pass (after the ParticleCount allreduce, so
! the global per-bin counts are valid). Purely local - no collectives.
!
! Accumulates, for THIS rank over the current calibration window:
!   Y(k) - measured seconds in timer channel k
!   X(k) - the workload-model basis term k summed over the bins this
!          rank owns under the current bin-to-rank map
! matched channel by channel to the terms of the per-bin model
!   W(bin) = C_L*Np + C_PP*sPP*Np^2 + C_MAP*sMAP*Np*Nc
!          + C_PROJ*sPROJ*Np*Nc + C_CELL*Nc  (unnormalized form).
!
! The host may zero the PERF accumulators at its own logging cadence
! (ppiclf_solve_LogPerformance); a negative timer delta means such a
! reset occurred, in which case the current value IS the accumulation
! since the reset (the sliver between the previous sample and the
! reset is lost to the fit, which is harmless).
!
! Active in every build: the channel timers it reads are
! unconditional source code.
!
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
      INTEGER*4 b, k
      REAL*8    sPP, sPROJ, sMAP, rn, rc, tch(5)
      REAL*8    tcal0, tcal1

      IF(.NOT. ppiclf_LB_docal) RETURN
      tcal0 = MPI_WTIME()

      ! ---- timer channels matched to the model terms ----
      ! ch1: linear-in-Np work: interpolation, RK update, real-
      !      particle sub-bin map, and the user SetYdot bracket
      !      (which CONTAINS the force/heat models) MINUS the P2P
      !      search time nested inside it. Listing the force leaves
      !      separately here would double-count them, and leaving the
      !      P2P time in would contaminate the linear channel with
      !      quadratic work.
      ! ch2: P2P nearest-neighbor search, plus the fine-grid particle
      !      sub-bin map that serves it (subbinFineParticleMap loops
      !      over real+ghost particles, so it is particle-, not
      !      cell-proportional; it is the P2P search structure build).
      ! ch3: P2C map build + search
      ! ch4: projection (compute) -- capped stencil, linear in Np
      ! ch5: per-overlap-cell map + per-step per-cell communication:
      !      fluid-field transfer IN (TMPI_moveInt) and projection
      !      fields OUT (TMPI_movePro), plus the overlap-cell map
      !      transfer on rebuild stages. TMPI_moveInt is the first
      !      synchronizing collective of the stage and absorbs any
      !      host-solver stagger; set ppiclf_perf_sync=.TRUE. to move
      !      that wait into TEntrySync so the channel is clean.
      tch(1) = PPICLF_TInterp + PPICLF_TIntegrate
     >       + PPICLF_TsubbinRealMap
     >       + PPICLF_TUserYdot - PPICLF_TPPNNSearch
      tch(2) = PPICLF_TPPNNSearch + PPICLF_TsubbinFineMap
      tch(3) = PPICLF_TPCNNSearch + PPICLF_TsubbinCellMap
      tch(4) = PPICLF_TProject
      tch(5) = PPICLF_TMapOverlap + PPICLF_TMPI_moveOvlp
     >       + PPICLF_TMPI_moveInt + PPICLF_TMPI_movePro

      ! The PERF accumulators are zeroed at every IntegrateParticle
      ! entry, so at this end-of-stage sample point tch(k) IS the
      ! full-stage channel time; accumulate it directly. (The older
      ! snapshot-delta scheme sampled mid-stage, before the P2C map,
      ! projection, and cell communication had run, and systematically
      ! undercounted channels 3-5.)
      DO k = 1,5
        ppiclf_LB_Y(k) = ppiclf_LB_Y(k) + tch(k)
      END DO

      ! ---- basis sums over the bins THIS rank owns ----
      ! Accumulated only on stages where FindParticlePartLB just
      ! refreshed the global counts (ppiclf_LB_countsfresh), which
      ! keeps this O(N_b) loop on the same ppiclf_LB_checkfreq
      ! cadence as the weighted scan. The Y channels above still
      ! accumulate every stage; both sides scale by the same uniform
      ! factor (checkfreq) across ranks and channels, so the fitted
      ! coefficient RATIOS - the only thing the partition uses - are
      ! unaffected. sMAP comes from the per-bin cache.
      IF(ppiclf_LB_countsfresh) THEN
        IF (.NOT. ppiclf_LB_started) THEN
          ! First-ever sample: discard it (startup transients would
          ! contaminate the basis), open the short warmup window, and
          ! stamp the window start so the walltime adaptation is live
          ! at the FIRST refit.
          ppiclf_LB_started = .TRUE.
          ppiclf_LB_X = 0.0D0
          ppiclf_LB_Y = 0.0D0
          ppiclf_LB_calsteps = 0
          IF (ppiclf_LB_caltime .GT. 0.0D0) THEN
            ppiclf_LB_calfreq =
     >        MIN(ppiclf_LB_calfreq, ppiclf_LB_calfirst)
            IF (ppiclf_nid .EQ. 0) ppiclf_LB_lastfit = MPI_WTIME()
          END IF
          tcal1 = MPI_WTIME() - tcal0
          PPICLF_TLBCalib = PPICLF_TLBCalib + tcal1
          RETURN
        END IF
        CALL ppiclf_comm_LBStencilRatios(sPP, sPROJ)
        DO b = 0,ppiclf_totalBins-1
          IF(ppiclf_BinToRankMap(b) .EQ. ppiclf_nid) THEN
            rn = DBLE(ppiclf_ParticleCount(b))
            rc = DBLE(ppiclf_CellCount(b))
            sMAP = ppiclf_binSmap(b)
            ppiclf_LB_X(1) = ppiclf_LB_X(1) + rn
            ppiclf_LB_X(2) = ppiclf_LB_X(2) + sPP*rn*rn
            ppiclf_LB_X(3) = ppiclf_LB_X(3) + sMAP*rn*rc
            ! Projection basis is Np (capped stencil), not sPROJ*n*c.
            ppiclf_LB_X(4) = ppiclf_LB_X(4) + rn
            ppiclf_LB_X(5) = ppiclf_LB_X(5) + rc
          END IF
        END DO
        ppiclf_LB_calsteps = ppiclf_LB_calsteps + 1
      END IF
      tcal1 = MPI_WTIME() - tcal0
      PPICLF_TLBCalib = PPICLF_TLBCalib + tcal1
      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_LBCalibrate
!
! Online cost-model calibration, part 2 of 2 (refit). COLLECTIVE -
! must be called by all ranks at a globally consistent point (it is
! invoked from FindParticlePartLB, whose control flow is identical on
! all ranks). No-op until ppiclf_LB_calfreq samples have accumulated.
!
! Each rank contributes one observation per channel: (X_r, Y_r) with
! Y_r ~ c_k * X_r. The through-origin least-squares fit across ranks,
!   c_k = SUM_r(Y_r X_r) / SUM_r(X_r^2),
! is computed from two 5-vector allreduce sums. Because only the
! RATIOS of the coefficients influence the partition, the raw fits
! (seconds per basis unit) are rescaled so channel 1 retains the
! magnitude of its analytic prior, then shrunk 25% toward the priors
! (ridge regularization against noisy windows) and blended 50/50 with
! the previous values (EWMA, prevents thrash). Channels with no
! activity in the window (X = 0, e.g. P2P search disabled) keep
! their current value. The futility baseline of the rebalance
! criterion is invalidated since the weight scale changed.
!
! Active in every build: the channel timers it reads are
! unconditional source code.
!
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
      INTEGER*4 k, ierr, newfreq
      REAL*8    sloc(10), sglb(10), raw(5), scale, cmix
      REAL*8    tcal0, tcal1, tnow, elapsed
      REAL*8    WPRIOR
      PARAMETER(WPRIOR = 0.25D0)

      IF(.NOT. ppiclf_LB_docal) RETURN
      IF(ppiclf_LB_calsteps .LT. ppiclf_LB_calfreq) RETURN
      tcal0 = MPI_WTIME()

      DO k = 1,5
        sloc(k)   = ppiclf_LB_Y(k)*ppiclf_LB_X(k)
        sloc(5+k) = ppiclf_LB_X(k)*ppiclf_LB_X(k)
      END DO
      CALL MPI_ALLREDUCE(sloc, sglb, 10, MPI_DOUBLE_PRECISION,
     >                   MPI_SUM, ppiclf_comm, ierr)

      ! Raw through-origin fits (seconds per basis unit)
      DO k = 1,5
        IF(sglb(5+k) .GT. 1.0D-300) THEN
          raw(k) = MAX(0.0D0, sglb(k)/sglb(5+k))
        ELSE
          raw(k) = -1.0D0     ! channel inactive this window
        END IF
      END DO

      ! The linear channel anchors the normalization; without a
      ! usable fit for it, skip this window entirely.
      IF(raw(1) .GT. 0.0D0) THEN
        scale = ppiclf_LB_prior(1)/raw(1)
        DO k = 1,5
          IF(raw(k) .GE. 0.0D0) THEN
            cmix = (1.0D0 - WPRIOR)*raw(k)*scale
     >           + WPRIOR*ppiclf_LB_prior(k)
            IF (ppiclf_LB_nrefits .EQ. 0) THEN
              ! First refit: replace the prior-initialized values
              ! outright instead of EWMA-blending with them, so the
              ! first measured coefficients take effect immediately.
              ppiclf_LB_C(k) = cmix
            ELSE
              ppiclf_LB_C(k) = 0.5D0*ppiclf_LB_C(k) + 0.5D0*cmix
            END IF
          END IF
        END DO
        ppiclf_LB_nrefits = ppiclf_LB_nrefits + 1
        ! Weight scale changed: re-baseline the futility criterion.
        ppiclf_lbBaseMax = -1.0D0
        IF(ppiclf_nid .EQ. 0) THEN
          PRINT*, 'PPICLF: LB coefficients recalibrated (C_L, C_PP,',
     >            ' C_MAP, C_PROJ, C_CELL):'
          PRINT*, '        ', ppiclf_LB_C(1), ppiclf_LB_C(2),
     >            ppiclf_LB_C(3), ppiclf_LB_C(4), ppiclf_LB_C(5)
        END IF
      END IF

      ppiclf_LB_X = 0.0D0
      ppiclf_LB_Y = 0.0D0
      ppiclf_LB_calsteps = 0

      ! ---- walltime-adaptive refit window ----
      ! Rescale the sample count so the NEXT window spans about
      ! ppiclf_LB_caltime seconds of walltime: freshness on expensive
      ! problems, sample depth on cheap ones. Rank 0 measures the
      ! elapsed window on its own clock (MPI clocks need not be
      ! synchronized across ranks) and broadcasts the new window
      ! length, so every rank gates the next refit identically.
      IF (ppiclf_LB_caltime .GT. 0.0D0) THEN
        newfreq = ppiclf_LB_calfreq
        IF (ppiclf_nid .EQ. 0) THEN
          tnow = MPI_WTIME()
          IF (ppiclf_LB_lastfit .GT. 0.0D0) THEN
            elapsed = tnow - ppiclf_LB_lastfit
            IF (elapsed .GT. 1.0D-3) THEN
              newfreq = NINT(DBLE(ppiclf_LB_calfreq)
     >                       *ppiclf_LB_caltime/elapsed)
            ELSE
              newfreq = ppiclf_LB_calmax
            END IF
            newfreq = MAX(ppiclf_LB_calmin,
     >                    MIN(ppiclf_LB_calmax, newfreq))
            IF (newfreq .NE. ppiclf_LB_calfreq) THEN
              WRITE(6,'(A,I6,A,F8.2,A)')
     >          ' ppiclF LB calibration window -> ', newfreq,
     >          ' samples (last window ', elapsed, ' s)'
            END IF
          END IF
          ppiclf_LB_lastfit = tnow
        END IF
        CALL MPI_BCAST(newfreq, 1, MPI_INTEGER4, 0,
     >                 ppiclf_comm, ierr)
        ppiclf_LB_calfreq = newfreq
      END IF

      tcal1 = MPI_WTIME() - tcal0
      PPICLF_TLBCalib = PPICLF_TLBCalib + tcal1
      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION ppiclf_comm_BinWeight(bin, aPP, bMAP,
     >                                                bPROJ, cCell)
!
! W(bin) = (1 + bPROJ)*Np + aPP*Np^2 + bMAP*sMAP*Np*Nc + cCell*Nc
! with the per-bin P2C stencil ratio sMAP from the bin's largest
! cell edges (see LBWeightCoef). Bins without cells: Nc terms vanish.
! NOTE (model revision, validated on the 80-rank blast dataset): the
! projection stencil is CAPPED at nnearest <= 27 cells per particle,
! so projection work is linear in Np, not proportional to Np*Nc
! (measured through-origin R^2: 0.996 for ~Np vs 0.467 for ~Np*Nc,
! with only a 1.2x per-particle rate spread across ranks). bPROJ is
! therefore the per-particle projection cost relative to C_L.
!
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INTEGER*4, INTENT(IN) :: bin
      REAL*8, INTENT(IN)    :: aPP, bMAP, bPROJ, cCell
      REAL*8    rNp, rNc, sMAP

      rNp = DBLE(ppiclf_ParticleCount(bin))
      rNc = DBLE(ppiclf_CellCount(bin))

      ! Table lookup: sMAP is cached per bin by FindCellPartLB
      ! whenever its inputs change (see ppiclf_binSmap).
      sMAP = ppiclf_binSmap(bin)

      ppiclf_comm_BinWeight = (1.0D0 + bPROJ)*rNp + aPP*rNp*rNp
     >                      + bMAP*sMAP*rNp*rNc
     >                      + cCell*rNc

      RETURN
      END FUNCTION

!-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION ppiclf_comm_BinMapStencil(bin, rNc)
!
! Per-bin P2C stencil ratio sMAP, mirroring SBParticleToCellMap.
! The kernel walks each bin at its OWN sub-cell resolution
! (nsf_d from the bin's own coarsest cell) but over a search reach
! set by the COARSEST cell in the +-1 bin neighborhood (see
! binReachDilate), capped at the bin size. The fraction of a bin's
! cells a particle's box spans, counted in own sub-cells, is
!   frac_d = 2*reach_d/bins_dx_d + 1/nsf_d,  reach_d = MIN(1.5*Ldil_d, bins_dx_d)
! capped at 3 (a full-bin reach touches at most 3 bins per dim).
! Coarse path (no fine map): 3x3x3 BIN stencil, sMAP = 27.
!
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INTEGER*4, INTENT(IN) :: bin
      REAL*8, INTENT(IN)    :: rNc
      INTEGER*4 d, nsf, ii, jj, kk, i2, j2, k2, nb1, nb2, nb3, src
      REAL*8    sMAP, reach, frac, ldil(3)
      LOGICAL   fineOK

      sMAP = 1.0D0
      IF(rNc .GT. 0.0D0) THEN
        nb1 = ppiclf_n_bins(1)
        nb2 = ppiclf_n_bins(2)
        nb3 = ppiclf_n_bins(3)
        kk  = bin/(nb1*nb2)
        jj  = (bin - kk*nb1*nb2)/nb1
        ii  = bin - kk*nb1*nb2 - jj*nb1

        ! +-1 neighborhood max of the per-bin max cell edge
        ! (global field, MPI_MAX-reduced in FindCellPartLB)
        ldil = 0.0D0
        DO k2 = MAX(0,kk-1), MIN(nb3-1,kk+1)
         DO j2 = MAX(0,jj-1), MIN(nb2-1,jj+1)
          DO i2 = MAX(0,ii-1), MIN(nb1-1,ii+1)
            src = i2 + nb1*j2 + nb1*nb2*k2
            DO d = 1,3
              ldil(d) = MAX(ldil(d), ppiclf_CellMaxLen(d,src))
            END DO
          END DO
         END DO
        END DO

        fineOK = .FALSE.
        DO d = 1,3
          IF(ppiclf_CellMaxLen(d,bin) .GT. 0.0D0) THEN
            nsf = MAX(1, FLOOR(ppiclf_bins_dx(d)
     >                    /(1.5D0*ppiclf_CellMaxLen(d,bin))))
          ELSE
            nsf = 1
          END IF
          IF(ppiclf_bins_dx(d)/DBLE(nsf) .LT.
     >       0.5D0*ppiclf_filter(d)) fineOK = .TRUE.

          reach = MIN(1.5D0*ldil(d), ppiclf_bins_dx(d))
          IF(reach .LE. 0.0D0) reach = ppiclf_bins_dx(d)
          frac  = 2.0D0*reach/ppiclf_bins_dx(d) + 1.0D0/DBLE(nsf)
          sMAP  = sMAP * MIN(3.0D0, frac)
        END DO
        IF(.NOT. fineOK) sMAP = 2.7D1
      END IF

      ppiclf_comm_BinMapStencil = sMAP

      RETURN
      END FUNCTION
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_PartLoadBalance

      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
      REAL*8 ppiclf_pt0, tstart, tfinal


      INTEGER*4   ierr, i, bin, irank
     >           ,nb1, nb2, nb1xnb2, j, k
     >           ,iloop, jloop, kloop
     >           ,partSum, prevPartSum, capPart, maxCnt
      INTEGER*4, ALLOCATABLE :: rankCnt(:)
      LOGICAL     forceCut
      REAL*8      binWeight, weightSum, prevWeightSum
     >           ,targetWeight, remainingWeight
     >           ,aPP, bMAP, bPROJ, cCell
      REAL*8      ppiclf_comm_BinWeight
      EXTERNAL    ppiclf_comm_BinWeight
      INTEGER*4   stride(3), stride_L, stride_M, stride_S
      INTEGER*4   bin_L, bin_M

      ! Calculate BTRM on root processor and broadcast
      ! to all others. This ensures all processors have same
      ! key global mapping.
      ppiclf_pt0 = MPI_WTIME()
      ! Refresh per-bin overlap cell counts if bin geometry changed.
      ! Collective call - must execute on ALL ranks.
      CALL ppiclf_comm_FindCellPartLB

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

        ! Per-bin work estimate (derivation in LBWeightCoef):
        !   W(bin) = Np + aPP*Np^2 + bPC*Np*Nc + cCell*Nc
        ! normalized by the linear per-particle cost; Nc counted
        ! centroid-in-bin by ppiclf_comm_FindCellPartLB.
        ! Analytic cost model coefficients (see LBWeightCoef)
        CALL ppiclf_comm_LBWeightCoef(aPP, bMAP, bPROJ, cCell)
        remainingWeight = 0.0D0
        DO i = 0,ppiclf_totalBins-1
          remainingWeight = remainingWeight
     >        + ppiclf_comm_BinWeight(i, aPP, bMAP, bPROJ, cCell)
        END DO
        targetWeight = remainingWeight/DBLE(ppiclf_np)
        weightSum = 0.0D0
        irank = 0
        ! HARD CAPACITY CONSTRAINT: buffers in MoveParticlePartLB and
        ! the iprop/y common blocks hold at most PPICLF_LPART
        ! particles per rank. Weight balance is allowed to skew
        ! particle COUNTS, so counts must be capped explicitly: force
        ! a cut before any bin that would push the running rank past
        ! 90% of PPICLF_LPART (margin covers drift until the next
        ! rebalance). A single bin above the cap cannot be split; it
        ! is isolated on its own rank and the validation below aborts
        ! cleanly if it exceeds PPICLF_LPART itself.
        capPart = PPICLF_LPART - PPICLF_LPART/10
        partSum = 0
        !Iterate through all bins using the fast, branchless stride logic
        DO iloop = 0,(ppiclf_n_bins(ppiclf_dL) - 1)
          bin_L = iloop * stride_L
          DO jloop = 0,(ppiclf_n_bins(ppiclf_dM) - 1)
            bin_M = bin_L + jloop * stride_M
            DO kloop = 0,(ppiclf_n_bins(ppiclf_dS) - 1)
              bin = bin_M + kloop * stride_S
              binWeight =
     >          ppiclf_comm_BinWeight(bin, aPP, bMAP, bPROJ, cCell)
              prevWeightSum = weightSum
              weightSum = weightSum + binWeight
              prevPartSum = partSum
              partSum = partSum + ppiclf_ParticleCount(bin)
              forceCut = (partSum .GT. capPart) .AND.
     >                   (prevPartSum .GT. 0)
              IF(weightSum .GE. targetWeight .OR. forceCut) THEN
                IF(irank .LT. ppiclf_np - 1) THEN
                  ! Check if the previous state (undershoot) was better
                  ! than the current state (overshoot). A capacity
                  ! forceCut always cuts BEFORE the current bin.
                  IF(forceCut .OR.
     >               (targetWeight - prevWeightSum) .LE.
     >               (weightSum - targetWeight)         ) THEN

                    remainingWeight = remainingWeight
     >                                - prevWeightSum
                    ! Assign the current bin to the NEXT rank.
                    irank = irank + 1
                    ppiclf_BinToRankMap(bin) = irank
                    ! Start the new rank's tally with the current bin
                    weightSum = binWeight
                    partSum = ppiclf_ParticleCount(bin)
                  ELSE
                    ! Keep this bin on the current rank.
                    remainingWeight = remainingWeight
     >                                - weightSum
                    ppiclf_BinToRankMap(bin) = irank
                    irank = irank + 1
                    ! Reset weight tally for the next rank
                    weightSum = 0.0D0
                    partSum = 0
                  END IF
                  targetWeight = remainingWeight
     >                             / DBLE(ppiclf_np - irank)
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

        ! --- Capacity validation -----------------------------------
        ! A rank assigned more than PPICLF_LPART particles overruns
        ! the crystal-router receive buffers in MoveParticlePartLB
        ! (segfault inside gslib memcpy) BEFORE the post-transfer
        ! LPART check can fire. Fail here, at the source, with an
        ! actionable message instead.
        ALLOCATE(rankCnt(0:ppiclf_np-1))
        rankCnt = 0
        DO i = 0,ppiclf_totalBins-1
          j = ppiclf_BinToRankMap(i)
          rankCnt(j) = rankCnt(j) + ppiclf_ParticleCount(i)
        END DO
        maxCnt = 0
        DO i = 0,ppiclf_np-1
          maxCnt = MAX(maxCnt, rankCnt(i))
        END DO
        DEALLOCATE(rankCnt)
        IF(maxCnt .GT. PPICLF_LPART) THEN
          PRINT*, '***ERROR*** Load balance assigns', maxCnt,
     >            ' particles to one rank; PPICLF_LPART =',
     >            PPICLF_LPART
          PRINT*, '***ERROR*** Increase PPICLF_LPART, add ranks,',
     >            ' or reduce weight concentration.'
          CALL ppiclf_exittr('LB exceeds PPICLF_LPART$',0.0D0,maxCnt)
        END IF

        ! Inform the futility guard in FindParticlePartLB that the
        ! next check should re-baseline against the achieved balance.
      END IF ! root Processor
      ppiclf_lbJustBal = .TRUE.

      ! Share BTRM to all processors 
      tstart = MPI_WTIME()
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
      ppiclf_pt0 = MPI_WTIME()
      CALL ppiclf_comm_setRankBoundaries
#ifdef PERF
      PPICLF_TRankBounds = PPICLF_TRankBounds
     >     + (MPI_WTIME() - ppiclf_pt0)
#endif
      ppiclf_pt0 = MPI_WTIME()
      CALL ppiclf_comm_setEmptyIndicator
#ifdef PERF
      PPICLF_TEmptyInd = PPICLF_TEmptyInd
     >     + (MPI_WTIME() - ppiclf_pt0)
#endif
      ppiclf_pt0 = MPI_WTIME()
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
      ! Heap, not automatic: at production LPART this buffer is
      ! O(100 MB). As an automatic array the columns beyond the local
      ! npart are first TOUCHED by gslib's receive memcpy, so a rank
      ! that receives many more particles than it sent walks into
      ! unmapped stack pages and segfaults inside the crystal router.
      REAL*8, ALLOCATABLE :: rtemp(:,:)
      INTEGER*4 i, icount, j0, ierr
      REAL*8    tstart, tfinal
!
      ALLOCATE(rtemp(rtempLim,PPICLF_LPART))

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

      tstart = MPI_WTIME()
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
        
      DEALLOCATE(rtemp)
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
     >          centeri(3), inv_dx(3)
      LOGICAL   partl, ErrorFound, MapCell 
      LOGICAL, ALLOCATABLE :: rank_is_mapped(:)
      REAL*8    tstart, tfinal
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

      PPICLF_T_OverlapCells_received  = ppiclf_nCells_FV2PICL
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_moveOvlp = PPICLF_TMPI_moveOvlp + tfinal
      PPICLF_TMapOverlap = PPICLF_TMapOverlap - tfinal

      ! Find distance check for interpolation.
      ! This is 1.5*MaxCellLength to ensure that at least
      ! 27 neighboring cells are mapped. ppiclf_MaxCellLen is the
      ! per-rank (local) max cell length, kept in common so the
      ! sub-bin refinement decision can read it.
      ppiclf_MaxCellLen(1) = 0.0D0
      ppiclf_MaxCellLen(2) = 0.0D0
      ppiclf_MaxCellLen(3) = 0.0D0
      ! Loop through overlapcells mapped to bin
      DO ie = 1,ppiclf_nCells_FV2PICL 
        DO l = 1,3
          ! Find max cell lengths in all dimensions
          IF(ppiclf_picl_grid(3+l,ie) .GT. ppiclf_MaxCellLen(l))
     >      ppiclf_MaxCellLen(l) = ppiclf_picl_grid(3+l,ie)
        END DO !l
      END DO !ie
      ! Update ppiclf filter, which is important as bin 
      ! boundaries grow and captures larger cells
      DO l = 1,3
        ppiclf_filter(l) = ppiclf_MaxCellLen(l)*1.51D0
      END DO
      ! Find max filter across processors (all overlap cells considered)
      tstart = MPI_WTIME()
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_filter
     >                   ,3 ,MPI_DOUBLE_PRECISION
     >                   ,MPI_MAX ,ppiclf_comm, ierr)
      tfinal = MPI_WTIME() - tstart
      PPICLF_TMPI_allreduces = PPICLF_TMPI_allreduces + tfinal
      PPICLF_TMapOverlap = PPICLF_TMapOverlap - tfinal

      DO l = 1,3
        ! Multiply by 1.5 so particle near face will
        ! find center one cell over in farthest direction
        ! Only need cells on this processor - not all overlap cells
        ppiclf_interp_dchk(l) = ppiclf_MaxCellLen(l)*1.5D0
      END DO

      ! New overlap-cell enumeration is live: bump the epoch so any
      ! CSR built against the previous enumeration reads as stale.
      CALL ppiclf_bumpMapEpoch
      RETURN
      END SUBROUTINE 
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_subbinRealParticleMap
      
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"


      INTEGER*4  i, tempSBin, i_SBin(3) 
     
      ! Quadruple it for ghost particles in the future. The global
      ! counts may be up to (checkfreq-1) stages stale; the 4x
      ! headroom plus the monotone MAX make that irrelevant here.
      ppiclf_maxParticlePerBin = MAX(ppiclf_maxParticlePerBin,
     >                               4*MAXVAL(ppiclf_ParticleCount))

      CALL ppiclf_allocate_BTP(ppiclf_total_SBin,
     >                         ppiclf_maxParticlePerBin)

      ! ROOT-CAUSE FIX: always invalidate the per-subbin particle
      ! lists, even when this rank currently has no particles.
      ! The old early return skipped the reset, leaving last
      ! stage's list (INCLUDING appended negative ghost entries)
      ! visible to consumers such as CreateGhostPartLB once
      ! particles return to this rank -> ppiclf_y(:,negative).
      ppiclf_binPartCount = 0

      IF(ppiclf_npart .LT. 1) RETURN

      ! This creates a list of local particles contained within
      ! each bin that resides on this rank.
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
      ! {-1,0,+1}. 13 == no wrap.
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

            ! Instant bailout for empty bins. Test the LOCAL per-bin
            ! count (rebuilt by subbinRealParticleMap every stage):
            ! it is fresh even when the GLOBAL counts are on the
            ! ppiclf_LB_checkfreq cadence, and it also skips bins
            ! whose particles are all remote (they generate no ghosts
            ! from this rank anyway).
            IF(ppiclf_binPartCount(tempSBin) .EQ. 0) CYCLE
 
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
              ! Negative entries are ghost particles (see
              ! subbinGhostParticleMap). Ghosts never seed ghosts
              ! and their coordinates are not in ppiclf_y.
              IF (ip .LE. 0) CYCLE
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
      REAL*8    tstart, tfinal
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
!     it by POSITION, using a compressed (CSR) layout: a flat candidate
!     list ppiclf_fineFlat indexed by 1-based row pointers
!     ppiclf_fineOffset. Fine bin b owns
!       ppiclf_fineFlat( fineOffset(b) : fineOffset(b+1)-1 ).
!     This is a separate, finer grid than the coarse iprop sub-bins
!     (sized by ppiclf_filter) used for comm and ghost creation.
!     Storage convention (matches ppiclf_solve_NearestNeighborSB):
!       real  particle i -> stored as +i (position ppiclf_y(1:3,i))
!       ghost particle i -> stored as -i (ppiclf_rprop_gp(1:3,i))
!
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
!
! Internal:
!
      INTEGER*4  i, d, tempBin, ng, ntot, pos, nretry
      REAL*8     xlo(3), xhi(3), boxLen(3), fineCut, fineLen(3), pad,
     >           xlo0(3), xhi0(3)
      INTEGER*8  nfb8, nbytes8
      ! 250 MBytes
      INTEGER*8, PARAMETER :: FINEBIN_BYTE_BUDGET = 250000000_8

      IF(ppiclf_npart .LT. 1) THEN
        ! Nothing to build: veto the fine path so the NN search
        ! falls back to the (freshly reset) coarse sub-bin list
        ! instead of walking stale fine-grid arrays.
        ppiclf_useFineGrid = .FALSE.
        RETURN
      END IF
      ng = ppiclf_npart_gp
      IF(ng .LT. 0) ng = 0
      fineCut = ppiclf_nndist

      ! Bounding box over real + ghost positions. The pad below keeps
      ! boxLen >= fineCut even for a single coincident point, so every
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

      ! Size the grid; coarsen (double fineCut) if the bin count would
      ! overflow INTEGER*4 or the CSR footprint exceeds the budget. The
      ! cost is now ~(total_fineBin+1 + npart+ng) INTEGER*4 words -- the
      ! old maxPartPerFineBin multiplier is gone. fineLen only grows
      ! under coarsening, so fineLen >= nndist still holds and the
      ! consumer's +/-1 stencil stays complete.
      DO d = 1,3
        xlo0(d) = xlo(d)
        xhi0(d) = xhi(d)
      END DO
      nretry = 0
      DO
        IF(nretry .EQ. 1) PRINT*,
     >    'fine subbin coarsened: nid',ppiclf_nid,
     >    ' nndist',ppiclf_nndist
        pad  = 0.5D0*fineCut
        nfb8 = 1_8
        DO d = 1,3
          xlo(d)               = xlo0(d) - pad
          xhi(d)               = xhi0(d) + pad
          boxLen(d)            = xhi(d) - xlo(d)
          ppiclf_nFine(d)      = MAX(1, FLOOR(boxLen(d)/fineCut))
          fineLen(d)           = boxLen(d)/DBLE(ppiclf_nFine(d))
          ppiclf_fineLo(d)     = xlo(d)
          ppiclf_fineInvLen(d) = 1.0D0/fineLen(d)
          nfb8 = nfb8 * INT(ppiclf_nFine(d),8)
        END DO
        nbytes8 = (nfb8 + 1_8 + INT(ppiclf_npart+ng,8))*4_8
        IF(nfb8    .LE. INT(HUGE(ppiclf_total_fineBin),8) .AND.
     >     nbytes8 .LE. FINEBIN_BYTE_BUDGET) EXIT
        nretry  = nretry + 1
        fineCut = fineCut*2.0D0
        IF(nretry .GT. 64) CALL ppiclf_exittr(
     >    'Fine grid keeps exceeding INT*4/budget$',
     >     fineCut/2.0D0, ppiclf_nid)
      END DO
      ppiclf_total_fineBin = INT(nfb8,4)

      ntot = ppiclf_npart + ng
      CALL ppiclf_allocate_FineCSR(ppiclf_total_fineBin, ntot)

      ! ---- Pass 1: count items per fine bin --------------------
      DO i = 0,ppiclf_total_fineBin-1
        ppiclf_finePartCount(i) = 0
      END DO
      DO i = 1,ppiclf_npart
        tempBin = 
     >   MAX(0,MIN(FLOOR((ppiclf_y(1,i)-ppiclf_fineLo(1))
     >   *ppiclf_fineInvLen(1)), ppiclf_nFine(1)-1))+ ppiclf_nFine(1)
     >   *MAX(0,MIN(FLOOR((ppiclf_y(2,i)-ppiclf_fineLo(2))
     >   *ppiclf_fineInvLen(2)), ppiclf_nFine(2)-1))
     >   + ppiclf_nFine(1)*ppiclf_nFine(2)
     >   *MAX(0,MIN(FLOOR((ppiclf_y(3,i)-ppiclf_fineLo(3))
     >   *ppiclf_fineInvLen(3)), ppiclf_nFine(3)-1))

        ppiclf_finePartCount(tempBin) = ppiclf_finePartCount(tempBin)+1
      END DO
      DO i = 1,ng
        tempBin = 
     >   MAX(0,MIN(FLOOR((ppiclf_rprop_gp(1,i)-ppiclf_fineLo(1))
     >   *ppiclf_fineInvLen(1)), ppiclf_nFine(1)-1)) + ppiclf_nFine(1)
     >   *MAX(0,MIN(FLOOR((ppiclf_rprop_gp(2,i)-ppiclf_fineLo(2))
     >   *ppiclf_fineInvLen(2)), ppiclf_nFine(2)-1))
     >   + ppiclf_nFine(1)*ppiclf_nFine(2)
     >   *MAX(0,MIN(FLOOR((ppiclf_rprop_gp(3,i)-ppiclf_fineLo(3))
     >   *ppiclf_fineInvLen(3)), ppiclf_nFine(3)-1))

        ppiclf_finePartCount(tempBin) = ppiclf_finePartCount(tempBin)+1
      END DO

      ! ---- Prefix sum -> row pointers, then seed cursors -------
      ppiclf_fineOffset(0) = 1
      DO i = 0,ppiclf_total_fineBin-1
        ppiclf_fineOffset(i+1) = ppiclf_fineOffset(i)
     >                         + ppiclf_finePartCount(i)
      END DO
      DO i = 0,ppiclf_total_fineBin-1
        ppiclf_finePartCount(i) = ppiclf_fineOffset(i)  ! cursor
      END DO

      ! ---- Pass 2: scatter (+i reals, -i ghosts) ---------------
      DO i = 1,ppiclf_npart
        tempBin =  
     >    MAX(0,MIN(FLOOR((ppiclf_y(1,i)-ppiclf_fineLo(1))
     >   *ppiclf_fineInvLen(1)), ppiclf_nFine(1)-1))+ ppiclf_nFine(1)
     >   *MAX(0,MIN(FLOOR((ppiclf_y(2,i)-ppiclf_fineLo(2))
     >   *ppiclf_fineInvLen(2)), ppiclf_nFine(2)-1))
     >   + ppiclf_nFine(1)*ppiclf_nFine(2)
     >   *MAX(0,MIN(FLOOR((ppiclf_y(3,i)-ppiclf_fineLo(3))
     >   *ppiclf_fineInvLen(3)), ppiclf_nFine(3)-1))

        pos     = ppiclf_finePartCount(tempBin)
        ppiclf_fineFlat(pos)          = i
        ppiclf_finePartCount(tempBin) = pos + 1
      END DO
      DO i = 1,ng
        tempBin =  
     >   MAX(0,MIN(FLOOR((ppiclf_rprop_gp(1,i)-ppiclf_fineLo(1))
     >   *ppiclf_fineInvLen(1)), ppiclf_nFine(1)-1)) + ppiclf_nFine(1)
     >   *MAX(0,MIN(FLOOR((ppiclf_rprop_gp(2,i)-ppiclf_fineLo(2))
     >   *ppiclf_fineInvLen(2)), ppiclf_nFine(2)-1))
     >   + ppiclf_nFine(1)*ppiclf_nFine(2)
     >   *MAX(0,MIN(FLOOR((ppiclf_rprop_gp(3,i)-ppiclf_fineLo(3))
     >   *ppiclf_fineInvLen(3)), ppiclf_nFine(3)-1))

        pos     = ppiclf_finePartCount(tempBin)
        ppiclf_fineFlat(pos)          = -i
        ppiclf_finePartCount(tempBin) = pos + 1
      END DO

      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------
!
      SUBROUTINE ppiclf_comm_cellBinImages(ie, wbin, iside, nplace)
!
!     Enumerate the WINDOW-BIN placements of overlap cell ie: its
!     direct bin plus single-rank periodic images (up to 8 total).
!     iside(d,m) records, per dimension, how placement m arrived:
!       0 = direct (use the cell's true position for sub-indexing)
!       1 = periodic image: cell near domain LO, placed in the TOP
!           bin -> top boundary sub-cell layer
!       2 = periodic image: cell near domain HI, placed in bin 0
!           -> bottom boundary sub-cell layer (index 0)
!     Mirrors the coarse map's periodic-image logic.
!
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INTEGER*4, INTENT(IN)  :: ie
      INTEGER*4, INTENT(OUT) :: wbin(8), iside(3,8), nplace
      INTEGER*4 i, j, k, i_SBin(3), iT(3), isd(3)

      DO i = 1,3
        i_SBin(i) = ppiclf_cell_map(4+i,ie) - ppiclf_binOffset(i)
      END DO

      nplace = 0
      DO i = 0,1
        IF(i .EQ. 0) THEN
          iT(1) = i_SBin(1)
          isd(1) = 0
          IF(iT(1).LT.0 .OR. iT(1).GT.ppiclf_nSBin(1)-1) CYCLE
        ELSE
          IF(ppiclf_linperiodic(1).AND.ppiclf_EqualDomain(1)) THEN
            IF(i_SBin(1).LE.0) THEN
              iT(1) = ppiclf_nSBin(1)-1
              isd(1) = 1
              IF(iT(1).EQ.i_SBin(1)) CYCLE
            ELSE IF(i_SBin(1).GE.ppiclf_nSBin(1)-1) THEN
              iT(1) = 0
              isd(1) = 2
              IF(iT(1).EQ.i_SBin(1)) CYCLE
            ELSE
              CYCLE
            END IF
          ELSE
            CYCLE
          END IF
        END IF
        DO j = 0,1
          IF(j .EQ. 0) THEN
            iT(2) = i_SBin(2)
            isd(2) = 0
            IF(iT(2).LT.0 .OR. iT(2).GT.ppiclf_nSBin(2)-1) CYCLE
          ELSE
            IF(ppiclf_linperiodic(2).AND.ppiclf_EqualDomain(2)) THEN
              IF(i_SBin(2).LE.0) THEN
                iT(2) = ppiclf_nSBin(2)-1
                isd(2) = 1
                IF(iT(2).EQ.i_SBin(2)) CYCLE
              ELSE IF(i_SBin(2).GE.ppiclf_nSBin(2)-1) THEN
                iT(2) = 0
                isd(2) = 2
                IF(iT(2).EQ.i_SBin(2)) CYCLE
              ELSE
                CYCLE
              END IF
            ELSE
              CYCLE
            END IF
          END IF
          DO k = 0,1
            IF(k .EQ. 0) THEN
              iT(3) = i_SBin(3)
              isd(3) = 0
              IF(iT(3).LT.0 .OR. iT(3).GT.ppiclf_nSBin(3)-1) CYCLE
            ELSE
              IF(ppiclf_linperiodic(3).AND.ppiclf_EqualDomain(3)) THEN
                IF(i_SBin(3).LE.0) THEN
                  iT(3) = ppiclf_nSBin(3)-1
                  isd(3) = 1
                  IF(iT(3).EQ.i_SBin(3)) CYCLE
                ELSE IF(i_SBin(3).GE.ppiclf_nSBin(3)-1) THEN
                  iT(3) = 0
                  isd(3) = 2
                  IF(iT(3).EQ.i_SBin(3)) CYCLE
                ELSE
                  CYCLE
                END IF
              ELSE
                CYCLE
              END IF
            END IF
            nplace = nplace + 1
            wbin(nplace) = iT(1) + ppiclf_nSBin(1)*iT(2)
     >                   + ppiclf_nSBin(1)*ppiclf_nSBin(2)*iT(3)
            iside(1,nplace) = isd(1)
            iside(2,nplace) = isd(2)
            iside(3,nplace) = isd(3)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE

      SUBROUTINE ppiclf_comm_fluidCellSubImages(ie, place, nplace)
!
!     Return the PACKED per-bin fine sub-cell index/indices for overlap
!     cell ie, including single-rank periodic images. Direct placements
!     sub-index by the cell's position in the bin's OWN grid; periodic
!     images land in the image bin's boundary sub-cell layer (the true
!     image position is within one cell length of the domain face).
!
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INTEGER*4, INTENT(IN)  :: ie
      INTEGER*4, INTENT(OUT) :: place(8), nplace
      INTEGER*4 wbin(8), iside(3,8), m, d, b, fw(3), ig(3)
      REAL*8    binLow, h

      CALL ppiclf_comm_cellBinImages(ie, wbin, iside, nplace)

      DO m = 1,nplace
        b = wbin(m)
        ! global (i,j,k) of window bin b, for its low corner
        ig(1) = MOD(b, ppiclf_nSBin(1))
        ig(2) = MOD(b/ppiclf_nSBin(1), ppiclf_nSBin(2))
        ig(3) = b/(ppiclf_nSBin(1)*ppiclf_nSBin(2))
        DO d = 1,3
          IF(iside(d,m) .EQ. 0) THEN
            h = ppiclf_bins_dx(d)/DBLE(ppiclf_binNsf(d,b))
            binLow = ppiclf_binb(2*d-1)
     >             + (ig(d)+ppiclf_binOffset(d))*ppiclf_bins_dx(d)
            fw(d) = FLOOR((ppiclf_picl_grid(d,ie)-binLow)/h)
            fw(d) = MAX(0, MIN(fw(d), ppiclf_binNsf(d,b)-1))
          ELSE IF(iside(d,m) .EQ. 1) THEN
            fw(d) = ppiclf_binNsf(d,b)-1
          ELSE
            fw(d) = 0
          END IF
        END DO
        place(m) = ppiclf_binSubOff(b) + fw(1)
     >           + ppiclf_binNsf(1,b)*fw(2)
     >           + ppiclf_binNsf(1,b)*ppiclf_binNsf(2,b)*fw(3)
      END DO

      RETURN
      END SUBROUTINE

      SUBROUTINE ppiclf_comm_binReachDilate
!
!     Convert the per-bin max cell edges (held in ppiclf_binReach as
!     scratch) into the per-bin search reach: 1.5x the largest cell
!     edge over the bin and its +-1 window neighbors, capped at the
!     bin size. The neighborhood max lets a particle near a bin face
!     reach the coarser cells of the adjacent bin; the cap keeps the
!     search interval within +-1 bin (guaranteed sufficient by the
!     Lfp bin-size constraint). In-place separable max filter.
!
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE
      INCLUDE "PPICLF"
      INTEGER*4 d, l, ii, jj, kk, m, mlo, mhi, b, src
      REAL*8, ALLOCATABLE :: rtmp(:,:)

      ALLOCATE(rtmp(3,0:ppiclf_total_SBin-1))
      DO kk = 0,ppiclf_nSBin(3)-1
       DO jj = 0,ppiclf_nSBin(2)-1
        DO ii = 0,ppiclf_nSBin(1)-1
          b = ii + ppiclf_nSBin(1)*jj
     >      + ppiclf_nSBin(1)*ppiclf_nSBin(2)*kk
          mlo = MAX(0, ii-1)
          mhi = MIN(ppiclf_nSBin(1)-1, ii+1)
          DO l = 1,3
            rtmp(l,b) = 0.0D0
            DO m = mlo,mhi
              src = m + ppiclf_nSBin(1)*jj
     >            + ppiclf_nSBin(1)*ppiclf_nSBin(2)*kk
              rtmp(l,b) = MAX(rtmp(l,b), ppiclf_binReach(l,src))
            END DO
          END DO
        END DO
       END DO
      END DO
      DO kk = 0,ppiclf_nSBin(3)-1
       DO jj = 0,ppiclf_nSBin(2)-1
        mlo = MAX(0, jj-1)
        mhi = MIN(ppiclf_nSBin(2)-1, jj+1)
        DO ii = 0,ppiclf_nSBin(1)-1
          b = ii + ppiclf_nSBin(1)*jj
     >      + ppiclf_nSBin(1)*ppiclf_nSBin(2)*kk
          DO l = 1,3
            ppiclf_binReach(l,b) = 0.0D0
            DO m = mlo,mhi
              src = ii + ppiclf_nSBin(1)*m
     >            + ppiclf_nSBin(1)*ppiclf_nSBin(2)*kk
              ppiclf_binReach(l,b) = MAX(ppiclf_binReach(l,b),
     >                                   rtmp(l,src))
            END DO
          END DO
        END DO
       END DO
      END DO
      DO kk = 0,ppiclf_nSBin(3)-1
       mlo = MAX(0, kk-1)
       mhi = MIN(ppiclf_nSBin(3)-1, kk+1)
       DO jj = 0,ppiclf_nSBin(2)-1
        DO ii = 0,ppiclf_nSBin(1)-1
          b = ii + ppiclf_nSBin(1)*jj
     >      + ppiclf_nSBin(1)*ppiclf_nSBin(2)*kk
          DO l = 1,3
            rtmp(l,b) = 0.0D0
            DO m = mlo,mhi
              src = ii + ppiclf_nSBin(1)*jj
     >            + ppiclf_nSBin(1)*ppiclf_nSBin(2)*m
              rtmp(l,b) = MAX(rtmp(l,b), ppiclf_binReach(l,src))
            END DO
          END DO
        END DO
       END DO
      END DO
      DO b = 0,ppiclf_total_SBin-1
        DO l = 1,3
          ppiclf_binReach(l,b) = MIN(1.5D0*rtmp(l,b),
     >                               ppiclf_bins_dx(l))
          ! Bins with no cells anywhere nearby: harmless default
          IF(ppiclf_binReach(l,b) .LE. 0.0D0)
     >       ppiclf_binReach(l,b) = ppiclf_bins_dx(l)
        END DO
      END DO
      DEALLOCATE(rtmp)

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
      INTEGER*4  ie, i, j, k, d, iTemp_SBin(3)
     >          ,tempSBin, i_SBin(3), icount, newMax
     >          ,nplace, place(8), iside(3,8), pos
      INTEGER*8  nfb8, nbytes8
      ! 250 MBytes
      INTEGER*8, PARAMETER :: FLUIDBIN_BYTE_BUDGET = 250000000_8

      ! Build-skip optimization only: the fluid-cell maps depend on
      ! geometry, not particles, so an empty rank may skip the
      ! build. The stamp is deliberately NOT touched here - it
      ! must keep describing the geometry of the last FULL build,
      ! so that when particles later arrive under a newer geometry
      ! the state check in PostTimeStepPartLB forces the rebuild.
      ! (Stamping here falsely marked stale arrays current and
      ! caused an out-of-bounds binReach read; the reader-side
      ! check is gated on npart>0 to avoid rebuild storms from
      ! perpetually empty ranks.)
      IF(ppiclf_npart .LT. 1) RETURN

      IF(ppiclf_nCells_FV2PICL .LE. 0 .AND. ppiclf_npart .GT. 0) THEN
        PRINT*,'ERROR: ',ppiclf_npart, 'Particles mapped to bin:'
     >         ,ppiclf_nid
        PRINT*,'No cells mapped to bin for Interpolation/Projection.'
        CALL ppiclf_exittr('Failure in particle to cell mapping$',
     >                      0.D0,0)
      END IF

      ! ---- Decide whether to refine the fluid cell map ----------------
      ! PER-BIN sizing: each window bin subdivides by its OWN coarsest
      ! cell, nsf_d(b) = max(1, floor(bins_dx_d/(1.5*binMax_d(b)))),
      ! so particles among fine cells search a fine grid while bins
      ! containing a refinement boundary keep the coarse resolution
      ! their particles genuinely need. This bounds the candidate
      ! count by LOCAL mesh contrast instead of the rank-wide
      ! (partition-dependent) contrast of the former uniform grid.
      ! The fine path is enabled when any bin refines meaningfully
      ! (sub-cell < 0.5*filter), matching the LB model's path check.
      CALL ppiclf_allocate_BinSubGrid(ppiclf_total_SBin)
      ! Per-bin max cell edges over direct + periodic-image placements
      ppiclf_binReach = 0.0D0   ! reuse as scratch for binMax first
      DO ie = 1,ppiclf_nCells_FV2PICL
        CALL ppiclf_comm_cellBinImages(ie, place, iside, nplace)
        DO i = 1,nplace
          DO d = 1,3
            IF(ppiclf_picl_grid(3+d,ie) .GT.
     >         ppiclf_binReach(d,place(i)))
     >         ppiclf_binReach(d,place(i)) = ppiclf_picl_grid(3+d,ie)
          END DO
        END DO
      END DO
      ppiclf_useFineFluid = .FALSE.
      DO i = 0,ppiclf_total_SBin-1
        DO d = 1,3
          IF(ppiclf_binReach(d,i) .GT. 0.0D0) THEN
            ppiclf_binNsf(d,i) = MAX(1, FLOOR(ppiclf_bins_dx(d)
     >                          /(1.5D0*ppiclf_binReach(d,i))))
          ELSE
            ppiclf_binNsf(d,i) = 1
          END IF
          IF(ppiclf_bins_dx(d)/DBLE(ppiclf_binNsf(d,i)) .LT.
     >       0.5D0*ppiclf_filter(d)) ppiclf_useFineFluid = .TRUE.
        END DO
      END DO

      IF(ppiclf_useFineFluid) THEN
        ! Memory guard: shrink the most-refined bins until the packed
        ! grid fits the byte budget (or refinement is exhausted).
        DO
          nfb8 = 0_8
          DO i = 0,ppiclf_total_SBin-1
            nfb8 = nfb8 + INT(ppiclf_binNsf(1,i),8)
     >                   *INT(ppiclf_binNsf(2,i),8)
     >                   *INT(ppiclf_binNsf(3,i),8)
          END DO
          nbytes8 = (2_8*nfb8 + 1_8
     >             + 8_8*INT(ppiclf_nCells_FV2PICL,8))*4_8
          IF(nfb8    .LE. INT(HUGE(ppiclf_total_fluidSBin),8) .AND.
     >       nbytes8 .LE. FLUIDBIN_BYTE_BUDGET) EXIT
          newMax = 1
          DO i = 0,ppiclf_total_SBin-1
            DO d = 1,3
              IF(ppiclf_binNsf(d,i) .GT. 1) THEN
                ppiclf_binNsf(d,i) = MAX(1, ppiclf_binNsf(d,i)/2)
                newMax = MAX(newMax, ppiclf_binNsf(d,i))
              END IF
            END DO
          END DO
          IF(newMax .EQ. 1) THEN
            ppiclf_useFineFluid = .FALSE.
            EXIT
          END IF
        END DO
      END IF

      IF(ppiclf_useFineFluid .AND. ppiclf_nid .EQ. 0) THEN
        newMax = 0
        DO i = 0,ppiclf_total_SBin-1
          newMax = MAX(newMax, ppiclf_binNsf(1,i)*ppiclf_binNsf(2,i)
     >                         *ppiclf_binNsf(3,i))
        END DO
        PRINT*, 'FineFluid per-bin map: packed subcells=',INT(nfb8,4),
     >          ' max nsf product=', newMax
      END IF

      IF(ppiclf_useFineFluid) THEN
        ! Per-bin search reach: 1.5x the largest cell edge over the bin
        ! and its +-1 window neighbors, so near-face particles reach
        ! coarser neighbor cells; capped at the bin size (which the
        ! Lfp bin constraint guarantees is sufficient).
        CALL ppiclf_comm_binReachDilate
        ! Packed sub-cell offsets per bin
        ppiclf_binSubOff(0) = 0
        DO i = 0,ppiclf_total_SBin-1
          ppiclf_binSubOff(i+1) = ppiclf_binSubOff(i)
     >        + ppiclf_binNsf(1,i)*ppiclf_binNsf(2,i)
     >          *ppiclf_binNsf(3,i)
        END DO
      END IF

      ! ================= COARSE (filter-sized) cell map ================
      IF(.NOT. ppiclf_useFineFluid) THEN
        IF(ppiclf_maxCellsPerBin .LT. 1) ppiclf_maxCellsPerBin = 92
        CALL ppiclf_allocate_BTC(ppiclf_total_SBin,
     >                           ppiclf_maxCellsPerBin)
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
     >           iTemp_SBin(1) .GT. ppiclf_nSBin(1) - 1) CYCLE 
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
     >             iTemp_SBin(2) .GT. ppiclf_nSBin(2) - 1) CYCLE
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
     >               iTemp_SBin(3) .GT. ppiclf_nSBin(3) - 1) CYCLE
                ELSE ! k .EQ. 1
                  ! This takes care of periodicity for single processor
                  IF(ppiclf_linperiodic(3) .AND. 
     >                                      ppiclf_EqualDomain(3)) THEN 
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
     >                   + ppiclf_nSBin(1)*ppiclf_nSBin(2)*iTemp_SBin(3)
                IF(tempSBin .LT. 0 .OR. 
     >             tempSBin .GT. ppiclf_total_SBin-1) THEN
                 PRINT*,'ERROR:Bad Subbin Index in Overlap Cell Mapping'
     >                   , tempSBin
                  CALL ppiclf_exittr('',0.0D0,0)
                END IF

                icount = ppiclf_binCellCount(tempSBin) + 1
                ppiclf_binCellCount(tempSBin) = icount

                ! Grow on demand. icount rises by 1 per hit, so a single
                ! doubling always restores capacity (newMax >= icount).
                IF(icount .GT. ppiclf_maxCellsPerBin) THEN
                  newMax = 2*ppiclf_maxCellsPerBin + 1
                  CALL ppiclf_reallocate_BTC(ppiclf_total_SBin,
     >                                       ppiclf_maxCellsPerBin,
     >                                       newMax)
                  ppiclf_maxCellsPerBin = newMax
                END IF
                ppiclf_binCellList(tempSBin, icount) = ie
              END DO !k
            END DO !j 
          END DO !i
        END DO !ie

        ! Record the geometry these maps were built for
        CALL ppiclf_setBTCStamp(ppiclf_nSBin(1),ppiclf_nSBin(2),
     >       ppiclf_nSBin(3),ppiclf_binOffset(1),
     >       ppiclf_binOffset(2),ppiclf_binOffset(3),
     >       ppiclf_total_SBin)
        RETURN
      END IF

      ! ================= FINE (CSR) fluid cell map =====================
      ! Rows are the PACKED per-bin sub-cells: bin b occupies rows
      ! [binSubOff(b), binSubOff(b+1)), local ordering fast-x.
      ppiclf_total_fluidSBin = INT(nfb8,4)
      CALL ppiclf_allocate_FluidCSR(ppiclf_total_fluidSBin,
     >                              8*ppiclf_nCells_FV2PICL)

      ! Pass 1: count placements per fine sub-bin (incl. periodic images)
      DO i = 0,ppiclf_total_fluidSBin-1
        ppiclf_fluidCellCount(i) = 0
      END DO
      DO ie = 1,ppiclf_nCells_FV2PICL
        CALL ppiclf_comm_fluidCellSubImages(ie, place, nplace)
        DO i = 1,nplace
          ppiclf_fluidCellCount(place(i)) =
     >                  ppiclf_fluidCellCount(place(i)) + 1
        END DO
      END DO

      ! Prefix sum -> row pointers, then seed cursors
      ppiclf_fluidCellOffset(0) = 1
      DO i = 0,ppiclf_total_fluidSBin-1
        ppiclf_fluidCellOffset(i+1) = ppiclf_fluidCellOffset(i)
     >                              + ppiclf_fluidCellCount(i)
      END DO
      DO i = 0,ppiclf_total_fluidSBin-1
        ppiclf_fluidCellCount(i) = ppiclf_fluidCellOffset(i)  ! cursor
      END DO

      ! Pass 2: scatter cell IDs into the flat list
      DO ie = 1,ppiclf_nCells_FV2PICL
        CALL ppiclf_comm_fluidCellSubImages(ie, place, nplace)
        DO i = 1,nplace
          pos = ppiclf_fluidCellCount(place(i))
          ppiclf_fluidCellFlat(pos)       = ie
          ppiclf_fluidCellCount(place(i)) = pos + 1
        END DO
      END DO

      ! Record the geometry these maps were built for
        CALL ppiclf_setBTCStamp(ppiclf_nSBin(1),ppiclf_nSBin(2),
     >       ppiclf_nSBin(3),ppiclf_binOffset(1),
     >       ppiclf_binOffset(2),ppiclf_binOffset(3),
     >       ppiclf_total_SBin)
      RETURN
      END SUBROUTINE
!-----------------------------------------------------------------------
! This is for statistics gathering only
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
