MODULE ppiclf_DynamicAllocation

  IMPLICIT NONE

  INTEGER*4 :: ppiclf_dL, ppiclf_dM, ppiclf_dS
  INTEGER*4 :: ppiclf_total_SBin,ppiclf_nSBin(3)
  INTEGER*4 :: ppiclf_maxParticlePerBin
  INTEGER*4 :: ppiclf_binOffset(3)
  INTEGER*4, ALLOCATABLE :: ppiclf_ParticleCount(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_CellCount(:)
  REAL*8, ALLOCATABLE :: ppiclf_CellMaxLen(:,:)
  ! Cached per-bin P2C stencil ratio sMAP (see BinMapStencil). Its
  ! inputs change only when FindCellPartLB reruns, so consumers on the
  ! per-stage path (weighted scan, calibration) use this table instead
  ! of re-deriving it for every bin every stage.
  REAL*8, ALLOCATABLE :: ppiclf_binSmap(:)
  LOGICAL :: ppiclf_CellCountValid = .FALSE.
  REAL*8  :: ppiclf_lbBaseMax     = -1.0D0  
  LOGICAL :: ppiclf_lbJustBal     = .FALSE. 
  INTEGER*4, ALLOCATABLE :: ppiclf_BinToRankMap(:)
  LOGICAL,   ALLOCATABLE :: ppiclf_LRankBoundary(:,:)
  LOGICAL,   ALLOCATABLE :: ppiclf_LMapFluid(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_binCellCount(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_binCellList(:,:)
  INTEGER*4 :: ppiclf_maxCellsPerBin = 0
  INTEGER*4, ALLOCATABLE :: ppiclf_binPartCount(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_binPartList(:,:)
  LOGICAL   :: ppiclf_useFineGrid = .FALSE.
  INTEGER*4 :: ppiclf_nFine(3), ppiclf_total_fineBin
  REAL*8    :: ppiclf_fineLo(3), ppiclf_fineInvLen(3)
  INTEGER*4, ALLOCATABLE :: ppiclf_finePartCount(:)  
  INTEGER*4, ALLOCATABLE :: ppiclf_fineOffset(:)     
  INTEGER*4, ALLOCATABLE :: ppiclf_fineFlat(:)       
  LOGICAL   :: ppiclf_useFineFluid = .FALSE.
  INTEGER*4 :: ppiclf_nSubFluid(3) = 1
  INTEGER*4 :: ppiclf_total_fluidSBin = 1
  INTEGER*4, ALLOCATABLE :: ppiclf_binNsf(:,:)    ! (3,0:nwin-1)
  INTEGER*4, ALLOCATABLE :: ppiclf_binSubOff(:)   ! (0:nwin)
  REAL*8,    ALLOCATABLE :: ppiclf_binReach(:,:)  ! (3,0:nwin-1)
  INTEGER*4, ALLOCATABLE :: ppiclf_fluidCellCount(:)  ! scratch: count/cursor
  INTEGER*4, ALLOCATABLE :: ppiclf_fluidCellOffset(:) ! (0:nbin) row pointers
  INTEGER*4, ALLOCATABLE :: ppiclf_fluidCellFlat(:)   ! (1:cap) flat list

  ! ===================================================================
  ! Online load-balance coefficient calibration (see
  ! ppiclf_comm_LBCalibAccum / ppiclf_comm_LBCalibrate in
  ! ppiclf_PBloading.f). ppiclf_LB_C holds the EFFECTIVE per-candidate
  ! costs (C_L, C_PP, C_MAP, C_PROJ, C_CELL) used by
  ! ppiclf_comm_LBWeightCoef. They are initialized to the analytic
  ! flop-count priors and, in PERF builds, periodically refit against
  ! the measured per-rank kernel timers. Only ratios matter to the
  ! partition. The five channel timers are unconditional source code,
  ! so calibration is active in EVERY build (PERF is only needed for
  ! the full instrumentation and CSV logging). Set
  ! ppiclf_LB_docal = .FALSE. to freeze the priors at run time.
  ! ===================================================================
  ! C_PROJ (4th entry) is on the per-PARTICLE basis (capped projection
  ! stencil); its prior is the measured projection/linear rate ratio
  ! (~0.31) times the C_L prior. See BinWeight for the model note.
  REAL*8    :: ppiclf_LB_prior(5) = &
               (/3.5D3, 1.2D1, 2.5D1, 1.1D3, 1.5D3/)
  REAL*8    :: ppiclf_LB_C(5) = &
               (/3.5D3, 1.2D1, 2.5D1, 1.1D3, 1.5D3/)
  LOGICAL   :: ppiclf_LB_docal    = .TRUE.  ! master switch (PERF only)
  INTEGER*4 :: ppiclf_LB_calfreq  = 100     ! steps per refit window
  INTEGER*4 :: ppiclf_LB_calsteps = 0       ! samples in current window
  REAL*8    :: ppiclf_LB_X(5)  = 0.0D0      ! basis sums, this rank
  REAL*8    :: ppiclf_LB_Y(5)  = 0.0D0      ! channel seconds, this rank
  REAL*8    :: ppiclf_LB_T0(5) = 0.0D0      ! last timer-channel sample
  LOGICAL   :: ppiclf_LB_snapinit = .FALSE. ! (deprecated; kept for ABI)
  ! Profiling-only entry barrier (see ppiclf_solve_IntegrateParticle):
  ! when .TRUE. in a PERF build, an MPI_BARRIER is taken at ppiclF
  ! entry and its wait time recorded as TEntrySync, separating host
  ! (fluid-solver) imbalance from ppiclF's own timers. Default off.
  LOGICAL   :: ppiclf_perf_sync = .FALSE.
  ! Cadence of the global-count reduction + rebalance criterion in
  ! FindParticlePartLB (stages between checks; 3 = once per RK3 step,
  ! 1 = legacy per-stage behavior). Per-particle bin assignment and
  ! migration remain per-stage regardless. ppiclf_LB_countsfresh is
  ! TRUE only during stages on which the counts were just reduced; the
  ! EIB refresh and the calibration basis accumulation key off it.
  INTEGER*4 :: ppiclf_LB_checkfreq  = 3
  INTEGER*4 :: ppiclf_LB_stagectr   = 0
  LOGICAL   :: ppiclf_LB_countsfresh = .FALSE.
  ! Walltime-adaptive refit window. The refit itself is one 10-word
  ! reduction (microseconds), so its cost is negligible at any
  ! cadence; what matters is fit quality and coefficient freshness.
  ! Per-sample signal scales with step cost, so a fixed-WALLTIME
  ! window gives roughly constant fit quality across problem sizes:
  ! large (expensive) problems refit every few steps, small (cheap)
  ! problems accumulate many samples per refit. After each refit,
  ! rank 0 rescales ppiclf_LB_calfreq so the next window targets
  ! ppiclf_LB_caltime seconds, clamped to [calmin, calmax] samples,
  ! and broadcasts it (collectively consistent). Set caltime <= 0 to
  ! disable adaptation and keep a fixed calfreq.
  REAL*8    :: ppiclf_LB_caltime = 3.0D1
  INTEGER*4 :: ppiclf_LB_calmin  = 10
  INTEGER*4 :: ppiclf_LB_calmax  = 1000
  REAL*8    :: ppiclf_LB_lastfit = -1.0D0
  ! Warmup: the very first sample (startup stage: allocations,
  ! first-touch paging) is discarded, then a SHORT first window of
  ! ppiclf_LB_calfirst samples puts the first refit within the first
  ! few time steps -- per-sample signal scales with step cost, so a
  ! handful of samples is already a usable fit for exactly the large
  ! problems where early coefficients matter, and the ridge blend
  ! damps the noise. The first refit is taken at full weight (no EWMA
  ! with the prior-initialized values), and the walltime targeting
  ! rescales the window from the measured first-window duration
  ! onward. calfirst applies only in adaptive mode (caltime > 0).
  INTEGER*4 :: ppiclf_LB_calfirst = 3
  LOGICAL   :: ppiclf_LB_started  = .FALSE.
  INTEGER*4 :: ppiclf_LB_nrefits  = 0

  ! Build stamp for the particle->cell map arrays: subbin geometry
  ! (nSBin, binOffset, total_SBin) the maps were last built for.
  ! Lets readers detect ANY missed rebuild trigger (state-based
  ! check, in addition to the event flags).
  INTEGER*4, SAVE :: ppiclf_BTC_stamp(8) = -1
  ! Monotone build counter for the fluid-overlap cell enumeration.
  ! Bumped every time MapOverlapGrid(PartLB) completes; captured in
  ! stamp slot 8 so the CSR is keyed to the EXACT cell-set build it
  ! was made from. This catches bin-box translations that change the
  ! overlap cells while leaving all window indices numerically equal.
  INTEGER*4, SAVE :: ppiclf_mapEpoch = 0

CONTAINS

  ! =====================================================================
  ! Global Bin Data Allocation
  ! =====================================================================
  SUBROUTINE ppiclf_dyn_alloc(nbins_in, np_in)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nbins_in, np_in
    INTEGER*4 current_size

    IF(ALLOCATED(ppiclf_ParticleCount)) THEN
      current_size = SIZE(ppiclf_ParticleCount)
      IF(current_size .NE. nbins_in) THEN
        DEALLOCATE(ppiclf_ParticleCount)
        ALLOCATE(ppiclf_ParticleCount(0:nbins_in-1))
      END IF
    ELSE
      ALLOCATE(ppiclf_ParticleCount(0:nbins_in-1))
    END IF

    IF(ALLOCATED(ppiclf_CellCount)) THEN
      IF(SIZE(ppiclf_CellCount) .NE. nbins_in) THEN
        DEALLOCATE(ppiclf_CellCount)
        ALLOCATE(ppiclf_CellCount(0:nbins_in-1))
        ppiclf_CellCount = 0
        ppiclf_CellCountValid = .FALSE.
      END IF
    ELSE
      ALLOCATE(ppiclf_CellCount(0:nbins_in-1))
      ppiclf_CellCount = 0
      ppiclf_CellCountValid = .FALSE.
    END IF

    IF(ALLOCATED(ppiclf_CellMaxLen)) THEN
      IF(SIZE(ppiclf_CellMaxLen,2) .NE. nbins_in) THEN
        DEALLOCATE(ppiclf_CellMaxLen)
        ALLOCATE(ppiclf_CellMaxLen(3,0:nbins_in-1))
        ppiclf_CellMaxLen = 0.0D0
        ppiclf_CellCountValid = .FALSE.
      END IF
    ELSE
      ALLOCATE(ppiclf_CellMaxLen(3,0:nbins_in-1))
      ppiclf_CellMaxLen = 0.0D0
      ppiclf_CellCountValid = .FALSE.
    END IF

    IF(ALLOCATED(ppiclf_binSmap)) THEN
      IF(SIZE(ppiclf_binSmap) .NE. nbins_in) THEN
        DEALLOCATE(ppiclf_binSmap)
        ALLOCATE(ppiclf_binSmap(0:nbins_in-1))
        ppiclf_binSmap = 2.7D1
      END IF
    ELSE
      ALLOCATE(ppiclf_binSmap(0:nbins_in-1))
      ppiclf_binSmap = 2.7D1
    END IF

    IF(ALLOCATED(ppiclf_BinToRankMap)) THEN
      IF(SIZE(ppiclf_BinToRankMap) .NE. nbins_in) THEN
        DEALLOCATE(ppiclf_BinToRankMap)
        ALLOCATE(ppiclf_BinToRankMap(0:nbins_in-1))
        ppiclf_BinToRankMap = -1
      END IF
    ELSE
      ALLOCATE(ppiclf_BinToRankMap(0:nbins_in-1))
      ppiclf_BinToRankMap = -1
    END IF
    
    IF(ALLOCATED(ppiclf_LRankBoundary)) THEN
      IF(SIZE(ppiclf_LRankBoundary, DIM=1) .NE. nbins_in) THEN
        DEALLOCATE(ppiclf_LRankBoundary)
        ALLOCATE(ppiclf_LRankBoundary(0:nbins_in-1, 6))
      END IF
    ELSE
      ALLOCATE(ppiclf_LRankBoundary(0:nbins_in-1, 6))
    END IF

    IF(ALLOCATED(ppiclf_LMapFluid)) THEN
      IF(SIZE(ppiclf_LMapFluid) .NE. nbins_in) THEN
        DEALLOCATE(ppiclf_LMapFluid)
        ALLOCATE(ppiclf_LMapFluid(0:nbins_in-1))
      END IF
    ELSE
      ALLOCATE(ppiclf_LMapFluid(0:nbins_in-1))
    END IF
    
  END SUBROUTINE

  ! =====================================================================
  ! Sub-Bin Particle Mapping Allocation
  ! =====================================================================
  SUBROUTINE ppiclf_allocate_BTP(nbins_in, maxPartPerBin)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nbins_in, maxPartPerBin
    INTEGER*4 size_dim1, size_dim2

    IF (ALLOCATED(ppiclf_binPartCount)) THEN
        IF (SIZE(ppiclf_binPartCount) .NE. nbins_in) THEN
            DEALLOCATE(ppiclf_binPartCount)
            ALLOCATE(ppiclf_binPartCount(0:nbins_in-1))
        END IF
    ELSE
        ALLOCATE(ppiclf_binPartCount(0:nbins_in-1))
    END IF

    IF (ALLOCATED(ppiclf_binPartList)) THEN
      size_dim1 = SIZE(ppiclf_binPartList, 1) 
      size_dim2 = SIZE(ppiclf_binPartList, 2)
      
      ! Reallocate ONLY if dimensions have actually grown
      IF (size_dim1 .NE. nbins_in .OR. size_dim2 .LT. maxPartPerBin) THEN
        DEALLOCATE(ppiclf_binPartList)
        ALLOCATE(ppiclf_binPartList(0:nbins_in-1, maxPartPerBin))
      END IF
    ELSE
      ALLOCATE(ppiclf_binPartList(0:nbins_in-1, maxPartPerBin))
    END IF
  END SUBROUTINE

  ! =====================================================================
  ! Sub-Bin Particle Mapping RE-Allocation (For Ghost Injection)
  ! =====================================================================
  SUBROUTINE ppiclf_reallocate_BTP(nbins_in, currentMax, newMax)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nbins_in, currentMax, newMax
    INTEGER*4, ALLOCATABLE :: temp_list(:,:)

    ALLOCATE(temp_list(0:nbins_in-1, newMax))
    temp_list = 0
    
    ! Copy the old data into the new, larger array
    temp_list(:, 1:currentMax) = ppiclf_binPartList(:,1:currentMax)

    ! Move pointers. 'temp_list' is automatically deallocated.
    CALL MOVE_ALLOC(FROM=temp_list, TO=ppiclf_binPartList)
    
  END SUBROUTINE

  ! =====================================================================
  ! Fine Sub-Bin Particle Mapping Allocation (P2P search, by nndist)
  ! CSR layout: row-pointer array (0:nbins), per-bin scratch (0:nbins-1),
  ! and a flat candidate list sized to the actual entry count. The flat
  ! list grows on demand and never shrinks.
  ! =====================================================================
  SUBROUTINE ppiclf_allocate_FineCSR(nbins_in, nentries_in)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nbins_in, nentries_in

    ! Row pointers: one extra slot so fineOffset(nbins_in) is valid.
    IF (ALLOCATED(ppiclf_fineOffset)) THEN
        IF (SIZE(ppiclf_fineOffset) .NE. nbins_in+1) THEN
            DEALLOCATE(ppiclf_fineOffset)
            ALLOCATE(ppiclf_fineOffset(0:nbins_in))
        END IF
    ELSE
        ALLOCATE(ppiclf_fineOffset(0:nbins_in))
    END IF

    ! Per-bin scratch (count, then reused as the write cursor).
    IF (ALLOCATED(ppiclf_finePartCount)) THEN
        IF (SIZE(ppiclf_finePartCount) .NE. nbins_in) THEN
            DEALLOCATE(ppiclf_finePartCount)
            ALLOCATE(ppiclf_finePartCount(0:nbins_in-1))
        END IF
    ELSE
        ALLOCATE(ppiclf_finePartCount(0:nbins_in-1))
    END IF

    ! Flat candidate list: grow-on-demand, never shrink.
    IF (ALLOCATED(ppiclf_fineFlat)) THEN
        IF (SIZE(ppiclf_fineFlat) .LT. nentries_in) THEN
            DEALLOCATE(ppiclf_fineFlat)
            ALLOCATE(ppiclf_fineFlat(MAX(nentries_in,1)))
        END IF
    ELSE
        ALLOCATE(ppiclf_fineFlat(MAX(nentries_in,1)))
    END IF
  END SUBROUTINE

  ! =====================================================================
  ! Fine FLUID Cell Mapping Allocation (CSR). Mirrors FineCSR: row
  ! pointers (0:nbins), per-bin scratch (0:nbins-1), and a flat cell-ID
  ! list sized to the actual placement count (grows, never shrinks).
  ! =====================================================================
  SUBROUTINE ppiclf_allocate_BinSubGrid(nwin_in)
    INTEGER*4, INTENT(IN) :: nwin_in
    IF (ALLOCATED(ppiclf_binSubOff)) THEN
      IF (SIZE(ppiclf_binSubOff) .NE. nwin_in+1) THEN
        DEALLOCATE(ppiclf_binNsf, ppiclf_binSubOff, ppiclf_binReach)
        ALLOCATE(ppiclf_binNsf(3,0:nwin_in-1))
        ALLOCATE(ppiclf_binSubOff(0:nwin_in))
        ALLOCATE(ppiclf_binReach(3,0:nwin_in-1))
      END IF
    ELSE
      ALLOCATE(ppiclf_binNsf(3,0:nwin_in-1))
      ALLOCATE(ppiclf_binSubOff(0:nwin_in))
      ALLOCATE(ppiclf_binReach(3,0:nwin_in-1))
    END IF
  END SUBROUTINE

  SUBROUTINE ppiclf_allocate_FluidCSR(nbins_in, nentries_in)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nbins_in, nentries_in

    IF (ALLOCATED(ppiclf_fluidCellOffset)) THEN
        IF (SIZE(ppiclf_fluidCellOffset) .NE. nbins_in+1) THEN
            DEALLOCATE(ppiclf_fluidCellOffset)
            ALLOCATE(ppiclf_fluidCellOffset(0:nbins_in))
        END IF
    ELSE
        ALLOCATE(ppiclf_fluidCellOffset(0:nbins_in))
    END IF

    IF (ALLOCATED(ppiclf_fluidCellCount)) THEN
        IF (SIZE(ppiclf_fluidCellCount) .NE. nbins_in) THEN
            DEALLOCATE(ppiclf_fluidCellCount)
            ALLOCATE(ppiclf_fluidCellCount(0:nbins_in-1))
        END IF
    ELSE
        ALLOCATE(ppiclf_fluidCellCount(0:nbins_in-1))
    END IF

    IF (ALLOCATED(ppiclf_fluidCellFlat)) THEN
        IF (SIZE(ppiclf_fluidCellFlat) .LT. nentries_in) THEN
            DEALLOCATE(ppiclf_fluidCellFlat)
            ALLOCATE(ppiclf_fluidCellFlat(MAX(nentries_in,1)))
        END IF
    ELSE
        ALLOCATE(ppiclf_fluidCellFlat(MAX(nentries_in,1)))
    END IF
  END SUBROUTINE

  ! =====================================================================
  ! Sub-Bin Cell Mapping Allocation
  ! =====================================================================
  SUBROUTINE ppiclf_allocate_BTC(nbins_in, maxCellsPerBin)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nbins_in, maxCellsPerBin
    INTEGER*4 size_dim1, size_dim2

    IF (ALLOCATED(ppiclf_binCellCount)) THEN
        IF (SIZE(ppiclf_binCellCount) .NE. nbins_in) THEN
            DEALLOCATE(ppiclf_binCellCount)
            ALLOCATE(ppiclf_binCellCount(0:nbins_in-1))
        END IF
    ELSE
        ALLOCATE(ppiclf_binCellCount(0:nbins_in-1))
    END IF

    IF (ALLOCATED(ppiclf_binCellList)) THEN
      size_dim1 = SIZE(ppiclf_binCellList, 1) 
      size_dim2 = SIZE(ppiclf_binCellList, 2)
      
      IF (size_dim1 .NE. nbins_in .OR. size_dim2 .LT. maxCellsPerBin) THEN
        DEALLOCATE(ppiclf_binCellList)
        ALLOCATE(ppiclf_binCellList(0:nbins_in-1, maxCellsPerBin))
      END IF
    ELSE
      ALLOCATE(ppiclf_binCellList(0:nbins_in-1, maxCellsPerBin))
    END IF
  END SUBROUTINE

  ! =====================================================================
  ! Sub-Bin Cell Mapping RE-Allocation (grow-on-demand, content safe)
  ! =====================================================================
  SUBROUTINE ppiclf_reallocate_BTC(nbins_in, currentMax, newMax)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nbins_in, currentMax, newMax
    INTEGER*4, ALLOCATABLE :: temp_list(:,:)

    ALLOCATE(temp_list(0:nbins_in-1, newMax))
    ! Preserve everything written so far; new columns are filled by the
    ! caller before they are ever read (count-controlled).
    temp_list(:, 1:currentMax) = ppiclf_binCellList(:, 1:currentMax)
    CALL MOVE_ALLOC(FROM=temp_list, TO=ppiclf_binCellList)
  END SUBROUTINE

  ! =====================================================================
  ! Particle->cell map build stamp (set by subbinCellMap, checked in
  ! PostTimeStepPartLB before SBParticleToCellMap).
  ! =====================================================================
  SUBROUTINE ppiclf_setBTCStamp(n1,n2,n3,o1,o2,o3,tot)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: n1,n2,n3,o1,o2,o3,tot
    ppiclf_BTC_stamp(1) = n1
    ppiclf_BTC_stamp(2) = n2
    ppiclf_BTC_stamp(3) = n3
    ppiclf_BTC_stamp(4) = o1
    ppiclf_BTC_stamp(5) = o2
    ppiclf_BTC_stamp(6) = o3
    ppiclf_BTC_stamp(7) = tot
    ppiclf_BTC_stamp(8) = ppiclf_mapEpoch
  END SUBROUTINE

  SUBROUTINE ppiclf_bumpMapEpoch
    IMPLICIT NONE
    ppiclf_mapEpoch = ppiclf_mapEpoch + 1
  END SUBROUTINE

  LOGICAL FUNCTION ppiclf_BTCStale(n1,n2,n3,o1,o2,o3,tot)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: n1,n2,n3,o1,o2,o3,tot
    ppiclf_BTCStale = (ppiclf_BTC_stamp(1) .NE. n1)  .OR. &
                      (ppiclf_BTC_stamp(2) .NE. n2)  .OR. &
                      (ppiclf_BTC_stamp(3) .NE. n3)  .OR. &
                      (ppiclf_BTC_stamp(4) .NE. o1)  .OR. &
                      (ppiclf_BTC_stamp(5) .NE. o2)  .OR. &
                      (ppiclf_BTC_stamp(6) .NE. o3)  .OR. &
                      (ppiclf_BTC_stamp(7) .NE. tot)  .OR. &
                      (ppiclf_BTC_stamp(8) .NE. ppiclf_mapEpoch)
  END FUNCTION

END MODULE ppiclf_DynamicAllocation
