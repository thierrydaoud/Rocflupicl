MODULE ppiclf_DynamicAllocation

  IMPLICIT NONE

  INTEGER*4 :: ppiclf_dL, ppiclf_dM, ppiclf_dS
  INTEGER*4 :: ppiclf_total_SBin,ppiclf_nSBin(3)
  INTEGER*4 :: ppiclf_maxParticlePerBin
  INTEGER*4 :: ppiclf_binOffset(3)

  INTEGER*4, ALLOCATABLE :: ppiclf_ParticleCount(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_BinToRankMap(:)
  LOGICAL,   ALLOCATABLE :: ppiclf_LRankBoundary(:,:)
  LOGICAL,   ALLOCATABLE :: ppiclf_LMapFluid(:)
  
  INTEGER*4, ALLOCATABLE :: ppiclf_binCellCount(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_binCellList(:,:)
  ! Persisted per-bin capacity for the overlap-cell list. Carried across
  ! calls so the steady state needs no reallocation; only ever grows.
  INTEGER*4 :: ppiclf_maxCellsPerBin = 0
  
  INTEGER*4, ALLOCATABLE :: ppiclf_binPartCount(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_binPartList(:,:)

  ! --- Fine sub-bin grid (sized by ppiclf_nndist) for the P2P search ---
  ! Separate from the coarse iprop sub-bins (sized by ppiclf_filter).
  ! ppiclf_useFineGrid is set by ppiclf_comm_subbinFineParticleMap; when
  ! .FALSE. the nearest-neighbor search should fall back to the coarse
  ! ppiclf_binPartCount / ppiclf_binPartList arrays.
  !
  ! CSR (compressed sparse) layout: ppiclf_fineOffset(0:total_fineBin)
  ! holds 1-based row pointers into the flat candidate list
  ! ppiclf_fineFlat. Fine bin b owns
  !   ppiclf_fineFlat( fineOffset(b) : fineOffset(b+1)-1 ).
  ! ppiclf_finePartCount is build-time scratch (per-bin count, then
  ! reused as the per-bin write cursor). Memory now scales with actual
  ! occupancy (npart+ng) rather than total_fineBin*maxPartPerFineBin.
  LOGICAL   :: ppiclf_useFineGrid = .FALSE.
  INTEGER*4 :: ppiclf_nFine(3), ppiclf_total_fineBin
  REAL*8    :: ppiclf_fineLo(3), ppiclf_fineInvLen(3)
  INTEGER*4, ALLOCATABLE :: ppiclf_finePartCount(:)  ! scratch: count/cursor
  INTEGER*4, ALLOCATABLE :: ppiclf_fineOffset(:)     ! (0:nbin) row pointers
  INTEGER*4, ALLOCATABLE :: ppiclf_fineFlat(:)       ! (1:cap) flat list

  ! --- Optional FINE FLUID cell sub-bin grid (CSR) -----------------------
  ! When the coarse (filter-sized) sub-bins hold too many cells (cells are
  ! small relative to ppiclf_filter), subdivide each coarse sub-bin by
  ! ppiclf_nSubFluid(d) per dimension and map overlap cells into the finer
  ! grid by POSITION. Enabled per-rank by ppiclf_useFineFluid. Same CSR
  ! layout as the P2P fine grid: bin b owns
  !   ppiclf_fluidCellFlat( fluidCellOffset(b) : fluidCellOffset(b+1)-1 ).
  ! Consumed by ppiclf_solve_SBParticleToCellMap.
  LOGICAL   :: ppiclf_useFineFluid = .FALSE.
  INTEGER*4 :: ppiclf_nSubFluid(3) = 1
  INTEGER*4 :: ppiclf_total_fluidSBin = 1
  INTEGER*4, ALLOCATABLE :: ppiclf_fluidCellCount(:)  ! scratch: count/cursor
  INTEGER*4, ALLOCATABLE :: ppiclf_fluidCellOffset(:) ! (0:nbin) row pointers
  INTEGER*4, ALLOCATABLE :: ppiclf_fluidCellFlat(:)   ! (1:cap) flat list

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

    IF(ALLOCATED(ppiclf_BinToRankMap)) THEN
      IF(SIZE(ppiclf_BinToRankMap) .NE. nbins_in) THEN
        DEALLOCATE(ppiclf_BinToRankMap)
        ALLOCATE(ppiclf_BinToRankMap(0:nbins_in-1))
        ! BUGFIX(defensive): a freshly (re)allocated map is otherwise
        ! undefined, yet FindParticlePartLB reads it before
        ! PartLoadBalance fills it. Initialize to an invalid rank so an
        ! early/stale read is deterministic rather than undefined. The
        ! primary fix forces a rebuild on any binchanged step (see
        ! ppiclf_comm_FindParticlePartLB).
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

END MODULE ppiclf_DynamicAllocation
