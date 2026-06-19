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
  LOGICAL   :: ppiclf_useFineGrid = .FALSE.
  INTEGER*4 :: ppiclf_nFine(3), ppiclf_total_fineBin
  INTEGER*4 :: ppiclf_maxPartPerFineBin = 0
  REAL*8    :: ppiclf_fineLo(3), ppiclf_fineInvLen(3)
  INTEGER*4, ALLOCATABLE :: ppiclf_finePartCount(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_finePartList(:,:)

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
  ! =====================================================================
  SUBROUTINE ppiclf_allocate_FineBTP(nbins_in, maxPartPerBin)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nbins_in, maxPartPerBin
    INTEGER*4 size_dim1, size_dim2

    IF (ALLOCATED(ppiclf_finePartCount)) THEN
        IF (SIZE(ppiclf_finePartCount) .NE. nbins_in) THEN
            DEALLOCATE(ppiclf_finePartCount)
            ALLOCATE(ppiclf_finePartCount(0:nbins_in-1))
        END IF
    ELSE
        ALLOCATE(ppiclf_finePartCount(0:nbins_in-1))
    END IF

    IF (ALLOCATED(ppiclf_finePartList)) THEN
      size_dim1 = SIZE(ppiclf_finePartList, 1)
      size_dim2 = SIZE(ppiclf_finePartList, 2)
      IF (size_dim1 .NE. nbins_in .OR. size_dim2 .LT. maxPartPerBin) THEN
        DEALLOCATE(ppiclf_finePartList)
        ALLOCATE(ppiclf_finePartList(0:nbins_in-1, maxPartPerBin))
      END IF
    ELSE
      ALLOCATE(ppiclf_finePartList(0:nbins_in-1, maxPartPerBin))
    END IF
  END SUBROUTINE

  ! =====================================================================
  ! Fine Sub-Bin Particle Mapping RE-Allocation (grow-on-demand, safe)
  ! =====================================================================
  SUBROUTINE ppiclf_reallocate_FineBTP(nbins_in, currentMax, newMax)
    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nbins_in, currentMax, newMax
    INTEGER*4, ALLOCATABLE :: temp_list(:,:)

    ALLOCATE(temp_list(0:nbins_in-1, newMax))
    ! Preserve everything written so far; new columns are filled by the
    ! caller before they are ever read (count-controlled).
    temp_list(:, 1:currentMax) = ppiclf_finePartList(:, 1:currentMax)
    CALL MOVE_ALLOC(FROM=temp_list, TO=ppiclf_finePartList)
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
