MODULE ppiclf_DynamicAllocation

  IMPLICIT NONE

  INTEGER*4, ALLOCATABLE :: ppiclf_ParticleCount(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_BinToRankMap(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_IRankBoundary(:,:)
  INTEGER*4, ALLOCATABLE :: ppiclf_binCellCount(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_binCellList(:,:)
  INTEGER*4, ALLOCATABLE :: ppiclf_binPartCount(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_binPartList(:,:)
  LOGICAL,   ALLOCATABLE :: ppiclf_LRankBoundary(:,:)
  LOGICAL,   ALLOCATABLE :: ppiclf_LMapFluid(:)
  INTEGER*4 :: ppiclf_dL, ppiclf_dM, ppiclf_dS
 
CONTAINS

  SUBROUTINE ppiclf_dyn_alloc(nbins_in,np_in)

    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nbins_in
    INTEGER*4, INTENT(IN) :: np_in

    IF(ALLOCATED(ppiclf_ParticleCount)) THEN
      DEALLOCATE(ppiclf_ParticleCount)
    END IF

    ALLOCATE(ppiclf_ParticleCount(0:nbins_in-1))

    IF(ALLOCATED(ppiclf_BinToRankMap)) THEN
      DEALLOCATE(ppiclf_BinToRankMap)
    END IF

    ALLOCATE(ppiclf_BinToRankMap(0:nbins_in-1))

#ifdef TEST
    IF(ALLOCATED(ppiclf_IRankBoundary)) THEN
      DEALLOCATE(ppiclf_IRankBoundary)
    END IF

    ALLOCATE(ppiclf_IRankBoundary(0:np_in-1,6))
#else
    IF(.NOT. ALLOCATED(ppiclf_IRankBoundary)) THEN
      ALLOCATE(ppiclf_IRankBoundary(0:np_in-1,6))
    END IF
#endif

    IF(ALLOCATED(ppiclf_LRankBoundary)) THEN
      DEALLOCATE(ppiclf_LRankBoundary)
    END IF

    ALLOCATE(ppiclf_LRankBoundary(0:nbins_in-1,6))

    IF(ALLOCATED(ppiclf_LMapFluid)) THEN
      DEALLOCATE(ppiclf_LMapFluid)
    END IF

    ALLOCATE(ppiclf_LMapFluid(0:nbins_in-1))


  END SUBROUTINE

  SUBROUTINE ppiclf_allocate_BTP(nbins_in,maxPartPerBin)

    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nbins_in
    INTEGER*4, INTENT(IN) :: maxPartPerBin

    ! In future, maybe don't deallocate. Doing now for safety
    IF(ALLOCATED(ppiclf_binPartCount)) THEN
      DEALLOCATE(ppiclf_binPartCount)
    END IF
    IF(ALLOCATED(ppiclf_binPartList)) THEN
      DEALLOCATE(ppiclf_binPartList)
    END IF
    
    ALLOCATE(ppiclf_binPartList(0:nbins_in-1,maxPartPerBin)) 
    ALLOCATE(ppiclf_binPartCount(0:nbins_in-1)) 

  END SUBROUTINE

  SUBROUTINE ppiclf_allocate_BTC(nbins_in, nCells)

    IMPLICIT NONE
    INTEGER*4, INTENT(IN) :: nbins_in
    INTEGER*4, INTENT(IN) :: nCells

    ! In future, maybe don't deallocate. Doing now for safety
    IF(ALLOCATED(ppiclf_binCellCount)) THEN
      DEALLOCATE(ppiclf_binCellCount)
    END IF

    IF(ALLOCATED(ppiclf_binCellList)) THEN
      DEALLOCATE(ppiclf_binCellList)
    END IF

    ALLOCATE(ppiclf_binCellCount(0:nbins_in-1))
    ALLOCATE(ppiclf_binCellList(0:nbins_in-1,nCells))

  END SUBROUTINE


END MODULE ppiclf_DynamicAllocation
