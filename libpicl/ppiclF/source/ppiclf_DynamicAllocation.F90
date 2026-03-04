MODULE ppiclf_DynamicAllocation

  IMPLICIT NONE

  INTEGER*4, ALLOCATABLE :: ppiclf_ParticleCount(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_BinToRankMap(:)
  INTEGER*4, ALLOCATABLE :: ppiclf_IRankBoundary(:,:)
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

    IF(ALLOCATED(ppiclf_IRankBoundary)) THEN
      DEALLOCATE(ppiclf_IRankBoundary)
    END IF

    ALLOCATE(ppiclf_IRankBoundary(0:np_in-1,6))

    IF(ALLOCATED(ppiclf_LRankBoundary)) THEN
      DEALLOCATE(ppiclf_LRankBoundary)
    END IF

    ALLOCATE(ppiclf_LRankBoundary(0:nbins_in-1,6))

    IF(ALLOCATED(ppiclf_LMapFluid)) THEN
      DEALLOCATE(ppiclf_LMapFluid)
    END IF

    ALLOCATE(ppiclf_LMapFluid(0:nbins_in-1))


  END SUBROUTINE

END MODULE ppiclf_DynamicAllocation
