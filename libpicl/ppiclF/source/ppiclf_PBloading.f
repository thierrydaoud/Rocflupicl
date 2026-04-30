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
        BinBuffer(i) = MAX(ppiclf_filter(i),ppiclf_nndist/2.0D0)
      END DO

      ! Looping through particles on this processor
      ! to find bin boundary locations
      DO i = 1,3
        local_extremes(2*i)   = -1.0D10 ! Dummy Max Val
        local_extremes(2*i-1) =  1.0D10 ! Dummy Min Val
      END DO
      DO i=1,ppiclf_npart
        DO j = 1,3
          ! Finding min/max particle extremes.
          ! Add buffer so that layer of outer cells 
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

      ppiclf_binchanged = .FALSE.
      ! If all particles are within last RK Stage binboundaries,
      ! do not calculate bins again and do not remap overlap grid
      DO i = 1,3
        IF((ppiclf_binb(2*i-1) + BinBuffer(i)) .LT.
     >             ppiclf_previousbinb(2*i-1)) THEN
          ppiclf_binchanged = .TRUE.
        END IF
        IF((ppiclf_binb(2*i)   - BinBuffer(i)) .GT.
     >             ppiclf_previousbinb(2*i))   THEN
          ppiclf_binchanged = .TRUE.
          EXIT
        END IF
      END DO

      CALL MPI_ALLREDUCE(MPI_IN_PLACE,ppiclf_binchanged,1
     >                   ,MPI_LOGICAL, MPI_LOR
     >                   ,ppiclf_comm, ierr)
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
        ppiclf_printbinvtu = .FALSE.
        RETURN
      ELSE
        ppiclf_printbinvtu = .TRUE.
      END IF

#ifdef TEST
       ! Set p-domain equal to f-domain for localized approach
       ! comparison
      IF(.FALSE.) THEN ! Need to think of correct conditional
        DO i = 1,3
          ppiclf_binb(2*i)    = ppiclf_xdrange(2,i) 
          ppiclf_binb(2*i-1)  = ppiclf_xdrange(1,i)
          ppiclf_BinDomLen(i) = ppiclf_binb(2*i) - ppiclf_binb(2*i-1) 
        END DO
      END IF
#endif
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
        ppiclf_BinDomLen(i) = ppiclf_binb(2*i) - ppiclf_binb(2*i-1) 
      END DO

      ! Set previous bin boundaries for next RK Stage check
      DO i = 1,6
        ppiclf_previousbinb(i) = ppiclf_binb(i)
      END DO

      ppiclf_totalBins = 1
      DO i = 1,3
        ppiclf_n_bins(i)  = MAX( 1,
     >                      FLOOR(ppiclf_BinDomLen(i)/BinMinLen(i)))
        ppiclf_bins_dx(i) = ppiclf_BinDomLen(i)/DBLE(ppiclf_n_bins(i))
        ppiclf_totalBins  = ppiclf_totalBins*ppiclf_n_bins(i)
      END DO

      ! Allocate all arrays dependant on number of bins or processors
      CALL ppiclf_dyn_alloc(ppiclf_totalBins,ppiclf_np)

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE ppiclf_comm_FindParticlePartLB

      USE ppiclf_DynamicAllocation
      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      INTEGER*4  i, j, ii, jj, kk, nbin, ierr, partcheck, NumBins, nRank
      INTEGER*4  maxPartPerBin
      REAL*8     LB_local, LB_criteria, LB_target

      ! If NumPart > LB_criteria*TargetNumPart -> Reassign BTRM
      ! (Bin to Rank Map)
      LB_criteria = 1.2D0

      DO i = 0,(ppiclf_totalBins - 1)
        ppiclf_ParticleCount(i) = 0
      END DO
      ppiclf_particleMoved = .FALSE.
#ifdef TEST
      ppiclf_particleMoved = .TRUE.
#endif
      DO i = 1,ppiclf_npart
        ! Calculates particle's bin index
        ii = FLOOR((ppiclf_y(1,i)-ppiclf_binb(1))/ppiclf_bins_dx(1))
        jj = FLOOR((ppiclf_y(2,i)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
        kk = FLOOR((ppiclf_y(3,i)-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
        ii = MAX(0, MIN(ii, ppiclf_n_bins(1)-1))
        jj = MAX(0, MIN(jj, ppiclf_n_bins(2)-1))
        kk = MAX(0, MIN(kk, ppiclf_n_bins(3)-1))
        ! Calculates particle's bin
        nbin = ii + ppiclf_n_bins(1)*jj + 
     >         ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk

        ! Maps particle to correct processor based on active bin number
        !ppiclf_iprop(4,i) = nrank ! Processor to send to
        ppiclf_iprop(5,i) = ii    ! x bin #
        ppiclf_iprop(6,i) = jj    ! y bin #
        ppiclf_iprop(7,i) = kk    ! z bin #
        ppiclf_iprop(8,i) = nbin ! total bin number
        nRank = ppiclf_BinToRankMap(nbin) ! This is based on old BTRM
        IF(nRank .NE. ppiclf_iprop(4,i)) THEN
          ppiclf_particleMoved = .TRUE.
          ppiclf_iprop(4,i) = ppiclf_BinToRankMap(nbin)
        END IF
        ppiclf_ParticleCount(nbin) = ppiclf_ParticleCount(nbin) + 1
      END DO
      ! Now sum particles per bin across MPI Ranks
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_ParticleCount
     >                   ,ppiclf_totalBins ,MPI_INTEGER4, MPI_SUM
     >                   ,ppiclf_comm, ierr)

      ! Logical OR comparison across MPI Ranks
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_particleMoved
     >                   ,1, MPI_LOGICAL, MPI_LOR
     >                   ,ppiclf_comm, ierr)

      ! Check if new BTRM required
      LB_target   = CEILING(DBLE(ppiclf_glnpart)/DBLE(ppiclf_np))
      LB_criteria = LB_criteria*LB_target
      ppiclf_rebalance = .FALSE.
      LB_local = ppiclf_ParticleCount(0) ! First bin count
      DO i = 1,(ppiclf_totalBins - 1)
        IF(ppiclf_BinToRankMap(i-1) .NE. ppiclf_BinToRankMap(i)) THEN
          LB_local = 0.0D0 !new Rank
        END IF
        LB_local = LB_local + ppiclf_ParticleCount(i)
        IF(LB_local .GT. LB_criteria) THEN
          ppiclf_rebalance = .TRUE.
          EXIT
        END IF
      END DO
      
      ! Logical OR comparison across MPI Ranks
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_rebalance
     >                   ,1, MPI_LOGICAL, MPI_LOR
     >                   ,ppiclf_comm, ierr)

      RETURN

      END SUBROUTINE
!----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_PartLoadBalance
   
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      INTEGER*4   ierr, i, j, targetParticleCnt, particleSum
     >           ,irank, bin, ii, jj, kk, nb1, nb2, nb3, nb1xnb2
     >           ,iloop, jloop, kloop, temp_dSize(3), itemp

      !CALL MPI_BARRIER(ppiclf_comm,ierr)

      ! Find the order to loop through dimensions. 
      ! Using largest dimension minimizes surface area between 
      ! processors, which likely minimizes ghost particle and
      ! overlap cell communication

      ! Sorting the domain lengths: ppiclf_dL = largest dimension,
      ! ppiclf_dM = medium dimension, ppiclf_dS = smallest dimension
      IF(ppiclf_nid .EQ. 0) THEN
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
          END DO
        END DO
      END IF

      CALL MPI_BCAST(temp_dSize,3,MPI_INTEGER4,0,ppiclf_comm,ierr)

      ppiclf_dL = temp_dSize(1) 
      ppiclf_dM = temp_dSize(2) 
      ppiclf_dS = temp_dSize(3) 

      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)
      nb1xnb2 = nb1*nb2

      targetParticleCnt = CEILING(DBLE(ppiclf_glnpart)/DBLE(ppiclf_np))

      particleSum   = 0
      irank = 0

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
          PRINT*, 'ppiclf_dL:',ppiclf_dL
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
            PRINT*, 'ppiclf_dM:',ppiclf_dM
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
              PRINT*, 'ppiclf_dS:',ppiclf_dS
              RETURN
            END IF
            bin  = ii + nb1*jj + nb1xnb2*kk
            particleSum = particleSum + ppiclf_ParticleCount(bin)
            ! This maps the bin to the MPI-rank
            ppiclf_BinToRankMap(bin) = irank
            IF(particleSum .GE. targetParticleCnt) THEN
              ! Reset Particle counter and advance MPI rank iteration
              particleSum = 0
              irank = irank + 1
              IF(irank .GT. ppiclf_np - 1) THEN
                ! Last rank will hold more than target number of
                ! particles
                irank = ppiclf_np - 1
              END IF
            END IF
          END DO !kloop (shortest bin dimension)
        END DO !jloop (medium bin dimension)
      END DO !iloop (longest bin dimension)

      ! Assign correct rank to each particle
      DO i =1,ppiclf_npart
        ! Now map particle to MPI Rank since we have a bin->rank map
        ii = ppiclf_iprop(8,i)
        ppiclf_iprop(4,i) = ppiclf_BinToRankMap(ii) ! Owning MPI Rank
      END DO

      CALL ppiclf_comm_setRankBoundaries
      CALL ppiclf_comm_setEmptyIndicator
      CALL ppiclf_comm_setInterfaceIndicator

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
        DO jj = 0, (nb2 - 1)
          DO ii = 0, (nb1 - 1)
            ! If THIS rank owns this bin, expand the bounding box
            IF(ppiclf_BinToRankMap(iBin) .EQ. ppiclf_nid) THEN
              min_i = MIN(min_i, ii)
              max_i = MAX(max_i, ii)
              
              min_j = MIN(min_j, jj)
              max_j = MAX(max_j, jj)
              
              min_k = MIN(min_k, kk)
              max_k = MAX(max_k, kk)
            END IF
            iBin = iBin + 1
          END DO
        END DO
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
        END DO
      END IF
      
      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------

      SUBROUTINE ppiclf_comm_setEmptyIndicator
  
      USE ppiclf_DynamicAllocation
      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      INTEGER*4   ierr, iBin, LbinCheck, i, j, k
     >           ,imin, imax, jmin, jmax, kmin, kmax
     >           ,ii, jj, kk, nb1, nb2, nb3, nb1xnb2, jInc, kInc
     >           ,di, dj, dk, ni, nj, nk, neighborBin

      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)
      nb1xnb2 = nb1*nb2

      iBin = 0
      DO kk = 0, nb3-1
        DO jj = 0, nb2-1
          DO ii = 0, nb1-1

            iBin = ii + nb1*jj + nb1xnb2*kk

            ppiclf_LMapFluid(iBin) = .FALSE.

            DO dk = -1,1
              DO dj = -1,1
                DO di = -1,1

                  ni = ii + di
                  nj = jj + dj
                  nk = kk + dk

                  CALL ppiclf_wrap_bin_index(ni, nb1, 1)
                  CALL ppiclf_wrap_bin_index(nj, nb2, 2)
                  CALL ppiclf_wrap_bin_index(nk, nb3, 3)

                  IF(ni .LT. 0 .OR. ni .GT. nb1-1) CYCLE
                  IF(nj .LT. 0 .OR. nj .GT. nb2-1) CYCLE
                  IF(nk .LT. 0 .OR. nk .GT. nb3-1) CYCLE

                  neighborBin = ni + nb1*nj + nb1xnb2*nk

                  IF(ppiclf_ParticleCount(neighborBin) .GT. 0) THEN
                    ppiclf_LMapFluid(iBin) = .TRUE.
                    EXIT
                  END IF

                END DO
                IF(ppiclf_LMapFluid(iBin)) EXIT
              END DO
              IF(ppiclf_LMapFluid(iBin)) EXIT
            END DO

          END DO
        END DO
      END DO

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
      INTEGER*4 side, di, dj, dk, ni, nj, nk, neighborBin, nRank
      LOGICAL sideHit

      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)
      nb1xnb2 = nb1*nb2

      ppiclf_LRankBoundary = .FALSE.

      iBin = 0
      DO kk = 0, nb3-1
        DO jj = 0, nb2-1
          DO ii = 0, nb1-1

            IF(.NOT. ppiclf_LMapFluid(iBin)) THEN
              iBin = iBin + 1
              CYCLE
            END IF

            iRank = ppiclf_BinToRankMap(iBin)

            DO side = 1,6
              sideHit = .FALSE.

              DO dk = -1,1
                DO dj = -1,1
                  DO di = -1,1

                    ! Only check the 9 bins on the requested side
                    IF(side .EQ. 1 .AND. di .NE. -1) CYCLE   ! -x
                    IF(side .EQ. 2 .AND. di .NE.  1) CYCLE   ! +x
                    IF(side .EQ. 3 .AND. dj .NE. -1) CYCLE   ! -y
                    IF(side .EQ. 4 .AND. dj .NE.  1) CYCLE   ! +y
                    IF(side .EQ. 5 .AND. dk .NE. -1) CYCLE   ! -z
                    IF(side .EQ. 6 .AND. dk .NE.  1) CYCLE   ! +z

                    ni = ii + di
                    nj = jj + dj
                    nk = kk + dk

                    CALL ppiclf_wrap_bin_index(ni, nb1, 1)
                    CALL ppiclf_wrap_bin_index(nj, nb2, 2)
                    CALL ppiclf_wrap_bin_index(nk, nb3, 3)

                    IF(ni .LT. 0 .OR. ni .GT. nb1-1) CYCLE
                    IF(nj .LT. 0 .OR. nj .GT. nb2-1) CYCLE
                    IF(nk .LT. 0 .OR. nk .GT. nb3-1) CYCLE

                    neighborBin = ni + nb1*nj + nb1xnb2*nk

                    IF(.NOT. ppiclf_LMapFluid(neighborBin)) CYCLE

                    nRank = ppiclf_BinToRankMap(neighborBin)
                    IF(nRank .NE. iRank) THEN
                      ppiclf_LRankBoundary(iBin,side) = .TRUE.
                      sideHit = .TRUE.
                      EXIT
                    END IF

                  END DO
                  IF(sideHit) EXIT
                END DO
                IF(sideHit) EXIT
              END DO
            END DO

            iBin = iBin + 1
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------

      SUBROUTINE ppiclf_wrap_bin_index(idx, nbin, jj)

      IMPLICIT NONE
      INCLUDE "PPICLF"

      INTEGER*4 idx, nbin, jj
      LOGICAL lperiodic, lequal

      lperiodic = ppiclf_linperiodic(jj)
      lequal    = ppiclf_EqualDomain(jj)

      IF(lperiodic .AND. lequal .AND. nbin .GT. 1) THEN
        IF(idx .EQ. -1)   idx = nbin - 1
        IF(idx .EQ. nbin) idx = 0
      END IF

      RETURN
      END SUBROUTINE

!----------------------------------------------------------------------

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
     
      RETURN

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
     >          nl, nii, njj, nrr, iic, jjc, kkc, ierr, nb1, nb2,
     >          nb3, nb1xnb2 
      INTEGER*4 ix, iy, iz, ixLow, ixHigh, iyLow,
     >          iyHigh, izLow, izHigh, ibin, jbin, kbin, nbin, bin,
     >          nRankMaps, RankMaps(27,5)
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

      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)
      nb1xnb2 = nb1*nb2

      ! Loops through number of fluid FV cells on this processor
      DO ie=1,ppiclf_nFVCells  
        centeri(1) = ppiclf_fluid_grid(1,ie)
        centeri(2) = ppiclf_fluid_grid(2,ie)
        centeri(3) = ppiclf_fluid_grid(3,ie)
        ! Cycles if fluid cell center is outside of any bin boundaries 
        IF (centeri(1) .GT. ppiclf_binb(2)) CYCLE
        IF (centeri(2) .GT. ppiclf_binb(4)) CYCLE
        IF (centeri(3) .GT. ppiclf_binb(6)) CYCLE
        IF (centeri(1) .LT. ppiclf_binb(1)) CYCLE
        IF (centeri(2) .LT. ppiclf_binb(3)) CYCLE
        IF (centeri(3) .LT. ppiclf_binb(5)) CYCLE
        ! Determines what bin the fluid cell is nominally mapped to
        ibin = FLOOR((centeri(1)-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
        jbin = FLOOR((centeri(2)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
        kbin = FLOOR((centeri(3)-ppiclf_binb(5))/ppiclf_bins_dx(3))
        ibin = MAX(0, MIN(ibin, nb1-1))
        jbin = MAX(0, MIN(jbin, nb2-1))
        kbin = MAX(0, MIN(kbin, nb3-1))
         ! Calculates processor rank
        bin  = ibin + nb1*jbin + nb1xnb2*kbin
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
     >                               .AND. nb1 .GT. 1) THEN
                IF(ii .EQ. nb1) ii = 0
                IF(ii .EQ. -1)  ii = nb1 - 1
              END IF
              IF(ppiclf_linperiodic(2) .AND. ppiclf_EqualDomain(2)
     >                               .AND. nb2 .GT. 1) THEN
                IF(jj .EQ. nb2) jj = 0
                IF(jj .EQ. -1)  jj = nb2 - 1
              END IF
              IF(ppiclf_linperiodic(3) .AND. ppiclf_EqualDomain(3)
     >                               .AND. nb3 .GT. 1) THEN
                IF(kk .EQ. nb3) kk = 0
                IF(kk .EQ. -1)  kk = nb3 - 1
              END IF
              
              IF (ii .LT. 0 .OR. ii .GT. nb1-1) CYCLE
              IF (jj .LT. 0 .OR. jj .GT. nb2-1) CYCLE
              IF (kk .LT. 0 .OR. kk .GT. nb3-1) CYCLE
              nbin  = ii + nb1*jj + nb1xnb2*kk
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
                  RankMaps(nRankMaps,2) = nbin  !original index
                  RankMaps(nRankMaps,3) = ibin !original index
                  RankMaps(nRankMaps,4) = jbin !original index
                  RankMaps(nRankMaps,5) = kbin !original index
                 END IF !MapCell
              END IF !ppiclf_LMapFluid
            END DO !iz
          END DO !iy
        END DO !ix
        IF(nRankMaps .GT. 0) THEN
          DO i = 1,nRankMaps
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
            ppiclf_cell_map(7,ppiclf_nCells_FV2PICL) = RankMaps(i,5)
          END DO !nRankMaps
        END IF
      END DO !ie

      DO ie=1,ppiclf_nCells_FV2PICL 
        ! These copy all indicies since Fortran is column-major
        iee = ppiclf_cell_map(1,ie)
        CALL ppiclf_copy(ppiclf_picl_grid(1,ie)
     >                 ,ppiclf_fluid_grid(1,iee),7)
      END DO

      ! Copy mapping since it is need to send fluid properties in interp
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
        ppiclf_filter(l) = Max_CellLen(l)*2.0D0
      END DO
      ! Find max filter across processors
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, ppiclf_filter
     >                   ,3 ,MPI_DOUBLE_PRECISION
     >                   ,MPI_MAX ,ppiclf_comm, ierr)

      DO l = 1,3
        ! Multiply by 1.5 so particle near face will
        ! find center one cell over in farthest direction
        ppiclf_interp_dchk(l) = ppiclf_filter(l)/2.0D0*1.5D0
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

      USE ppiclf_DynamicAllocation
      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      REAL*8     GhostPos(3), PeriodicShift(3), buffer, distchk
      REAL*8     xlo(3), xhi(3)
      REAL*8     createdPos(3,27)
      INTEGER*4  ip, idum, iip, jjp, kkp, iig, jjg, kkg
      INTEGER*4  nrank, j, k, l, nbin
      INTEGER*4  nb1, nb2, nb3, nb1xnb2, iBin
      INTEGER*4  ix, iy, iz, ixLow, ixHigh, iyLow, iyHigh
      INTEGER*4  izLow, izHigh

      DO l = 1,3
        IF(ppiclf_linperiodic(l)) THEN
          PeriodicShift(l) = ppiclf_binb(2*l)
     >                     - ppiclf_binb(2*l-1)
        ELSE
          PeriodicShift(l) = 0.0D0
        END IF
      END DO

      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)
      nb1xnb2 = nb1*nb2

      buffer  = 1.1D0
      distchk = ppiclf_nndist*buffer

      ppiclf_npart_gp = 0

      DO ip = 1,ppiclf_npart

        idum = 0
        DO j = 1,PPICLF_LRS
          idum = idum + 1
          ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip)
        END DO

        DO j = 1,PPICLF_LRP
          idum = idum + 1
          ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip)
        END DO

        iip  = ppiclf_iprop(5,ip)
        jjp  = ppiclf_iprop(6,ip)
        kkp  = ppiclf_iprop(7,ip)
        iBin = ppiclf_iprop(8,ip)

        ixLow  = 0
        ixHigh = 0
        iyLow  = 0
        iyHigh = 0
        izLow  = 0
        izHigh = 0

        xlo(1) = ppiclf_binb(1) + DBLE(iip)*ppiclf_bins_dx(1)
        xhi(1) = xlo(1) + ppiclf_bins_dx(1)

        xlo(2) = ppiclf_binb(3) + DBLE(jjp)*ppiclf_bins_dx(2)
        xhi(2) = xlo(2) + ppiclf_bins_dx(2)

        xlo(3) = ppiclf_binb(5) + DBLE(kkp)*ppiclf_bins_dx(3)
        xhi(3) = xlo(3) + ppiclf_bins_dx(3)

        IF(ppiclf_LRankBoundary(iBin,1) .OR.
     >     (ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)
     >     .AND. iip .EQ. 0)) THEN
          IF(ABS(ppiclf_cp_map(1,ip)-xlo(1)) .LT. distchk)
     >      ixLow = -1
        END IF

        IF(ppiclf_LRankBoundary(iBin,2) .OR.
     >     (ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)
     >     .AND. iip .EQ. nb1-1)) THEN
          IF(ABS(ppiclf_cp_map(1,ip)-xhi(1)) .LT. distchk)
     >      ixHigh = 1
        END IF

        IF(ppiclf_LRankBoundary(iBin,3) .OR.
     >     (ppiclf_linperiodic(2) .AND. ppiclf_EqualDomain(2)
     >     .AND. jjp .EQ. 0)) THEN
          IF(ABS(ppiclf_cp_map(2,ip)-xlo(2)) .LT. distchk)
     >      iyLow = -1
        END IF

        IF(ppiclf_LRankBoundary(iBin,4) .OR.
     >     (ppiclf_linperiodic(2) .AND. ppiclf_EqualDomain(2)
     >     .AND. jjp .EQ. nb2-1)) THEN
          IF(ABS(ppiclf_cp_map(2,ip)-xhi(2)) .LT. distchk)
     >      iyHigh = 1
        END IF

        IF(ppiclf_LRankBoundary(iBin,5) .OR.
     >     (ppiclf_linperiodic(3) .AND. ppiclf_EqualDomain(3)
     >     .AND. kkp .EQ. 0)) THEN
          IF(ABS(ppiclf_cp_map(3,ip)-xlo(3)) .LT. distchk)
     >      izLow = -1
        END IF

        IF(ppiclf_LRankBoundary(iBin,6) .OR.
     >     (ppiclf_linperiodic(3) .AND. ppiclf_EqualDomain(3)
     >     .AND. kkp .EQ. nb3-1)) THEN
          IF(ABS(ppiclf_cp_map(3,ip)-xhi(3)) .LT. distchk)
     >      izHigh = 1
        END IF

        DO ix = ixLow,ixHigh
          DO iy = iyLow,iyHigh
            DO iz = izLow,izHigh

              IF(ix .EQ. 0 .AND. iy .EQ. 0 .AND.
     >           iz .EQ. 0) CYCLE

              iig = iip + ix
              jjg = jjp + iy
              kkg = kkp + iz

              GhostPos(1) = ppiclf_cp_map(1,ip)
              GhostPos(2) = ppiclf_cp_map(2,ip)
              GhostPos(3) = ppiclf_cp_map(3,ip)

              IF(iig .LT. 0 .OR. iig .GT. nb1-1) THEN
                IF(ppiclf_linperiodic(1) .AND.
     >             ppiclf_EqualDomain(1)) THEN
                  IF(iig .LT. 0) THEN
                    iig = nb1 - 1
                    GhostPos(1) = GhostPos(1) + PeriodicShift(1)
                  ELSE
                    iig = 0
                    GhostPos(1) = GhostPos(1) - PeriodicShift(1)
                  END IF
                ELSE
                  CYCLE
                END IF
              END IF

              IF(jjg .LT. 0 .OR. jjg .GT. nb2-1) THEN
                IF(ppiclf_linperiodic(2) .AND.
     >             ppiclf_EqualDomain(2)) THEN
                  IF(jjg .LT. 0) THEN
                    jjg = nb2 - 1
                    GhostPos(2) = GhostPos(2) + PeriodicShift(2)
                  ELSE
                    jjg = 0
                    GhostPos(2) = GhostPos(2) - PeriodicShift(2)
                  END IF
                ELSE
                  CYCLE
                END IF
              END IF

              IF(kkg .LT. 0 .OR. kkg .GT. nb3-1) THEN
                IF(ppiclf_linperiodic(3) .AND.
     >             ppiclf_EqualDomain(3)) THEN
                  IF(kkg .LT. 0) THEN
                    kkg = nb3 - 1
                    GhostPos(3) = GhostPos(3) + PeriodicShift(3)
                  ELSE
                    kkg = 0
                    GhostPos(3) = GhostPos(3) - PeriodicShift(3)
                  END IF
                ELSE
                  CYCLE
                END IF
              END IF

              nbin  = iig + nb1*jjg + nb1xnb2*kkg
              nrank = ppiclf_BinToRankMap(nbin)
              IF(nrank .EQ. ppiclf_nid) THEN
                IF(ppiclf_linperiodic(1) .AND.
     >             ppiclf_EqualDomain(1)
     >             .OR.
     >             ppiclf_linperiodic(2) .AND.
     >             ppiclf_EqualDomain(2)
     >             .OR.
     >             ppiclf_linperiodic(3) .AND.
     >             ppiclf_EqualDomain(3)) THEN
                  ! Keep going
                ELSE
                  CYCLE
                END IF
              END IF

              ppiclf_npart_gp = ppiclf_npart_gp + 1

              ppiclf_iprop_gp(1,ppiclf_npart_gp) = ppiclf_iprop(1,ip)
              ppiclf_iprop_gp(2,ppiclf_npart_gp) = ppiclf_iprop(2,ip)
              ppiclf_iprop_gp(3,ppiclf_npart_gp) = ppiclf_iprop(3,ip)

              ppiclf_iprop_gp(4,ppiclf_npart_gp) = nrank

              ppiclf_iprop_gp(5,ppiclf_npart_gp) = iig
              ppiclf_iprop_gp(6,ppiclf_npart_gp) = jjg
              ppiclf_iprop_gp(7,ppiclf_npart_gp) = kkg
              ppiclf_iprop_gp(8,ppiclf_npart_gp) = nbin

              ppiclf_rprop_gp(1,ppiclf_npart_gp) = GhostPos(1)
              ppiclf_rprop_gp(2,ppiclf_npart_gp) = GhostPos(2)
              ppiclf_rprop_gp(3,ppiclf_npart_gp) = GhostPos(3)

              DO k = 4,PPICLF_LRP_GP
                ppiclf_rprop_gp(k,ppiclf_npart_gp) =
     >             ppiclf_cp_map(PPICLF_LRS + k, ip)
              END DO

            END DO
          END DO
        END DO

      END DO

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
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_subbinParticleMap

      USE ppiclf_DynamicAllocation

      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
!
! Input:
!
      INTEGER*4  i, n_SBin(3), total_SBin
     >          ,maxParticlePerBin, binOffset(3),tempSBin, i_SBin(3) 
     >          ,iterationLoop
      
      IF(ppiclf_npart .LT. 1) RETURN

      total_SBin = 1
      DO i = 1,3
        ! +3 to account for one layer of overlap cells or ghost
        ! particles in adjacent bins on either side of rank boundary
        n_SBin(i) = ppiclf_binBIndex(2*i)-ppiclf_binBIndex(2*i-1) + 3 
        binOffset(i) = ppiclf_BinBIndex(2*i-1) - 1
        ! Following IF corrects when full particle domain is used
        IF(n_SBin(i) .GT. ppiclf_n_bins(i)) THEN
          n_SBin(i) = ppiclf_n_bins(i)
          binOffset(i) = 0
        END IF
        total_SBin = total_SBin*n_SBin(i)    
      END DO

      maxParticlePerBin = 1
      DO iterationLoop = 1,2
        IF(iterationLoop .EQ. 2) THEN
          maxParticlePerBin = MAX(1,MAXVAL(ppiclf_binPartCount))
        END IF
        ! Allocates binPartList and binPartCount Arrays
        CALL ppiclf_allocate_BTP(total_SBin,maxParticlePerBin)

        ! Assign Subbin counters to 0
        ppiclf_binPartList  = 0
        ppiclf_binPartCount = 0

        ! Map each real particle to a subbin
        DO i = 1,ppiclf_npart
          i_SBin(1) = ppiclf_iprop(5,i) - binOffset(1)
          i_SBin(2) = ppiclf_iprop(6,i) - binOffset(2)
          i_SBin(3) = ppiclf_iprop(7,i) - binOffset(3)
          tempSBin = i_SBin(1) + n_SBin(1)*i_SBin(2) +
     >               n_SBin(1)*n_SBin(2)*i_SBin(3)
          IF(tempSBin .LT. 0 .OR. tempSBin .GT. total_SBin-1) THEN
            PRINT*, 'ERROR: Bad Subbin Index in Real Particle Mapping',
     >              tempSBin
            CALL ppiclf_exittr('',0.0D0,0)
          END IF
          ppiclf_binPartCount(tempSBin) = 
     >                       ppiclf_binPartCount(tempSBin)+1
          IF(iterationLoop .EQ. 2) THEN
            IF(ppiclf_binPartCount(tempSBin).GT. maxParticlePerBin) THEN
              PRINT*, 'ERROR: More Particles per Subbin than expected'
              CALL ppiclf_exittr('',0.0D0,0)
            END IF
            ppiclf_binPartList(tempSBin,
     >                  ppiclf_binPartCount(tempSBin)) = i
          END IF
        END DO ! real particle loop
        ! Map each ghost particle to a subbin
        DO i = 1,ppiclf_npart_gp
          i_SBin(1) = ppiclf_iprop_gp(5,i) - binOffset(1)
          i_SBin(2) = ppiclf_iprop_gp(6,i) - binOffset(2)
          i_SBin(3) = ppiclf_iprop_gp(7,i) - binOffset(3)
          tempSBin = i_SBin(1) + n_SBin(1)*i_SBin(2) +
     >               n_SBin(1)*n_SBin(2)*i_SBin(3)
          IF(tempSBin .LT. 0 .OR. tempSBin .GT. total_SBin-1) THEN
            PRINT*, 'ERROR: Bad Subbin Index in Ghost Particle Mapping',
     >              tempSBin
            CALL ppiclf_exittr('',0.0D0,0)
          END IF
          ppiclf_binPartCount(tempSBin)=ppiclf_binPartCount(tempSBin)+1
          IF(iterationLoop .EQ. 2) THEN
            IF(ppiclf_binPartCount(tempSBin) .GT.maxParticlePerBin) THEN
              PRINT*, 'ERROR: More Particles per Subbin than expected'
              CALL ppiclf_exittr('',0.0D0,0)
            END IF
            ! negative in subbin map means it is ghost particle
            ppiclf_binPartList(tempSBin,
     >                   ppiclf_binPartCount(tempSBin)) = -i
          END IF
        END DO ! gp loop
      END DO !iterationLoop
      RETURN

      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE ppiclf_comm_subbinCellMap

      USE ppiclf_DynamicAllocation

      IMPLICIT NONE
!
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
!
! Input:
!
      INTEGER*4  ie, i, j, k, n_SBin(3), total_SBin, iTemp_SBin(3)
     >          ,maxCellsPerBin, binOffset(3),tempSBin, i_SBin(3)
     >          , iterationLoop 
      
      IF(ppiclf_npart .LT. 1) RETURN

      IF(ppiclf_nCells_Interp .LE. 0 .AND. ppiclf_npart .GT. 0) THEN
        PRINT*,'ERROR: ',ppiclf_npart, 'Particles mapped to bin:'
     >         ,ppiclf_nid
        PRINT*,'No cells mapped to bin for Interpolation/Projection.'
        CALL ppiclf_exittr('Failure in particle to cell mapping$',
     >                      0.D0,0)
      END IF
 
      total_SBin = 1
      DO i = 1,3
        ! +3 to account for one layer of overlap cells or ghost
        ! particles in adjacent bins on either side of rank boundary
        n_SBin(i) = ppiclf_binBIndex(2*i)-ppiclf_binBIndex(2*i-1) + 3 
        binOffset(i) = ppiclf_BinBIndex(2*i-1) - 1
        ! Following IF corrects when full particle domain is used
        IF(n_SBin(i) .GT. ppiclf_n_bins(i)) THEN
          n_SBin(i) = ppiclf_n_bins(i)
          binOffset(i) = 0
        END IF
        total_SBin = total_SBin*n_SBin(i)    
      END DO

      ! iterationLoop 1 -> Find the maximum cell count per bin for array
      ! allocation
      ! iterationLoop 2 -> Creates list of cells per bin
      maxCellsPerBin = 1
      DO iterationLoop = 1,2
        IF(iterationLoop .EQ. 2) THEN
          maxCellsPerBin = MAX(1,MAXVAL(ppiclf_binCellCount))
        END IF
        ! Allocated binCellList and binCellCount Arrays
        CALL ppiclf_allocate_BTC(total_SBin,maxCellsPerBin)
        ! Assign Subbin counters to 0
        ppiclf_binCellList  = 0
        ppiclf_binCellCount = 0

        ! NOTE: Overlap cells have not been shifted for periodicity
        !       Both the cell centroid and cell bin indices are of the
        !       original fluid cell.
        ! Map each overlap cells to subbins
        ! In the i,j,k loops below, 0 takes care of non-periodic mapping
        ! and 1 takes care of periodic mapping.  If a cell is in corner,
        ! it'll be mapped to 2*2*2=8 subbins.
        DO ie = 1,ppiclf_nCells_Interp
          i_SBin(1) = ppiclf_cell_map(5,ie) - binOffset(1)
          i_SBin(2) = ppiclf_cell_map(6,ie) - binOffset(2)
          i_SBin(3) = ppiclf_cell_map(7,ie) - binOffset(3)
          DO i = 0,1
            IF(i .EQ. 0) THEN
              iTemp_SBin(1) = i_SBin(1)
              IF(iTemp_SBin(1) .LT. 0 .OR. 
     >           iTemp_SBin(1) .GT. n_SBin(1) - 1) CYCLE 
            ELSE ! i .EQ. 1
              IF(ppiclf_linperiodic(1) .AND. ppiclf_EqualDomain(1)) THEN
                IF(i_SBin(1) .LE. 0) THEN
                  iTemp_SBin(1) = n_SBin(1) - 1
                  IF(iTemp_SBin(1) .EQ. i_SBin(1)) CYCLE 
                ELSE IF(i_SBin(1) .GE. n_SBin(1) - 1) THEN
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
     >             iTemp_SBin(2) .GT. n_SBin(2) - 1) CYCLE
              ELSE ! j .EQ. 1
                ! This takes care of periodicity for single processor
                IF(ppiclf_linperiodic(2).AND.ppiclf_EqualDomain(2)) THEN
                  IF(i_SBin(2) .LE. 0) THEN
                    iTemp_SBin(2) = n_SBin(2) - 1
                    IF(iTemp_SBin(2) .EQ. i_SBin(2)) CYCLE
                  ELSE IF(i_SBin(2) .GE. n_SBin(2) - 1) THEN
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
     >               iTemp_SBin(3) .GT. n_SBin(3) - 1) CYCLE
                ELSE ! k .EQ. 1
                  ! This takes care of periodicity for single processor
                  IF(ppiclf_linperiodic(3) .AND. 
     >                                      ppiclf_EqualDomain(3)) THEN 
                    IF(i_SBin(3) .LE. 0) THEN
                      iTemp_SBin(3) = n_SBin(3) - 1
                      IF(iTemp_SBin(3) .EQ. i_SBin(3)) CYCLE
                    ELSE IF(i_SBin(3) .GE. n_SBin(3) - 1) THEN
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
                tempSBin = iTemp_SBin(1) + n_SBin(1)*iTemp_SBin(2) +
     >                    n_SBin(1)*n_SBin(2)*iTemp_SBin(3)
                IF(tempSBin .LT. 0 .OR. tempSBin .GT. total_SBin-1) THEN
                PRINT*, 'ERROR:Bad Subbin Index in Overlap Cell Mapping'
     >                   , tempSBin
                  CALL ppiclf_exittr('',0.0D0,0)
                END IF
                ppiclf_binCellCount(tempSBin) = 
     >                        ppiclf_binCellCount(tempSBin) + 1
                IF(iterationLoop .EQ. 2) THEN
                  IF(ppiclf_binCellCount(tempSBin) 
     >                                    .GT. maxCellsPerBin) THEN
                    PRINT*, 'ERROR: More Cells per Subbin than expected'
                    CALL ppiclf_exittr('',0.0D0,0)
                  END IF
                  ppiclf_binCellList(tempSBin,
     >                               ppiclf_binCellCount(tempSBin)) = ie
                END IF
              END DO !k
            END DO !j 
          END DO !i
        END DO !ie
      END DO !iterationLoop

      RETURN

      END SUBROUTINE
!-----------------------------------------------------------------------
