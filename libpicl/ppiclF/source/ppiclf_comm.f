!-----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_comm_InitMPI(comm,id,np)
     > bind(C, name="ppiclc_comm_InitMPI")
#else
      SUBROUTINE ppiclf_comm_InitMPI(comm,id,np)
#endif
!
!     This subroutine is called from rocflu/RFLU_InitFlowSolver.F90
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input: 
!
      INTEGER*4 comm
      INTEGER*4 id
      INTEGER*4 np
!
! Code:
!
      ! Ensures a later subroutine init wasn't called out of order
      IF (PPICLF_LINIT .OR. PPICLF_OVERLAP)
     >   CALL ppiclf_exittr('InitMPI must be called first$',0.0d0,0)

      ! set ppiclf_processor information
      ppiclf_comm = comm
      ppiclf_nid  = id
      ppiclf_np   = np

      ! GSlib call
      CALL ppiclf_prints('   *Begin InitCrystal$')
         CALL ppiclf_comm_InitCrystal
      CALL ppiclf_prints('    End InitCrystal$')

      ! check to make sure subroutine is called in correct order later
      ! on in the code sequence
      PPICLF_LCOMM = .TRUE.

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_InitCrystal
!
!     This subroutine is called form ppiclf_comm_InitMPI
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"

!
! Code:
!
      ! GSlib call to setup crystal router for communication across
      ! processors
      CALL pfgslib_crystal_setup(ppiclf_cr_hndl,ppiclf_comm,ppiclf_np)

      RETURN
      END
!-----------------------------------------------------------------------
#ifdef PPICLC
      SUBROUTINE ppiclf_comm_InitOverlapGrid(ncell,fluidGrid)
     > bind(C, name="ppiclc_comm_InitOverlapGrid")
#else
      SUBROUTINE ppiclf_comm_InitOverlapGrid(ncell,fluidGrid)
#endif
!
! This subroutine is called from rocpicl/PICL_TEMP_InitSolver.F90
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4 ncell
      REAL*8    fluidGrid(7,ncell)
      ! Expected size: fluidGrid(7,ncell)
      ! Indicies 1-3: Centroid x,y,z position
      ! Indicies 4-6: Max cell dx,dy,dz based on any vertex combination
      ! Index      7: Cell Volume 
!
! External:
!
      INTEGER*4 ie, i, ierr
!
      ppiclf_overlap = .TRUE.

      IF(.NOT. PPICLF_LCOMM)
     > CALL ppiclf_exittr('InitMPI must be before InitOverlap$',0.0d0,0)
      IF(.NOT. PPICLF_LINIT)
     > CALL ppiclf_exittr('InitParticle must be before InitOverlap$'
     >                  ,0.0d0,0)

      IF(ncell .GT. PPICLF_LEE .OR. ncell .LT. 0) THEN
        PRINT*, '***ERROR*** PPICLF_LEE', PPICLF_LEE, 'in', 
     >   'InitMapOverlapGrid must be greater than', ncell
        CALL MPI_BARRIER(ppiclf_comm,ierr) 
        CALL ppiclf_exittr('Increase LEE in InitOverlap$',0.0d0,ncell)
      END IF

      ! Number of finite volume cells from fluid solver
      ppiclf_nFVCells = ncell

      DO ie=1,ppiclf_nFVCells
        DO i = 1,7
          ppiclf_fluid_grid(i,ie) = fluidGrid(i,ie)
        END DO
      END DO

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_CreateBin

      IMPLICIT NONE
      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
!
! Internal:
!
      INTEGER*4 ix, iy, iz, npt_total, i, j, k, idum, jdum, kdum, 
     >          total_bin, targetTotBin, idealBin(3), Temp_iBin(3),
     >          iBin(3), iBinTot, tempi, ideal_bin_index(3), largeBin, 
     >          medBin, smallBin,
     >          ppiclf_iglsum, NBMax, ierr, MaxPotentialBins(3), 
     >          maxbincount1, maxbincount2, BinCheck, minbin(3),
     >          ppiclf_iglmax, binNegBound, binPosBound, binIterations
     >          
      REAL*8    xmin, ymin, zmin, xmax, ymax, zmax, temp1, temp2,
     >          binb_length(3), BinMinLen(3), ppiclf_glmin, 
     >          ppiclf_glmax, ppiclf_glsum, periodicDistCheck,
     >          BinBuffer(3), binsReal(3), binError, minBinError,
     >          increaseRatio
      EXTERNAL  ppiclf_iglsum, ppiclf_glmin, ppiclf_glmax, ppiclf_glsum,
     >          ppiclf_iglmax
      LOGICAL   MaxBinsAchieved(3), TwoSmallBins
!

      ix = 1
      iy = 2
      iz = 3
        
      ! iglsum is integer global sum across MPI ranks.
      npt_total = ppiclf_iglsum(ppiclf_npart,1)

      ! Bin must be larger than nearest neighbor search distance
      ! and the ppiclf_filter(1:3).  This makes a buffer around the bin
      ! domain. Increase if you desire to bin less frequently.
      ! ppiclf_filter is 2x cell length in each direction
      DO i = 1,3
        BinMinLen(i) = MAX(ppiclf_filter(i),ppiclf_nndist)
        ! Need ppiclf_filter to make sure you have 1 layer
        ! of outer fluid cells
        ! Need ppiclf_nndist/2 to ensure BinMinLen is never violated
        BinBuffer(i) = MAX(ppiclf_filter(i)/2,ppiclf_nndist/2)
        MaxBinsAchieved(i) = .FALSE.
      END DO

      xmin =  1D10
      ymin =  1D10
      zmin =  1D10
      xmax = -1D10
      ymax = -1D10
      zmax = -1D10

      ! Looping through particles on this processor
      ! to find bin boundary locations
      DO i=1,ppiclf_npart
         ! Finding min/max particle extremes.
         ! Add buffer so that layers of outer cells 
         ! are available for interpolation/projection.
         temp1 = ppiclf_y(ix,i) - BinBuffer(ix)
         temp2 = ppiclf_y(ix,i) + BinBuffer(ix)
         IF(temp1 .LT. xmin) xmin = temp1
         IF(temp2 .GT. xmax) xmax = temp2

         temp1 = ppiclf_y(iy,i) - BinBuffer(iy)
         temp2 = ppiclf_y(iy,i) + BinBuffer(iy)
         IF(temp1 .LT. ymin) ymin = temp1
         IF(temp2 .GT. ymax) ymax = temp2

         temp1 = ppiclf_y(iz,i) - BinBuffer(iz)
         temp2 = ppiclf_y(iz,i) + BinBuffer(iz)
         IF(temp1 .LT. zmin) zmin = temp1
         IF(temp2 .GT. zmax) zmax = temp2
      END DO

      ! Finds global bin domain boundaries across MPI ranks
      ppiclf_binb(1) = ppiclf_glmin(xmin,1)
      ppiclf_binb(2) = ppiclf_glmax(xmax,1)
      ppiclf_binb(3) = ppiclf_glmin(ymin,1)
      ppiclf_binb(4) = ppiclf_glmax(ymax,1)
      ppiclf_binb(5) = ppiclf_glmin(zmin,1)
      ppiclf_binb(6) = ppiclf_glmax(zmax,1)

      CALL MPI_BARRIER(ppiclf_comm,ierr)

      ! If all particles within last RK Stage binbound, do not calculate
      ! bins again and do not remap overlap grid.
      BinCheck = 0
      DO i = 1,3
        IF((ppiclf_binb(2*i-1) + BinBuffer(i)) .LT.
     >                             ppiclf_previousbinb(2*i-1)) THEN
          BinCheck = 1
          EXIT
        END IF
        IF((ppiclf_binb(2*i)   - BinBuffer(i)) .GT.
     >                             ppiclf_previousbinb(2*i)) THEN
          BinCheck = 1
          EXIT
        END IF
      END DO

      CALL MPI_BARRIER(ppiclf_comm,ierr)

      BinCheck = ppiclf_iglmax(BinCheck,1)

#ifdef TEST
      ! So that CreateBin can be called repeatedly
      ! for different number of processors 
      BinCheck = 1
#endif

      IF(BinCheck .EQ. 0) THEN
        DO i = 1,3
          ppiclf_binb(2*i-1) = ppiclf_previousbinb(2*i-1)
          ppiclf_binb(2*i)   = ppiclf_previousbinb(2*i)
        END DO
        ppiclf_binchanged = .FALSE.
        RETURN
      ELSE
        ppiclf_binchanged  = .TRUE.
        ppiclf_printbinvtu = .TRUE.
      END IF

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

      ! Set previous bin bound for next RK Stage check
      DO i = 1,6
        ppiclf_previousbinb(i) = ppiclf_binb(i)
      END DO

      ! Return from subroutine if no particles present      
      IF(npt_total .LT. 1) RETURN

      ! Find ppiclf bin domain lengths
      DO i = 1,3
        binb_length(i) = ppiclf_binb(2*i) -
     >                         ppiclf_binb(2*i-1)
      END DO

      targetTotBin = ppiclf_np
      ! This sets the number of bin permutations to try to find optimal
      ! balance
      binPosBound = 1
      binNegBound = 1
      binIterations = binPosBound + binNegBound 
      DO i = 1,3
        MaxPotentialBins(i) = FLOOR(binb_length(i)/BinMinLen(i))
        IF(MaxPotentialBins(i) .LT. 1) THEN
          CALL ppiclf_exittr('BinMinLen() criteria violated.',0.0D0,0)
        END IF
      END DO

      ! Number of bins calculated based on bin surface
      ! area minimization and bin aspect ratio close to 1
      binsReal(1) = (targetTotBin**(1.0D0/3.0D0))*
     >              (binb_length(1)**(2.0D0/3.0D0))/ 
     >              ((binb_length(2)**(1.0D0/3.0D0))*
     >              (binb_length(3))**(1.0D0/3.0D0)) 
      
      binsReal(2) = (targetTotBin**(1.0D0/3.0D0))*
     >              (binb_length(2)**(2.0D0/3.0D0))/ 
     >              ((binb_length(1)**(1.0D0/3.0D0))*
     >              (binb_length(3))**(1.0D0/3.0D0)) 
     
      binsReal(3) = (targetTotBin**(1.0D0/3.0D0))*
     >              (binb_length(3)**(2.0D0/3.0D0))/ 
     >              ((binb_length(2)**(1.0D0/3.0D0))*
     >              (binb_length(1))**(1.0D0/3.0D0)) 

      ! The loop below ensures bins don't exceed number of
      ! processors or minimum bin length criteria, while
      ! maximizing the other dimensions.
      maxbincount2 = 0
      tempi = 0
      DO
        tempi = tempi + 1
        IF(tempi .GT. 100) THEN
          PRINT*, 'CreateBin stuck in infinite loop'
          CALL ppiclf_exittr('Error in Createbins',0.0,0)
        END IF
        ! Ensures most bins don't exceed number of 
        ! processors in case where one direction is
        ! much longer than the others.
        ! Ranking the bins Real numbers
        largeBin = 1 
        smallBin = 1
        medBin   = 1
        DO i = 1,3
          ! INT(x+0.5) is equivalent to ROUND(x). 
          ! Round isn't a built in fortran function.
          IF(binsReal(i) .LT. 0.50D0) binsReal(i) = 0.51D0
          ppiclf_n_bins(i) = INT(binsReal(i)+0.5D0)
          minBin(i) = ppiclf_n_bins(i) - binNegBound
          IF(minBin(i) .LT. 1) minBin(i) = 1
          IF(binsReal(i) .GT. binsReal(largeBin)) largeBin = i
          IF(binsReal(i) .LT. binsReal(smallBin)) smallBin = i
        END DO
        DO i = 1,3
          IF(i .EQ. largeBin) CYCLE
          IF(i .EQ. smallBin) CYCLE
          medBin = i
        END DO
        ! This checks to make sure at least one combination is within
        ! criteria of the largest possible total bins
        IF(minBin(1)*minBin(2)*minBin(3) .GT. targetTotBin) THEN
          IF(binsReal(medBin) .LT. 1) THEN
            ! Make largest dimension equal to number of processors.
            ! Other two dimensions are small
            DO i = 1,3
              IF(i .EQ. largeBin) THEN
                ppiclf_n_bins(largeBin) = targetTotBin 
                binsReal(largeBin) = REAL(ppiclf_n_bins(largeBin))
              ELSE
                ppiclf_n_bins(i) = 1
                binsReal(i) = 0.51D0
              END IF
            END DO
          ELSE
            ! Set bins of two large dimensions proportionally by length
            ! Set small dimension with 1 bins.
            binsReal(smallBin) = 0.51D0
            ppiclf_n_bins(smallBin) = 1

            binsReal(medBin) = SQRT(REAL(targetTotBin)) *
     >                    binb_length(medBin)/binb_length(largeBin)
            ppiclf_n_bins(medBin) = INT(binsReal(medBin)+0.5D0)

            binsReal(largeBin) = SQRT(REAL(targetTotBin)) *
     >                     binb_length(largeBin)/binb_length(medBin)
            ppiclf_n_bins(largeBin) = INT(binsReal(largeBin)+0.5D0)
          END IF 
        END IF
        ! Ensures bins don't violate minimum length criteria
        maxbincount1  = 0
        increaseRatio = 1.0D0
        DO i = 1,3
          IF(ppiclf_n_bins(i) .GE. MaxPotentialBins(i)) THEN
            ppiclf_n_bins(i) = MaxPotentialBins(i)
            ! This will be used to increase non-maximized bins
            increaseRatio  = (binsReal(i)/REAL(MaxPotentialBins(i)))* 
     >                        increaseRatio
            binsReal(i) = REAL(MaxPotentialBins(i))
            MaxBinsAchieved(i) = .TRUE.
            maxbincount1 = maxbincount1 + 1
          END IF
        END DO
        ! Either they both equal 0 and exit on the first pass,
        ! or they iterate and are equal on second or third pass.
        IF(maxbincount1 .EQ. maxbincount2) EXIT

        ! Increases other bins if one or two dimensions 
        ! are at maximum size based on minimum length criteria
        IF(MaxBinsAchieved(1) .OR. MaxBinsAchieved(2) .OR.
     >                             MaxBinsAchieved(3)) THEN
          maxbincount2 = 0
          DO i = 1,3
            IF(MaxBinsAchieved(i)) THEN
              maxbincount2 = maxbincount2 + 1
            END IF
          END DO
          DO i = 1,3
            IF(MaxBinsAchieved(i)) CYCLE
            IF(maxbincount2 .EQ. 1) THEN
              binsReal(i) = SQRT(increaseRatio)*binsReal(i)
            ELSE
              binsReal(i) = binsReal(i)*increaseRatio
            END IF
          END DO
        END IF
      END DO

      ! Since bin must be an integer, check -1:+1 number of bins
      ! for each bin dimension. Ideal number of bins will be max value
      ! while less than number of total target of bins.
      ! Minimize total error of sum of real bin calc - bin integer
      ! This will ensure larger dimensions get more bins.
      total_bin = 0
      minBinError = 1.0D9
      ideal_bin_index = 0
      DO ix = 0,binIterations
        iBin(1) = ppiclf_n_bins(1) + (ix-binNegBound)
        ppiclf_bins_dx(1) = binb_length(1)/iBin(1)
        IF(ppiclf_bins_dx(1) .LT. BinMinLen(1) .OR.
     >                           iBin(1) .LT. 1) CYCLE
        DO iy = 0,binIterations
          iBin(2) = ppiclf_n_bins(2) + (iy-binNegBound)
          ppiclf_bins_dx(2) = binb_length(2)/iBin(2)
          IF(ppiclf_bins_dx(2) .LT. BinMinLen(2) .OR.
     >                             iBin(2) .LT. 1) CYCLE
          DO iz = 0,binIterations
            iBin(3) = ppiclf_n_bins(3) + (iz-binNegBound)
            ppiclf_bins_dx(3) = binb_length(3)/iBin(3)
            IF(ppiclf_bins_dx(3) .LT. BinMinLen(3) .OR.
     >                               iBin(3) .LT. 1) CYCLE
            iBinTot = iBin(1)*iBin(2)*iBin(3)
            IF(iBinTot .GE. total_bin .AND. iBinTot .GT. 0 .AND.
     >                     iBinTot .LE. targetTotBin) THEN
              ! This resets the error min when a larger amount of
              ! total bins are found.
              IF(iBinTot .GT. total_bin) THEN
                minBinError = 1.0D9
              END IF
              ! If this is equal maximum # bins, find the best
              ! combination of bins per dimension
              binError = 0.0D0
              DO i = 1,3
                binError = binError + ABS(binsReal(i) - REAL(iBin(i))) 
              END DO
              IF(binError .LT. minBinError) THEN
                  minBinError = binError
                  ideal_bin_index(1) = ix
                  ideal_bin_index(2) = iy
                  ideal_bin_index(3) = iz
                  total_bin = iBinTot
              END IF            
            END IF
          END DO !iz
        END DO !iy
      END DO !ix
      IF(total_bin .EQ. 0 .AND. ppiclf_nid .EQ. 0) THEN
        PRINT*, 'correct bin combination not found'
      END IF

      tempi = 0
      total_bin = 1      
      DO i = 1,3
        ! Set common ppiclf bin arrays based on above calculation
        ! same number of bin equation as in above loops with best
        ! indices
        ppiclf_n_bins(i) = INT(binsReal(i)+0.5D0)
     >                     + (ideal_bin_index(i) - binNegBound)
        ppiclf_bins_dx(i) = binb_length(i)/ppiclf_n_bins(i)
        total_bin = total_bin*ppiclf_n_bins(i)
        IF(total_bin .GT. ppiclf_np) THEN
          IF(ppiclf_nid .EQ. 0) THEN
            PRINT*, 'binERROR: Num Bins > NumProcessors',total_bin,
     >            ppiclf_np, targetTotBin
          END IF
          CALL ppiclf_exittr('Error in Createbins',0.0,0)
        END IF
        IF(ppiclf_n_bins(i) .LT. 1) THEN
          PRINT*, 'binERROR: this dimension has negative bins', i,
     >            'Bins per dimension:', ppiclf_bins_dx(1),
     >            ppiclf_bins_dx(2), ppiclf_bins_dx(3)
        END IF
        IF(ppiclf_n_bins(i) .GT. tempi) THEN
          NBMax = i ! Dimension with max number of bins
          tempi = ppiclf_n_bins(i)
        END IF
      END DO
     
      ! Loop to see if we can add one to dimension with largest number of bins
      ! Choose this dimension because it is smallest incremental increase to total bins 
      DO
        IF((total_bin/ppiclf_n_bins(NBMax))*
     >      (ppiclf_n_bins(NBMax)+1) .LT. targetTotBin) THEN
          ! Add a bin and set new bin dx length
          ppiclf_n_bins(NBMax) = ppiclf_n_bins(NBMax)+1
          ppiclf_bins_dx(NBMax) = binb_length(NBMax)/
     >                              ppiclf_n_bins(NBMax)
          IF(ppiclf_bins_dx(NBMax) .LT. BinMinLen(NBMax)) THEN
            ! If BinMinLen criteria violated, return to previous bin configuration
            ppiclf_n_bins(NBMax) = ppiclf_n_bins(NBMax)-1
            ppiclf_bins_dx(NBMax) = binb_length(NBMax)/
     >                                ppiclf_n_bins(NBMax)
            EXIT
          END IF
          total_bin = 1
          DO i = 1,3
            total_bin = total_bin*ppiclf_n_bins(i)
          END DO
        ELSE
          EXIT
        END IF
      END DO

      ! Find this processor's x,y,z bin indicies
      idum = modulo(ppiclf_nid,ppiclf_n_bins(1))
      jdum = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
      kdum = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))

      ! Calculate this processor's bin min/max position in each dimension.
      ! Note that this value stays with this MPI rank.
      ppiclf_bin_pos(1,1) = ppiclf_binb(1) + idum    *ppiclf_bins_dx(1)
      ppiclf_bin_pos(2,1) = ppiclf_binb(1) + (idum+1)*ppiclf_bins_dx(1)
      ppiclf_bin_pos(1,2) = ppiclf_binb(3) + jdum    *ppiclf_bins_dx(2)
      ppiclf_bin_pos(2,2) = ppiclf_binb(3) + (jdum+1)*ppiclf_bins_dx(2)
      ppiclf_bin_pos(1,3) = ppiclf_binb(5) + kdum    *ppiclf_bins_dx(3)
      ppiclf_bin_pos(2,3) = ppiclf_binb(5) + (kdum+1)*ppiclf_bins_dx(3)

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_FindParticle
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
      INTEGER*4  i, ii, jj, kk, nrank, ierr, partcheck
      EXTERNAL   ppiclf_iglmax
      INTEGER*4  ppiclf_iglmax
!
      partcheck = 0
      DO i=1,ppiclf_npart
         ! Calculates particle's bin index
         ii  = FLOOR((ppiclf_y(1,i)-ppiclf_binb(1))/ppiclf_bins_dx(1))
         jj  = FLOOR((ppiclf_y(2,i)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
         kk  = FLOOR((ppiclf_y(3,i)-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
        
         ! Calculates particle's bin
         nrank  = ii + ppiclf_n_bins(1)*jj + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
         IF(nrank .NE. ppiclf_iprop(4,i)) partcheck = 1

         ! Maps particle to correct processor based on active bin number
         ! ***Use BinToProcMap for active/inactive bin***
         ppiclf_iprop(4,i) = nrank ! Processor to send to
         ppiclf_iprop(5,i) = ii    ! x bin #
         ppiclf_iprop(6,i) = jj    ! y bin #
         ppiclf_iprop(7,i) = kk    ! z bin #
         ppiclf_iprop(8,i) = nrank ! total bin number
      END DO
      ppiclf_particleMoved = ppiclf_iglmax(partcheck,1)
      CALL mpi_barrier(ppiclf_comm,ierr)

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_MoveParticle
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
      END

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_MapOverlapGrid
!
! This subroutine is called from ppiclf_solve_InitSolve
!
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
      INTEGER*4 nkey(2), i, j, k, l, ie, iee, ii, jj, kk, nrank,
     >          nl, nii, njj, nrr, iic, jjc, kkc il, ierr
      INTEGER*4 ix, iy, iz, ixLow, ixHigh, iyLow,
     >          iyHigh, izLow, izHigh 
      REAL*8    rxval, ryval, rzval, EleSizei(3), MaxPoint(3),
     >          MinPoint(3), ppiclf_vlmin, ppiclf_vlmax,
     >          centeri(3), exchCellMultiplier, Max_CellLen(3)
      LOGICAL   partl, ErrorFound
      EXTERNAL  ppiclf_vlmin, ppiclf_vlmax
#ifdef PERF
      REAL*8    tstart, tfinal
#endif
!
! Code Start:
!
      ! Number of fluid finite volume cells that map to particle bins
      ppiclf_nCells_FV2PICL = 0 
      
      ! Multiplies by x.4999 the cell length to ensure that x layers of cells
      ! outside of the ppiclf bin are mapped for interpolation.
      ! This could be changed based on the desired frequency of
      ! ppiclf bin creation and mapping.
      exchCellMultiplier  = 2.499D0

      ! Loops through number of fluid FV cells on this processor
      DO ie=1,ppiclf_nFVCells  
        DO l=1,3
         !indicies 1:3
         centeri(l) = ppiclf_fluid_grid(l,ie)
         !indicies 4:6
         EleSizei(l) =  exchCellMultiplier*ppiclf_fluid_grid(3+l,ie)
         IF(EleSizei(l) .GT. ppiclf_bins_dx(l)/2) THEN
           EleSizei(l) =  1.499D0*ppiclf_fluid_grid(3+l,ie)
         END IF
        END DO !l 

        ! Fluid Cell vertex position without additional length
        rxval = centeri(1)
        ryval = centeri(2)
        rzval = centeri(3)
      
        ! Exits if fluid cell center is outside of any bin boundaries 
        IF (rxval .GT. ppiclf_binb(2)) CYCLE
        IF (rxval .LT. ppiclf_binb(1)) CYCLE
        IF (ryval .GT. ppiclf_binb(4)) CYCLE
        IF (ryval .LT. ppiclf_binb(3)) CYCLE
        IF (rzval .GT. ppiclf_binb(6)) CYCLE
        IF (rzval .LT. ppiclf_binb(5)) CYCLE
 
        ! Determines what bin the fluid cell is nominally mapped to
        ii    = FLOOR((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
        jj    = FLOOR((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
        kk    = FLOOR((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3))

        ! Default is Do loop with ix=iy=iz=2 for fluid cells not near
        ! bin boundary

        ixLow =2
        ixHigh=2
        iyLow =2
        iyHigh=2
        izLow =2
        izHigh=2

        ! These series of if statements check if bin mapping changes
        ! when adding/subtracting multiple of fluid cell length defined
        ! by EleSizei(l). 
        ! This is used to map fluid cells slightly outside of the ppiclf
        ! bin boundary.  If any .NE. 2, then fluid cell is mapped to
        ! multiple ppiclf bins. 
        
        IF(FLOOR((rxval + EleSizei(1) - ppiclf_binb(1))
     >       /ppiclf_bins_dx(1)) .NE. ii)  ixHigh = 3

        IF(FLOOR((rxval - EleSizei(1) - ppiclf_binb(1))
     >       /ppiclf_bins_dx(1)) .NE. ii)  ixLow = 1
        IF(FLOOR((ryval + EleSizei(2) - ppiclf_binb(3))
     >       /ppiclf_bins_dx(2)) .NE. jj)  iyHigh = 3

        IF(FLOOR((ryval - EleSizei(2) - ppiclf_binb(3))
     >       /ppiclf_bins_dx(2)) .NE. jj)  iyLow = 1

        IF(ppiclf_ndim .GT. 2 .AND. FLOOR((rzval + EleSizei(3)
     >    - ppiclf_binb(5))/ppiclf_bins_dx(3)) .NE. kk)  izHigh = 3

        IF(ppiclf_ndim .GT. 2 .AND. FLOOR((rzval - EleSizei(3)
     >    - ppiclf_binb(5))/ppiclf_bins_dx(3)) .NE. kk)  izLow = 1

        DO ix=ixLow,ixHigh
          DO iy=iyLow,iyHigh
            DO iz=izLow,izHigh
              ! Change cell position by EleSizei if ix,iy,or iz NE 2
              rxval = centeri(1) + (ix-2)*EleSizei(1)
              ryval = centeri(2) + (iy-2)*EleSizei(2)
              rzval = centeri(3) + (iz-2)*EleSizei(3)
              ! Find bin for adjusted rval
              ii    = FLOOR((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
              jj    = FLOOR((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
              kk    = FLOOR((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
              

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


              ! Calculates processor rank
              nrank  = ii + ppiclf_n_bins(1)*jj + 
     >                     ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
              ppiclf_nCells_FV2PICL = ppiclf_nCells_FV2PICL + 1
              IF(ppiclf_nCells_FV2PICL .GT. PPICLF_LEE) THEN
                PRINT*, '***ERROR*** PPICLF_LEE',PPICLF_LEE, 'in', 
     >           'MapOverlapGrid must be greater than',
     >            ppiclf_nCells_FV2PICL 
                CALL ppiclf_exittr('Increase PPICLF_LEE$ (MapOverlap)',0.0D0
     >               ,ppiclf_nCells_FV2PICL)
              END IF

              ! make sure it is mapped to active nrank and map rank to
              ! processor. *** FOR ACTIVE BINNING ***

              ! Stores cells to rank mapping.
              ! Fluid solver cell ID
              ppiclf_cell_map(1,ppiclf_nCells_FV2PICL) = ie
              ! Fluid solver cell rank
              ppiclf_cell_map(2,ppiclf_nCells_FV2PICL) = ppiclf_nid
              ! Particle solver cell rank and bin indicies
              ppiclf_cell_map(3,ppiclf_nCells_FV2PICL) = nrank
              ppiclf_cell_map(4,ppiclf_nCells_FV2PICL) = ii
              ppiclf_cell_map(5,ppiclf_nCells_FV2PICL) = jj
              ppiclf_cell_map(6,ppiclf_nCells_FV2PICL) = kk
            END DO !iz
          END DO !iy
        END DO !ix
      END DO !ie

      DO ie=1,ppiclf_nCells_FV2PICL 
        ! These copy all indicies since Fortran is column-major
        iee = ppiclf_cell_map(1,ie)
        CALL ppiclf_copy(ppiclf_picl_grid(1,ie)
     >                 ,ppiclf_fluid_grid(1,iee),7)
 
        ! ppiclf_filter initially set in PICL_TEMP_InitSolver
        ! Want to only consider cells that reside in the particle domain
        ! Update ppiclf_filter for next binning cycle 2.1*dx since next
        ! layer of cells in a growing particle domain may be slightly larger
        ! and we want filter equal a minimum of 2 cells.
        DO l = 1,3
          IF(ppiclf_filter(l) .LT. ppiclf_fluid_grid(3+l,iee)) THEN
            ppiclf_filter(l) = 2.1D0*ppiclf_fluid_grid(3+l,iee)
          END IF
        END DO
      END DO

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
      END 
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_CreateGhost
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
      END

!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_GhostDistCheck(ix,Pos,distchk,
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
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_LinearPeriodicityGhost(iig,l,Pos,PerShift)
      
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
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_MoveGhost
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
      END

!-----------------------------------------------------------------------
! The following subroutines are for the new particle-based load balance
! approach. The goal is to make ppiclf_npart approximately equal per
! processor.
!-----------------------------------------------------------------------

      SUBROUTINE ppiclf_comm_CreateBinPartLB

      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"
 
      INTEGER*4 i, ix, iy, iz, ierr, ppiclf_iglsum
      EXTERNAL  ppiclf_iglsum

      REAL*8    BinMinLen(3), xmin, ymin, zmin, xmax, ymax, zmax
     >          ,BinBuffer(3), temp1, temp2, periodicDistCheck
     >          ,idum, jdum, kdum, ppiclf_glmin, ppiclf_glmax 
      EXTERNAL  ppiclf_glmin, ppiclf_glmax 


      ix = 1
      iy = 2
      iz = 3
        
      ! iglsum is integer global sum across MPI ranks.
      ppiclf_glnpart = ppiclf_iglsum(ppiclf_npart,1)

      ! Bin must be larger than nearest neighbor search distance
      ! and the ppiclf_filter(1:3).  This makes a buffer around the bin
      ! domain. Increase if you desire to bin less frequently.
      ! ppiclf_filter is 2x cell length in each direction
      DO i = 1,3
        BinMinLen(i) = MAX(ppiclf_filter(i),ppiclf_nndist)
        ! Need ppiclf_filter to make sure you have 1 layer
        ! of outer fluid cells
        ! Need ppiclf_nndist/2 to ensure BinMinLen is never violated
        BinBuffer(i) = MAX(ppiclf_filter(i)/2,ppiclf_nndist/2)
      END DO

      xmin =  1D10
      ymin =  1D10
      zmin =  1D10
      xmax = -1D10
      ymax = -1D10
      zmax = -1D10

      ! Looping through particles on this processor
      ! to find bin boundary locations
      DO i=1,ppiclf_npart
         ! Finding min/max particle extremes.
         ! Add buffer so that layers of outer cells 
         ! are available for interpolation/projection.
         temp1 = ppiclf_y(ix,i) - BinBuffer(ix)
         temp2 = ppiclf_y(ix,i) + BinBuffer(ix)
         IF(temp1 .LT. xmin) xmin = temp1
         IF(temp2 .GT. xmax) xmax = temp2

         temp1 = ppiclf_y(iy,i) - BinBuffer(iy)
         temp2 = ppiclf_y(iy,i) + BinBuffer(iy)
         IF(temp1 .LT. ymin) ymin = temp1
         IF(temp2 .GT. ymax) ymax = temp2

         temp1 = ppiclf_y(iz,i) - BinBuffer(iz)
         temp2 = ppiclf_y(iz,i) + BinBuffer(iz)
         IF(temp1 .LT. zmin) zmin = temp1
         IF(temp2 .GT. zmax) zmax = temp2
      END DO

      ! Finds global bin domain boundaries across MPI ranks
      ppiclf_binb(1) = ppiclf_glmin(xmin,1)
      ppiclf_binb(2) = ppiclf_glmax(xmax,1)
      ppiclf_binb(3) = ppiclf_glmin(ymin,1)
      ppiclf_binb(4) = ppiclf_glmax(ymax,1)
      ppiclf_binb(5) = ppiclf_glmin(zmin,1)
      ppiclf_binb(6) = ppiclf_glmax(zmax,1)

      CALL MPI_BARRIER(ppiclf_comm,ierr)

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
        ppiclf_n_bins(i) = INT( (ppiclf_binb(2)-ppiclf_binb(1)) /
     >                          BinMinLen(i) )
        ppiclf_bins_dx(i) = (ppiclf_binb(2) - ppiclf_binb(1)) / 
     >                      ppiclf_n_bins(i)
        ppiclf_totalBins = ppiclf_totalBins*ppiclf_n_bins(i)
      END DO

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_FindParticlePartLB

      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      INTEGER*4  i, ii, jj, kk, nbin, ierr, partcheck, NumBins
     >           ,ParticleCount(0:ppiclf_totalBins-1)
     >           ,BinToRankMapping(0:ppiclf_totalBins-1)
      EXTERNAL   ppiclf_iglmax
      INTEGER*4  ppiclf_iglmax

      DO i = 0,(ppiclf_totalBins - 1)
        ParticleCount(i) = 0
      END DO

      DO i=1,ppiclf_npart
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
        ParticleCount(nbin) = ParticleCount(nbin) + 1
      END DO

      CALL ppiclf_comm_partLoadBalance(ParticleCount,BinToRankMapping)

      RETURN
      END
!----------------------------------------------------------------------
      SUBROUTINE ppiclf_comm_partLoadBalance(PC,BTRM)
   
      IMPLICIT NONE

      INCLUDE "PPICLF"
      INCLUDE "mpif.h"

      INTEGER*4  PC(0:ppiclf_totalBins-1), ierr, i
     >           ,BTRM(0:ppiclf_totalBins-1)
     >           ,targetParticleCnt, particleSum
     >           ,irank, dL, bin
     >           ,d2, d3, ii, jj, kk, nb1, nb2, nb3
     >           ,islice
      REAL*8     largestDomainDim

      ! Find the dimension to slice along. Using largest dimension
      ! minimizes surface area between processors, which likely
      ! minimizes ghost particle and overlap cell communication

      !dL = largest dimension
      dL = 1
      largestDomainDim = ppiclf_binb(2) - ppiclf_binb(1)
      DO i = 2,3
        IF(ppiclf_binb(2*i)-ppiclf_binb(2*i-1) .GT. 
     >                                        largestDomainDim) THEN
          dL = i
          largestDomainDim = ppiclf_binb(2*i) - ppiclf_binb(2*i-1)
        END IF
      END DO

      nb1 = ppiclf_n_bins(1)
      nb2 = ppiclf_n_bins(2)
      nb3 = ppiclf_n_bins(3)

      ! The following code now determines 
      ! Now sum particles per bin across MPI Ranks
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, PC, ppiclf_totalBins
     >                   ,MPI_INTEGER4, MPI_SUM
     >                   ,ppiclf_comm, ierr)
      targetParticleCnt = INT(REAL(ppiclf_glnpart)/REAL(ppiclf_np)) + 1

      particleSum   = 0
      irank = 0
      ! Iterate through slices of largest dimension.
      ! Increment slice by one and loop through all other dimensions.
      DO islice = 0,(ppiclf_n_bins(dL) - 1)
        IF(dL .EQ. 1) THEN
          ii = islice
          DO jj = 0, nb2-1
            DO kk = 0, nb3-1
              bin  = ii + nb1*jj + nb1*nb2*kk
              particleSum = particleSum + PC(bin)
              ! This maps the bin to the MPI-rank
              BTRM(bin) = irank 
            END DO
          END DO
        ELSEIF(dL .EQ. 2) THEN
          jj = islice
          DO ii = 0, (nb1 - 1)
            DO kk = 0, (nb3 - 1)
              bin  = ii + nb1*jj + nb1*nb2*kk
              particleSum = particleSum + PC(bin)
              ! This maps the bin to the MPI-rank
              BTRM(bin) = irank 
            END DO
          END DO
        ELSEIF(dL .EQ. 3) THEN
          kk = islice
          DO ii = 0, (nb1 - 1)
            DO jj = 0, (nb2 - 1)
              bin  = ii + nb1*jj + nb1*nb2*kk
              particleSum = particleSum + PC(bin)
              ! This maps the bin to the MPI-rank
              BTRM(bin) = irank 
            END DO
          END DO
        ELSE
          PRINT*, 'ERROR in ppiclf_comm_partLoadBalance'
          RETURN
        END IF

        IF(particleSum .GE. targetParticleCnt) THEN
          ! Can also write MPI-rank domain boundary here
          particleSum = 0
          irank = irank + 1
          irank = MIN(irank,ppiclf_np-1)
        END IF
      END DO

      RETURN
      END
!----------------------------------------------------------------------
!      subroutine ppiclf_comm_AngularCreateGhost
!!
!      implicit none
!!
!      include "PPICLF"
!!
!! Internal:
!!
!      real*8 xdlen,ydlen,zdlen,rxdrng(3),rxnew(3), rfac, rxval, ryval,
!     >       rzval, rxl, ryl, rzl, rxr, ryr, rzr, distchk, dist, gFilt
!      integer*4 iadd(3),gpsave(27)
!      real*8 map(PPICLF_LRP_PRO)
!      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
!     >           nfacegp, nedgegp, ncornergp, iperiodicx, iperiodicy,
!     >           iperiodicz, jx, jy, jz, ip, idum, iip, jjp, kkp, ii1,
!     >           jj1, kk1, iig, jjg, kkg, iflgx, iflgy, iflgz,
!     >           isave, iflgsum, ndumn, nrank, ibctype, i, ifc, ist, j,
!     >           k
!      ! 08/27/24 - Thierry - added for angular periodicty starts here
!      real*8 alpha
!      integer*4 xrank, yrank, zrank
!      ! 08/27/24 - Thierry - added for angular periodicty ends here
!      ! 09/26/24 - Thierry - added for angular periodicty starts here
!      real*8 dist1, dist2
!      ! 09/26/24 - Thierry - added for angular periodicty ends here
!!
!
!c     face, edge, and corner number, x,y,z are all inline, so stride=3
!      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
!      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
!     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
!     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
!      el_corner_num = (/ -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
!     >                   -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)
!
!      nfacegp   = 4  ! number of faces
!      nedgegp   = 4  ! number of edges
!      ncornergp = 0  ! number of corners
!
!      if (ppiclf_ndim .gt. 2) then
!         nfacegp   = 6  ! number of faces
!         nedgegp   = 12 ! number of edges
!         ncornergp = 8  ! number of corners
!      endif
!
!      iperiodicx = ppiclf_iperiodic(1)
!      iperiodicy = ppiclf_iperiodic(2)
!      iperiodicz = ppiclf_iperiodic(3)
!
!! ------------------------
!c CREATING GHOST PARTICLES
!! ------------------------
!      jx    = 1
!      jy    = 2
!      jz    = 3
!
!      ! Thierry - we dont use xdlen and ydlen in this algorithm. no need to modify them.
!      xdlen = ppiclf_binb(2) - ppiclf_binb(1)
!      ydlen = ppiclf_binb(4) - ppiclf_binb(3)
!      zdlen = -1.
!      if (ppiclf_ndim .gt. 2) 
!!     >   zdlen = ppiclf_binb(6) - ppiclf_binb(5)
!     >   zdlen = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
!      if (iperiodicx .ne. 0) xdlen = -1
!      if (iperiodicy .ne. 0) ydlen = -1
!      if (iperiodicz .ne. 0) zdlen = -1
!
!      rxdrng(1) = xdlen
!      rxdrng(2) = ydlen
!      rxdrng(3) = zdlen
!
!      ppiclf_npart_gp = 0
!
!      rfac = 1.0d0
!      gFilt = MAX(ppiclf_nndist,ppiclf_filter(1),
!     >            ppiclf_filter(2),ppiclf_filter(3))
!
!
!      do ip=1,ppiclf_npart
!
!         call ppiclf_user_MapProjPart(map,ppiclf_y(1,ip)
!     >         ,ppiclf_ydot(1,ip),ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip))
!
!c        idum = 1
!c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
!c        idum = 2
!c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
!c        idum = 3
!c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
!
!         idum = 0
!         do j=1,PPICLF_LRS
!            idum = idum + 1
!            ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip) ! ppiclf_y(PPICLF_JX/ JY/ JZ/ JVX/ JVY/ JVZ/ JT, ip)
!         enddo
!         idum = PPICLF_LRS
!         do j=1,PPICLF_LRP
!            idum = idum + 1
!            ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip) ! ppiclf_rprop(PPICLF_R_JRHOP/ R_JRHOF/ .../ R_WDOTZ, ip)
!         enddo
!         idum = PPICLF_LRS+PPICLF_LRP
!         do j=1,PPICLF_LRP_PRO
!            idum = idum + 1
!            ppiclf_cp_map(idum,ip) = map(j) ! map(PPICLF_P_JPHIP/ JFX/ .../ JPHIPW) - these are found in ppiclf_user_MapProjPart
!         enddo
!
!         rxval = ppiclf_cp_map(1,ip) ! ppiclf_y(PPICLF_JX,ip)
!         ryval = ppiclf_cp_map(2,ip) ! ppiclf_y(PPICLF_JY,ip)
!         rzval = 0.0d0
!         if (ppiclf_ndim .gt. 2) rzval = ppiclf_cp_map(3,ip) ! ppiclf_y(PPICLF_JZ,ip)
!
!         iip    = ppiclf_iprop(4,ip) ! ith coordinate of bin
!         jjp    = ppiclf_iprop(5,ip) ! jth coordinate of bin
!         kkp    = ppiclf_iprop(6,ip) ! kth coordinate of bin
!
!         rxl = ppiclf_binb(1) + ppiclf_bins_dx(1)*iip ! min x of bin
!         rxr = rxl + ppiclf_bins_dx(1)                ! max x of bin
!         ryl = ppiclf_binb(3) + ppiclf_bins_dx(2)*jjp
!         ryr = ryl + ppiclf_bins_dx(2)
!         rzl = 0.0d0
!         rzr = 0.0d0
!         if (ppiclf_ndim .gt. 2) then
!            rzl = ppiclf_binb(5) + ppiclf_bins_dx(3)*kkp
!            rzr = rzl + ppiclf_bins_dx(3)
!         endif
!
!         isave = 0
!
!         ! faces
!         do ifc=1,nfacegp
!            ist = (ifc-1)*3
!            ii1 = iip + el_face_num(ist+1) 
!            jj1 = jjp + el_face_num(ist+2)
!            kk1 = kkp + el_face_num(ist+3)
!
!            iig = ii1
!            jjg = jj1
!            kkg = kk1
!
!            distchk = 0.0d0
!            dist = 0.0d0
!            if (ii1-iip .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
!               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
!            endif
!            if (jj1-jjp .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
!               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
!            endif
!            if (ppiclf_ndim .gt. 2) then
!            if (kk1-kkp .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
!               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
!            endif
!            endif
!            distchk = sqrt(distchk)
!            dist = sqrt(dist)
!
!            if (ang_case==1) then  ! for wedge geometry
!
!               ! Thierry - I dont think it's code efficient to call this subroutine
!               !           for every particle, every ghost face, at every time step
!               !           I'm wondering if it's better if we make the plane values 
!               !           as global values that are initialized in the beginning 
!            
!               call ppiclf_solve_InitAngularPlane(ip,
!     >                                 ang_per_rin  , ang_per_rout  ,
!     >                                 ang_per_angle, ang_per_xangle,
!     >                                 dist1, dist2)
!               if ((dist .gt. distchk).and.(dist1.gt.distchk)
!     >           .and.(dist2.gt.distchk)) cycle
!            else
!               if (dist .gt. distchk) cycle
!            endif
!
!            iflgx = 0
!            iflgy = 0
!            iflgz = 0
!!-----------------------------------------------------------------------
!            ! 08/27/24 - Thierry - modification for angular periodicty starts here
!
!               ! angle between particle and x-axis
!                alpha = atan2(ppiclf_y(PPICLF_JY,ip), 
!     >                        ppiclf_y(PPICLF_JX,ip))
!                
!
!                call ppiclf_solve_InvokeAngularPeriodic(ip, 
!     >                                                  ang_per_flag,
!     >                                                  alpha,         
!     >                                                  ang_per_angle,  
!     >                                                  ang_per_xangle, 
!     >                                                  0)
!
!              ! Thierry - this is how FindParticle implements it
!              ! need to find a way to make the code deal with negative xrot values
!
!            xrank = iig ; yrank=jjg; zrank = kkg
!            ! Thierry - previously placed before the CheckPeriodicBC call, had to move them for the periodic check
!            iadd(1) = ii1
!            iadd(2) = jj1
!            iadd(3) = kk1
!            rxnew(1) = rxval
!            rxnew(2) = ryval
!            rxnew(3) = rzval ! z-coordinate does not change when angular periodicity is invoked
!            
!            ! Angular periodicity check in x- and y-directions
!            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
!              iflgx = 1
!              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
!              if (iperiodicx .ne. 0) cycle
!              iig = xrank
!              jjg = yrank
!            end if
!            
!            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
!              iflgy = 1
!              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
!              if (iperiodicy .ne. 0) cycle
!              iig = xrank
!              jjg = yrank
!            end if
!            
!            ! Linear periodicity check in z-direction
!            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
!              iflgz = 1
!              kkg =modulo(kkg,ppiclf_n_bins(3))
!              if (iperiodicz .ne. 0) cycle
!              ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
!              ! rxdrng(3) = -1.0  if not periodic in Z
!              if (rxdrng(3) .gt. 0) then 
!                if (iadd(3) .ge. ppiclf_n_bins(3)) then ! particle leaving from max z-face
!                  rxnew(3) = rxnew(3) - rxdrng(3)
!                elseif (iadd(3) .lt. 0) then ! particle leaving from min z-face
!                  rxnew(3) = rxnew(3) + rxdrng(3)
!                end if ! iadd
!              end if ! rxrdrng
!            else ! z-linear periodicity not applicable
!              kkg = zrank
!            end if ! kkg
!            
!            iflgsum = iflgx + iflgy + iflgz
!            ndumn  = iig + ppiclf_n_bins(1)*jjg + 
!     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
!             nrank = ndumn
!
!            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle
!
!            ! 08/27/24 - Thierry - modification for angular periodicty ends here
!!-----------------------------------------------------------------------
!
!            do i=1,isave
!               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 111
!            enddo
!            isave = isave + 1
!            gpsave(isave) = nrank
!
!            ibctype = iflgx+iflgy+iflgz
!            
!            rxnew(1) = xrot(1)
!            rxnew(2) = xrot(2)
!            ppiclf_cp_map(4,ip) = vrot(1)
!            ppiclf_cp_map(5,ip) = vrot(2)
!                 
!            ppiclf_npart_gp = ppiclf_npart_gp + 1
!            ppiclf_iprop_gp(1,ppiclf_npart_gp) = ppiclf_iprop(1,ip)
!            ppiclf_iprop_gp(2,ppiclf_npart_gp) = ppiclf_iprop(2,ip)
!            ppiclf_iprop_gp(3,ppiclf_npart_gp) = nrank
!            ppiclf_iprop_gp(4,ppiclf_npart_gp) = iig
!            ppiclf_iprop_gp(5,ppiclf_npart_gp) = jjg
!            ppiclf_iprop_gp(6,ppiclf_npart_gp) = kkg
!            ppiclf_iprop_gp(7,ppiclf_npart_gp) = nrank
!
!            ! Thierry - we don't need ppiclf_comm_CheckPeriodicBC anymore for the angular periodic ghost algorithm
!            !           as this is now taken care of when anticipating where the particle might be when calling
!            !           ppiclf_comm_InvokeAngularPeriodic
!            !           we only need to assign xr and vr to ppiclf_rprop_gp
!
!            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1) ! ppiclf_y(PPICLF_JX, ip) for the periodic ghost particle
!            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2) ! JY
!            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3) ! JZ
!            
!            do k=4,PPICLF_LRP_GP
!               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
!            enddo
!  111 continue
!         enddo
!
!         ! edges
!         do ifc=1,nedgegp
!            ist = (ifc-1)*3
!            ii1 = iip + el_edge_num(ist+1) 
!            jj1 = jjp + el_edge_num(ist+2)
!            kk1 = kkp + el_edge_num(ist+3)
!
!            iig = ii1
!            jjg = jj1
!            kkg = kk1
!
!            distchk = 0.0d0
!            dist = 0.0d0
!            if (ii1-iip .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
!               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
!            endif
!            if (jj1-jjp .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
!               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
!            endif
!            if (ppiclf_ndim .gt. 2) then
!            if (kk1-kkp .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
!               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
!            endif
!            endif
!            distchk = sqrt(distchk)
!            dist = sqrt(dist)
!
!            if (ang_case==1) then  ! for wedge geometry
!
!               call ppiclf_solve_InitAngularPlane(ip,
!     >                                 ang_per_rin  , ang_per_rout  ,
!     >                                 ang_per_angle, ang_per_xangle,
!     >                                 dist1, dist2)
!               if ((dist .gt. distchk).and.(dist1.gt.distchk)
!     >           .and.(dist2.gt.distchk)) cycle
!            else
!               if (dist .gt. distchk) cycle
!            endif
!
!            iflgx = 0
!            iflgy = 0
!            iflgz = 0
!            ! periodic if out of domain - add some ifsss
!!-----------------------------------------------------------------------
!            ! 08/27/24 - Thierry - modification for angular periodicty starts here
!
!               ! angle between particle and x-axis
!                alpha = atan2(ppiclf_y(PPICLF_JY,ip), 
!     >                        ppiclf_y(PPICLF_JX,ip))
!                
!
!                call ppiclf_solve_InvokeAngularPeriodic(ip, 
!     >                                                  ang_per_flag,
!     >                                                  alpha,         
!     >                                                  ang_per_angle,  
!     >                                                  ang_per_xangle, 
!     >                                                  0)
!
!              ! Thierry - this is how FindParticle implements it
!              ! need to find a way to make the code deal with negative xrot values
!
!            xrank = iig ; yrank=jjg; zrank = kkg
!            ! Thierry - previously placed before the CheckPeriodicBC call, had to move them for the periodic check
!            iadd(1) = ii1
!            iadd(2) = jj1
!            iadd(3) = kk1
!            rxnew(1) = rxval
!            rxnew(2) = ryval
!            rxnew(3) = rzval ! z-coordinate does not change when angular periodicity is invoked
!            
!            ! Angular periodicity check in x- and y-directions
!            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
!              iflgx = 1
!              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
!              if (iperiodicx .ne. 0) cycle
!              iig = xrank
!              jjg = yrank
!            end if
!            
!            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
!              iflgy = 1
!              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
!              if (iperiodicy .ne. 0) cycle
!              iig = xrank
!              jjg = yrank
!            end if
!            
!            ! Linear periodicity check in z-direction
!            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
!              iflgz = 1
!              kkg =modulo(kkg,ppiclf_n_bins(3))
!              if (iperiodicz .ne. 0) cycle
!              ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
!              ! rxdrng(3) = -1.0  if not periodic in Z
!              if (rxdrng(3) .gt. 0) then ! particle leaving from max z-face
!                if (iadd(3) .ge. ppiclf_n_bins(3)) then
!                  rxnew(3) = rxnew(3) - rxdrng(3)
!                elseif (iadd(3) .lt. 0) then
!                  rxnew(3) = rxnew(3) + rxdrng(3)
!                end if ! iadd
!              end if ! rxrdrng
!            else ! z-linear periodicity not applicable
!              kkg = zrank
!            end if ! kkg
!            
!            iflgsum = iflgx + iflgy + iflgz
!            ndumn  = iig + ppiclf_n_bins(1)*jjg + 
!     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
!             nrank = ndumn
!
!            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle
!
!            ! 08/27/24 - Thierry - modification for angular periodicty ends here
!!-----------------------------------------------------------------------
!
!            do i=1,isave
!               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 222
!            enddo
!            isave = isave + 1
!            gpsave(isave) = nrank
!
!            ibctype = iflgx+iflgy+iflgz
!
!            rxnew(1) = xrot(1)
!            rxnew(2) = xrot(2)
!            ppiclf_cp_map(4,ip) = vrot(1)
!            ppiclf_cp_map(5,ip) = vrot(2)
!                 
!            ppiclf_npart_gp = ppiclf_npart_gp + 1
!
!            ppiclf_iprop_gp(1,ppiclf_npart_gp) = ppiclf_iprop(1,ip)
!            ppiclf_iprop_gp(2,ppiclf_npart_gp) = ppiclf_iprop(2,ip)
!            ppiclf_iprop_gp(3,ppiclf_npart_gp) = nrank
!            ppiclf_iprop_gp(4,ppiclf_npart_gp) = iig
!            ppiclf_iprop_gp(5,ppiclf_npart_gp) = jjg
!            ppiclf_iprop_gp(6,ppiclf_npart_gp) = kkg
!            ppiclf_iprop_gp(7,ppiclf_npart_gp) = nrank
!
!            ! Thierry - we don't need ppiclf_comm_CheckPeriodicBC anymore for the angular periodic ghost algorithm
!            !           as this is now taken care of when anticipating where the particle might be when calling
!            !           ppiclf_comm_InvokeAngularPeriodic
!            !           we only need to assign xr and vr to ppiclf_rprop_gp
!
!            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1) ! ppiclf_y(PPICLF_JX, ip) for the periodic ghost particle
!            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2) ! JY
!            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3) ! JZ
!            
!            do k=4,PPICLF_LRP_GP
!               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
!            enddo
!  222 continue
!         enddo
!
!         ! corners
!         do ifc=1,ncornergp
!            ist = (ifc-1)*3
!            ii1 = iip + el_corner_num(ist+1) 
!            jj1 = jjp + el_corner_num(ist+2)
!            kk1 = kkp + el_corner_num(ist+3)
!
!            iig = ii1
!            jjg = jj1
!            kkg = kk1
!
!            distchk = 0.0d0
!            dist = 0.0d0
!            if (ii1-iip .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
!               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
!            endif
!            if (jj1-jjp .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
!               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
!            endif
!            if (ppiclf_ndim .gt. 2) then
!            if (kk1-kkp .ne. 0) then
!               distchk = distchk + (rfac*gFilt)**2
!               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
!               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
!            endif
!            endif
!            distchk = sqrt(distchk)
!            dist = sqrt(dist)
!
!            if (ang_case==1) then  ! for wedge geometry
!            
!               call ppiclf_solve_InitAngularPlane(ip,
!     >                                 ang_per_rin  , ang_per_rout  ,
!     >                                 ang_per_angle, ang_per_xangle,
!     >                                 dist1, dist2)
!               if ((dist .gt. distchk).and.(dist1.gt.distchk)
!     >           .and.(dist2.gt.distchk)) cycle
!            else
!               if (dist .gt. distchk) cycle
!            endif
!
!            iflgx = 0
!            iflgy = 0
!            iflgz = 0
!
!!-----------------------------------------------------------------------
!            ! 08/27/24 - Thierry - modification for angular periodicty starts here
!
!               ! angle between particle and x-axis
!                alpha = atan2(ppiclf_y(PPICLF_JY,ip), 
!     >                        ppiclf_y(PPICLF_JX,ip))
!                
!
!                call ppiclf_solve_InvokeAngularPeriodic(ip, 
!     >                                                  ang_per_flag,
!     >                                                  alpha,         
!     >                                                  ang_per_angle,  
!     >                                                  ang_per_xangle, 
!     >                                                  0)
!
!              ! Thierry - this is how FindParticle implements it
!              ! need to find a way to make the code deal with negative xrot values
!
!            xrank = iig ; yrank=jjg; zrank = kkg
!            ! Thierry - previously placed before the CheckPeriodicBC call, had to move them for the periodic check
!            iadd(1) = ii1
!            iadd(2) = jj1
!            iadd(3) = kk1
!            rxnew(1) = rxval
!            rxnew(2) = ryval
!            rxnew(3) = rzval ! z-coordinate does not change when angular periodicity is invoked
!            
!            ! Angular periodicity check in x- and y-directions
!            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
!              iflgx = 1
!              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
!              if (iperiodicx .ne. 0) cycle
!              iig = xrank
!              jjg = yrank
!            end if
!            
!            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
!              iflgy = 1
!              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
!              if (iperiodicy .ne. 0) cycle
!              iig = xrank
!              jjg = yrank
!            end if
!            
!            ! Linear periodicity check in z-direction
!            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
!              iflgz = 1
!              kkg =modulo(kkg,ppiclf_n_bins(3))
!              if (iperiodicz .ne. 0) cycle
!              ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
!              ! rxdrng(3) = -1.0  if not periodic in Z
!              if (rxdrng(3) .gt. 0) then ! particle leaving from max z-face
!                if (iadd(3) .ge. ppiclf_n_bins(3)) then
!                  rxnew(3) = rxnew(3) - rxdrng(3)
!                elseif (iadd(3) .lt. 0) then
!                  rxnew(3) = rxnew(3) + rxdrng(3)
!                end if ! iadd
!              end if ! rxrdrng
!            else ! z-linear periodicity not applicable
!              kkg = zrank
!            end if ! kkg
!            
!            iflgsum = iflgx + iflgy + iflgz
!            ndumn  = iig + ppiclf_n_bins(1)*jjg + 
!     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
!             nrank = ndumn
!
!            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle
!
!            ! 08/27/24 - Thierry - modification for angular periodicty ends here
!!-----------------------------------------------------------------------
!            do i=1,isave
!               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 333
!            enddo
!            isave = isave + 1
!            gpsave(isave) = nrank
!
!            ibctype = iflgx+iflgy+iflgz
!
!            rxnew(1) = xrot(1)
!            rxnew(2) = xrot(2)
!            ppiclf_cp_map(4,ip) = vrot(1)
!            ppiclf_cp_map(5,ip) = vrot(2)
!
!            ppiclf_npart_gp = ppiclf_npart_gp + 1
!
!            ppiclf_iprop_gp(1,ppiclf_npart_gp) = ppiclf_iprop(1,ip)
!            ppiclf_iprop_gp(2,ppiclf_npart_gp) = ppiclf_iprop(2,ip)
!            ppiclf_iprop_gp(3,ppiclf_npart_gp) = nrank
!            ppiclf_iprop_gp(4,ppiclf_npart_gp) = iig
!            ppiclf_iprop_gp(5,ppiclf_npart_gp) = jjg
!            ppiclf_iprop_gp(6,ppiclf_npart_gp) = kkg
!            ppiclf_iprop_gp(7,ppiclf_npart_gp) = nrank
!
!            ! Thierry - we don't need ppiclf_comm_CheckPeriodicBC anymore for the angular periodic ghost algorithm
!            !           as this is now taken care of when anticipating where the particle might be when calling
!            !           ppiclf_comm_InvokeAngularPeriodic
!            !           we only need to assign xr and vr to ppiclf_rprop_gp
!
!            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1) ! ppiclf_y(PPICLF_JX, ip) for the periodic ghost particle
!            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2) ! JY
!            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3) ! JZ
!
!            do k=4,PPICLF_LRP_GP
!               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
!            enddo
!  333 continue
!         enddo
!
!      enddo ! ip 
!
!      return
!      end
!!----------------------------------------------------------------------
!      subroutine ppiclf_comm_CheckAngularBC(xrank, yrank, zrank)
!!
!      implicit none
!!
!      include "PPICLF"
!!
!! Local:
!!
!      integer*4 xrank, yrank, zrank
!!
!! Output:
!!
!
!      SELECT CASE (ang_case)
!        CASE(1) ! general wedge ; 0 <= angle < 90
!!          print*, "Wedge CheckAngularBC"
!          xrank  = FLOOR((xrot(1)-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
!          yrank  = FLOOR((xrot(2)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
!          zrank  = FLOOR((xrot(3)-ppiclf_binb(5))/ppiclf_bins_dx(3))
!
!        CASE(2) ! quarter cylinder ; angle = 90
!!          print*, "Quarter Cylinder CheckAngularBC"
!          xrank  = FLOOR((abs(xrot(1))-ppiclf_binb(1))
!     >                    /ppiclf_bins_dx(1)) 
!          yrank  = FLOOR((abs(xrot(2))-ppiclf_binb(3))
!     >                   /ppiclf_bins_dx(2)) 
!          zrank  = FLOOR((xrot(3)-ppiclf_binb(5))/ppiclf_bins_dx(3))
!
!        CASE(3) ! half cylinder ; angle = 180
!          print*, "Half Cylinder CheckAngularBC"
!
!        CASE DEFAULT
!            call ppiclf_exittr('Invalid Ghost Rotational Case!$',
!     >       0.0d0 ,ppiclf_nid)
!          END SELECT
!
!      return
!      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! David's old binning method left below for now
!      finished(1) = 0
!      finished(2) = 0
!      finished(3) = 0
!      total_bin = 1 
!
!      BinMinLen(1) = MAX(BinMinLen(1),BinMinLen(2),BinMinLen(3))
!      do i=1,ppiclf_ndim
!         finished(i) = 0
!         exit_1_array(i) = ppiclf_bins_set(i)
!         exit_2_array(i) = 0
!         if (ppiclf_bins_set(i) .ne. 1) ppiclf_n_bins(i) = 1
!         ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                        ppiclf_binb(2*(i-1)+1)  ) / 
!     >                       ppiclf_n_bins(i)
!         ! Make sure exit_2 is not violated by user input
!         if (ppiclf_bins_dx(i) .lt. BinMinLen(1)) then
!            do while (ppiclf_bins_dx(i) .lt. BinMinLen(1))
!               ppiclf_n_bins(i) = max(1, ppiclf_n_bins(i)-1)
!               ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                              ppiclf_binb(2*(i-1)+1)  ) / 
!     >                             ppiclf_n_bins(i)
!         WRITE(*,*) "Inf. loop in CreateBin", i, 
!     >              ppiclf_bins_dx(i), BinMinLen(1)
!         call ppiclf_exittr('Inf. loop in CreateBin$',0.0,0)
!            enddo
!         endif
!         total_bin = total_bin*ppiclf_n_bins(i)
!      enddo
!
!      ! Make sure exit_1 is not violated by user input
!      count = 0
!      do while (total_bin > ppiclf_np)
!          count = count + 1;
!          i = modulo((ppiclf_ndim-1)+count,ppiclf_ndim)+1
!          ppiclf_n_bins(i) = max(ppiclf_n_bins(i)-1,1)
!          ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                         ppiclf_binb(2*(i-1)+1)  ) / 
!     >                        ppiclf_n_bins(i)
!          total_bin = 1
!          do j=1,ppiclf_ndim
!             total_bin = total_bin*ppiclf_n_bins(j)
!          enddo
!          if (total_bin .le. ppiclf_np) exit
!       enddo
!
!       exit_1 = .false.
!       exit_2 = .false.
!
!       do while (.not. exit_1 .and. .not. exit_2)
!          do i=1,ppiclf_ndim
!             if (exit_1_array(i) .eq. 0) then
!                ppiclf_n_bins(i) = ppiclf_n_bins(i) + 1
!                ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                               ppiclf_binb(2*(i-1)+1)  ) / 
!     >                              ppiclf_n_bins(i)
!
!                ! Check conditions
!                ! exit_1
!                total_bin = 1
!                do j=1,ppiclf_ndim
!                   total_bin = total_bin*ppiclf_n_bins(j)
!                enddo
!                if (total_bin .gt. ppiclf_np) then
!                   ! two exit arrays aren't necessary for now, but
!                   ! to make sure exit_2 doesn't slip through, we
!                   ! set both for now
!                   exit_1_array(i) = 1
!                   exit_2_array(i) = 1
!                   ppiclf_n_bins(i) = ppiclf_n_bins(i) - 1
!                   ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                                  ppiclf_binb(2*(i-1)+1)  ) / 
!     >                                  ppiclf_n_bins(i)
!                   exit
!                endif
!                
!                ! exit_2
!                if (ppiclf_bins_dx(i) .lt. BinMinLen(1)) then
!                   ! two exit arrays aren't necessary for now, but
!                   ! to make sure exit_2 doesn't slip through, we
!                   ! set both for now
!                   exit_1_array(i) = 1
!                   exit_2_array(i) = 1
!                   ppiclf_n_bins(i) = ppiclf_n_bins(i) - 1
!                   ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
!     >                                  ppiclf_binb(2*(i-1)+1)  ) / 
!     >                                  ppiclf_n_bins(i)
!                   exit
!                endif
!             endif
!          enddo
!
!          ! full exit_1
!          sum_value = 0
!          do i=1,ppiclf_ndim
!             sum_value = sum_value + exit_1_array(i)
!          enddo
!          if (sum_value .eq. ppiclf_ndim) then
!             exit_1 = .true.
!          endif
!
!          ! full exit_2
!          sum_value = 0
!          do i=1,ppiclf_ndim
!             sum_value = sum_value + exit_2_array(i)
!          enddo
!          if (sum_value .eq. ppiclf_ndim) then
!             exit_2 = .true.
!          endif
!       enddo
!      ! Check for too small bins 
!      rthresh = 1E-12
!      total_bin = 1
!      do i=1,ppiclf_ndim
!         total_bin = total_bin*ppiclf_n_bins(i)
!         if (ppiclf_bins_dx(i) .lt. rthresh) ppiclf_bins_dx(i) = 1.0
!      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


