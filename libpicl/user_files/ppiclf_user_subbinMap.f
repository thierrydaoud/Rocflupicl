!-----------------------------------------------------------------------
!
! Created Oct. 18, 2024
!
! Subroutine to map both real and ghost particles to subbins
! for nearest neighbor search
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ppiclf_user_subbinMap(i_Bin, n_SBin, tot_SBin,
     >                                  SBin_counter, SBin_map)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4  SBin_map( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_nndist)/ppiclf_nndist) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_nndist)/ppiclf_nndist)
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_nndist)/ppiclf_nndist) 
     >       + 1) - 1), (ppiclf_npart+ppiclf_npart_gp))
      INTEGER*4  SBin_counter( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_nndist)/ppiclf_nndist) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_nndist)/ppiclf_nndist)
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_nndist)/ppiclf_nndist) 
     >       + 1) - 1))
      INTEGER*4  i_Bin(3), n_SBin(3), tot_SBin
     >          
!
! Internal:
!
      REAL*8    xp(3), bin_xMin(3)
      INTEGER*4 temp_SBin, i_SBin(3),ibinTemp, i, j, k, l 

!
! Code:
!
        ! Determine ppiclf bin in each dimension for this processor
        ! All real particles are in the same bin.  Look at 1st r particle
        DO l = 1,3
            i_Bin(l) = FLOOR((ppiclf_y(l,1) - ppiclf_binb(2*l-1))
     >                 /ppiclf_bins_dx(l))
            bin_xMin(l) = ppiclf_binb(2*l-1)+i_Bin(l)*ppiclf_bins_dx(l) 
        END DO 
        ! Determine the number of subbins in each dimension
        DO l = 1,3
          IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
            n_SBin(l) = FLOOR((ppiclf_bins_dx(l)+2*ppiclf_nndist)
     >                       /ppiclf_nndist)
          ELSE
            n_SBin(l) = 0
          END IF
        END DO

        ! Determine total number of subbins
        tot_SBin = (n_SBin(1)+1)*(n_SBin(2)+1)
        IF (ppiclf_ndim .EQ. 3) THEN
          tot_SBin = tot_SBin*(n_SBin(3)+1)
        END IF
        ! Assign Subbin counters to 0
        SBin_counter = 0

        ! Map each real particle to a subbin
        DO i = 1 , ppiclf_npart
           DO l = 1,3
              IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
                 xp(l) = ppiclf_y(l,i)
              ELSE
                 xp(l) = 0.0
              END IF
           END DO 
           ! Determine subbin
           DO l = 1,3
              IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
                  i_SBin(l) = FLOOR((xp(l) - (bin_xMin(l) 
     >            - ppiclf_nndist))/ppiclf_nndist) 
              ELSE
                 i_SBin(l) = 0
              END IF
           END DO
           temp_SBin = i_SBin(1) + n_SBin(1)*i_SBin(2) +
     >                n_SBin(1)*n_SBin(2)*i_SBin(3)
           SBin_counter(temp_SBin) = SBin_counter(temp_SBin) + 1
           SBin_map(temp_SBin,SBin_counter(temp_SBin)) = i
        END DO ! real particle loop


        ! Map each ghost particle to a subbin
        DO i = 1 , ppiclf_npart_gp
          DO l = 1,3
            IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
              xp(l) = ppiclf_rprop_gp(l,i)
            ELSE
              xp(l) = 0.0
            END IF
          END DO
          ! Only map ghost particles within one neighborwidth
          ! from bin edge to subbins. All others are outside
          ! of collision search distance.
          IF (xp(1) .GT. (bin_xMin(1)-ppiclf_nndist)
     >  .AND. xp(2) .GT. (bin_xMin(2)-ppiclf_nndist)
     >  .AND. xp(3) .GT. (bin_xMin(3)-ppiclf_nndist)
     >  .AND. xp(1) .LT. (bin_xMin(1)+ppiclf_bins_dx(1)+ppiclf_nndist)
     >  .AND. xp(2) .LT. (bin_xMin(2)+ppiclf_bins_dx(2)+ppiclf_nndist)
     >  .AND. xp(3) .LT. (bin_xMin(3)+ppiclf_bins_dx(3)+ppiclf_nndist)
     >        ) THEN
            ! Determine subbin
            DO l = 1,3
              IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
                i_SBin(l) = FLOOR((xp(l) - (bin_xMin(l) 
     >          - ppiclf_nndist))/ppiclf_nndist) 
              ELSE
                i_SBin(l) = 0
              END IF
            END DO
            temp_SBin = i_SBin(1) + n_SBin(1)*i_SBin(2) +
     >                 n_SBin(1)*n_SBin(2)*i_SBin(3)
            SBin_counter(temp_SBin) = SBin_counter(temp_SBin) + 1
            ! negative in subbin map means it is ghost particle
            SBin_map(temp_SBin,SBin_counter(temp_SBin)) = -i
          END IF 
        END DO ! gp loop

      RETURN
      END 
