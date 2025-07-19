      SUBROUTINE test_CreateParticles(gDom,dxr,dp,dxp,ncell,part_y)
      
      IMPLICIT NONE
      ! Input/Output
      ! gDom(1:2,1:3) - fluid domain min/max in x,y,z
      ! dxr - particle size to grid dx size ratio
      ! dp  - particle diameter
      ! dxp - particle spacing
      ! ncell - number of fluid cells to fill with particles
      ! part_y(1:3,1:30,000) - particle centroid x,y,z coordinates 
      REAL*8    dxr, dp, dxp, gDom(2,3), part_y(3,30000)
      INTEGER*4 ncell, Pcount, PartDomain, i, j, k, ii

      PartDomain = dxr*ncell
      Pcount = 0
      DO i = 1,PartDomain
        DO j = 1,PartDomain
          DO k = 1,PartDomain
            Pcount = Pcount + 1
            part_y(1,Pcount) = gDom(1,1) + dxp*(i-1)
            part_y(2,Pcount) = gDom(1,2) + dxp*(j-1)
            part_y(3,Pcount) = gDom(1,3) + dxp*(k-1)
          END DO
        END DO
      END DO

      ! Add a single particle in all domain corners besides origin for
      ! periodicity testing
      DO i = 1,2
        DO j = 1,2
          DO k = 1,2
            IF(i .EQ. 1 .AND. j. EQ. 1 .AND. k .EQ. 1) CYCLE
            Pcount = Pcount + 1
            part_y(1,Pcount) = gDom(2,1)*(i-1) - dp/4
            part_y(2,Pcount) = gDom(2,2)*(j-1) - dp/4
            part_y(3,Pcount) = gDom(2,3)*(k-1) - dp/4
            DO ii = 1,3
              IF(part_y(ii,Pcount) .LT. 0.0D0) part_y(ii,Pcount) = 0.0D0
            END DO
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE
