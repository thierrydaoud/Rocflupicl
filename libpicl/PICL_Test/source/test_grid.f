      SUBROUTINE test_CreateGrid(gIn,NCells,ProcID,NProc,gOut)
      
      IMPLICIT NONE

      ! Input/Output
      ! gridIn(1:2,1:3) - fluid domain 1:min,2:max faces in x,y,&z
      ! Cells - per processor
      ! Num Proc - number of processors
      ! gridOut 1:3 - centroids, 4:6 - cell lengths, 7 - cell volume
      INTEGER*4 NCells, NProc, ProcID
      REAL*8    gIn(2,3),gOut(7,NCells)

      ! Local
      INTEGER*4 i, j, k, cellCount, NCellsPerDim
      REAL*8    delta(3), xMin(3), dx(3)

      NCellsPerDim = INT(NCells**(1.0/3.0))
      ! Create uniform rectangular grid
      DO i = 1,3
        !processors fluid domain length
        delta(i) = (gIn(2,1) - gIn(1,1))/NProc
        xMin(i)  = delta(i)*NProc ! This processors Min x,y,z
        dx(i)    = delta(i)/NCellsPerDim
      END DO

      cellCount = 0
      DO i = 1,NCellsPerDim 
        DO j = 1,NCellsPerDim
          DO k = 1,NCellsPerDim
            cellCount = cellCount + 1
            gOut(1,cellCount) = xMin(1) + (i+0.5-1)*dx(1) !x centroid
            gOut(2,cellCount) = xMin(2) + (j+0.5-1)*dx(2) !y centroid
            gOut(3,cellCount) = xMin(3) + (k+0.5-1)*dx(3) !z centroid
            gOut(4,cellCount) = dx(1)
            gOut(5,cellCount) = dx(2)
            gOut(6,cellCount) = dx(3)
            gOut(7,cellCount) = dx(1)*dx(2)*dx(3)
          END DO !k
        END DO !j
      END DO !i

      IF(cellCount .NE. NCells) THEN
        PRINT*,'ERROR, wrong NCells generated on proc', NProc
      END IF

      RETURN
      END SUBROUTINE
