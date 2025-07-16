      SUBROUTINE test_grid_CreateGrid(x,y,z,dx,dy,dz)
      
      IMPLICIT NONE
      ! x(1:2): overall grid 1:min,2:max length in x,y,&z
      ! dx, dy, dz: cell lengths 

      REAL :: x(2),y(2),z(2),dx,dy,dz

      PRINT*, 'Create Grid called successfully!'

      RETURN
      END SUBROUTINE
