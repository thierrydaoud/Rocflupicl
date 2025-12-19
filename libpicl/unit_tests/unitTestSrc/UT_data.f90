#include "PPICLF_STD.h"
module UT_data
! General variables
      INTEGER*4 rootProc, nid, nproc, ierr, icomm, test


      REAL*8    PI, randNum


! Grid variables
      INTEGER*4 nCells(3), proc_ncells, numCells, iCend, iCstart, cellsPerProc


      REAL*8    p_grid(7,PPICLF_LEE), grid(7,PPICLF_LEE), gridDomain(2,3), gridDX(3), filter(3), nFilterCells, dx_min(3) 


      LOGICAL   gridBoundsDefined



! Particle variables
      REAL*8    part_y(PPICLF_LRS,PPICLF_LPART), pdia, p_part_y(PPICLF_LRS,PPICLF_LPART), part_dx(3), p_part_r(PPICLF_LRP,PPICLF_LPART), nndist


      INTEGER*4 particlesPerProc, npart_local, totalParticles, iPend, iPstart



end module