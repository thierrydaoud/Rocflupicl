#include "../../ppiclF/source/PPICLF_USER.h"
#include "../../ppiclF/source/PPICLF_STD.h"
! General variables
      INTEGER*4 rootProc, nid, nproc, ierr, icomm, test

      COMMON /Gen_Int/ 
     >          rootProc, nid, nproc, ierr, icomm, test

      REAL*8    PI, randNum

      COMMON /Gen_Real/ 
     >          PI

! Grid variables
      INTEGER*4 nCells(3), proc_ncells, numCells,
     >          iCend, iCstart, cellsPerProc

      COMMON /Grid_INT/ 
     >          nCells, proc_ncells, numCells,
     >          iCend, iCstart, cellsPerProc

      REAL*8    p_grid(7,PPICLF_LEE), grid(7,PPICLF_LEE),
     >          gridDomain(2,3), gridDX(3), filter(3), 
     >          nFilterCells, dx_min(3) 

      COMMON /Grid_Real/ 
     >          p_grid, grid,
     >          gridDomain, gridDX, filter, 
     >          nFilterCells, dx_min

      LOGICAL   gridBoundsDefined
      COMMON /Grid_Log/ gridBoundsDefined


! Particle variables
      REAL*8    part_y(PPICLF_LRS,PPICLF_LPART), pdia, 
     >          p_part_y(PPICLF_LRS,PPICLF_LPART), part_dx(3),
     >          p_part_r(PPICLF_LRP,PPICLF_LPART), nndist

      COMMON /Part_Real/ 
     >          part_y, pdia, 
     >          p_part_y, part_dx,
     >          p_part_r, nndist

      INTEGER*4 particlesPerProc, npart_local,
     >          totalParticles, iPend, iPstart

      COMMON /Part_Int/ 
     >          particlesPerProc, npart_local,
     >          totalParticles, iPend, iPstart

