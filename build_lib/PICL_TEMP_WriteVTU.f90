










!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
!******************************************************************************
!
! Purpose: 
!
! Description: none.
!
! Input: 
!
! Output:
!
! Notes: 
!
!******************************************************************************
!
! $Id: PICL_F90,v 1.0 2022/05/08 bdurant Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PICL_TEMP_WriteVTU( pRegion)

!  USE 

 USE ModDataTypes
  USE ModDataStruct, ONLY : t_level,t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt
USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                               RFLU_ConvertCvPrim2Cons

 USE ModInterfaces, ONLY: RFLU_DecideWrite !BRAD added for picl
 



!DEC$ NOFREEFORM

! number of timesteps kept in history kernels
! maximum number of triangular patch boundaries

! y, y1, ydot, ydotc: 12

! rprop: 48

! rprop4: PPICLF_LRP4 - Reynolds Subgrid Stress Components

! rprop5: PPICLF_LRP5 - Storing Force Models
!#define PPICLF_R_FQSX 1
!#define PPICLF_R_FQSY 2
!#define PPICLF_R_FQSZ 3
!#define PPICLF_R_FAMX 4
!#define PPICLF_R_FAMY 5
!#define PPICLF_R_FAMZ 6
!#define PPICLF_R_FAMBX 7
!#define PPICLF_R_FAMBY 8
!#define PPICLF_R_FAMBZ 9
!#define PPICLF_R_FCX 10
!#define PPICLF_R_FCY 11
!#define PPICLF_R_FCZ 12
!#define PPICLF_R_FVUX 13
!#define PPICLF_R_FVUY 14
!#define PPICLF_R_FVUZ 15
!#define PPICLF_R_QQ 16
!#define PPICLF_R_FPGX 17 
!#define PPICLF_R_FPGY 18 
!#define PPICLF_R_FPGZ 19 

! map: 22
!--- x,y,z Forces Fedback to Rocflu
!---
!--- Add comment about these terms 
!--- Reynolds Subgrid Stress Tensor
!--- Pseudo Turbulent Kinetic Energy






















!DEC$ FREEFORM


  IMPLICIT NONE


! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString


TYPE(t_global), POINTER :: global
TYPE(t_level), POINTER :: levels(:)
TYPE(t_region), POINTER :: pRegion
TYPE(t_grid), POINTER :: pGrid
!INTEGER :: errorFlag

  INTEGER :: errorFlag,icg      
  REAL(KIND=8) :: piclDtMin,piclCurrentTime, &
          temp_drudtMixt,temp_drvdtMixt,temp_drwdtMixt,energydotg


   
!******************************************************************************

  RCSIdentString = '$RCSfile: PICL_TEMP_WriteVTU.F90,v $ $Revision: 1.0 $'
 
  global => pRegion%global
  
  CALL RegisterFunction( global, 'PICL_TEMP_WriteVTU',"../rocpicl/PICL_TEMP_WriteVTU.F90" )



!PPICLF Write

   piclCurrentTime = REAL(global%currentTime,8)

   call ppiclf_solve_WriteVTU(piclCurrentTime)


! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PICL_TEMP_WriteVTU

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PICL_.F90,v $
!
!
!******************************************************************************

