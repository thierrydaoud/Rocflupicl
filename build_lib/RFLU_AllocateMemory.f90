










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
! Purpose: Allocate memory for Rocflu.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_AllocateMemory.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_AllocateMemory(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE RFLU_ModAllocateMemory, ONLY: RFLU_AllocateMemoryGSpeeds, &
                                    RFLU_AllocateMemorySol, & 
                                    RFLU_AllocateMemoryTStep, & 
                                    RFLU_AllocateMemoryTStep_C, & 
                                    RFLU_AllocateMemoryTStep_I 
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_AllocateMemory.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AllocateMemory',"../rocflu/RFLU_AllocateMemory.F90")

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Allocating memory for mixture...'
  END IF ! global%verbLevel

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel )

! ==============================================================================
!   Compressible fluid 
! ==============================================================================

    CASE ( FLUID_MODEL_COMP ) 
      CALL RFLU_AllocateMemorySol(pRegion)
      CALL RFLU_AllocateMemoryTStep(pRegion)
      CALL RFLU_AllocateMemoryTStep_C(pRegion)
      CALL RFLU_AllocateMemoryGSpeeds(pRegion)

! ==============================================================================
!   Incompressible fluid 
! ==============================================================================

    CASE ( FLUID_MODEL_INCOMP )
      CALL RFLU_AllocateMemorySol(pRegion)
      CALL RFLU_AllocateMemoryTStep(pRegion)      
      CALL RFLU_AllocateMemoryTStep_I(pRegion)
      CALL RFLU_AllocateMemoryGSpeeds(pRegion)

! ==============================================================================
!   Default
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,154)
  END SELECT ! pRegion%mixtInput%fluidModel

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Allocating memory for mixture done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AllocateMemory

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_AllocateMemory.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:47  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:59  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:56  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:00  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.23  2006/03/26 20:22:12  haselbac
! Removed error trap for GL model
!
! Revision 1.22  2004/12/19 15:48:55  haselbac
! Modified so can select different fluid models
!
! Revision 1.21  2004/03/19 21:21:06  haselbac
! Complete rewrite
!
! Revision 1.20  2004/03/03 23:55:40  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.19  2004/01/29 22:59:15  haselbac
! Added/deleted allocation for mfMixt/vfMixt arrays, clean-up
!
! Revision 1.18  2003/12/04 03:29:57  haselbac
! Added memory allocation for gradients, cleaned up
!
! Revision 1.17  2003/11/25 21:04:34  haselbac
! Added allocation for mass and volume fluxes on patches
!
! Revision 1.16  2003/11/03 03:50:15  haselbac
! Removed allocation of bf2bg list
!
! Revision 1.15  2003/05/13 23:48:36  haselbac
! Changed allocation of gs for GENX
!
! Revision 1.14  2003/04/18 20:00:42  haselbac
! Added explicit initialization (prevent Frost-problem)
!
! Revision 1.13  2003/03/31 16:16:26  haselbac
! Added disp array, some cosmetics
!
! Revision 1.12  2003/03/15 18:24:04  haselbac
! Some changes for parallel computations
!
! Revision 1.11  2003/01/28 14:20:33  haselbac
! Consolidation of allocation, some clean-up
!
! Revision 1.10  2002/11/02 02:03:09  wasistho
! Added TURB statistics
!
! Revision 1.9  2002/10/27 19:11:01  haselbac
! Proper allocation of grid motion variables
!
! Revision 1.8  2002/10/08 15:49:29  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.7  2002/09/09 15:26:31  haselbac
! global and mixtInput now under regions, new allocation statements
!
! Revision 1.6  2002/07/25 14:26:42  haselbac
! Added allocation for cell gradients
!
! Revision 1.5  2002/06/27 15:25:54  haselbac
! Changed allocation to *Tot for parallelization
!
! Revision 1.4  2002/06/17 13:34:12  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.3  2002/06/14 21:54:35  wasistho
! Added time avg statistics
!
! Revision 1.2  2002/06/14 20:21:19  haselbac
! Added grid speed stuff
!
! Revision 1.1  2002/05/04 17:01:59  haselbac
! Initial revision
!
!******************************************************************************

