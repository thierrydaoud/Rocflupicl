










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
! ******************************************************************************
!
! Purpose: Define data types related to chemical species.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ModSpecies.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

MODULE ModSpecies

  USE ModDataTypes
  USE ModMaterials, ONLY: t_material

  IMPLICIT NONE

! ******************************************************************************
! Species input
! ******************************************************************************

  TYPE t_spec_input
    LOGICAL :: sourceFlag,usedFlag
    INTEGER :: nSpecies,nSpeciesEE
    TYPE(t_spec_type), DIMENSION(:), POINTER :: specType
  END TYPE t_spec_input

! ******************************************************************************
! Species type
! ******************************************************************************

  TYPE t_spec_type

! ==============================================================================
!   Common to all species
! ==============================================================================

    LOGICAL :: frozenFlag,discreteFlag,settlingFlag
    INTEGER :: iCont,sourceType,velocityMethod
    REAL(RFREAL) :: initVal,schmidtNumber
    TYPE(t_material), POINTER :: pMaterial

! ==============================================================================
!   Specific to species representing particles
! ==============================================================================

    INTEGER :: iSpecEEv2iSpec,iSpec2iSpecEEv
    REAL(RFREAL) :: diameter,puffFactor
    REAL(RFREAL) :: effectiveDensity,effectiveVolume,materialVolume ! derived
    REAL(RFREAL) :: tauCoefficient
  END TYPE t_spec_type

! ******************************************************************************
! Species data
! ******************************************************************************

  TYPE t_spec
    INTEGER :: cvState,nSpecEqs
    INTEGER, DIMENSION(:), POINTER :: cvInfo
    REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,cvOld,cvVert
    REAL(RFREAL), DIMENSION(:,:), POINTER :: dv,tv
    REAL(RFREAL), DIMENSION(:,:), POINTER :: diss,fterm,rhs,rhsSum
    REAL(RFREAL), DIMENSION(:,:), POINTER :: lim
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: eev,eevVert, & 
                                               gradCell,gradFace
  END TYPE t_spec

END MODULE ModSpecies

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModSpecies.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:38  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:11  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:17  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.14  2006/08/19 15:38:57  mparmar
! Removed bGradFace
!
! Revision 1.13  2005/11/27 02:08:12  haselbac
! Added EEv support
!
! Revision 1.12  2005/11/14 16:57:33  haselbac
! Added iCont flag
!
! Revision 1.11  2005/11/10 02:22:09  haselbac
! Added settlingFlag
!
! Revision 1.10  2005/07/11 19:28:18  mparmar
! Added lim
!
! Revision 1.9  2005/04/20 14:41:12  haselbac
! Removed unifSpec, cosmetics
!
! Revision 1.8  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.7  2004/07/28 15:40:30  jferry
! added USED field to SPECIES input section
!
! Revision 1.6  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.5  2004/04/01 21:29:19  haselbac
! Added sourceFlag
!
! Revision 1.4  2004/01/29 22:57:29  haselbac
! Added arrays for second-order and viscous fluxes, clean-up
!
! Revision 1.3  2003/11/25 21:03:16  haselbac
! Added definition of variables
!
! Revision 1.2  2001/12/11 21:59:29  jblazek
! memory allocation added.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
! ******************************************************************************

