










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
! Purpose: Compute variables other than conserved variables at vertices.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region data
!
! Output: None.
!
! Notes: 
!   1. This routine should only be used to compute variables at vertices 
!      which can be computed exclusively from variables already interpolated 
!      to vertices.
!
! ******************************************************************************
!
! $Id: RFLU_ComputeVertexVariables.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ComputeVertexVariables(pRegion)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  
  USE ModInterfaces, ONLY: MixtGasLiq_C, &
                           MixtGasLiq_P, &
                           MixtLiq_C2_Bp, &
                           MixtLiq_D_DoBpPPoBtTTo, &
                           MixtPerf_C_DGP, & 
                           MixtPerf_C_GRT, &
                           MixtPerf_C2_GRT, &
                           MixtPerf_Cv_CpR, &
                           MixtPerf_D_PRT, &
                           MixtPerf_G_CpR, & 
                           MixtPerf_P_DEoGVm2, & 
                           MixtPerf_R_M, &
                           MixtPerf_T_CvEoVm2, & 
                           MixtPerf_T_DPR
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: indCp,indMol,ivg
  REAL(RFREAL) :: Bg2,Bl2,Bp,Bt,Bv2,Cg2,Cl2,Cv2,cp,cvg,cvl,Cvm,cvv,Eo,g,gc, &
                  ir,mw,p,po,r,rEo,rg,rgas,rl,ro,ru,rv,rvap,rw,rYg,rYl,rYv, &
                  to,u,v,vfg,vfl,vfv,vm2,w
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ComputeVertexVariables.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeVertexVariables', &
                        "../../utilities/post/RFLU_ComputeVertexVariables.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing vertex variables...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  
  pCv   => pRegion%mixt%cvVert ! NOTE vertex variables
  pDv   => pRegion%mixt%dvVert
  pGv   => pRegion%mixt%gvVert

  pCvSpec => pRegion%spec%cvVert

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol
  
! ******************************************************************************
! Compute dependent variables
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel )

! ==============================================================================  
!   Incompressible fluid model
! ==============================================================================  
   
    CASE ( FLUID_MODEL_INCOMP ) 

! ==============================================================================  
!   Compressible fluid model. NOTE check state of solution vector.
! ==============================================================================  
    
    CASE ( FLUID_MODEL_COMP )   
      IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
        CALL ErrorStop(global,ERR_CV_STATE_INVALID,179)
      END IF ! pRegion%mixt%cvState    

      SELECT CASE ( pRegion%mixtInput%gasModel ) 

! ------------------------------------------------------------------------------
!       Gas models without liquid phase
! ------------------------------------------------------------------------------      

        CASE ( GAS_MODEL_TCPERF, & 
               GAS_MODEL_TPERF, & 
               GAS_MODEL_MIXT_TCPERF, & 
               GAS_MODEL_MIXT_TPERF, & 
               GAS_MODEL_MIXT_PSEUDO )
          DO ivg = 1,pGrid%nVertTot
            r   = pCv(CV_MIXT_DENS,ivg)
            ru  = pCv(CV_MIXT_XMOM,ivg)
            rv  = pCv(CV_MIXT_YMOM,ivg)
            rw  = pCv(CV_MIXT_ZMOM,ivg)
            rEo = pCv(CV_MIXT_ENER,ivg)                

            ir  = 1.0_RFREAL/r
            u   = ir*ru
            v   = ir*rv
            w   = ir*rw    
            Eo  = ir*rEo
            vm2 = u*u + v*v + w*w

            cp = pGv(GV_MIXT_CP ,indCp *ivg)
            mw = pGv(GV_MIXT_MOL,indMol*ivg)

            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            p = MixtPerf_P_DEoGVm2(r,Eo,g,vm2)

            pDv(DV_MIXT_PRES,ivg) = p
            pDv(DV_MIXT_TEMP,ivg) = MixtPerf_T_DPR(r,p,gc)
            pDv(DV_MIXT_SOUN,ivg) = MixtPerf_C_DGP(r,g,p)     
          END DO ! ivg

! ------------------------------------------------------------------------------
!       Mixture of gas, liquid, and vapor
! ------------------------------------------------------------------------------      

        CASE ( GAS_MODEL_MIXT_GASLIQ ) 
          IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,227)
          END IF ! pRegion%mixt%cvState

          ro  = global%refDensityLiq
          po  = global%refPressLiq
          to  = global%refTempLiq
          Bp  = global%refBetaPLiq
          Bt  = global%refBetaTLiq
          cvl = global%refCvLiq

          rgas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rgas)

          rvap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rvap)

          DO ivg = 1,pGrid%nVertTot
            r   = pCv(CV_MIXT_DENS,ivg)
            ru  = pCv(CV_MIXT_XMOM,ivg)
            rv  = pCv(CV_MIXT_YMOM,ivg)
            rw  = pCv(CV_MIXT_ZMOM,ivg)
            rEo = pCv(CV_MIXT_ENER,ivg)

            ir  = 1.0_RFREAL/r
            u   = ir*ru
            v   = ir*rv
            w   = ir*rw
            Eo  = ir*rEo
            vm2 = u*u + v*v + w*w

            rYg  = pCvSpec(1,ivg)
            rYv  = pCvSpec(2,ivg)
            rYl  = r - rYg - rYv

            Cvm  = (rYl*cvl + rYv*cvv + rYg*cvg)/r

            pDv(DV_MIXT_TEMP,ivg) = MixtPerf_T_CvEoVm2(Cvm,Eo,vm2)

            Cl2  = MixtLiq_C2_Bp(Bp)
            Cv2  = MixtPerf_C2_GRT(1.0_RFREAL,rvap,pDv(DV_MIXT_TEMP,ivg))
            Cg2  = MixtPerf_C2_GRT(1.0_RFREAL,rgas,pDv(DV_MIXT_TEMP,ivg))

            pDv(DV_MIXT_PRES,ivg) = MixtGasLiq_P(rYl,rYv,rYg,Cl2,Cv2,Cg2,r,ro, &
                                                 po,to,Bp,Bt, &
                                                 pDv(DV_MIXT_TEMP,ivg))
            rl = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,pDv(DV_MIXT_PRES,ivg),po, &
                                        pDv(DV_MIXT_TEMP,ivg),to)
            rv = MixtPerf_D_PRT(pDv(DV_MIXT_PRES,ivg),rvap,pDv(DV_MIXT_TEMP,ivg))
            rg = MixtPerf_D_PRT(pDv(DV_MIXT_PRES,ivg),rgas,pDv(DV_MIXT_TEMP,ivg))

            vfg = rYg/rg
            vfv = rYv/rv
            vfl = 1.0_RFREAL - vfg - vfv

            Bl2 = -Bt/Bp
            Bv2 = rv*rvap
            Bg2 = rg*rgas

            pDv(DV_MIXT_SOUN,ivg) = MixtGasLiq_C(Cvm,r,pDv(DV_MIXT_PRES,ivg), &
                                                 rl,rv,rg,vfl,vfv,vfg,Cl2,Cv2, &
                                                 Cg2,Bl2,Bv2,Bg2) 
          END DO ! ivg       
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,295)
      END SELECT ! pRegion%mixt%gasModel

! ==============================================================================  
!   Default
! ==============================================================================  

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,303)
  END SELECT ! pRegion%mixtInput%fluidModel 

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing vertex variables done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeVertexVariables

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeVertexVariables.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:57  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:12  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:58:08  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.9  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.8  2006/03/26 20:28:45  haselbac
! Added support for GL model
!
! Revision 1.7  2005/11/10 02:47:15  haselbac
! Added support for variable properties
!
! Revision 1.6  2005/10/31 21:09:39  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.5  2004/11/14 19:56:57  haselbac
! Added setting of vertex variables for various fluid models
!
! Revision 1.4  2003/12/04 03:32:59  haselbac
! Added missing description of input argument
!
! Revision 1.3  2003/11/25 21:03:54  haselbac
! Added clarifying comment under Notes
!
! Revision 1.2  2003/03/20 20:07:19  haselbac
! Modified RegFun call to avoid probs with long "../../utilities/post/RFLU_ComputeVertexVariables.F90" names
!
! Revision 1.1  2003/03/15 19:16:54  haselbac
! Initial revision
!
! ******************************************************************************

