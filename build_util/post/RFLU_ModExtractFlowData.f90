










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
! Purpose: Collect routines to extract data from flow solution.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes:
!   1. These routines are hardcoded to extract data for particular flows on 
!      particular grids. This means that one CANNOT use these routines for any 
!      grid. 
!
! ******************************************************************************
!
! $Id: RFLU_ModExtractFlowData.F90,v 1.4 2016/05/05 22:25:13 rahul Exp $
!
! Copyright: (c) 2004-2007 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModExtractFlowData

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region


  USE ModSpecies, ONLY: t_spec_input, &
                        t_spec_type
  USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex


  USE RFLU_ModExactFlow
  USE RFLU_ModExtractFlowDataUtils

  USE ModInterfaces, ONLY: MixtPerf_G_CPR, &
                           MixtPerf_P_DEoGVm2, &
                           MixtPerf_R_M, &
                           MixtPerf_T_DPR

  !BBR
  USE RFLU_ModPlottingVars
  !BBR
 
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PARAMETER, PRIVATE :: &
    RCSIdentString = '$RCSfile: RFLU_ModExtractFlowData.F90,v $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_ExtractFlowData

! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: RFLU_ExtractFlowDataAcoustic, &
             RFLU_ExtractFlowDataBlasius, & 
             RFLU_ExtractFlowDataCyldet,&
             RFLU_ExtractFlowDataLineFarf, &
             RFLU_ExtractFlowDataNSCBC, &
             RFLU_ExtractFlowDataProudman, &
             RFLU_ExtractFlowDataSod, &
             RFLU_ExtractFlowDataSTG2D, &   
             RFLU_ExtractFlowDataSkews, &
             RFLU_ExtractFlowDataSphdet,&
             RFLU_ExtractFlowDataSurfCylds, &
             RFLU_ExtractFlowDataSurfCylWall, &
             RFLU_ExtractFlowDataSurfSphds, &
             RFLU_ExtractFlowDataTaylorVortex, &
             RFLU_ExtractFlowDataVolumeCylds, &
             RFLU_WriteMeshBump, &
             RFLU_ExtractFlowDataShktb 

! ******************************************************************************
! Subroutines and functions
! ******************************************************************************

  CONTAINS








! ******************************************************************************
!
! Purpose: Extract data from flow solution.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowData(pRegion)
                 
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
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
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowData', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extracting data from flow solution...'
                             
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME,'Case: ',TRIM(global%casename)
    END IF ! global%verbLevel                         
  END IF ! global%verbLevel
  
! ******************************************************************************
! Initialize flow field based on user input
! ******************************************************************************

  SELECT CASE ( global%casename )
  
! ==============================================================================
!   Acoustic flow case 
! ==============================================================================

    CASE ( "acoustic" )
      CALL RFLU_ExtractFlowDataAcoustic(pRegion)

! ==============================================================================
!   Bump test case 'bumpq10'
! ==============================================================================

    CASE ( "bumpq10" )
      CALL RFLU_WriteMeshBump(pRegion)
 
! ==============================================================================
!   NSCBC farfield 
! ==============================================================================

    CASE ( "farf" )
      CALL RFLU_ExtractFlowDataLineFarf(pRegion)
 
! ==============================================================================
!   Incompressible laminar flat plate
! ==============================================================================  
  
    CASE ( "lfpbli-16.48x8.8x3"   ,"lfpbli-32.96x16.16x3", & 
           "lfpbli-64.192x32.32x3","lfpbli-128.384x32.32x3"  )
      CALL RFLU_ExtractFlowDataBlasius(pRegion)
    CASE( "lfpblim-64x16x1" , "lfpblim-128x32x1", & 
          "lfpblim-256x64x1", "lfpblim-512x128x1" ) 
      CALL RFLU_ExtractFlowDataBlasius(pRegion)

! ==============================================================================
!   NSCBC test cases 'nscbcX'
! ==============================================================================

    CASE ( "nscbc1","nscbc2","nscbc3","nscbc4","nscbc5","nscbc6","nscbc7", & 
           "nscbc8" )
      CALL RFLU_ExtractFlowDataNSCBC(pRegion)

! ==============================================================================
!   ONERA C0
! ==============================================================================

    CASE ( "onera_c0_2d_100x50" )
      CALL RFLU_ExtractFlowDataProudman(pRegion)

! ==============================================================================
!   Skews diffracting shock
! ==============================================================================

    CASE ( "skews_ms2p0","skews_ms3p0","skews_ms4p0" )
      CALL RFLU_ExtractFlowDataSkews(pRegion)
 
! ==============================================================================
!   Shock tubes
! ==============================================================================

    CASE ( "st_sod1","st_sod1_mp2","st_sod2","st_sod2_mp2" )
      CALL RFLU_ExtractFlowDataSod(pRegion)
    CASE ( "stg1d" )
      CALL RFLU_ExtractFlowDataSTG1D(pRegion)      
    CASE ( "stg2d" )
      CALL RFLU_ExtractFlowDataSTG2D(pRegion)      

! ==============================================================================
!   Sommerfeld shock-particle interaction
! ==============================================================================

    CASE ( "somm_spi" )
      CALL RFLU_ExtractFlowDataSommSPI(pRegion)

! ==============================================================================
!   Taylor vortex problem 
! ==============================================================================

    CASE ( "taylorvortex" )
      CALL RFLU_ExtractFlowDataTaylorVortex(pRegion)

! ==============================================================================
!   Cyldet case
! ==============================================================================

    CASE ( "cyldet" )
      CALL RFLU_ExtractFlowDataCyldet(pRegion)

! ==============================================================================
!   Cylinder diffracting shock
! ==============================================================================

    !CASE ( "cylds" )
    CASE ( "cylpotential" )
      CALL RFLU_ExtractFlowDataSurfCylds(pRegion)
!      CALL RFLU_ExtractFlowDataVolumeCylds(pRegion)


! Saptarshi - 01/14/2015 - begin
! ==============================================================================
!   Square ShockTube problems
! ==============================================================================

    CASE ( "shktb" )  
      CALL RFLU_ExtractFlowDataShktb(pRegion)

! Saptarshi - 01/14/2015 - end

! ==============================================================================
!   Cylinder near wall 
! ==============================================================================

    CASE ( "cylwall" )
      CALL RFLU_ExtractFlowDataSurfCylWall(pRegion)

! ==============================================================================
!   Sphdet case
! ==============================================================================

    CASE ( "sphdet" )
      CALL RFLU_ExtractFlowDataSphdet(pRegion)

! ==============================================================================
!   Sphere case
! ==============================================================================

    CASE ( "sphds" )
      CALL RFLU_ExtractFlowDataSurfSphds(pRegion)
 
! ==============================================================================
!   Default - due to input error or missing CALL in this routine
! ==============================================================================  
            
    CASE DEFAULT 
      global%warnCounter = global%warnCounter + 1
    
      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A,2(1X,A))') SOLVER_NAME,'*** WARNING ***', & 
              'Extraction of data not available.', & 
              'Returning to calling procedure.'                                   
      END IF ! global%verbLevel
  END SELECT ! global%casename

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extracting data from flow solution done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowData








! ******************************************************************************
!
! Purpose: Extract data for acoustic case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataAcoustic(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,indCp,indMol
  REAL(RFREAL) :: A1,A2,cp,d,dc,de,idc,g,gc,Mo,mw,p,pc,pe,po,ro,t,u,uc,ue,v, &
                  vc,ve,w,wc,we,x
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataAcoustic',"../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract acoustic flow data:'
  END IF ! global%verbLevel

! ******************************************************************************
! Open file for data
! ******************************************************************************

  WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.xy_',global%currentTime,'.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,447,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Extract 1D data
! ******************************************************************************

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  A1 = pRegion%mixtInput%prepRealVal1
  A2 = pRegion%mixtInput%prepRealVal2
  Mo = pRegion%mixtInput%prepRealVal3
  ro = pRegion%mixtInput%prepRealVal4
  po = pRegion%mixtInput%prepRealVal5
  
  DO icg = 1,pRegion%grid%nCells
    x = pRegion%grid%cofg(XCOORD,icg)
                           
    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
      g  = global%refGamma
    ELSE
      mw = pRegion%mixt%gv(GV_MIXT_MOL,indMol*icg)
      cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *icg)
                   
      gc = MixtPerf_R_M(mw)
      g  = MixtPerf_G_CpR(cp,gc)
    END IF ! solverType

    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
      t = global%currentTime
      CALL RFLU_ComputeExactFlowAcoustic(global,x,t,ro,po,Mo,g, &
                                         A1,A2,d,u,v,w,p)

      ue = u

      t = global%currentTime + global%dtImposed/2.0_RFREAL
      CALL RFLU_ComputeExactFlowAcoustic(global,x,t,ro,po,Mo,g, &
                                         A1,A2,d,u,v,w,p)

      de = d
      pe = p

      dc = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      uc = pRegion%mixt%cv(CV_MIXT_XVEL,icg)
      pc = pRegion%mixt%cv(CV_MIXT_PRES,icg)
    ELSE
      t = global%currentTime
      CALL RFLU_ComputeExactFlowAcoustic(global,x,t,ro,po,Mo,g, &
                                         A1,A2,d,u,v,w,p)

      ue = u
      de = d
      pe = p

      dc  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      idc = 1.0_RFREAL/dc

      uc = idc*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
      pc = pRegion%mixt%dv(DV_MIXT_PRES,icg)
    END IF ! solverType

    WRITE(IF_EXTR_DATA1,'(10(1X,E23.16))') x,dc,uc,pc,de,ue,pe,dc-de,uc-ue,pc-pe
  END DO ! icg

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,524,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract data for acoustic case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataAcoustic










! ******************************************************************************
!
! Purpose: Extract data for Blasius flow
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Many assumptions are made, most of which are documented directly in the 
!      code below. Others are listed in the following: 
!   2. Plate is assumed to lie on x-z plane, with the plate leading edge at x=0,
!      and plate being defined by plane y=0.
!   3. Plate length is assumed to be unity, so that x-stations at which 
!      velocity profiles are extracted are equally spaced in [0,1].
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataBlasius(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1,iFileName2
  INTEGER, PARAMETER :: NPOINTS = 18
  INTEGER :: errorFlag,i,icg,ifl,ile,ite,ivp,j,jbl,jfs,k,nvp
  INTEGER, DIMENSION(:), ALLOCATABLE :: iloc
  REAL(RFREAL) :: cf,cfTheory,delta0,delta0Theory,delta1,delta1Theory,delta2, &
                  delta2Theory,delta3,delta3Theory,dist1,dist2,dy,eta,inter, &
                  ir,lRef,r,ReRef,Rex,ru,rv,slope,term1,term2,u,uNorm,uRef, &
                  v,vNorm,x,y
  REAL(RFREAL), DIMENSION(NPOINTS), PARAMETER :: etaBlas = & 
    (/0.0000_RFREAL,0.1000_RFREAL,0.2000_RFREAL,0.3000_RFREAL,0.4000_RFREAL, &
      0.5000_RFREAL,0.6000_RFREAL,0.8000_RFREAL,1.0000_RFREAL,1.2000_RFREAL, &
      1.4000_RFREAL,1.6000_RFREAL,1.8000_RFREAL,2.0000_RFREAL,2.5000_RFREAl, &
      3.0000_RFREAL,3.5000_RFREAL,4.0000_RFREAL/)                  
  REAL(RFREAL), DIMENSION(NPOINTS), PARAMETER :: uNormBlas = & 
    (/0.0000_RFREAL,0.0664_RFREAL,0.1328_RFREAL,0.1989_RFREAL,0.2647_RFREAL, &
      0.3298_RFREAL,0.3938_RFREAL,0.5168_RFREAL,0.6298_RFREAL,0.7290_RFREAL, &
      0.8115_RFREAL,0.8761_RFREAL,0.9233_RFREAL,0.9555_RFREAL,0.9916_RFREAl, &
      0.9990_RFREAL,0.9999_RFREAL,1.0000_RFREAL/)
  REAL(RFREAL), DIMENSION(NPOINTS), PARAMETER :: vNormBlas = & 
    (/0.0000_RFREAL,0.0033_RFREAL,0.0133_RFREAL,0.0298_RFREAL,0.0528_RFREAL, &
      0.0821_RFREAL,0.1173_RFREAL,0.2033_RFREAL,0.3048_RFREAL,0.4136_RFREAL, &
      0.5206_RFREAL,0.6172_RFREAL,0.6972_RFREAL,0.7581_RFREAL,0.8373_RFREAl, &
      0.8482_RFREAL,0.8600_RFREAL,0.8604_RFREAL/)      
  REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: yu
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataBlasius', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract data for Blasius flow...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

  nvp = 9 ! Number of equally-spaced stations at which profiles are extracted

! ******************************************************************************
! Get data about freestream conditions
! ******************************************************************************

  uRef  = global%refVelocity
  ReRef = global%refReNum
  lRef  = global%refLength 

! ******************************************************************************
! Get data about dimensions
!   ile         i-index of first cell on plate (leading edge)
!   ite         i-index of last cell on plate (trailing edge)
!   jbl         j-index of (nominally) last cell in boundary layer
!   jfs         j-index of cell abutting freestream boundary
!   k           k-index of middle layer of cells
! ******************************************************************************

  SELECT CASE ( global%casename ) 
    CASE ( "lfpbli-16.48x8.8x3" ) 
      ile =  17 
      ite =  64
      jbl =   8
      jfs =  16
      k   =   2
    CASE ( "lfpbli-32.96x16.16x3" ) 
      ile =  33 
      ite = 128
      jbl =  16
      jfs =  32
      k   =   2  
    CASE ( "lfpbli-64.192x32.32x3" ) 
      ile =  65 
      ite = 256
      jbl =  32
      jfs =  64
      k   =   2 
    CASE ( "lfpbli-128.384x64.64x3" ) 
      ile = 129 
      ite = 512
      jbl =  64
      jfs = 128
      k   =   2   
    CASE ( "lfpblim-64x16x1" ) 
      ile =  17
      ite =  64
      jbl =   8
      jfs =  16
      k   =   1 
    CASE ( "lfpblim-128x32x1" ) 
      ile =  33
      ite = 128
      jbl =  16
      jfs =  32
      k   =   1      
    CASE ( "lfpblim-256x64x1" ) 
      ile =  65
      ite = 256
      jbl =  32
      jfs =  64
      k   =   1
    CASE ( "lfpblim-512x128x1" ) 
      ile = 129
      ite = 512
      jbl =  64
      jfs = 128
      k   =   1                               
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,704)
  END SELECT ! global%casename

! ******************************************************************************
! Allocate temporary memory
! ******************************************************************************

  ALLOCATE(yu(2,0:jfs),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,714,'yu')
  END IF ! global%error

  ALLOCATE(iloc(nvp),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,720,'iloc')
  END IF ! global%error  

! ******************************************************************************
! Find i-indices of x-stations at which data about velocity profiles is to be
! extracted. NOTE assume that the cell preceding a given cell in the stream-
! wise direction has i-index smaller by one. 
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
      'Determining stations at which velocity data is extracted...'
  END IF ! global%verbLevel 
   
  i = ile 

  DO ivp = 1,nvp
    x = ivp*0.1_RFREAL
  
    emptyLoop: DO 
      icg = i + (k-1)*ite*jfs ! NOTE no j-term because always zero here 

      dist1 = pGrid%cofg(XCOORD,icg  ) - x
      dist2 = pGrid%cofg(XCOORD,icg-1) - x 
            
      IF ( (dist1 > 0.0_RFREAL) .AND. (dist2 < 0.0_RFREAL) ) THEN 
        IF ( ABS(dist1) < ABS(dist2) ) THEN 
          iloc(ivp) = i
          icg = i + (k-1)*ite*jfs
        ELSE 
          iloc(ivp) = i-1  
          icg = i - 1 + (k-1)*ite*jfs            
        END IF ! ABS

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(STDOUT,'(A,5X,A,1X,I2,1X,A,1X,I3,A,1X,E13.6)') SOLVER_NAME, &
            'Station',ivp,'located at i=',iloc(ivp),', x=',pGrid%cofg(XCOORD,icg)
        END IF ! global%verbLevel
        
        EXIT emptyLoop  
      ELSE 
        i = i + 1
      END IF ! dist1      
    END DO emptyLoop
  END DO ! ivp

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
      'Determining stations at which velocity data is extracted done.'
  END IF ! global%verbLevel

! ******************************************************************************
! Extract velocity profiles
! ******************************************************************************

  DO ivp = 1,nvp
    i = iloc(ivp)

    yu(1,0) = 0.0_RFREAL
    yu(2,0) = 0.0_RFREAL

! ==============================================================================
!   Open file
! ==============================================================================
      
    WRITE(iFileName2,'(A,I2.2,A)') 'blasius-vel',ivp,'.dat'

    OPEN(IF_EXTR_DATA2,FILE=iFileName2,FORM='FORMATTED',STATUS='UNKNOWN', &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,791,'File: '//TRIM(iFileName2))
    END IF ! global%error 

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Writing velocity-profile data to file: '// & 
                               TRIM(iFileName2)
    END IF ! global%verbLevel              

! ==============================================================================
!   Normalize velocity and write to file at selected stations
! ==============================================================================
            
    DO j = 1,jfs
      icg = i + (j-1)*ite + (k-1)*ite*jfs
      
      x = pGrid%cofg(XCOORD,icg)
      y = pGrid%cofg(YCOORD,icg)
      
      Rex = ReRef*x/lRef
      eta = y/(2.0_RFREAL*x)*SQRT(Rex)
      
      r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      ru = pRegion%mixt%cv(CV_MIXT_XMOM,icg)
      rv = pRegion%mixt%cv(CV_MIXT_YMOM,icg)      
      
      ir = 1.0_RFREAL/r
      u  = ir*ru
      v  = ir*rv
      
      uNorm = u/uRef
      vNorm = v/uRef
      
      yu(1,j) = y
      yu(2,j) = uNorm
      
      WRITE(IF_EXTR_DATA2,'(3(1X,E13.6))') eta,uNorm,vNorm*SQRT(Rex)
    END DO ! j

! ==============================================================================
!   Close file
! ==============================================================================
      
    CLOSE(IF_EXTR_DATA2,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,837,'File: '// & 
                     TRIM(iFileName2))
    END IF ! global%error       
  END DO ! ivp

! ******************************************************************************
! Extract thicknesses
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  iFileName1 = 'blasius-thick.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,856,'File: '//TRIM(iFileName1))
  END IF ! global%error  

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                             'Writing thickness data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  

! ==============================================================================
! Extract data
! ==============================================================================

  DO i = ile,ite
    icg = i + (k-1)*ite*jfs ! NOTE no j-term because always zero here 
      
    x = pGrid%cofg(XCOORD,icg)
      
    Rex = ReRef*x/lRef
 
! ------------------------------------------------------------------------------
!   Compute normalized velocity profile
! ------------------------------------------------------------------------------  
  
    yu(1,0) = 0.0_RFREAL
    yu(2,0) = 0.0_RFREAL  
  
    DO j = 1,jfs
      icg = i + (j-1)*ite + (k-1)*ite*jfs
      
      y = pGrid%cofg(YCOORD,icg)
      
      r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      ru = pRegion%mixt%cv(CV_MIXT_XMOM,icg)
      
      ir = 1.0_RFREAL/r
      u  = ir*ru
      
      uNorm = u/uRef
      
      yu(1,j) = y
      yu(2,j) = uNorm
    END DO ! j  
  
! ------------------------------------------------------------------------------  
!   Determine boundary layer thickness (indicated by normalized streamwise 
!   velocity of 0.99).
! ------------------------------------------------------------------------------  

    delta0       = CRAZY_VALUE_INT
    delta0Theory = 5.0_RFREAL*x/SQRT(Rex)
        
    DO j = 1,jfs
      IF ( yu(2,j-1) < 0.99_RFREAL .AND. yu(2,j) >= 0.99_RFREAL ) THEN 
        slope = (yu(2,j  )         - yu(2,j-1)          )/(yu(1,j) - yu(1,j-1)) 
        inter = (yu(2,j-1)*yu(1,j) - yu(2,j  )*yu(1,j-1))/(yu(1,j) - yu(1,j-1)) 
        delta0 = (0.99_RFREAL - inter)/slope
      END IF ! yu            
    END DO ! j    
        
! ------------------------------------------------------------------------------  
!   Integrate normalized velocity to get displacement, momentum, and energy 
!   thicknesses and write thicknesses to file. For the moment, use simple 
!   trapezoidal rule. 
! ------------------------------------------------------------------------------  

    delta1 = 0.0_RFREAL
    delta2 = 0.0_RFREAL
    delta3 = 0.0_RFREAL
    
    delta1Theory = 1.720_RFREAL*x/SQRT(Rex)
    delta2Theory = 0.664_RFREAL*x/SQRT(Rex)
    delta3Theory = 1.044_RFREAL*x/SQRT(Rex)            
    
    DO j = 1,jbl
      dy = yu(1,j) - yu(1,j-1)
                       
      term1 = yu(2,j-1)                 
      term2 = yu(2,j  )                       
                                    
      delta1 = delta1 & 
             + 0.5_RFREAL*dy*(        (1.0_RFREAL-      term1) & 
                              +       (1.0_RFREAL-      term2))
      delta2 = delta2 & 
             + 0.5_RFREAL*dy*(  term1*(1.0_RFREAL-      term1) & 
                              + term2*(1.0_RFREAL-      term2))
      delta3 = delta3 & 
             + 0.5_RFREAL*dy*(  term1*(1.0_RFREAL-term1*term1) & 
                              + term2*(1.0_RFREAL-term2*term2))                                             
    END DO ! j

    WRITE(IF_EXTR_DATA1,'(9(1X,E13.6))') Rex,delta0/x, &
                                         delta1/delta0, &
                                         delta2/delta0, &
                                         delta3/delta0, & 
                                         delta0Theory/x, &
                                         delta1Theory/delta0Theory, & 
                                         delta2Theory/delta0Theory, &
                                         delta3Theory/delta0Theory                                         
  END DO ! i

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,964,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! Extract density profiles
! ******************************************************************************

  DO ivp = 1,nvp
    i = iloc(ivp)

    yu(1,0) = 0.0_RFREAL
    yu(2,0) = 0.0_RFREAL

! ==============================================================================
!   Open file
! ==============================================================================
      
    WRITE(iFileName2,'(A,I2.2,A)') 'blasius-rho',ivp,'.dat'

    OPEN(IF_EXTR_DATA2,FILE=iFileName2,FORM='FORMATTED',STATUS='UNKNOWN', &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,987,'File: '//TRIM(iFileName2))
    END IF ! global%error 

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Writing density-profile data to file: '// & 
                               TRIM(iFileName2)
    END IF ! global%verbLevel              

! ==============================================================================
!   write density to file at selected stations
! ==============================================================================
            
    DO j = 1,jfs
      icg = i + (j-1)*ite + (k-1)*ite*jfs
      
      x = pGrid%cofg(XCOORD,icg)
      y = pGrid%cofg(YCOORD,icg)
      
      Rex = ReRef*x/lRef
      eta = y/(2.0_RFREAL*x)*SQRT(Rex)
      
      r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      ru = pRegion%mixt%cv(CV_MIXT_XMOM,icg)
      rv = pRegion%mixt%cv(CV_MIXT_YMOM,icg)      
      
      ir = 1.0_RFREAL/r
      u  = ir*ru
      v  = ir*rv
      
      uNorm = u/uRef
      vNorm = v/uRef
      
      yu(1,j) = y
      yu(2,j) = uNorm
      
      WRITE(IF_EXTR_DATA2,'(4(1X,E13.6))') eta,r,uNorm,vNorm*SQRT(Rex)
    END DO ! j

! ==============================================================================
!   Close file
! ==============================================================================
      
    CLOSE(IF_EXTR_DATA2,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,1033,'File: '// & 
                     TRIM(iFileName2))
    END IF ! global%error       
  END DO ! ivp

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(yu,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,1045,'yu')
  END IF ! global%error

  DEALLOCATE(iloc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,1051,'iloc')
  END IF ! global%error 

! ******************************************************************************
! Write file with skin-friction data. NOTE assume that flat plate is always on
! patch 2. 
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  iFileName1 = 'blasius-cf.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,1069,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                             'Writing skin-friction data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  

! ==============================================================================
! Write data
! ==============================================================================

  pPatch => pRegion%patches(2) ! NOTE assumption

  DO ifl = 1+(k-1)*(ite-ile+1),k*(ite-ile+1)
    x = pPatch%fc(XCOORD,ifl)
        
    Rex = ReRef*x/lRef
  
    cfTheory = 0.664_RFREAL/SQRT(Rex)  
    cf       = pPatch%cf(XCOORD,ifl)
  
    WRITE(IF_EXTR_DATA1,'(3(1X,E13.6))') Rex,cfTheory,cf
  END DO ! ifl

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,1102,'File: '//TRIM(iFileName1))
  END IF ! global%error    

! ******************************************************************************
! Write file with exact velocity profile data for comparison purposes
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  iFileName1 = 'blasius-vel-exact.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,1119,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ==============================================================================
! Write data
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                             'Writing exact velocity-profile data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  

  DO j = 1,NPOINTS
    WRITE(IF_EXTR_DATA1,'(3(1X,E11.4))') etaBlas(j),uNormBlas(j),vNormBlas(j)
  END DO ! j

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,1143,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract data for Blasius flow done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataBlasius









! ******************************************************************************
!
! Purpose: Extract Line data for NSCBC case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Assume domain to lie in z=constant plane.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataNSCBC(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: nExtract
  INTEGER :: errorFlag,icg,ix
  REAL(RFREAL) :: a,p,r,u,v,w,M,xx,yy,angle
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataNSCBC', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract Line data:'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  pPatch => pRegion%patches(1)

! ******************************************************************************
! Open file for data
! ******************************************************************************

  WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.line_',global%currentTime,'.plt'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,1240,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to  '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Compute number of cells
! ******************************************************************************

  nExtract = pPatch%nBFaces

! ******************************************************************************
! Extract along cylinder surface
! ******************************************************************************
  DO ix = 1,nExtract
    icg  = pPatch%bf2c(ix)

    xx = pGrid%cofg(XCOORD,icg)
    yy = pGrid%cofg(YCOORD,icg)

    r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg)/r
    v = pRegion%mixt%cv(CV_MIXT_YMOM,icg)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg)

    WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') xx,yy,r,u,v,p
  END DO ! ix

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,1278,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract Line data for NSCBC case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataNSCBC








! ******************************************************************************
!
! Purpose: Extract Surf data for Farfield 
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Assume domain to lie in z=constant plane.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataLineFarf(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1,iFileName2,iFileName3,iFileName4
  INTEGER :: nExtract
  INTEGER :: errorFlag,icg,ix,icg1,icg2,icg3,icg4
  REAL(RFREAL) :: a,p,r,u,v,w,M,xx,yy,angle,radius,distance
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataLineFarf', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract Surf data:'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  pPatch => pRegion%patches(1)

! ******************************************************************************
! Open file for data
! ******************************************************************************

  WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.line1_',global%currentTime,'.plt'

  WRITE(iFileName2,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.line2_',global%currentTime,'.plt'

  WRITE(iFileName3,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.line3_',global%currentTime,'.plt'

  WRITE(iFileName4,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.line4_',global%currentTime,'.plt'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,1386,'File: '//TRIM(iFileName1))
  END IF ! global%error

  OPEN(IF_EXTR_DATA2,FILE=iFileName2,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,1393,'File: '//TRIM(iFileName2))
  END IF ! global%error

  OPEN(50,FILE=iFileName3,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,1400,'File: '//TRIM(iFileName3))
  END IF ! global%error

  OPEN(51,FILE=iFileName4,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,1407,'File: '//TRIM(iFileName4))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Compute number of cells
! ******************************************************************************

  nExtract = 95 
  radius   = 0.006125_RFREAL

! ******************************************************************************
! Extract along cylinder surface
! ******************************************************************************

  DO ix = 1,nExtract
    icg   = pPatch%bf2c(ix)
    icg1  =   1 + (ix-1)*258 
    icg2  =  65 + (ix-1)*258 
    icg3  = 130 + (ix-1)*258 
    icg4  = 195 + (ix-1)*258 

    xx = pGrid%cofg(XCOORD,icg1)
    yy = pGrid%cofg(YCOORD,icg1)
    distance = SQRT(xx*xx + yy*yy)

    r = pRegion%mixt%cv(CV_MIXT_DENS,icg1)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg1)/r
    v = pRegion%mixt%cv(CV_MIXT_YMOM,icg1)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg1)
    WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') distance-radius,r,u,v,p 

    xx = pGrid%cofg(XCOORD,icg2)
    yy = pGrid%cofg(YCOORD,icg2)
    distance = SQRT(xx*xx + yy*yy)

    r = pRegion%mixt%cv(CV_MIXT_DENS,icg2)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg2)/r
    v = pRegion%mixt%cv(CV_MIXT_YMOM,icg2)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg2)
    WRITE(IF_EXTR_DATA2,'(6(1X,E23.16))') distance-radius,r,u,v,p 

    xx = pGrid%cofg(XCOORD,icg3)
    yy = pGrid%cofg(YCOORD,icg3)
    distance = SQRT(xx*xx + yy*yy)

    r = pRegion%mixt%cv(CV_MIXT_DENS,icg3)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg3)/r
    v = pRegion%mixt%cv(CV_MIXT_YMOM,icg3)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg3)
    WRITE(50,'(6(1X,E23.16))') distance-radius,r,u,v,p 

    xx = pGrid%cofg(XCOORD,icg4)
    yy = pGrid%cofg(YCOORD,icg4)
    distance = SQRT(xx*xx + yy*yy)

    r = pRegion%mixt%cv(CV_MIXT_DENS,icg4)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg4)/r
    v = pRegion%mixt%cv(CV_MIXT_YMOM,icg4)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg4)
    WRITE(51,'(6(1X,E23.16))') distance-radius,r,u,v,p 
  END DO ! ix

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,1481,'File: '//TRIM(iFileName1))
  END IF ! global%error

  CLOSE(IF_EXTR_DATA2,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,1487,'File: '//TRIM(iFileName2))
  END IF ! global%error

  CLOSE(50,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,1493,'File: '//TRIM(iFileName3))
  END IF ! global%error

  CLOSE(51,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,1499,'File: '//TRIM(iFileName4))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract Surface data for Farf case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataLineFarf










! ******************************************************************************
!
! Purpose: Extract data for Proudman flow
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Many assumptions are made, most of which are documented directly in the 
!      code below, otherwise also assume that case satisfies same restrictions
!      as those for defining exact solution.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataProudman(pRegion)

  USE RFLU_ModExactFlow, ONLY: RFLU_ComputeExactFlowProudman
  USE RFLU_ModFlowHardCode, ONLY: RFLU_GetParamsHardCodeProudman

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,i,icg,ivp,j,nvp,nx,ny
  INTEGER, DIMENSION(:), ALLOCATABLE :: iloc
  REAL(RFREAL) :: dInc,height,ir,mInj,vInj,p,pTot,r,ru,rv,u,uNorm,v,vNorm,w,x,y
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataProudman', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract data for Blasius flow...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

  nvp = 9

! ******************************************************************************
! Get data about flow and geometry
! ******************************************************************************

  CALL RFLU_GetParamsHardCodeProudman(dInc,mInj,vInj,pTot)

  height = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))

! ******************************************************************************
! Get data about dimensions
! ******************************************************************************

  SELECT CASE ( global%casename ) 
    CASE ( "onera_c0_2d_100x50" ) 
      nx = 100 
      ny =  50            
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,1609)
  END SELECT ! global%casename

! ******************************************************************************
! Extract data at selected stations
! ******************************************************************************

  DO ivp = 1,nvp
    i = nx/(nvp+1)*ivp

! ==============================================================================
!   Open file
! ==============================================================================
      
    WRITE(iFileName1,'(A,I2.2,A)') 'onera_c0-vel',ivp,'.dat'

    OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,1629,'File: '//TRIM(iFileName1))
    END IF ! global%error 

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Writing velocity-profile data to file: '// & 
                               TRIM(iFileName1)
    END IF ! global%verbLevel              

! ==============================================================================
!   Normalize velocity and write to file at selected stations. NOTE assume grid 
!   is uniform in x-direction so that stations are equally spaced.
! ==============================================================================
            
    DO j = 1,ny
      icg = i + (j-1)*nx ! NOTE assumption about cell numbering

      x = pGrid%cofg(XCOORD,icg)
      y = pGrid%cofg(YCOORD,icg)

      IF ( global%verbLevel > VERBOSE_NONE .AND. j == 1 ) THEN
        WRITE(STDOUT,'(A,5X,A,1X,I2,1X,A,1X,E13.6)') SOLVER_NAME,'Station', &
                                                     ivp,'located at x=',x
      END IF ! global%verbLevel      
      
      r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      ru = pRegion%mixt%cv(CV_MIXT_XMOM,icg)
      rv = pRegion%mixt%cv(CV_MIXT_YMOM,icg)      
      
      ir = 1.0_RFREAL/r
      u  = ir*ru
      v  = ir*rv
      
      uNorm = -u/(0.5_RFREAL*global%pi*x/height*vInj)
      vNorm =  v/vInj
      
      WRITE(IF_EXTR_DATA1,'(3(1X,E13.6))') (height-y)/height,uNorm,vNorm
    END DO ! j

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,1675,'File: '// & 
                     TRIM(iFileName1))
    END IF ! global%error
  END DO ! ivp       

! ******************************************************************************
! Write file with exact velocity profile data for comparison purposes
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  iFileName1 = 'onera_c0-vel-exact.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,1694,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ==============================================================================
! Write data
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                             'Writing exact velocity-profile data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  

  DO j = 1,ny
    icg = 1 + (j-1)*nx ! NOTE set i to 1 as only need y-coordinate variation

    x = pGrid%cofg(XCOORD,icg)
    y = pGrid%cofg(YCOORD,icg)

    CALL RFLU_ComputeExactFlowProudman(global,x,y,height,dInc,vInj,pTot, &
                                       r,u,v,w,p)

    uNorm = -u/(0.5_RFREAL*global%pi*x/height*vInj)
    vNorm =  v/vInj

    WRITE(IF_EXTR_DATA1,'(3(1X,E13.6))') (height-y)/height,uNorm,vNorm
  END DO ! j

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,1729,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract data for Proudman flow done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataProudman






! ******************************************************************************
!
! Purpose: Extract data for Sod shocktube case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Assume domain to lie along x-axis.
!   2. Assume cross-section to have 3x3 cells, so can find number of cells
!      along x-axis by dividing total number of cells by 9.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSod(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,icgBeg,icgEnd,nCellsX
  REAL(RFREAL) :: a,cp,mw,p,r,T,u,v
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSod', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract data for Sod shocktube case...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Open file for data 
! ******************************************************************************

  iFileName1 = 'sodst.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,1821,'File: '//TRIM(iFileName1))
  END IF ! global%error  

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  

! ******************************************************************************
! Find cell indices for extraction
! ******************************************************************************

  nCellsX = pGrid%nCellsTot/9 ! NOTE integer division

  icgBeg = 4*nCellsX + 1
  icgEnd = 5*nCellsX

! ******************************************************************************
! Extract along line of cells 
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%gasModel )
    CASE ( GAS_MODEL_TCPERF ) 
      DO icg = icgBeg,icgEnd
        r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        u = pRegion%mixt%cv(CV_MIXT_XMOM,icg)/r
        p = pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T = pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a = pRegion%mixt%dv(DV_MIXT_SOUN,icg)

        WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                              r,u,p,T,a
      END DO ! icg    
    CASE ( GAS_MODEL_MIXT_TCPERF,GAS_MODEL_MIXT_PSEUDO ) 
      DO icg = icgBeg,icgEnd
        r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        u  = pRegion%mixt%cv(CV_MIXT_XMOM,icg)/r
        p  = pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T  = pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a  = pRegion%mixt%dv(DV_MIXT_SOUN,icg)
        mw = pRegion%mixt%gv(GV_MIXT_MOL ,icg)
        cp = pRegion%mixt%gv(GV_MIXT_CP  ,icg)

        WRITE(IF_EXTR_DATA1,'(8(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                              r,u,p,T,a,mw,cp
      END DO ! icg        
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,1868)
  END SELECT ! pRegion%mixtInput%gasModel

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,1878,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extracting data for shocktube done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSod





! ******************************************************************************
!
! Purpose: Extract data for Sommerfeld shock-particle-interaction case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Assume domain to lie along x-axis.
!   2. Assume cross-section to have 1 cell layer in z.
!   3. Assume patches 1 and 2 to be for x=xlow and x=xhigh boundaries.
!   4. Assumptions 2 and 3 mean that number of y-layers can be determined by 
!      number of faces on patch 1 or 2.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSommSPI(pRegion)

  USE ModTools, ONLY: FloatEqual

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,icgBeg,icgEnd,nCellLayersY,nCellsX
  REAL(RFREAL) :: a,cp,idx,ir,lx,mw,p,r,T,u,v,xMax,xMin,xs
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pPv
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSommSPI', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extracting data for Sommerfeld SPI case...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  pCv   => pRegion%mixt%cv
  pDv   => pRegion%mixt%dv

  IF ( ASSOCIATED(pRegion%plot%pv) .EQV. .TRUE. ) THEN 
    pPv => pRegion%plot%pv
  END IF ! ASSOCIATED

! ******************************************************************************
! Find number of cell layers in y, cell indices, and grid spacing for extraction
! ******************************************************************************

  nCellLayersY = pRegion%patches(1)%nBFaces
  nCellsX      = pGrid%nCellsTot/nCellLayersY ! NOTE integer division

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A,1X,I5)') SOLVER_NAME,'Number of cell layers in y:', &
                                   nCellLayersY
    WRITE(STDOUT,'(A,3X,A,1X,I5)') SOLVER_NAME,'Number of cells in x:      ', &
                                   nCellsX
  END IF ! global%verbLevel
   
  icgBeg = (nCellLayersY-1)*nCellsX + 1
  icgEnd = nCellLayersY*nCellsX

  xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nCellsTot))
  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nCellsTot))
  idx  = REAL(nCellsX,KIND=RFREAL)/(xMax-xMin)


! ******************************************************************************
! Write usual data
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  WRITE(iFileName1,'(A,1PE11.5)') 'somm_spi.dat_',global%currentTime  

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,2028,'File: '//TRIM(iFileName1))
  END IF ! global%error  
      
! ==============================================================================
! Extract usual data along line of cells 
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing solution data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel

  SELECT CASE ( pRegion%mixtInput%gasModel )
    CASE ( GAS_MODEL_TCPERF )             
      DO icg = icgBeg,icgEnd
        r = pCv(CV_MIXT_DENS,icg)
        ir = 1.0_RFREAL/r

        u = ir*pCv(CV_MIXT_XMOM,icg)
        p =    pDv(DV_MIXT_PRES,icg)
        T =    pDv(DV_MIXT_TEMP,icg)
        a =    pDv(DV_MIXT_SOUN,icg)
        WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                              r,u,p,T,a                                                                                                    
      END DO ! icg    
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,2068)
  END SELECT ! pRegion%mixtInput%gasModel
  
! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,2078,'File: '//TRIM(iFileName1))
  END IF ! global%error 
    
! ******************************************************************************
! Find shock position and write to file
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  WRITE(iFileName1,'(A,1PE11.5)') 'somm_spi_xs.dat_',global%currentTime

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,2095,'File: '//TRIM(iFileName1))
  END IF ! global%error   

! ==============================================================================
! Extract data and write to file
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing shock position to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel
  
  CALL RFLU_ExtractDiscontLocation1D(pRegion,icgBeg,icgEnd,nCellsX,pCv, & 
                                     CV_MIXT_DENS,xs)

  WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') global%currentTime,xs

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,2119,'File: '//TRIM(iFileName1))
  END IF ! global%error                                                                         


! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
          'Extracting data for Sommerfeld SPI case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSommSPI






! ******************************************************************************
!
! Purpose: Extract data for generic 1d shocktube case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Assume domain to lie along x-axis.
!   2. Assume have 2 species if running with mixture gas models.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSTG1D(pRegion)

  USE ModSpecies, ONLY: t_spec_type

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,icgShock
  REAL(RFREAL) :: a,cp,ir,mw,p,r,T,u,v,xp,xs
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
  TYPE(t_grid), POINTER :: pGrid
  INTEGER :: iSpec,iSpecEEv
  REAL(RFREAL) :: Y1,Y2  
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pEEv   
  TYPE(t_spec_type), POINTER :: pSpecType  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSTG1D', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract data for generic 1d shocktube case...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  pCv   => pRegion%mixt%cv

! ******************************************************************************
! Write usual data
! ******************************************************************************

! ==============================================================================
! Open file 
! ==============================================================================

  WRITE(iFileName1,'(A,1PE11.5)') 'stg1d.dat_',global%currentTime

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,2352,'File: '//TRIM(iFileName1))
  END IF ! global%error  

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  
  
! ==============================================================================
! Extract usual data along line of cells 
! ==============================================================================

  SELECT CASE ( pRegion%mixtInput%gasModel )
    CASE ( GAS_MODEL_TCPERF ) 
      DO icg = 1,pGrid%nCellsTot
        r = pCv(CV_MIXT_DENS,icg)
        ir = 1.0_RFREAL/r

        u = ir*pCv(CV_MIXT_XMOM,icg)
        p =    pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T =    pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a =    pRegion%mixt%dv(DV_MIXT_SOUN,icg)

        WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                              r,u,p,T,a
      END DO ! icg    
    CASE ( GAS_MODEL_MIXT_TCPERF,GAS_MODEL_MIXT_PSEUDO ) 
      DO icg = 1,pGrid%nCellsTot
        r  = pCv(CV_MIXT_DENS,icg)
        ir = 1.0_RFREAL/r
 
        u  = ir*pCv(CV_MIXT_XMOM,icg)
        p  =    pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T  =    pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a  =    pRegion%mixt%dv(DV_MIXT_SOUN,icg)
        mw =    pRegion%mixt%gv(GV_MIXT_MOL ,icg)
        cp =    pRegion%mixt%gv(GV_MIXT_CP  ,icg)
        Y1 = ir*pRegion%spec%cv(1,icg)
        Y2 = ir*pRegion%spec%cv(2,icg)

        WRITE(IF_EXTR_DATA1,'(10(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                               r,u,p,T,a,Y1,Y2,mw,cp
      END DO ! icg        
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,2398)
  END SELECT ! pRegion%mixtInput%gasModel

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,2408,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! Find shock position and write to file
! ******************************************************************************

! ==============================================================================
! Open file
! ==============================================================================

  WRITE(iFileName1,'(A,1PE11.5)') 'stg1d_xs.dat_',global%currentTime

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,2425,'File: '//TRIM(iFileName1))
  END IF ! global%error   

! ==============================================================================
! Extract data and write to file
! ==============================================================================

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing shock position to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel
  
  CALL RFLU_ExtractDiscontLocation1D(pRegion,1,pGrid%nCells,pGrid%nCells,pCv, &
                                     CV_MIXT_DENS,xs)

  WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') global%currentTime,xs 

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,2449,'File: '//TRIM(iFileName1))
  END IF ! global%error         


! ******************************************************************************
! Write EE data
! ******************************************************************************

  IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 

! ==============================================================================
!   Loop over species
! ==============================================================================

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)
     
! ------------------------------------------------------------------------------
!     Write EE data if species evolved with EE approach
! ------------------------------------------------------------------------------

      IF ( pSpecType%velocityMethod == SPEC_METHV_EQEUL ) THEN
        iSpecEEv = pSpecType%iSpec2iSpecEEv
  
        pEEv => pRegion%spec%eev

! ----- Open file --------------------------------------------------------------

        WRITE(iFileName1,'(A,I2.2,A)') 'stg1d',iSpecEEv,'.dat'

        OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
             IOSTAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_OPEN,2529,'File: '// &
                         TRIM(iFileName1))
        END IF ! global%error  

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I2,1X,A)') SOLVER_NAME, & 
                'Writing eev data for species',iSpec,'to file: '// &
                TRIM(iFileName1)
        END IF ! global%verbLevel          
        
! ----- Write data -------------------------------------------------------------        
        
        DO icg = 1,pGrid%nCellsTot
          WRITE(IF_EXTR_DATA1,'(5(1X,E23.16))') & 
                pGrid%cofg(XCOORD,icg), & 
!                pEEv(EEV_SPEC_XVEL,iSpecEEv,icg), & 
!                pEEv(EEV_SPEC_YVEL,iSpecEEv,icg), &
!                pEEv(EEV_SPEC_ZVEL,iSpecEEv,icg), &  
!                pEEv(EEV_SPEC_TEMP,iSpecEEv,icg)
                pEEv(1,iSpecEEv,icg), & 
                pEEv(2,iSpecEEv,icg), &
                pEEv(3,iSpecEEv,icg), &  
                pEEv(4,iSpecEEv,icg)                
        END DO ! icg 

! ----- Close file -------------------------------------------------------------
        
        CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_CLOSE,2559,'File: '// & 
                         TRIM(iFileName1))
        END IF ! global%error                    
      END IF ! pSpecType%velocityMethod
    END DO ! iSpecEE
  END IF ! pRegion%specInput%nSpeciesEE

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
          'Extracting data for generic 1d shocktube case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSTG1D




! ******************************************************************************
!
! Purpose: Extract data for generic 2d shocktube case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Assume domain to lie along x-axis.
!   2. Assume cross-section to have 3 cells, so can find number of cells
!      along x-axis by dividing total number of cells by 3.
!   3. Assume have 2 species if running with mixture gas models.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSTG2D(pRegion)

  USE ModSpecies, ONLY: t_spec_type

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,icgBeg,icgEnd,nCellsX
  REAL(RFREAL) :: a,cp,ir,mw,p,r,T,u,v
  TYPE(t_grid), POINTER :: pGrid
  INTEGER :: iSpec,iSpecEEv
  REAL(RFREAL) :: Y1,Y2  
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pEEv   
  TYPE(t_spec_type), POINTER :: pSpecType  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSTG2D', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract data for generic 2d shocktube case...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Find cell indices for extraction
! ******************************************************************************

  nCellsX = pGrid%nCellsTot/3 ! NOTE integer division

  icgBeg =   nCellsX + 1
  icgEnd = 2*nCellsX

! ******************************************************************************
! Write usual data
! ******************************************************************************

! ==============================================================================
! Open file 
! ==============================================================================

  iFileName1 = 'stg2d.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,2679,'File: '//TRIM(iFileName1))
  END IF ! global%error  

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  
  
! ==============================================================================
! Extract usual data along line of cells 
! ==============================================================================

  SELECT CASE ( pRegion%mixtInput%gasModel )
    CASE ( GAS_MODEL_TCPERF ) 
      DO icg = icgBeg,icgEnd
        r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        ir = 1.0_RFREAL/r

        u = ir*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
        p =    pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T =    pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a =    pRegion%mixt%dv(DV_MIXT_SOUN,icg)

        WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                              r,u,p,T,a
      END DO ! icg    
    CASE ( GAS_MODEL_MIXT_TCPERF,GAS_MODEL_MIXT_PSEUDO ) 
      DO icg = icgBeg,icgEnd
        r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        ir = 1.0_RFREAL/r
 
        u  = ir*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
        p  =    pRegion%mixt%dv(DV_MIXT_PRES,icg)
        T  =    pRegion%mixt%dv(DV_MIXT_TEMP,icg)
        a  =    pRegion%mixt%dv(DV_MIXT_SOUN,icg)
        mw =    pRegion%mixt%gv(GV_MIXT_MOL ,icg)
        cp =    pRegion%mixt%gv(GV_MIXT_CP  ,icg)
        Y1 = ir*pRegion%spec%cv(1,icg)
        Y2 = ir*pRegion%spec%cv(2,icg)

        WRITE(IF_EXTR_DATA1,'(10(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                               r,u,p,T,a,Y1,Y2,mw,cp
      END DO ! icg        
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,2725)
  END SELECT ! pRegion%mixtInput%gasModel

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,2735,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! Write EE data
! ******************************************************************************

  IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 

! ==============================================================================
!   Loop over species
! ==============================================================================

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)
     
! ------------------------------------------------------------------------------
!     Write EE data if species evolved with EE approach
! ------------------------------------------------------------------------------

      IF ( pSpecType%velocityMethod == SPEC_METHV_EQEUL ) THEN
        iSpecEEv = pSpecType%iSpec2iSpecEEv
  
        pEEv => pRegion%spec%eev

! ----- Open file --------------------------------------------------------------

        WRITE(iFileName1,'(A,I2.2,A)') 'stg2d',iSpecEEv,'.dat'

        OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
             IOSTAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_OPEN,2769,'File: '// &
                         TRIM(iFileName1))
        END IF ! global%error  

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I2,1X,A)') SOLVER_NAME, & 
                'Writing eev data for species',iSpec,'to file: '// &
                TRIM(iFileName1)
        END IF ! global%verbLevel          
        
! ----- Write data -------------------------------------------------------------        
        
        DO icg = icgBeg,icgEnd
          WRITE(IF_EXTR_DATA1,'(5(1X,E23.16))') & 
                pGrid%cofg(XCOORD,icg), & 
!                pEEv(EEV_SPEC_XVEL,iSpecEEv,icg), & 
!                pEEv(EEV_SPEC_YVEL,iSpecEEv,icg), &
!                pEEv(EEV_SPEC_ZVEL,iSpecEEv,icg), &  
!                pEEv(EEV_SPEC_TEMP,iSpecEEv,icg)
                pEEv(1,iSpecEEv,icg), & 
                pEEv(2,iSpecEEv,icg), &
                pEEv(3,iSpecEEv,icg), &  
                pEEv(4,iSpecEEv,icg)                
        END DO ! icg 

! ----- Close file -------------------------------------------------------------
        
        CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_CLOSE,2799,'File: '// & 
                         TRIM(iFileName1))
        END IF ! global%error                    
      END IF ! pSpecType%velocityMethod
    END DO ! iSpecEE
  END IF ! pRegion%specInput%nSpeciesEE

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
          'Extracting data for generic 2d shocktube case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSTG2D







! ******************************************************************************
!
! Purpose: Extract data for Skews shock diffraction case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Assume domain to lie in z=constant plane.
!   2. Assume domain to be divided into two pieces, a chamber (into which the
!      shock diffracts) and a tube (from which the shock emanates).
!   3. Assume that the chamber is bounded by (xMin,xMax) x (yMin,yMax) and tube
!      by (xMin,0) x (0,yMax).
!   4. Assume that grid is uniform and that cells in tube are numbered in
!      y-direction fastest. 
!   5. All these assumptions are satisfied by grids generated by gg_skews.f90
!      grid generator...
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSkews(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,icgOffs,ix,iy,nx,ny
  REAL(RFREAL) :: a,h,p,r,T,u,v,xMax,xMin,yMax,yMin
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSkews', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extracting data for Skews case...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Open file for data 
! ******************************************************************************

  iFileName1 = 'skews.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,2904,'File: '//TRIM(iFileName1))
  END IF ! global%error  

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// & 
                             TRIM(iFileName1)
  END IF ! global%verbLevel  

! ******************************************************************************
! Find cell from which to start extracting
! ******************************************************************************

! ==============================================================================
! Compute spacing from coordinate extrema
! ==============================================================================

  xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  yMax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))

  h = SQRT((xMax*(yMax-yMin)-xMin*yMax)/DBLE(pGrid%nCells))

! ==============================================================================
! Compute number of cells in duct portion of domain 
! ==============================================================================

  nx = INT(-xMin/h+0.5_RFREAL)
  ny = INT( yMax/h+0.5_RFREAL)

! ==============================================================================
! Compute number of cells in chamber portion of domain and use as offset 
! ==============================================================================

  icgOffs = INT(xMax/h+0.5_RFREAL)*INT((yMax-yMin)/h+0.5_RFREAL)  

! ******************************************************************************
! Extract along line of cells halfway up duct portion of domain
! ******************************************************************************

  iy = nx/2 + 1

  DO ix = 1,nx
    icg = icgOffs + iy + (ix - 1)*ny

    r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg)
    T = pRegion%mixt%dv(DV_MIXT_TEMP,icg)
    a = pRegion%mixt%dv(DV_MIXT_SOUN,icg)

    WRITE(IF_EXTR_DATA1,'(7(1X,E23.16))') pGrid%cofg(XCOORD,icg), & 
                                          pGrid%cofg(YCOORD,icg), &
                                          r,u,p,T,a
  END DO ! ix

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,2967,'File: '//TRIM(iFileName1))
  END IF ! global%error  

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extracting data for Skews case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSkews








! ******************************************************************************
!
! Purpose: Extract total kinetic energy for taylor vortex problem.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataTaylorVortex(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: fileExists
  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,indCp,indMol
  REAL(RFREAL) :: cp,d,g,gc,ke,kee,L,mw,pe,pi,r,re,refD,refL,refNu,refP,refU, &
                  ru,rv,rw,t,u,ue,v,ve,vol,w,we,x,y
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start, Set pointers and variables
! ******************************************************************************

  global => pRegion%global
  pGrid  => pRegion%grid

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataTaylorVortex', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract KE data for taylor-vortex problem:'
  END IF ! global%verbLevel

! ******************************************************************************
! Extract ke from all cells 
! ******************************************************************************

  ke  = 0.0_RFREAL
  kee = 0.0_RFREAL

  refL  = global%refLength
  refNu = global%refVisc/global%refDensity
  refU  = global%refVelocity
  refD  = global%refDensity
  refP  = global%refPressure

  pi = global%pi
  d  = refD
  L  = pRegion%mixtInput%prepRealVal1

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
!    refL  = global%refLength
!    refNu = global%refVisc/global%refDensity
!    refU  = global%refVelocity
!    refD  = global%refDensity
!    refP  = global%refPressure

!    pi = global%pi
!    d  = refD
!    L  = pRegion%mixtInput%prepRealVal1

    DO icg = 1,pGrid%nCells
      vol = pRegion%grid%vol(icg)

      r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      u = pRegion%mixt%cv(CV_MIXT_XVEL,icg)
      v = pRegion%mixt%cv(CV_MIXT_YVEL,icg)
      w = pRegion%mixt%cv(CV_MIXT_ZVEL,icg)

      ke = ke + 0.5_RFREAL*(u*u + v*v + w*w)*r*vol

      x = pGrid%cofg(XCOORD,icg)
      y = pGrid%cofg(YCOORD,icg)

      g = global%refGamma

      re = refD 
      t  = global%currentTime
      CALL RFLU_ComputeExactFlowTaylorVortex(t,pi,x,y,L,refL,refNu, &
                                             refU,refD,refP,ue,ve,we,pe)

      kee = kee + 0.5_RFREAL*(ue*ue + ve*ve + we*we)*re*vol
    END DO ! icg
  ELSE
    DO icg = 1,pGrid%nCells
      vol = pRegion%grid%vol(icg)

      mw = pRegion%mixt%gv(GV_MIXT_MOL,indMol*icg)
      cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *icg)

      gc = MixtPerf_R_M(mw)
      g  = MixtPerf_G_CpR(cp,gc)

      r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      ru = pRegion%mixt%cv(CV_MIXT_XMOM,icg)
      rv = pRegion%mixt%cv(CV_MIXT_YMOM,icg)
      rw = pRegion%mixt%cv(CV_MIXT_ZMOM,icg)

      ke = ke + 0.5_RFREAL*(ru*ru + rv*rv + rw*rw)*vol/r

      x = pGrid%cofg(XCOORD,icg)
      y = pGrid%cofg(YCOORD,icg)

      re = refD 
      t  = global%currentTime
      CALL RFLU_ComputeExactFlowTaylorVortex(t,pi,x,y,L,refL,refNu, &
                                             refU,refD,refP,ue,ve,we,pe)

      kee = kee + 0.5_RFREAL*(ue*ue + ve*ve + we*we)*re*vol
    END DO ! icg
  END IF ! solverType

! ******************************************************************************
! Open file for data
! ******************************************************************************

! TEMPORARY: manoj: Writing only one file instead of many files for each time
!  WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
!                              TRIM(global%casename)// &
!                              '.ke_',global%currentTime,'.dat'
  WRITE(iFileName1,'(A)') TRIM(global%outDir)//TRIM(global%casename)// &
                              '_ke.dat'

!  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
!       IOSTAT=errorFlag)

  INQUIRE(FILE=iFileName1,EXIST=fileExists)

  IF ( fileExists .EQV. .TRUE. ) THEN
    OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='OLD', &
         POSITION='APPEND',IOSTAT=errorFlag) 
  ELSE
    OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='NEW', &
         IOSTAT=errorFlag)
  END IF ! file

  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,3152,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Write data
! ******************************************************************************

!  WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') ke,kee
  WRITE(IF_EXTR_DATA1,'(3(1X,E23.16))') global%currentTime,ke,kee
! END TEMPORARY

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,3175,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract KE data for TaylorVortex done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataTaylorVortex



! ******************************************************************************
!
! Purpose: Extract Surf data for Spherical detonation case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSphdet(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: cell,nExtract
  INTEGER :: errorFlag,icg
  REAL,PARAMETER :: PI = 3.141592653589793
  REAL(RFREAL) :: p,phi,r,rad,the,u,v,w,xx,yy,zz
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSphdet', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract Surf data:'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Open file for data
! ******************************************************************************

  WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.surf_',global%currentTime,'.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,3264,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Compute number of cells
! ******************************************************************************

  nExtract = pGrid%nCells

! ******************************************************************************
! Extract data from interior cells
! ******************************************************************************
 DO icg = 1,nExtract
    !icg  = pPatch%bf2c(ix)
    cell = pGrid%hex2CellGlob(icg)
    xx   = pGrid%cofg(XCOORD,icg)
    yy   = pGrid%cofg(YCOORD,icg)
    zz   = pGrid%cofg(ZCOORD,icg)
    rad  = DSQRT(xx**2 + yy**2 + zz**2) 
    the  = ACOS(zz/rad) * 180/PI
    phi  = ATAN(yy/xx)  * 180/PI
    !phi  = ASIN(yy/SQRT(xx**2.0_RFREAL + yy**2.0_RFREAL)) * 180/PI
    r    = pRegion%mixt%cv(CV_MIXT_DENS,icg)
    p    = PRegion%mixt%dv(DV_MIXT_PRES,icg)
    u    = pRegion%mixt%cv(CV_MIXT_XMOM,icg)/r
    v    = pRegion%mixt%cv(CV_MIXT_YMOM,icg)/r
    w    = pRegion%mixt%cv(CV_MIXT_ZMOM,icg)/r
    WRITE(IF_EXTR_DATA1,'(1(1X,I8),(8(1X,E23.16)))') cell,xx,r,u,w,v,p
  END DO ! ix

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,3306,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract Surface data for Sphdet case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSphdet



! ******************************************************************************
!
! Purpose: Extract Surf data for Cylindrical detonation case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataCyldet(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion
! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) ::iFileName1,iFileName2,iFileName3, &
                      iFileName4 !Fred - Fractal Calcs
  INTEGER :: errorFlag,icg,icgBeg,icgEnd,Imax,Jmax,Kmax,nCellsSecIPerLayer, &
             rLayer,rLayers0,tLayer,tLayers0,zLayer,nCell,k,TotCells,iPcl
  DOUBLE PRECISION :: PI
  REAL(RFREAL) :: VR !MGF
  REAL(RFREAL), DIMENSION(:) , ALLOCATABLE :: densGas,presGas,radius,Vol, &
					      uGas,vFracE,rho,theta,densGaspCell,PVF,zLoc,&
                                              pVelo,sos,Temp
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global



   INTEGER :: iCvSpecAir,iCvSpecExplosive,iCvSpecProducts,iSpec
   REAL(RFREAL) :: YProducts,YProducts_beg,YProducts_end
   REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
   REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: Yp,Inv_Yp,Mix_num
   TYPE(t_spec_input), POINTER :: pSpecInput
   TYPE(t_spec_type), POINTER :: pSpecType


! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataCyldet', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract Surf data:'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

  Jmax = pGrid%nCells/pRegion%patches(2)%nBQuads
  Imax = pRegion%patches(1)%nBQuads/Jmax
  nCellsSecIPerLayer = (Imax/4)*(Imax/4) ! No of cells within charge for every z-layer
  Kmax = (pRegion%patches(2)%nBQuads - nCellsSecIPerLayer)/ Imax 
  rLayers0 = Imax/8 ! No. of radial layers in Sec I
  TotCells = pGrid%nCells

  IF (global%specUsed .EQV. .TRUE.) THEN 
  pCvSpec    => pRegion%spec%cv
  pSpecInput => pRegion%specInput
  END IF


! ******************************************************************************
! Allocate memory
! ******************************************************************************

  ALLOCATE(densGas(Kmax+rLayers0),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3485,'densGas')
  END IF ! global%error 

  ALLOCATE(presGas(Kmax+rLayers0),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3491,'presGas')
  END IF ! global%error 

  ALLOCATE(radius(Kmax+rLayers0),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3497,'radius')
  END IF ! global%error 

  ALLOCATE(uGas(Kmax+rLayers0),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3503,'uGas')
  END IF ! global%error 

  ALLOCATE(vFracE(Kmax+rLayers0),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3509,'vFracE')
  END IF ! global%error 

  ALLOCATE(pVelo(Kmax+rLayers0),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3515,'vFracE')
  END IF ! global%error

  ALLOCATE(sos(Kmax+rLayers0),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3521,'sos')
  END IF ! global%error

  ALLOCATE(Temp(Kmax+rLayers0),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3527,'uGas')
  END IF ! global%error

  ALLOCATE(Vol(Kmax+rLayers0),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3533,'Vol')
  END IF ! global%error 

  !CellData
  ALLOCATE(rho(TotCells),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3540,'rho')
  END IF ! global%error

  ALLOCATE(theta(TotCells),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3546,'theta')
  END IF ! global%error

  ALLOCATE(zLoc(TotCells),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3552,'zLoc')
  END IF ! global%error

  !PartData
  ALLOCATE(densGaspCell(TotCells),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3559,'densGaspCell')
  END IF ! global%error

  ALLOCATE(PVF(TotCells),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3565,'PVF')
  END IF ! global%error

  IF (global%specUsed .EQV. .TRUE.) THEN
  ALLOCATE(Yp(Kmax+rLayers0),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3573,'presGas')
  END IF ! global%error

  ALLOCATE(Inv_Yp(Kmax+rLayers0),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3579,'presGas')
  END IF ! global%error

  ALLOCATE(Mix_num(Kmax+rLayers0),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,3585,'presGas')
  END IF ! global%error
  END IF !specUsed

! ********************************************************************************
! Initialize all variables to null
! ********************************************************************************

  DO rlayer = 1,Kmax+rLayers0
    densGas(rlayer) = 0.0_RFREAL
    presGas(rlayer) = 0.0_RFREAL
    radius(rlayer)  = 0.0_RFREAL
    uGas(rLayer)    = 0.0_RFREAL
    sos(rlayer)     = 0.0_RFREAL 
    Temp(rlayer)    = 0.0_RFREAL
    vFracE(rlayer)  = 0.0_RFREAL
    Vol(rlayer)     = 0.0_RFREAL
  END DO

  IF (global%specUsed .EQV. .TRUE.) THEN
  DO rlayer = 1,Kmax+rLayers0
    Yp(rlayer)       = 0.0_RFREAL 
    Inv_Yp(rlayer)   = 0.0_RFREAL
    Mix_num(rlayer)  = 0.0_RFREAL
  END DO

  iCvSpecAir = SPEC_GetSpeciesIndex(global,pSpecInput,'AIR')
  iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                   'PRODUCTS')
  END IF !specUsed



! ********************************************************************************
! Compute reqd properties - Sum over all cells loacted in a given radial location
! ********************************************************************************
 
! ============== Applicable to sections II and III - Polar Grid ==================

 IF (global%plagUsed .EQV. .TRUE.) THEN 
 END IF !plagUsed 

  DO rLayer = 1,Kmax
    icg = (rLayer-1)*(Imax) + 1
    radius(rLayer+rLayers0) = DSQRT(pGrid%cofg(XCOORD,icg)**2 &
                                  + pGrid%cofg(YCOORD,icg)**2)
    DO zLayer = 1,Jmax 
      icgBeg = (zLayer-1)*(Kmax*Imax + nCellsSecIPerLayer) + &
               (rLayer-1)*Imax + 1
      icgEnd = (zLayer-1)*(Kmax*Imax + nCellsSecIPerLayer) + &
               (rLayer)*Imax 
      DO icg = icgBeg,icgEnd
        densGas(rLayer+rLayers0) = densGas(rlayer+rLayers0) &
                                 + pRegion%mixt%cv(CV_MIXT_DENS,icg) * pGrid%vol(icg)
        presGas(rLayer+rlayers0) = presGas(rlayer+rLayers0) &
                                 + pRegion%mixt%dv(DV_MIXT_PRES,icg) * pGrid%vol(icg)
        Temp(rLayer+rlayers0)    = Temp(rlayer+rLayers0) &
                                 + pRegion%mixt%dv(DV_MIXT_TEMP,icg) * pGrid%vol(icg)
        sos(rLayer+rlayers0)     = sos(rlayer+rLayers0) &
                                 + pRegion%mixt%dv(DV_MIXT_SOUN,icg) * pGrid%vol(icg)

        uGas(rLayer+rlayers0)    = uGas(rlayer+rLayers0) &
                                 + DSQRT(pRegion%mixt%cv(CV_MIXT_XMOM,icg)**2 &
                                 + pRegion%mixt%cv(CV_MIXT_YMOM,icg)**2) &
			         *  pGrid%vol(icg)/pRegion%mixt%cv(CV_MIXT_DENS,icg)

 
        IF (global%specUsed .EQV. .TRUE.) THEN
          YProducts = pCvSpec(iCvSpecProducts,icg)
          Yp(rLayer+rLayers0) = Yp(rLayer+rLayers0) + (YProducts/pRegion%mixt%cv(CV_MIXT_DENS,icg)) &
                                                     * pGrid%vol(icg)
          Inv_Yp(rLayer+rLayers0) = Inv_Yp(rLayer+rLayers0) + &
                         (1.0_RFREAL-(YProducts/pRegion%mixt%cv(CV_MIXT_DENS,icg)))* pGrid%vol(icg)

          Mix_num(rLayer+rLayers0) = Mix_num(rLayer+rLayers0) + &
                  (YProducts/pRegion%mixt%cv(CV_MIXT_DENS,icg))* &
                  (1.0_RFREAL-(YProducts/pRegion%mixt%cv(CV_MIXT_DENS,icg))) * pGrid%vol(icg)   
         END IF !Spec 

        Vol(rLayer+rLayers0)     = Vol(rlayer+rLayers0) + pGrid%vol(icg)  
      END DO
    END DO
    densGas(rLayer+rLayers0) = densGas(rlayer+rLayers0)/Vol(rLayer+rLayers0) 
    presGas(rLayer+rLayers0) = presGas(rlayer+rLayers0)/Vol(rLayer+rLayers0) 
    uGas(rLayer+rLayers0)    = uGas(rlayer+rLayers0)/Vol(rLayer+rLayers0)
    Temp(rLayer+rLayers0)    = Temp(rlayer+rLayers0)/Vol(rLayer+rLayers0)
    sos(rLayer+rLayers0)     = sos(rlayer+rLayers0)/Vol(rLayer+rLayers0) 

    IF (global%plagUsed .EQV. .TRUE.) THEN
    Yp(rLayer+rLayers0)  = Yp(rlayer+rLayers0)/Vol(rLayer+rLayers0)
    Inv_Yp(rLayer+rLayers0)  = Inv_Yp(rlayer+rLayers0)/Vol(rLayer+rLayers0)
    Mix_num(rLayer+rLayers0)  = Mix_num(rlayer+rLayers0)/Vol(rLayer+rLayers0)
    END IF !Spec
  END DO

! ================== Applicable to section I - Cartesian Grid =================

  DO rLayer = 1,rLayers0
    icg = (Kmax*Imax) + (rLayer-1)*(Imax/4 + 1) + 1
    radius(rLayers0-rLayer+1) = DSQRT(pGrid%cofg(XCOORD,icg)**2 &
                              + pGrid%cofg(YCOORD,icg)**2)
    DO zLayer = 1,Jmax 
      tLayers0 = Imax/4 - 2*(rLayer-1)
      DO tLayer = 1,tLayers0
        icgBeg = zLayer*(Kmax*Imax) + (zLayer-1)*(nCellsSecIPerLayer) + &
                 (rLayer-1)*(Imax/4 + 1) + (tLayer-1)*(Imax/4) + 1
        icgEnd = icgBeg + (tLayers0-1)
        IF (tLayer == 1 .OR. tLayer == tLayers0) THEN
          DO icg = icgBeg,icgEnd
            densGas(rLayers0-rLayer+1) = densGas(rlayers0-rLayer+1) &
                                       + pRegion%mixt%cv(CV_MIXT_DENS,icg) &
                                       * pGrid%vol(icg)
            presGas(rLayers0-rLayer+1) = presGas(rlayers0-rLayer+1) &
                                       + pRegion%mixt%dv(DV_MIXT_PRES,icg) &
                                       * pGrid%vol(icg)
            Temp(rLayers0-rLayer+1)    = Temp(rlayers0-rLayer+1) &
                                       + pRegion%mixt%dv(DV_MIXT_TEMP,icg) &
                                       * pGrid%vol(icg)
            uGas(rLayers0-rlayer+1)    = uGas(rlayers0-rLayer+1) &
                                       + DSQRT(pRegion%mixt%cv(CV_MIXT_XMOM,icg)**2 &
                                       +       pRegion%mixt%cv(CV_MIXT_YMOM,icg)**2) &
		                       *  pGrid%vol(icg)/pRegion%mixt%cv(CV_MIXT_DENS,icg)
            sos(rLayers0-rLayer+1) = sos(rlayers0-rLayer+1) &
                                       + pRegion%mixt%dv(DV_MIXT_SOUN,icg) &
                                       * pGrid%vol(icg)

          IF (global%plagUsed .EQV. .TRUE.) THEN
          YProducts = pCvSpec(iCvSpecProducts,icg)
          Yp(rLayers0-rLayer+1) = Yp(rLayers0-rLayer+1)+(YProducts/pRegion%mixt%cv(CV_MIXT_DENS,icg)) &
                                                     * pGrid%vol(icg)
          Inv_Yp(rLayers0-rLayer+1) = Inv_Yp(rLayers0-rLayer+1) + &
                         (1.0_RFREAL-(YProducts/pRegion%mixt%cv(CV_MIXT_DENS,icg)))*pGrid%vol(icg)

          Mix_num(rLayers0-rLayer+1) = Mix_num(rLayers0-rLayer+1) + &
                  (YProducts/pRegion%mixt%cv(CV_MIXT_DENS,icg))* &
                  (1.0_RFREAL-(YProducts/pRegion%mixt%cv(CV_MIXT_DENS,icg)))*pGrid%vol(icg)

          END IF !Spec


            Vol(rLayers0-rLayer+1)     = Vol(rlayers0-rLayer+1) + pGrid%vol(icg)  
            END DO
        ELSE
          densGas(rLayers0-rLayer+1) = densGas(rlayers0-rLayer+1) &
                                     + ( pRegion%mixt%cv(CV_MIXT_DENS,icgBeg) &
                                     * pGrid%vol(icgBeg) ) &
                                     + ( pRegion%mixt%cv(CV_MIXT_DENS,icgEnd) &
                                     * pGrid%vol(icgEnd) )
          presGas(rLayers0-rLayer+1) = presGas(rlayers0-rLayer+1) &
                                     + ( pRegion%mixt%dv(DV_MIXT_PRES,icgBeg) &
                                     * pGrid%vol(icgBeg) ) &
                                     + ( pRegion%mixt%dv(DV_MIXT_PRES,icgEnd) &
                                     * pGrid%vol(icgEnd) )
          Temp(rLayers0-rLayer+1)    = Temp(rlayers0-rLayer+1) &
                                     + ( pRegion%mixt%dv(DV_MIXT_TEMP,icgBeg) &
                                     * pGrid%vol(icgBeg) ) &
                                     + ( pRegion%mixt%dv(DV_MIXT_TEMP,icgEnd) &
                                     * pGrid%vol(icgEnd) )
          uGas(rLayers0-rlayer+1)  = uGas(rlayers0-rLayer+1) &
                                     + ( SQRT(pRegion%mixt%cv(CV_MIXT_XMOM,icgBeg)**2 &
                                     +        pRegion%mixt%cv(CV_MIXT_YMOM,icgBeg)**2 ) &
		                     *  pGrid%vol(icgBeg)/pRegion%mixt%cv(CV_MIXT_DENS,icgBeg) ) &
                                     + ( SQRT(pRegion%mixt%cv(CV_MIXT_XMOM,icgEnd)**2 &
                                     +        pRegion%mixt%cv(CV_MIXT_YMOM,icgEnd)**2 ) &
		                     *  pGrid%vol(icgEnd)/pRegion%mixt%cv(CV_MIXT_DENS,icgEnd) )
          sos(rLayers0-rLayer+1) = sos(rlayers0-rLayer+1) &
                                     + ( pRegion%mixt%dv(DV_MIXT_SOUN,icgBeg) &
                                     * pGrid%vol(icgBeg) ) &
                                     + ( pRegion%mixt%dv(DV_MIXT_SOUN,icgEnd) &
                                     * pGrid%vol(icgEnd) )

          IF (global%plagUsed .EQV. .TRUE.) THEN
          YProducts_beg = pCvSpec(iCvSpecProducts,icgBeg)
          YProducts_end = pCvSpec(iCvSpecProducts,icgEnd)
          Yp(rLayers0-rLayer+1)    = Yp(rlayers0-rLayer+1) &
                                     + (YProducts_beg/pRegion%mixt%cv(CV_MIXT_DENS,icgBeg) &
                                     * pGrid%vol(icgBeg) ) &
                                     + (YProducts_end/pRegion%mixt%cv(CV_MIXT_DENS,icgEnd) &
                                     * pGrid%vol(icgEnd) )
          Inv_Yp(rLayers0-rLayer+1) = Inv_Yp(rLayers0-rLayer+1) + &
                         (1.0_RFREAL-(YProducts_beg/pRegion%mixt%cv(CV_MIXT_DENS,icgBeg)))*pGrid%vol(icgBeg)&
                       + (1.0_RFREAL-(YProducts_end/pRegion%mixt%cv(CV_MIXT_DENS,icgEnd)))*pGrid%vol(icgEnd)    

          Mix_num(rLayers0-rLayer+1) = Mix_num(rLayers0-rLayer+1) + &
                  (YProducts_beg/pRegion%mixt%cv(CV_MIXT_DENS,icgBeg))* &
                  (1.0_RFREAL-(YProducts_beg/pRegion%mixt%cv(CV_MIXT_DENS,icgBeg)))*pGrid%vol(icgBeg)+ &
                  (YProducts_end/pRegion%mixt%cv(CV_MIXT_DENS,icgEnd))* &
                  (1.0_RFREAL-(YProducts_end/pRegion%mixt%cv(CV_MIXT_DENS,icgEnd)))*pGrid%vol(icgEnd)
          END IF !specUsed 




          Vol(rLayers0-rLayer+1)     = Vol(rlayers0-rLayer+1) + pGrid%vol(icgBeg) &
                                     + pGrid%vol(icgEnd)
        END IF
      END DO
    END DO
    densGas(rLayers0-rLayer+1) = densGas(rlayers0-rLayer+1)/Vol(rLayers0-rLayer+1) 
    presGas(rLayers0-rLayer+1) = presGas(rlayers0-rLayer+1)/Vol(rLayers0-rLayer+1) 
    Temp(rLayers0-rLayer+1)    = Temp(rlayers0-rLayer+1)/Vol(rLayers0-rLayer+1)
    uGas(rLayers0-rLayer+1)    = uGas(rlayers0-rLayer+1)/Vol(rLayers0-rLayer+1)
    sos(rLayers0-rLayer+1)     = sos(rlayers0-rLayer+1)/Vol(rLayers0-rLayer+1) 

    IF (global%plagUsed .EQV. .TRUE.) THEN
    Yp(rLayers0-rLayer+1)     = Yp(rlayers0-rLayer+1)/Vol(rLayers0-rLayer+1)
    Inv_Yp(rLayers0-rLayer+1)     = Inv_Yp(rlayers0-rLayer+1)/Vol(rLayers0-rLayer+1)
    Mix_num(rLayers0-rLayer+1)     = Mix_num(rlayers0-rLayer+1)/Vol(rLayers0-rLayer+1)
    END IF !specUsed
  END DO

! ******************************************************************************
! Write volFracE & other ppts as a function of radius to data file
! ******************************************************************************

  WRITE(*,*) 'Writing data'
  WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.zavg_',global%currentTime,'.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)

  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,3921,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
   WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// &
                            TRIM(iFileName1)
  END IF

  DO rLayer = 1,Kmax+rLayers0
   IF (global%specUsed .EQV. .TRUE.) THEN
     WRITE(IF_EXTR_DATA1,'(11(1X,E23.16))') radius(rLayer),densGas(rLayer), &
                                          presGas(rLayer),Temp(rLayer),uGas(rLayer),&
                                          sos(rLayer),vFracE(rLayer),pVelo(rLayer),Yp(rLayer),&
                                          Inv_Yp(rLayer),Mix_num(rLayer)
   ELSE
    WRITE(IF_EXTR_DATA1,'(8(1X,E23.16))') radius(rLayer),densGas(rLayer), & 
                                          presGas(rLayer),Temp(rLayer),uGas(rLayer), &
                                          sos(rLayer),vFracE(rLayer),pVelo(rLayer)
  END IF !Spec
  END DO


!Fred - 2D Fractal Dimension Calculation: 5/29/19 

!End Fractal Dimension Calculator--------------------------------------------


! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(densGas,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4143,'densGas')
  END IF ! global%error

  DEALLOCATE(presGas,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4149,'presGas')
  END IF ! global%error

  DEALLOCATE(Temp,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4155,'presGas')
  END IF ! global%error

  DEALLOCATE(radius,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4161,'radius')
  END IF ! global%error

  DEALLOCATE(uGas,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4167,'uGas')
  END IF ! global%error

  DEALLOCATE(sos,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4173,'sos')
  END IF ! global%error

  DEALLOCATE(vFracE,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4179,'vFracE')
  END IF ! global%error

  DEALLOCATE(pVelo,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4185,'vFracE')
  END IF ! global%error


  DEALLOCATE(Vol,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4192,'Vol')
  END IF ! global%error

  DEALLOCATE(rho,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4198,'rho')
  END IF ! global%error

  DEALLOCATE(theta,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4204,'theta')
  END IF ! global%error

  DEALLOCATE(zLoc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4210,'zLoc')
  END IF ! global%error

  DEALLOCATE(densGaspCell,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4216,'densGaspCell')
  END IF ! global%error

  DEALLOCATE(PVF,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4222,'PVF')
  END IF ! global%error

  IF (global%specUsed .EQV. .TRUE.) THEN
  DEALLOCATE(Yp,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4230,'Vol')
  END IF ! global%error

  DEALLOCATE(Inv_Yp,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4236,'Vol')
  END IF ! global%error

  DEALLOCATE(Mix_num,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4242,'Vol')
  END IF ! global%error
  END IF !specUsed


! ******************************************************************************
! Close file
! ******************************************************************************

  
  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,4256,'File: '//TRIM(iFileName1))
  END IF ! global%error

  CLOSE(3245,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,4262,'File: '//TRIM(iFileName2))
  END IF ! global%error

  CLOSE(3246,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,4268,'File: '//TRIM(iFileName3))
  END IF ! global%error


  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataCyldet


! ******************************************************************************
! Saptarshi
!
! Purpose: Extract data for shock tube case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
!
! ******************************************************************************


SUBROUTINE RFLU_ExtractFlowDataShktb(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) ::iFileName1,iFileName2
  INTEGER :: errorFlag,Nx_g,Ny_g,Nz_g,i,j,k,icg

  REAL(RFREAL), DIMENSION(:) , ALLOCATABLE :: densGas,presGas,tempGas,veloGas
  REAL(RFREAL), DIMENSION(:) , ALLOCATABLE :: vFracE,pVelo

  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

  REAL(RFREAL) :: upf,dpf,coef,vflim
  INTEGER :: Np,Npe,iPcl,upfdone


    INTEGER :: iCvSpecAir,iCvSpecExplosive,iCvSpecProducts,iSpec
    REAL(RFREAL) :: YProducts,YProducts_beg,YProducts_end
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
    REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: Yp,Inv_Yp,Mix_num
    TYPE(t_spec_input), POINTER :: pSpecInput
    TYPE(t_spec_type), POINTER :: pSpecType
!Fred - Species quantities for averaging

! ******************************************************************************
! Start
! ******************************************************************************


  global => pRegion%global


  CALL RegisterFunction(global,'RFLU_ExtractFlowDataShktb', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract data:'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  Nx_g = pGrid%nCells/pRegion%patches(1)%nBQuads ! no of cells in x direction
  Ny_g = pGrid%nCells/pRegion%patches(4)%nBQuads ! no of cells in y direction
  Nz_g = pGrid%nCells/(Nx_g * Ny_g)              ! no of cells in z direction

      IF (global%specUsed .EQV. .TRUE.) THEN
        pCvSpec    => pRegion%spec%cv
        pSpecInput => pRegion%specInput
      END IF !specUsed

! ******************************************************************************
! Allocate memory
! ******************************************************************************

 ! Nx_g = 2000

  ALLOCATE(densGas(Nx_g),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,4431,'densGas')
  END IF ! global%error

  ALLOCATE(presGas(Nx_g),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,4437,'presGas')
  END IF ! global%error

  ALLOCATE(tempGas(Nx_g),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,4443,'tempGas')
  END IF ! global%error

  ALLOCATE(veloGas(Nx_g),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,4449,'veloGas')
  END IF ! global%error


       IF (global%specUsed .EQV. .TRUE.) THEN
        ALLOCATE(Yp(Nx_g),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
           CALL ErrorStop(global,ERR_ALLOCATE,4473,'presGas')
        END IF ! global%error

        ALLOCATE(Inv_Yp(Nx_g),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
           CALL ErrorStop(global,ERR_ALLOCATE,4479,'presGas')
        END IF ! global%error

        ALLOCATE(Mix_num(Nx_g),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,4485,'presGas')
        END IF ! global%error
       END IF !specUsed


! ********************************************************************************
! Initialize all variables to null
! ********************************************************************************

  densGas = 0.0_RFREAL
  presGas = 0.0_RFREAL
  tempGas = 0.0_RFREAL
  veloGas = 0.0_RFREAL
  vFracE  = 0.0_RFREAL
  pVelo   = 0.0_RFREAL

  IF (global%specUsed .EQV. .TRUE.) THEN
  Yp      = 0.0_RFREAL
  Inv_Yp  = 0.0_RFREAL
  Mix_num = 0.0_RFREAL

  iCvSpecAir = SPEC_GetSpeciesIndex(global,pSpecInput,'AIR')
  iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                 'PRODUCTS')
  END IF !specUsed
      
! ********************************************************************************
! Compute reqd properties - Sum over all spanwise cells at each x-location 
! ********************************************************************************

    DO k = 1,Nz_g
      DO j = 1,Ny_g
        DO i = 1,Nx_g
          icg  = i + (j-1)*Nx_g + (k-1)*Ny_g*Nx_g
          densGas(i) = densGas(i) + pRegion%mixt%cv(CV_MIXT_DENS,icg)&
                       /(Ny_g*Nz_g)
          presGas(i) = presGas(i) + pRegion%mixt%dv(DV_MIXT_PRES,icg)&
                       /(Ny_g*Nz_g)
          tempGas(i) = tempGas(i) + pRegion%mixt%dv(DV_MIXT_TEMP,icg)&
                       /(Ny_g*Nz_g)
          veloGas(i) = veloGas(i) + SQRT( &
                       pRegion%mixt%cv(CV_MIXT_XMOM,icg)**2 &
                     + pRegion%mixt%cv(CV_MIXT_YMOM,icg)**2 &
                     + pRegion%mixt%cv(CV_MIXT_ZMOM,icg)**2) &
                     / pRegion%mixt%cv(CV_MIXT_DENS,icg)/(Ny_g*Nz_g)


         IF (global%specUsed .EQV. .TRUE.) THEN
          YProducts = pCvSpec(iCvSpecProducts,icg)
          Yp(i) = Yp(i) + (YProducts/pRegion%mixt%cv(CV_MIXT_DENS,icg)) & 
                          /(Ny_g*Nz_g)
          Inv_Yp(i) = Inv_Yp(i) +(1.0_RFREAL- & 
                      (YProducts/pRegion%mixt%cv(CV_MIXT_DENS,icg))) & 
                         /(Ny_g*Nz_g)
          Mix_num(i) = Mix_num(i) + &
                  (YProducts/pRegion%mixt%cv(CV_MIXT_DENS,icg))* &   
                  (1.0_RFREAL-(YProducts/pRegion%mixt%cv(CV_MIXT_DENS,icg))) & 
                         /(Ny_g*Nz_g) 
         END IF !specUsed
        END DO
      END DO
    END DO


! ******************************************************************************
! Write parameters as function of x-location to data file
! ******************************************************************************

  WRITE(*,*) 'Writing data'
  WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.yzavg_',global%currentTime,'.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)

  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,4578,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
   WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// &
                            TRIM(iFileName1)
  END IF

  DO i = 1,Nx_g
    IF(global%plagUsed .EQV. .TRUE.)THEN
      IF (global%specUsed .EQV. .TRUE.) THEN !Fred - for species output 
  IF (global%specUsed .EQV. .TRUE.) THEN
      WRITE(IF_EXTR_DATA1,'(10(1X,E23.16))') pGrid%cofg(XCOORD,i),densGas(i), & 
                                          presGas(i),tempGas(i),veloGas(i), &
                                          vFracE(i),pVelo(i) &
                                          ,Yp(i),Inv_Yp(i),Mix_num(i)
  END IF !specUsed
      ELSE   
      WRITE(IF_EXTR_DATA1,'(7(1X,E23.16))') pGrid%cofg(XCOORD,i),densGas(i), &
                                          presGas(i),tempGas(i),veloGas(i), &
                                          vFracE(i),pVelo(i)
      END IF
    ELSE
      WRITE(IF_EXTR_DATA1,'(5(1X,E23.16))') pGrid%cofg(XCOORD,i),densGas(i), &
                                          presGas(i),tempGas(i),veloGas(i)
    END IF
  END DO

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,4611,'File: '//TRIM(iFileName1))
  END IF ! global%error



! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(densGas,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4700,'densGas')
  END IF ! global%error

  DEALLOCATE(presGas,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4706,'presGas')
  END IF ! global%error

  DEALLOCATE(tempGas,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4712,'tempGas')
  END IF ! global%error

  DEALLOCATE(veloGas,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4718,'veloGas')
  END IF ! global%error


 IF (global%specUsed .EQV. .TRUE.) THEN
  DEALLOCATE(Yp,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,4742,'Vol')  
  END IF ! global%error    

  DEALLOCATE(Inv_Yp,STAT=errorFlag) 
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4748,'Vol')
  END IF ! global%error 

  DEALLOCATE(Mix_num,STAT=errorFlag)   
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,4754,'Vol') 
  END IF ! global%error   
 END IF !specUsed
!Fred - End species additions  

! ******************************************************************************
! Finish up
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataShktb

! Saptarshi Ends
!*******************************************************************************


! ******************************************************************************
!
! Purpose: Extract Surf data for Cylds shock diffraction case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Assume domain to lie in z=constant plane.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSurfCylds(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: nExtract
  INTEGER :: errorFlag,icg,ix
  REAL(RFREAL) :: a,p,r,u,v,w,M,xx,yy,angle,ratio
  REAL(RFREAL) :: accX,cp,cp1,cp_acc,cp_potential,p_inf,q_inf,q_inf1,theta, &
                  Radius,u_rel,delTime
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSurfCylds', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract Surf data:'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  pPatch => pRegion%patches(1)

! ******************************************************************************
! Open file for data
! ******************************************************************************

  WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.surf_',global%currentTime,'.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,4848,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Compute number of cells
! ******************************************************************************

! Debug - Subbu
  nExtract = pPatch%nBFacesTot
! End Debug

! ******************************************************************************
! Extract along cylinder surface
! ******************************************************************************

  IF ( global%currentTime > global%mvfAccTs ) THEN
    delTime = global%currentTime - global%mvfAccTs
    accX    = global%mvfAccX
  ELSE
    delTime = 0.0_RFREAL
    accX    = 0.0_RFREAL
  END IF !

  IF ( global%currentTime > global%mvfAccTe ) THEN
    delTime = global%mvfAccTe - global%mvfAccTs
    accX    = 0.0_RFREAL
  END IF !

  u_rel  = global%refVelocity - (delTime)*(global%mvfAccX)
  q_inf  = 0.5_RFREAL*global%refDensity*u_rel*u_rel
  q_inf1 = 0.5_RFREAL*global%refDensity*global%refVelocity*global%refVelocity
  p_inf  = global%refPressure

  ratio  = q_inf1/q_inf 

  Radius = 1.000000_RFREAL

  DO ix = 1,nExtract
    icg  = pPatch%bf2c(ix)

    xx = pGrid%cofg(XCOORD,icg)
    yy = pGrid%cofg(YCOORD,icg)
! Debug - Subbu
    !angle = 180.0_RFREAL*ATAN2(yy,xx)/global%pi
    !theta = global%pi - global%pi*angle/180.0_RFREAL
    !angle = ATAN2((yy-4.0_RFREAL),xx)
    angle = ATAN2(yy,xx)
! End Debug 
    r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
    u = pRegion%mixt%cv(CV_MIXT_XMOM,icg)/r
    v = pRegion%mixt%cv(CV_MIXT_YMOM,icg)/r
    w = pRegion%mixt%cv(CV_MIXT_ZMOM,icg)/r
    p = pRegion%mixt%dv(DV_MIXT_PRES,icg)
    a = pRegion%mixt%dv(DV_MIXT_SOUN,icg)
    M = SQRT(u*u+v*v+w*w)/a
! Debug Subbu
    !cp_potential = 1.0_RFREAL &
    !             - 4.0_RFREAL*SIN(global%pi*angle/180.0_RFREAL)**2.0_RFREAL
    !cp_acc       = cp_potential + global%refDensity*Radius*(-accX)*COS(theta)/q_inf
    !cp           = (p-p_inf)/q_inf
! End Debug

    cp           =  pPatch%cp(ix)
! Debug Subbu
    !cp           =  pPatch%cp(ix)
    !cp1          = cp*ratio          ! ratio is defined earlier as q_inf/q_inf1
! End Debug 

    cp1           = (p-p_inf)/q_inf1
! END COMPUTATION OF Cp


! Debug Subbu
    !WRITE(IF_EXTR_DATA1,'(10(1X,E23.16))') xx,yy,angle,cp_potential,cp_acc,cp,cp1, &
    !                                       cp_acc-cp_potential, &
    !                                       cp-cp_potential, &
    !                                       cp1-cp_potential

! End Debug 
    WRITE(IF_EXTR_DATA1,'(3(1X,E23.16))') angle,cp,cp1
  END DO ! ix

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,4942,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract Surface data for Cylds case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSurfCylds








! ******************************************************************************
!
! Purpose: Extract Surf data for Cylinder near wall case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Assume domain to lie in z=constant plane.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSurfCylWall(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: nExtract
  INTEGER :: errorFlag,icg,ix
  REAL(RFREAL) :: a,p,r,u,v,w,M,xx,yy,angle,ratio
  REAL(RFREAL) :: accX,cp,cp1,cp_acc,cp_potential,p_inf,q_inf,q_inf1,theta, &
                  Radius,u_rel,delTime
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSurfCylWall', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract Surf data:'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  pPatch => pRegion%patches(1)
! ******************************************************************************
! Open file for data
! ******************************************************************************

  WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.surf_',global%currentTime,'.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,5039,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Compute number of cells
! ******************************************************************************

  nExtract = pPatch%nBFaces

! ******************************************************************************
! Extract along cylinder surface
! ******************************************************************************

  IF ( global%currentTime > global%mvfAccTs ) THEN
    delTime = global%currentTime - global%mvfAccTs
    accX    = global%mvfAccX
  ELSE
    delTime = 0.0_RFREAL
    accX    = 0.0_RFREAL
  END IF !

  IF ( global%currentTime > global%mvfAccTe ) THEN
    delTime = global%mvfAccTe - global%mvfAccTs
    accX    = 0.0_RFREAL
  END IF !

  u_rel  = global%refVelocity - (delTime)*(global%mvfAccX)
  q_inf  = 0.5_RFREAL*global%refDensity*u_rel*u_rel
  q_inf1 = 0.5_RFREAL*global%refDensity*global%refVelocity*global%refVelocity
  p_inf  = global%refPressure

  Radius = 1.000000_RFREAL

  DO ix = 1,nExtract
    icg  = pPatch%bf2c(ix)
    xx = pGrid%cofg(XCOORD,icg)
    yy = pGrid%cofg(YCOORD,icg)

    angle = ATAN2(yy,xx)

    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
     p   = pRegion%mixt%cv(CV_MIXT_PRES,icg)
    ELSE
     p   = pRegion%mixt%dv(DV_MIXT_PRES,icg)
    END IF

    cp  = pPatch%cp(ix)
    cp1 = (p-p_inf)/q_inf1
    WRITE(IF_EXTR_DATA1,'(3(1X,E23.16))') angle,cp,cp1
  END DO ! ix

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,5102,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract Surface data for CylWall case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSurfCylWall








! ******************************************************************************
!
! Purpose: Extract total kinetic energy for cylinder and axi-symmetric cases. 
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataVolumeCylds(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: errorFlag,icg,indCp,indMol,ir,nRadii
  REAL(RFREAL) :: cp,d,dh,dp,ds,Eo,factor,g,gc,h,keLocal,keTotal,L,momLocal, &
                  momTotal,mw,p,pe,peLocal,peTotal,pi,r,re,refD,refH,refL, &
                  refNu,refP,refU,refT,radialDist,radius,ru,rv,rw,t,u,ue,v,ve, &
                  vol,Vm2,w,we,x,y,xTransformed,yTransformed,radialDistTransformed,&
                  refA,PressCE
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: keRadii,momRadii,peRadii,radii
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  REAL(RFREAL) :: ksg

! ******************************************************************************
! Start, Set pointers and variables
! ******************************************************************************

  global => pRegion%global
  pGrid  => pRegion%grid

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataVolumeCylds', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract KE data for cylds problem:'
  END IF ! global%verbLevel

! ******************************************************************************
! Allocate memory for radial distances 
! ******************************************************************************

  nRadii = 150

  ALLOCATE(radii(nRadii),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,5195,'radii')
  END IF ! global%error 

  ALLOCATE(keRadii(nRadii),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,5201,'keRadii')
  END IF ! global%error 

  ALLOCATE(momRadii(nRadii),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,5207,'momRadii')
  END IF ! global%error 

  ALLOCATE(peRadii(nRadii),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,5213,'peRadii')
  END IF ! global%error 

  radius = 1.0_RFREAL !TEMPORARY - Subbu: radius = 1.0 - By Default radius = 0.5_RFREAL
  factor = 1.5_RFREAL/(150_RFREAL**3.0_RFREAL)

  radii(1) = radius + factor*(1)**3.0_RFREAL + .015_RFREAL
  DO ir=2,nRadii
    radii(ir) = radii(ir-1) + factor*(ir)**3.0_RFREAL + .015_RFREAL
    keRadii(ir)  = 0.0_RFREAL
    momRadii(ir) = 0.0_RFREAL
    peRadii(ir)  = 0.0_RFREAL
  END DO ! icg

! ******************************************************************************
! Extract ke from all cells 
! ******************************************************************************

  keTotal  = 0.0_RFREAL
  keLocal  = 0.0_RFREAL
  momTotal = 0.0_RFREAL
  momLocal = 0.0_RFREAL
  peTotal  = 0.0_RFREAL
  peLocal  = 0.0_RFREAL

  pi  = global%pi

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    refL  = global%refLength
    refNu = global%refVisc/global%refDensity
    refU  = global%refVelocity
    refD  = global%refDensity
    refP  = global%refPressure
    
    d  = refD
    L  = pRegion%mixtInput%prepRealVal1

    g = global%refGamma

    refH = (g/(g-1.0_RFREAL))*global%refPressure/global%refDensity

    DO icg = 1,pGrid%nCells
      vol = pRegion%grid%vol(icg)

      r = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      u = pRegion%mixt%cv(CV_MIXT_XVEL,icg) - pRegion%mixtInput%iniVelX
      v = pRegion%mixt%cv(CV_MIXT_YVEL,icg) - pRegion%mixtInput%iniVelY
      w = pRegion%mixt%cv(CV_MIXT_ZVEL,icg) - pRegion%mixtInput%iniVelZ

      x = pGrid%cofg(XCOORD,icg)
      y = pGrid%cofg(YCOORD,icg)
      radialDist = SQRT(x*x+y*y)

      IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
        keLocal = 0.5_RFREAL*(u*u + v*v + w*w)*r*vol*2.0_RFREAL*pi*y
        peLocal = 0.0_RFREAL
      ELSE
        keLocal = 0.5_RFREAL*(u*u + v*v + w*w)*r*vol
        peLocal = 0.0_RFREAL
      END IF ! pRegion%mixtInput%axiFlag 

      keTotal = keTotal + keLocal 
      peTotal = peTotal + peLocal 

      DO ir=1,nRadii
        IF ( radii(ir) > radialDist ) THEN
          keRadii(ir) = keRadii(ir) + keLocal
          peRadii(ir) = peRadii(ir) + peLocal
        END IF ! radii
      END DO ! ir
    END DO ! icg
  ELSE
    DO icg = 1,pGrid%nCells
      vol = pRegion%grid%vol(icg)

      mw = pRegion%mixt%gv(GV_MIXT_MOL,indMol*icg)
      cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *icg)

      gc = MixtPerf_R_M(mw)
      g  = MixtPerf_G_CpR(cp,gc)

      refH = (g/(g-1.0_RFREAL))*global%refPressure/global%refDensity
      refT = MixtPerf_T_DPR(global%refDensity,global%refPressure,gc)

      r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      ir = 1.0_RFREAL/r
      ru = pRegion%mixt%cv(CV_MIXT_XMOM,icg)
      rv = pRegion%mixt%cv(CV_MIXT_YMOM,icg)
      rw = pRegion%mixt%cv(CV_MIXT_ZMOM,icg)

      Eo = pRegion%mixt%cv(CV_MIXT_ENER,icg)/r

      Vm2 = (ru*ru + rv*rv + rw*rw)/(r*r)
      
      ksg = pRegion%mixt%piclKsg(icg)

      p = MixtPerf_P_DEoGVm2(r,Eo,g,Vm2,ksg)

      dp = p - global%refPressure
      dh = (g/(g-1.0_RFREAL))*p/r - refH
      ds = (log(p/global%refPressure) - g*log(r/global%refDensity))* &
           gc/(g-1.0_RFREAL)

      u = ru/r - pRegion%mixtInput%iniVelX
      v = rv/r - pRegion%mixtInput%iniVelY
      w = rw/r - pRegion%mixtInput%iniVelZ

      x = pGrid%cofg(XCOORD,icg)
      y = pGrid%cofg(YCOORD,icg)
      radialDist = SQRT(x*x+y*y)

      keLocal = 0.5_RFREAL*(u*u + v*v + w*w)*r*vol
      peLocal = (r*(dh - refT*ds) - dp)*vol

      momLocal = u*r*vol

      IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
        keLocal = keLocal*2.0_RFREAL*pi*y
        peLocal = peLocal*2.0_RFREAL*pi*y

        momLocal = momLocal*2.0_RFREAL*pi*y
      END IF ! pRegion%mixtInput%axiFlag 

      keTotal = keTotal + keLocal 
      peTotal = peTotal + peLocal 

      momTotal = momTotal + momLocal 

      DO ir=1,nRadii
        IF ( radii(ir) > radialDist ) THEN
          keRadii(ir) = keRadii(ir) + keLocal
          peRadii(ir) = peRadii(ir) + peLocal

          momRadii(ir) = momRadii(ir) + momLocal
        END IF ! radii
      END DO ! ir
    END DO ! icg
  END IF ! solverType

! ******************************************************************************
! Open file for data
! ******************************************************************************

  IF ( global%flowType == FLOW_STEADY ) THEN
    WRITE(iFileName1,'(A,I5.5,A)') TRIM(global%outDir)// &
                                TRIM(global%casename)// &
                                '.ke_',global%currentIter,'.dat'
  ELSE
    WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                                TRIM(global%casename)// &
                                '.ke_',global%currentTime,'.dat'
  END IF ! global%flowType  

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,5370,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Write data
! ******************************************************************************

  WRITE(*,*) 'Outer Boundary is:',50.0,'radii(nRadii):',radii(nRadii)
  WRITE(*,*) 'Total KE is:',keTotal,'keRadii(nRadii):',keRadii(nRadii)
  WRITE(*,*) 'Total PE is:',peTotal,'peRadii(nRadii):',peRadii(nRadii)
  WRITE(*,*) 'Total Mom is:',momTotal,'momRadii(nRadii):',momRadii(nRadii)
  DO ir=1,nRadii
    WRITE(IF_EXTR_DATA1,'(7(1X,E23.16))') radii(ir),keRadii(ir), peRadii(ir), &
                                          (keRadii(ir)+peRadii(ir)),momRadii(ir), &
                                          x,y
  END DO ! ir


! ******************************************************************************
! TEMPORARY: Subbu - To Extract Coefficient of pressure
! ******************************************************************************
 
  IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN
    refL  = global%refLength
    refU  = global%refVelocity
    refD  = global%refDensity
    refA  = global%forceRefArea
    refP  = global%refPressure
   
    WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                                TRIM(global%casename)// &
                                '.Cp_',global%currentTime,'.dat'
    
    OPEN(IF_EXTR_DATA1,FILE=iFileName1)
    
    DO icg=1,pGrid%nCells
       x = pGrid%cofg(XCOORD,icg)
       y = pGrid%cofg(YCOORD,icg)
       xTransformed = x
       yTransformed = y+4
       radialDistTransformed = SQRT(xTransformed*xTransformed +yTransformed*yTransformed)
       p = MixtPerf_P_DEoGVm2(r,Eo,g,Vm2,ksg)
       !IF( radialDistTransformedi == refL/2)THEN
          PressCE = (p-refP)/(0.5_RFREAL * refD * refU * refU * refA)
          WRITE(IF_EXTR_DATA1,'(4(1X,E23.16))') xTransformed,yTransformed,radialDistTransformed,& 
                                                PressCE
       !END IF
    END DO
    CLOSE(IF_EXTR_DATA1)
  
  END IF
! ******************************************************************************
! END TEMPORARY: Subbu
! ******************************************************************************



! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,5439,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract KE data for cylds done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataVolumeCylds 










! ******************************************************************************
!
! Purpose: Extract Surf data for Sphds case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Assume domain to lie in z=constant plane.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractFlowDataSurfSphds(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: nExtract
  INTEGER :: distrib,errorFlag,icg,ix
  REAL(RFREAL) :: a,p,r,u,v,w,M,M2,Mroot,xx,yy,zz,angle,ratio,cosTheta
  REAL(RFREAL) :: accX,cp,cp1,cp2,cp_acc,cp_potential,length,p_inf,q_inf, &
                  q_inf1,theta,Radius,u_rel,delTime
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractFlowDataSurfSphds', &
                        "../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract Surf data from sphere:'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  pPatch => pRegion%patches(1)

  distrib = pPatch%mixt%distrib

  Radius = 0.0125_RFREAL

  M = pRegion%patches(2)%mixt%vals(BCDAT_FARF_MACH,distrib*ix)
  M2 = M*M ;
  Mroot = sqrt(1-M2) ;

  nExtract = pPatch%nBFaces

! ******************************************************************************
! Open file for data
! ******************************************************************************

  WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.surf_',global%currentTime,'.plt'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,5549,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to file: '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Extract along cylinder surface
! ******************************************************************************

  IF ( global%currentTime > global%mvfAccTs ) THEN
    delTime = global%currentTime - global%mvfAccTs
    accX    = global%mvfAccX
  ELSE
    delTime = 0.0_RFREAL
    accX    = 0.0_RFREAL
  END IF !

  IF ( global%currentTime > global%mvfAccTe ) THEN
    delTime = global%mvfAccTe - global%mvfAccTs
    accX    = 0.0_RFREAL
  END IF !

  u_rel  = global%refVelocity - (delTime)*(global%mvfAccX)
  q_inf  = 0.5_RFREAL*global%refDensity*u_rel*u_rel
  q_inf1 = 0.5_RFREAL*global%refDensity*global%refVelocity*global%refVelocity
  p_inf  = global%refPressure

  ratio  = q_inf/q_inf1 

  DO ix = 1,nExtract
    icg  = pPatch%bf2c(ix)

    xx = pGrid%cofg(XCOORD,icg)
    yy = pGrid%cofg(YCOORD,icg)
    zz = pGrid%cofg(ZCOORD,icg)

    cosTheta = xx/(xx*xx+yy*yy+zz*zz)**0.5_RFREAL
    angle    = ACOS(cosTheta)*180.0_RFREAL/global%pi

    cp_potential = 0.25_RFREAL*(9.0_RFREAL*cosTheta**2 - 5.0_RFREAL)
    cp_acc       = cp_potential &
                 + 0.5_RFREAL*global%refDensity*Radius*(-accX)*cosTheta/q_inf
    cp           = pPatch%cp(ix)
    cp1 = cp_potential/Mroot ;
    cp2 = cp_potential/(Mroot + (M2/(1+Mroot))*cp_potential/2 ) 

! END COMPUTATION OF Cp

    WRITE(IF_EXTR_DATA1,'(9(1X,E23.16))') xx,yy,zz,angle,cp_potential,cp_acc,cp,cp1,cp2
  END DO ! ix

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,5610,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract Surface data for Sphds case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractFlowDataSurfSphds










! ******************************************************************************
!
! Purpose: Extract Mesh boundaries for Bump case.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!
! Output: None.
!
! Notes:
!   1. Assume domain to lie in z=constant plane.
!
! ******************************************************************************

SUBROUTINE RFLU_WriteMeshBump(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName1
  INTEGER :: iPatch,iv1,iv2
  INTEGER :: errorFlag,icg,ix
  REAL(RFREAL) :: a,p,r,u,v,w,M,xx,yy,angle
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_WriteMeshBump',"../../utilities/post/RFLU_ModExtractFlowData.F90")

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Extract boundary grid points data:'
  END IF ! global%verbLevel

! ******************************************************************************
! Open file for data
! ******************************************************************************

  WRITE(iFileName1,'(A,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '_bdry','.dat'

  OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,5700,'File: '//TRIM(iFileName1))
  END IF ! global%error

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Writing data to  '// &
                             TRIM(iFileName1)
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Extract along cylinder surface
! ******************************************************************************

!  DO iPatch = 1,pGrid%nPatches
  DO iPatch = 1,4
    pPatch => pRegion%patches(iPatch)
     
    WRITE(IF_EXTR_DATA1,*) pPatch%nBFaces+1 
     
    IF ( iPatch < 3 ) THEN
      ix = 1
      iv1 = pPatch%bv(pPatch%bf2v(1,ix))
      xx = pGrid%xyz(XCOORD,iv1) 
      yy = pGrid%xyz(YCOORD,iv1) 
      WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') xx,yy

      DO ix = 1,pPatch%nBFaces
        iv2 = pPatch%bv(pPatch%bf2v(2,ix))
        xx = pGrid%xyz(XCOORD,iv2) 
        yy = pGrid%xyz(YCOORD,iv2) 
        WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') xx,yy
      END DO ! ix
    ELSE
      ix = pPatch%nBFaces
      iv1 = pPatch%bv(pPatch%bf2v(1,ix))
      xx = pGrid%xyz(XCOORD,iv1) 
      yy = pGrid%xyz(YCOORD,iv1) 
      WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') xx,yy

      DO ix = pPatch%nBFaces,1,-1
        iv2 = pPatch%bv(pPatch%bf2v(2,ix))
        xx = pGrid%xyz(XCOORD,iv2) 
        yy = pGrid%xyz(YCOORD,iv2) 
        WRITE(IF_EXTR_DATA1,'(2(1X,E23.16))') xx,yy
      END DO ! ix
    END IF
  END DO ! iPatch

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,5760,'File: '//TRIM(iFileName1))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extract boundary grid points for Bump case done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_WriteMeshBump







END MODULE RFLU_ModExtractFlowData

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModExtractFlowData.F90,v $
! Revision 1.4  2016/05/05 22:25:13  rahul
! Bug fix in 'taylorvortex' case. Declared reference values outside the loop
! that checks type of flow solver.
!
! Revision 1.3  2015/08/12 19:41:43  brollin
! Updating module declaration in rfluinit.F90
!
! Revision 1.2  2015/07/27 04:45:42  brollin
! 1) Corrected bug in RFLUCONV where global%gridFormat was used instead of global%gridSrcFormat
! 2) Implemented new subroutine for shock tube problems (Shktb)
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.11  2009/07/08 19:12:29  mparmar
! Added RFLU_ExtractFlowDataVolumeCylds and modified RFLU_ExtractFlowDataTaylorVortex
!
! Revision 1.10  2008/12/06 08:43:58  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:12  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2008/05/29 01:35:33  mparmar
! Added extraction of density profiles for laminar flat plate case
!
! Revision 1.7  2008/01/20 22:56:14  haselbac
! Fix problem when disc cannot be located
!
! Revision 1.6  2008/01/17 13:12:52  haselbac
! Updated somm_spi case to extract reflected shock
!
! Revision 1.5  2008/01/10 18:47:32  haselbac
! Enabled extraction of particle contact for SommSPI cases
!
! Revision 1.4  2007/11/28 23:05:50  mparmar
! Added acoustic and taylorvortex, modified cylds cp data
!
! Revision 1.3  2007/06/18 18:14:23  mparmar
! Added data extractions cylinder and sphere data
!
! Revision 1.2  2007/04/12 12:13:03  haselbac
! Updated to take into account new shock extraction routine
!
! Revision 1.21  2007/04/05 01:14:28  haselbac
! Added USE of new extract module, added xs extraction for stg1d
!
! Revision 1.20  2007/03/27 00:47:27  haselbac
! Added extraction of particle data to Sommerfeld case
!
! Revision 1.19  2007/03/19 21:43:28  haselbac
! Fixed typo in somm_spi write statement
!
! Revision 1.18  2007/03/02 17:56:33  haselbac
! Updated SommSPI data extraction to allow for 1d/2d cases
!
! Revision 1.17  2007/02/27 13:23:40  haselbac
! Added stg1d case
!
! Revision 1.16  2007/02/17 20:57:20  haselbac
! Added shock-position extraction to Sommerfeld case
!
! Revision 1.15  2007/02/16 20:01:29  haselbac
! Added code for somm_spi case
!
! Revision 1.14  2006/08/19 15:41:17  mparmar
! Added data extractions for nscbc[1-8], farf, bumpq10
!
! Revision 1.13  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.12  2006/01/04 15:48:54  haselbac
! Changed computation of thicknesses for lfp
!
! Revision 1.11  2005/11/27 02:00:22  haselbac
! Added extraction for EEv, bug fix
!
! Revision 1.10  2005/11/17 14:50:48  haselbac
! Added extraction of mass fractions for STG
!
! Revision 1.9  2005/11/14 17:05:24  haselbac
! Added support for pseudo-gas model
!
! Revision 1.8  2005/11/11 17:21:04  haselbac
! Added extraction for generic 2d shocktube
!
! Revision 1.7  2005/11/10 02:51:30  haselbac
! Added support for MP sod cases, clean-up
!
! Revision 1.6  2005/10/09 15:11:40  haselbac
! Added ONERA C0 data extraction
!
! Revision 1.5  2005/07/19 19:19:25  haselbac
! Added new lfp cases, generalized and fixed extraction for lfp
!
! Revision 1.4  2005/06/14 01:09:50  haselbac
! Adapted to new case names for Skews case
!
! Revision 1.3  2005/03/18 23:09:34  haselbac
! Added routine to extract data for Skews case
!
! Revision 1.2  2004/10/27 18:08:51  haselbac
! Cosmetics only
!
! Revision 1.1  2004/10/26 15:19:11  haselbac
! Initial revision
!
! ******************************************************************************

