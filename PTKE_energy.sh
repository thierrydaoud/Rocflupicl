#!/bin/bash

# this script open all the files in which Energy term was modified to account for PTKE

f1=libfloflu/MixtPerf_E.F90
f2=libfloflu/MixtPerf_P.F90
f3=libfloflu/MixtGasLiq_E.F90
f4=libfloflu/MixtPerf_H.F90
f5=libfloflu/MixtPerf_C.F90
f6=libfloflu/MixtPerf_D.F90

f7=libflu/RFLU_EnforceBounds.F90
f8=libflu/RFLU_CheckPositivity.F90
f9=libflu/RFLU_SetDependentVars.F90

f10=modflu/RFLU_ModConvertCv.F90
f11=modflu/RFLU_ModBoundConvertCv.F90
f12=modflu/RFLU_ModNSCBC.F90

f13=utilities/init/RFLU_InitFlowHardCode.F90
f14=utilities/init/RFLU_InitFlowHardCodeLim.F90

f15=rocspecies/SPEC_RFLU_ModPBA.F90

f16=modfloflu/ModInterfacesMixt.F90 

vim -p $f1 $f2 $f3 $f4 $f5 $f6
#vim -p $f7 $f8 $f9
#vim -p $f10 $f11 $f12
#vim -p $f13 $f14 $f15

