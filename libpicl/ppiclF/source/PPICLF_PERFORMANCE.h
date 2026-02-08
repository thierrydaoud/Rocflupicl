#include "PPICLF_USER.h"
#include "PPICLF_STD.h"

! Time per operation per stage
      REAL*8  PPICLF_TCreateBin
     >       ,PPICLF_TSendParticles
     >       ,PPICLF_TSendGridOverlap
     >       ,PPICLF_TSendFluidFields
     >       ,PPICLF_TParticleParticleModels
     >       ,PPICLF_TFluidParticleModels
     >       ,PPICLF_TSendGhostParticles
     >       ,PPICLF_TMapParticlesCells
     >       ,PPICLF_TInterpolation 
     >       ,PPICLF_TProjection
     >       ,PPICLF_TWriteSolution
     >       ,PPICLF_TIntegration
     >       ,PPICLF_TPeriodicity
     >       ,PPICLF_TDataTransfers
     >       ,PPICLF_TTotal

      COMMON /PPICLF_RUNTIMES/ PPICLF_TCreateBin
     >       ,PPICLF_TSendParticles
     >       ,PPICLF_TSendGridOverlap
     >       ,PPICLF_TSendFluidFields
     >       ,PPICLF_TParticleParticleModels
     >       ,PPICLF_TFluidParticleModels
     >       ,PPICLF_TSendGhostParticles
     >       ,PPICLF_TMapParticlesCells
     >       ,PPICLF_TInterpolation 
     >       ,PPICLF_TProjection
     >       ,PPICLF_TWriteSolution
     >       ,PPICLF_TIntegration
     >       ,PPICLF_TPeriodicity
     >       ,PPICLF_TDataTransfers
     >       ,PPICLF_TTotal
