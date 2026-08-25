#include "PPICLF_USER.h"
#include "PPICLF_STD.h"

!=====================================================================
! PPICLF_PERFORMANCE.h
!
! Per-operation timing + count instrumentation for ppiclF.
!
! ALL instrumentation is compiled in only when the library is built
! with PERF=1 (see Makefile, which then adds -DPERF). With PERF unset
! every timing site is an empty #ifdef block, so the production build
! is byte-for-byte unaffected.
!
! TIMER SEMANTICS (all REAL*8, units = seconds, MPI_WTIME based)
! --------------------------------------------------------------
! Each timer ACCUMULATES the wall time spent in its operation over the
! current logging interval. ppiclf_solve_LogPerformance() reduces the
! timers across ranks, writes one CSV row, and zeroes the interval
! accumulators (call it once per step, or every N steps for amortized
! numbers; PPICLF_PERF_NSTEP records how many steps are in the row).
!
!   LEAF timers (mutually exclusive; sum <= TTotal):
!     TCreateBin TFindPart TLoadBalance TEmptyInd TInterfaceInd
!     TRankBounds TMapOverlap TCreateGhost TMoveGhost TsubbinRealMap
!     TsubbinGhostMap TsubbinFineMap TsubbinCellMap TPCNNSearch
!     TPPNNSearch TProject TInterp TIntegrate
!
!   COMMUNICATION timers (separate bucket, mutually exclusive w/
!   leaves) -- six channels, one per synchronization site:
!     TMPI_moveInt  - interp-field crystal transfer (InterpTupleTransfer)
!     TMPI_movePro  - projection crystal transfer (ProjectParticleGrid)
!     TMPI_moveGP   - ghost-particle crystal transfer (MoveGhostPartLB)
!     TMPI_moveRP   - real-particle crystal transfer (MoveParticlePartLB)
!     TMPI_moveOvlp - overlap-cell crystal transfer (MapOverlapGrid)
!     TMPI_allreduces - every timed MPI_ALLREDUCE / BCAST, incl. the
!                     per-stage remap_stale reduction in PostTimeStep
!   Each transfer is CARVED OUT of (subtracted from) the leaf that
!   contains it, so the leaves are compute-only. NOTE (measured,
!   112-rank blast case): TMPI_moveInt and TMPI_movePro are WAIT-
!   dominated (r = -0.91 vs the interp gather, r = -1.00 vs the P2C
!   search) -- they mirror upstream compute imbalance, they do not
!   measure per-cell communication work. Treat them as imbalance
!   indicators, not comm cost, and keep them out of the LB fit.
!
!   ADDITIONAL leaf timers (appended as CSV columns 32-35 so existing
!   column positions are unchanged):
!     TPeriodicShift - ppiclf_solve_PeriodicParticleShift body
!     TRemovePart    - ppiclf_solve_RemoveParticle body (carved out of
!                      TPCNNSearch, which brackets both P2C map calls)
!     TLBCalib       - online LB coefficient calibration
!                      (LBCalibAccum + LBCalibrate at end of stage)
!     TEntrySync     - column 35: OUTSIDE TTotal. Wait in the optional
!                      entry barrier (ppiclf_perf_sync=.TRUE.) = host
!                      fluid-solve imbalance absorbed before ppiclF
!                      timing starts. Not part of the Unaccounted sum.
!
!   SPLIT-INSTRUMENTATION leaves (CSV columns 36-38):
!     TInterpGather  - InitInterp + InterpField gather loop, closed
!                      BEFORE the interp crystal transfer. Overlap-
!                      cell proportional (r = +0.995 vs cells sent);
!                      feeds LB calibration channel 5. TInterp is now
!                      the IDW loop ONLY (linear in Np, channel 1).
!     TSortInt       - crystal_tuple_sort of the received interp
!                      tuples, carved out of TMPI_moveInt so that
!                      channel is crystal wait + wire only. Local
!                      cell-proportional compute; channel 5.
!     TMoveReal      - MoveParticlePartLB pack/unpack + rtemp
!                      alloc, transfer excluded. PERF-only
!                      diagnostic; previously the largest untimed
!                      contributor to Unaccounted.
!
!       -> Unaccounted = TTotal - sum(leaves 1..18)
!                        - (TUserYdot - TPPNNSearch)
!                        - sum(TMPI_*) - TPeriodicShift - TRemovePart
!                        - TLBCalib - TInterpGather - TSortInt
!                        - TMoveReal
!          (the force channels TQuasiSteady/TAddedMass/TPresGrad/
!          THeatTransfer are nested INSIDE the TUserYdot bracket and
!          are informational only -- never add them to the sum).
!
!   TTotal - wall time of the full per-step picl advance
!            (ppiclf_solve_IntegrateParticle body). The denominator.
!
! COUNT SEMANTICS (all INTEGER*4) - per-rank instantaneous snapshot
! taken in LogPerformance, reduced with MAX (and SUM for the mean):
!     T_RealPart      = PPICLF_NPART          (real particles / rank)
!     T_GhostPartRec  = PPICLF_NPART_GP       (ghosts held / rank)
!     T_GhostPartSent = (proxy) ghosts held   ! TODO set exact send cnt
!     T_OverlapCells  = PPICLF_NCELLS_INTERP  (overlap cells / rank)
!     T_FVCells       = (proxy) overlap cells ! TODO set exact FV cnt
!     T_GlbBins       = PPICLF_TOTALBINS       (global bin count)
!     T_LocalBins     = (proxy) 0             ! TODO set local bin cnt
!=====================================================================

!---------------------------------------------------------------------
! CALIBRATION CHANNEL TIMERS: the five timing channels that feed the
! online load-balance coefficient calibration (ch1: IDW interp +
! integrate + user ydot + real-particle map; ch2: P2P search + fine
! map; ch3: P2C map; ch4: projection; ch5: overlap-cell map/transfer
! + interp gather + tuple sort -- the moveInt/movePro transfers are
! EXCLUDED as wait-dominated) are UNCONDITIONAL source code - no build flag
! required - so the coefficients adapt at run time in every build.
! PERF gates only the full instrumentation and the CSV logging.
! Overhead of the always-on channels: a few dozen MPI_WTIME pairs per
! stage plus one pair per particle in the P2P search, well under 0.5%
! of a stage. Set ppiclf_LB_docal = .FALSE. to freeze the
! coefficients at run time.
!---------------------------------------------------------------------
! Time per operation per stage
      REAL*8  PPICLF_TCreateBin
     >       ,PPICLF_TFindPart
     >       ,PPICLF_TLoadBalance
     >       ,PPICLF_TEmptyInd
     >       ,PPICLF_TInterfaceInd
     >       ,PPICLF_TRankBounds
     >       ,PPICLF_TMapOverlap
     >       ,PPICLF_TCreateGhost
     >       ,PPICLF_TMoveGhost
     >       ,PPICLF_TsubbinRealMap
     >       ,PPICLF_TsubbinGhostMap
     >       ,PPICLF_TsubbinFineMap
     >       ,PPICLF_TsubbinCellMap
     >       ,PPICLF_TPCNNSearch
     >       ,PPICLF_TPPNNSearch
     >       ,PPICLF_TProject
     >       ,PPICLF_TInterp
     >       ,PPICLF_TMPI_allreduces
     >       ,PPICLF_TMPI_moveRP
     >       ,PPICLF_TMPI_moveGP
     >       ,PPICLF_TMPI_moveInt
     >       ,PPICLF_TMPI_movePro
     >       ,PPICLF_TMPI_moveOvlp
     >       ,PPICLF_TIntegrate
     >       ,PPICLF_TTotal
     >       ,PPICLF_TQuasiSteady
     >       ,PPICLF_TAddedMass
     >       ,PPICLF_TPresGrad
     >       ,PPICLF_THeatTransfer
     >       ,PPICLF_TUserYdot
     >       ,PPICLF_TIO
     >       ,PPICLF_TPeriodicShift
     >       ,PPICLF_TRemovePart
     >       ,PPICLF_TLBCalib
     >       ,PPICLF_TEntrySync
     >       ,PPICLF_TInterpGather
     >       ,PPICLF_TSortInt
     >       ,PPICLF_TMoveReal

      COMMON /PPICLF_RUNTIMES/ PPICLF_TCreateBin
     >       ,PPICLF_TFindPart
     >       ,PPICLF_TLoadBalance
     >       ,PPICLF_TEmptyInd
     >       ,PPICLF_TInterfaceInd
     >       ,PPICLF_TRankBounds
     >       ,PPICLF_TMapOverlap
     >       ,PPICLF_TCreateGhost
     >       ,PPICLF_TMoveGhost
     >       ,PPICLF_TsubbinRealMap
     >       ,PPICLF_TsubbinGhostMap
     >       ,PPICLF_TsubbinFineMap
     >       ,PPICLF_TsubbinCellMap
     >       ,PPICLF_TPCNNSearch
     >       ,PPICLF_TPPNNSearch
     >       ,PPICLF_TProject
     >       ,PPICLF_TInterp
     >       ,PPICLF_TMPI_allreduces
     >       ,PPICLF_TMPI_moveRP
     >       ,PPICLF_TMPI_moveGP
     >       ,PPICLF_TMPI_moveInt
     >       ,PPICLF_TMPI_movePro
     >       ,PPICLF_TMPI_moveOvlp
     >       ,PPICLF_TIntegrate
     >       ,PPICLF_TTotal
     >       ,PPICLF_TQuasiSteady
     >       ,PPICLF_TAddedMass
     >       ,PPICLF_TPresGrad
     >       ,PPICLF_THeatTransfer
     >       ,PPICLF_TUserYdot
     >       ,PPICLF_TIO
     >       ,PPICLF_TPeriodicShift
     >       ,PPICLF_TRemovePart
     >       ,PPICLF_TLBCalib
     >       ,PPICLF_TEntrySync
     >       ,PPICLF_TInterpGather
     >       ,PPICLF_TSortInt
     >       ,PPICLF_TMoveReal

      INTEGER*4  PPICLF_T_RealPart
     >          ,PPICLF_T_GhostPartSent
     >          ,PPICLF_T_GhostPartRec
     >          ,PPICLF_T_FVCells
     >          ,PPICLF_T_OverlapCells_sent
     >          ,PPICLF_T_OverlapCells_received
     >          ,PPICLF_T_GlbBins
     >          ,PPICLF_T_LocalBins

      COMMON /PPICLF_TIMERCOUNTS/ PPICLF_T_RealPart
     >          ,PPICLF_T_GhostPartSent
     >          ,PPICLF_T_GhostPartRec
     >          ,PPICLF_T_FVCells
     >          ,PPICLF_T_OverlapCells_sent
     >          ,PPICLF_T_OverlapCells_received
     >          ,PPICLF_T_GlbBins
     >          ,PPICLF_T_LocalBins

! Internal bookkeeping for the performance logger (not user data).
      INTEGER*4  PPICLF_PERF_NSTEP
     >          ,PPICLF_PERF_HDR
     >          ,PPICLF_PERF_UNIT
      COMMON /PPICLF_PERF_INTERNAL/ PPICLF_PERF_NSTEP
     >          ,PPICLF_PERF_HDR
     >          ,PPICLF_PERF_UNIT
