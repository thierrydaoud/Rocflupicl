#define PPICLF_LPART 20000
#define PPICLF_LRS 12
#define PPICLF_LRP 48
#define PPICLF_LEE 60000
#define PPICLF_LEX 2
#define PPICLF_LEY 2
#define PPICLF_LEZ 2
#define PPICLF_LRP_INT 30
#define PPICLF_LRP_PRO 22

! number of timesteps kept in history kernels
#define PPICLF_VU 50
! maximum number of triangular patch boundaries
#define PPICLF_LWALL 20

! y, y1, ydot, ydotc: PPICLF_LRS
#define PPICLF_JX 1
#define PPICLF_JY 2
#define PPICLF_JZ 3
#define PPICLF_JVX 4
#define PPICLF_JVY 5
#define PPICLF_JVZ 6
#define PPICLF_JT 7
#define PPICLF_JOX 8
#define PPICLF_JOY 9
#define PPICLF_JOZ 10
#define PPICLF_JMETAL 11
#define PPICLF_JOXIDE 12

! rprop: PPICLF_LRP
#define PPICLF_R_JRHOP 1
#define PPICLF_R_JRHOF 2
#define PPICLF_R_JDP 3
#define PPICLF_R_JVOLP 4
#define PPICLF_R_JPHIP 5
#define PPICLF_R_JUX 6
#define PPICLF_R_JUY 7
#define PPICLF_R_JUZ 8
#define PPICLF_R_JCS 9
#define PPICLF_R_JDPDX 10
#define PPICLF_R_JDPDY 11
#define PPICLF_R_JDPDZ 12
#define PPICLF_R_JSDRX 13
#define PPICLF_R_JSDRY 14
#define PPICLF_R_JSDRZ 15
#define PPICLF_R_JRHSR 16
#define PPICLF_R_JPGCX 17 
#define PPICLF_R_JPGCY 18
#define PPICLF_R_JPGCZ 19
#define PPICLF_R_JDPi  20
#define PPICLF_R_JDPe  21
#define PPICLF_R_JSPL 22
#define PPICLF_R_JSPT 23
#define PPICLF_R_JT   24
#define PPICLF_R_FLUCTFX 25
#define PPICLF_R_FLUCTFY 26
#define PPICLF_R_FLUCTFZ 27
#define PPICLF_R_WDOTX 28
#define PPICLF_R_WDOTY 29
#define PPICLF_R_WDOTZ 30
#define PPICLF_R_JXVOR 31
#define PPICLF_R_JYVOR 32
#define PPICLF_R_JZVOR 33
#define PPICLF_R_JIDP 34
#define PPICLF_R_JBRNT 35
#define PPICLF_R_JP 36
#define PPICLF_R_JRHOGX 37
#define PPICLF_R_JRHOGY 38
#define PPICLF_R_JRHOGZ 39
#define PPICLF_R_JDPVDX 40
#define PPICLF_R_JDPVDY 41
#define PPICLF_R_JDPVDZ 42
#define PPICLF_R_JSDOX 43
#define PPICLF_R_JSDOY 44
#define PPICLF_R_JSDOZ 45
#define PPICLF_R_XIPAR 46
#define PPICLF_R_XIPERP 47
#define PPICLF_R_XIT 48

! rprop4: PPICLF_LRP4 - Reynolds Subgrid Stress Components
#define PPICLF_R_JRSG11 1
#define PPICLF_R_JRSG12 2
#define PPICLF_R_JRSG13 3
#define PPICLF_R_JRSG21 4
#define PPICLF_R_JRSG22 5
#define PPICLF_R_JRSG23 6
#define PPICLF_R_JRSG31 7
#define PPICLF_R_JRSG32 8
#define PPICLF_R_JRSG33 9
#define PPICLF_R_JTSG1 10
#define PPICLF_R_JTSG2 11
#define PPICLF_R_JTSG3 12

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

! map: PPICLF_LRP_PRO
#define PPICLF_P_JPHIP 1
!--- x,y,z Forces Fedback to Rocflu
#define PPICLF_P_JFX 2
#define PPICLF_P_JFY 3
#define PPICLF_P_JFZ 4
!---
#define PPICLF_P_JE 5
!--- Add comment about these terms 
#define PPICLF_P_JPHIPD 6
#define PPICLF_P_JPHIPU 7
#define PPICLF_P_JPHIPV 8
#define PPICLF_P_JPHIPW 9
#define PPICLF_P_JPHIPT 10
!--- Reynolds Subgrid Stress Tensor
#define PPICLF_P_JRSG11 11
#define PPICLF_P_JRSG12 12
#define PPICLF_P_JRSG13 13
#define PPICLF_P_JRSG21 14
#define PPICLF_P_JRSG22 15
#define PPICLF_P_JRSG23 16
#define PPICLF_P_JRSG31 17
#define PPICLF_P_JRSG32 18
#define PPICLF_P_JRSG33 19
!--- Pseudo Turbulent Kinetic Energy
#define PPICLF_P_JTSG1 20
#define PPICLF_P_JTSG2 21
#define PPICLF_P_JTSG3 22
