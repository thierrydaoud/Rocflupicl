#include "PPICLF_STD.h"
#:include 'PPICLF_PARTMACROS.fypp'
module ppiclf_initsolve

    ! particle data
    use ppiclf_data, only: ppiclf_npart

    use ppiclf_m_particledata, only: @{USEMODVAR(PPICLF_t_particle, ppiclf_parts)}@
    ! grid data
    use ppiclf_data, only: ppiclf_int_fld
    ! particle options variables
    use ppiclf_data, only: ppiclf_lcomm, ppiclf_linit, ppiclf_ndim
    ! comm variables
    use ppiclf_data, only: 
    ! binning variables
    use ppiclf_data, only: 
    ! ghost particle variables
    use ppiclf_data, only: 
    ! wall support variables
    use ppiclf_data, only: ppiclf_nwall, ppiclf_wall_c, ppiclf_wall_n


    ! used functions/subroutines
    use ppiclf_op, only: ppiclf_exittr

    use ppiclf_user, only: ppiclf_user_InitZero

    implicit none
    private

    public :: ppiclf_solve_InitWall
    public :: ppiclf_solve_InitZero

    contains

    SUBROUTINE ppiclf_solve_InitZero
        !
        ! Internal:
        !
        INTEGER*4 i,j,ie
        !
        ! zero'ing real particle properties
        DO i=1, PPICLF_LPART
#:for particle, n in fyppmacros.Loop_All_Reals("ppiclf_parts(i)")
            DO j=1, ${n}$
                ${particle}$(j) = 0.0
            END DO
#:endfor
#:for particle, n in fyppmacros.Loop_All_Ints("ppiclf_parts(i)")
            DO j=1, ${n}$
                ${particle}$(j) = 0.0
            END DO
#:endfor
        END DO

        ppiclf_npart = 0
        ! zero'ing grid properties for interpolation
        DO ie=1,PPICLF_LEE
            DO j=1,PPICLF_LRP_INT
                ppiclf_int_fld(j,ie) = 0.0D0
            END DO
        END DO

        CALL ppiclf_user_InitZero

        RETURN
    END SUBROUTINE ppiclf_solve_InitZero

    SUBROUTINE ppiclf_solve_InitWall(xp1,xp2,xp3)
        ! 
        ! Input:
        ! 
        REAL*8 xp1(*)
        REAL*8 xp2(*)
        REAL*8 xp3(*)
        !
        ! Internal:
        !
        REAL*8 rpx1, rpy1, rpz1, rpx2, rpy2, rpz2, a_sum, theta, tri_area, ab_dot_ac, ab_mag, ac_mag, rise, run, rmag, rpx3, rpy3, rpz3
        INTEGER*4 istride, k, kmax, kp, kkp, kk
        REAL*8 A(3),B(3),C(3),AB(3),AC(3)
        !
        if (.not.PPICLF_LCOMM) CALL ppiclf_exittr('InitMPI must be before InitWall$',0.d0,0)
        if (.not.PPICLF_LINIT) CALL ppiclf_exittr('InitParticle must be before InitWall$',0.d0,0)

        ppiclf_nwall = ppiclf_nwall + 1 

        if (ppiclf_nwall .gt. PPICLF_LWALL) CALL ppiclf_exittr('Increase LWALL in user file$',0.d0,ppiclf_nwall)

        istride = ppiclf_ndim
        a_sum = 0.0d0
        kmax = 2
        if (ppiclf_ndim .EQ. 3) kmax = 3

        if (ppiclf_ndim .EQ. 3) then
            ppiclf_wall_c(1,ppiclf_nwall) = xp1(1)
            ppiclf_wall_c(2,ppiclf_nwall) = xp1(2)
            ppiclf_wall_c(3,ppiclf_nwall) = xp1(3)
            ppiclf_wall_c(4,ppiclf_nwall) = xp2(1)
            ppiclf_wall_c(5,ppiclf_nwall) = xp2(2)
            ppiclf_wall_c(6,ppiclf_nwall) = xp2(3)
            ppiclf_wall_c(7,ppiclf_nwall) = xp3(1)
            ppiclf_wall_c(8,ppiclf_nwall) = xp3(2)
            ppiclf_wall_c(9,ppiclf_nwall) = xp3(3)

            A(1) = (xp1(1) + xp2(1) + xp3(1))/3.0d0
            A(2) = (xp1(2) + xp2(2) + xp3(2))/3.0d0
            A(3) = (xp1(3) + xp2(3) + xp3(3))/3.0d0
        elseif (ppiclf_ndim .EQ. 2) then
            ppiclf_wall_c(1,ppiclf_nwall) = xp1(1)
            ppiclf_wall_c(2,ppiclf_nwall) = xp1(2)
            ppiclf_wall_c(3,ppiclf_nwall) = xp2(1)
            ppiclf_wall_c(4,ppiclf_nwall) = xp2(2)

            A(1) = (xp1(1) + xp2(1))/2.0d0
            A(2) = (xp1(2) + xp2(2))/2.0d0
            A(3) = 0.0d0
        ENDif

        ! compute area:
        do k=1,kmax 
            kp = k+1
            if (kp .gt. kmax) kp = kp-kmax ! cycle
            
            kk   = istride*(k-1)
            kkp  = istride*(kp-1)
            rpx1 = ppiclf_wall_c(kk+1,ppiclf_nwall)
            rpy1 = ppiclf_wall_c(kk+2,ppiclf_nwall)
            rpz1 = 0.0d0
            rpx2 = ppiclf_wall_c(kkp+1,ppiclf_nwall)
            rpy2 = ppiclf_wall_c(kkp+2,ppiclf_nwall)
            rpz2 = 0.0d0

            B(1) = rpx1
            B(2) = rpy1
            B(3) = 0.0d0
            
            C(1) = rpx2
            C(2) = rpy2
            C(3) = 0.0d0
            
            AB(1) = B(1) - A(1)
            AB(2) = B(2) - A(2)
            AB(3) = 0.0d0
            
            AC(1) = C(1) - A(1)
            AC(2) = C(2) - A(2)
            AC(3) = 0.0d0

            if (ppiclf_ndim .EQ. 3) then
                rpz1 = ppiclf_wall_c(kk+3,ppiclf_nwall)
                rpz2 = ppiclf_wall_c(kkp+3,ppiclf_nwall)
                B(3) = rpz1
                C(3) = rpz2
                AB(3) = B(3) - A(3)
                AC(3) = C(3) - A(3)
            
                AB_DOT_AC = AB(1)*AC(1) + AB(2)*AC(2) + AB(3)*AC(3)
                AB_MAG = sqrt(AB(1)**2 + AB(2)**2 + AB(3)**2)
                AC_MAG = sqrt(AC(1)**2 + AC(2)**2 + AC(3)**2)
                theta  = acos(AB_DOT_AC/(AB_MAG*AC_MAG))
                tri_area = 0.5d0*AB_MAG*AC_MAG*sin(theta)
            elseif (ppiclf_ndim .EQ. 2) then
                AB_MAG = sqrt(AB(1)**2 + AB(2)**2)
                tri_area = AB_MAG
            ENDif
            a_sum = a_sum + tri_area
        ENDdo
    
        ppiclf_wall_n(ppiclf_ndim+1,ppiclf_nwall) = a_sum

        ! wall normal:
        if (ppiclf_ndim .EQ. 2) then

            rise = xp2(2) - xp1(2)
            run  = xp2(1) - xp1(1)

            rmag = sqrt(rise**2 + run**2)
            rise = rise/rmag
            run  = run/rmag
            
            ppiclf_wall_n(1,ppiclf_nwall) = rise
            ppiclf_wall_n(2,ppiclf_nwall) = -run

        elseif (ppiclf_ndim .EQ. 3) then

            k  = 1
            kk = istride*(k-1)
            rpx1 = ppiclf_wall_c(kk+1,ppiclf_nwall)
            rpy1 = ppiclf_wall_c(kk+2,ppiclf_nwall)
            rpz1 = ppiclf_wall_c(kk+3,ppiclf_nwall)
            
            k  = 2
            kk = istride*(k-1)
            rpx2 = ppiclf_wall_c(kk+1,ppiclf_nwall)
            rpy2 = ppiclf_wall_c(kk+2,ppiclf_nwall)
            rpz2 = ppiclf_wall_c(kk+3,ppiclf_nwall)
            
            k  = 3
            kk = istride*(k-1)
            rpx3 = ppiclf_wall_c(kk+1,ppiclf_nwall)
            rpy3 = ppiclf_wall_c(kk+2,ppiclf_nwall)
            rpz3 = ppiclf_wall_c(kk+3,ppiclf_nwall)
        
            A(1) = rpx2 - rpx1
            A(2) = rpy2 - rpy1
            A(3) = rpz2 - rpz1

            B(1) = rpx3 - rpx2
            B(2) = rpy3 - rpy2
            B(3) = rpz3 - rpz2

            ppiclf_wall_n(1,ppiclf_nwall) = A(2)*B(3) - A(3)*B(2)
            ppiclf_wall_n(2,ppiclf_nwall) = A(3)*B(1) - A(1)*B(3)
            ppiclf_wall_n(3,ppiclf_nwall) = A(1)*B(2) - A(2)*B(1)

            rmag = sqrt(ppiclf_wall_n(1,ppiclf_nwall)**2 + ppiclf_wall_n(2,ppiclf_nwall)**2 + ppiclf_wall_n(3,ppiclf_nwall)**2)

            ppiclf_wall_n(1,ppiclf_nwall) = ppiclf_wall_n(1,ppiclf_nwall)/rmag
            ppiclf_wall_n(2,ppiclf_nwall) = ppiclf_wall_n(2,ppiclf_nwall)/rmag
            ppiclf_wall_n(3,ppiclf_nwall) = ppiclf_wall_n(3,ppiclf_nwall)/rmag

        END IF

        RETURN
    END SUBROUTINE ppiclf_solve_InitWall

end module ppiclf_initsolve