module ppiclf_op
    use mpi
    ! comm variables
    use ppiclf_data, only: ppiclf_comm, ppiclf_nid
    implicit none
    private

    public :: ppiclf_gop
    public :: ppiclf_igop
    public :: ppiclf_iglsum
    public :: ppiclf_glsum
    public :: ppiclf_glmax
    public :: ppiclf_iglmax
    public :: ppiclf_glmin
    public :: ppiclf_iglmin
    public :: ppiclf_vlmin
    public :: ppiclf_vlmax
    public :: ppiclf_copy
    public :: ppiclf_icopy
    public :: ppiclf_chcopy
    public :: ppiclf_exittr
    public :: ppiclf_printsri
    public :: ppiclf_printsi
    public :: ppiclf_printsr
    public :: ppiclf_prints
    public :: PPICLF_BLANK
    public :: PPICLF_INDX1
    public :: PPICLF_CHSTR
    public :: ppiclf_byte_open_mpi
    public :: ppiclf_byte_read_mpi
    public :: ppiclf_byte_write_mpi
    public :: ppiclf_byte_close_mpi
    public :: ppiclf_byte_set_view
    public :: ppiclf_bcast

    public fixed_str_len
    ! internal module variables
    integer, parameter :: fixed_str_len = 132
    contains

    ! group op ? gather op ?
    SUBROUTINE ppiclf_gop( x, w, op, n)
        !
        ! Input:
        !
        real*8 x(n), w(n)
        character*3 op
        integer*4 n
        !
        ! Internal:
        !
        integer*4 i, ie
        !
        if (op.eq.'+  ') then
            call mpi_allreduce(x,w,n,MPI_DOUBLE_PRECISION,mpi_sum,ppiclf_comm,ie)
        elseif (op.EQ.'M  ') then
            call mpi_allreduce(x,w,n,MPI_DOUBLE_PRECISION,mpi_max,ppiclf_comm,ie)
        elseif (op.EQ.'m  ') then
            call mpi_allreduce(x,w,n,MPI_DOUBLE_PRECISION,mpi_min,ppiclf_comm,ie)
        elseif (op.EQ.'*  ') then
            call mpi_allreduce(x,w,n,MPI_DOUBLE_PRECISION,mpi_prod,ppiclf_comm,ie)
        endif

        do i=1,n
            x(i) = w(i)
        enddo

        return
    END SUBROUTINE ppiclf_gop

    ! integer gop
    SUBROUTINE ppiclf_igop( x, w, op, n)

        !
        ! Input:
        !
        integer*4 x(n), w(n)
        character*3 op
        integer*4 n
        !
        ! Internal:
        !
        integer*4 i, ierr
        !
        if     (op.eq.'+  ') then
            call MPI_Allreduce (x,w,n,mpi_integer,mpi_sum ,ppiclf_comm,ierr)
        elseif (op.EQ.'M  ') then
            call MPI_Allreduce (x,w,n,mpi_integer,mpi_max ,ppiclf_comm,ierr)
        elseif (op.EQ.'m  ') then
            call MPI_Allreduce (x,w,n,mpi_integer,mpi_min ,ppiclf_comm,ierr)
        elseif (op.EQ.'*  ') then
            call MPI_Allreduce (x,w,n,mpi_integer,mpi_prod,ppiclf_comm,ierr)
        endif

        do i=1,n
            x(i) = w(i)
        enddo

        return
    END SUBROUTINE ppiclf_igop

    ! this method of having 
    integer*4 FUNCTION ppiclf_iglsum(a,n)

        ! 
        ! Input:
        ! 
        integer*4 a(n)
        integer*4 n
        ! 
        ! Internal:
        ! 
        integer*4 tsum
        integer*4 tmp(n),work(n)
        integer*4 i
        !
        tsum= 0
        do i=1,n
            tsum=tsum+a(i)
        enddo
        tmp(1)=tsum
        call ppiclf_igop(tmp,work,'+  ',1)
        ppiclf_iglsum=tmp(1)
        return
    END FUNCTION ppiclf_iglsum

    real*8 FUNCTION ppiclf_glsum (x,n)

        ! 
        ! Vars:
        ! 
        ! TLJ changed dimension of A
        real*8 x(n),tsum
        real*8 tmp(1),work(1)
        integer*4 i,n
        !
        TSUM = 0.0d0
        DO I=1,N
            TSUM = TSUM+X(I)
        END DO
        TMP(1)=TSUM
        CALL ppiclf_GOP(TMP,WORK,'+  ',1)
        ppiclf_GLSUM = TMP(1)
        return
    END FUNCTION ppiclf_glsum

    real*8 FUNCTION ppiclf_glmax(a,n)

        ! 
        ! Vars:
        ! 
        ! TLJ changed dimension of A
        REAL*8 A(n),tmax
        real*8 TMP(1),WORK(1)
        integer*4 i,n
        !
        TMAX=-99.0e20
        DO I=1,N
            TMAX=MAX(TMAX,A(I))
        END DO
        TMP(1)=TMAX
        CALL ppiclf_GOP(TMP,WORK,'M  ',1)
        ppiclf_GLMAX=TMP(1)
        return
    END FUNCTION ppiclf_glmax

    integer*4 FUNCTION ppiclf_iglmax(a,n)

        ! 
        ! Vars:
        ! 
        ! TLJ changed dimension of A
        integer*4 a(n),tmax
        integer*4 tmp(1),work(1)
        integer*4 i,n
        !
        tmax= -999999999
        do i=1,n
            tmax=max(tmax,a(i))
        enddo
        tmp(1)=tmax
        call ppiclf_igop(tmp,work,'M  ',1)
        ppiclf_iglmax=tmp(1)
        return
    END FUNCTION ppiclf_iglmax

    real*8 FUNCTION ppiclf_glmin(a,n)

        ! 
        ! Vars:
        ! 
        ! TLJ changed dimension of A
        REAL*8 A(n),tmin
        real*8 TMP(1),WORK(1)
        integer*4 i,n
        !
        TMIN=99.0e20
        DO I=1,N
            TMIN=MIN(TMIN,A(I))
        END DO
        TMP(1)=TMIN
        CALL ppiclf_GOP(TMP,WORK,'m  ',1)
        ppiclf_GLMIN = TMP(1)
        return
    END FUNCTION ppiclf_glmin

    integer*4 FUNCTION ppiclf_iglmin(a,n)
        ! 
        ! Vars:
        ! 
        ! TLJ changed dimension of a
        integer*4 a(n),tmin
        integer*4 tmp(1),work(1)
        integer*4 i, n
        !
        tmin=  999999999
        do i=1,n
            tmin=min(tmin,a(i))
        enddo
        tmp(1)=tmin
        call ppiclf_igop(tmp,work,'m  ',1)
        ppiclf_iglmin=tmp(1)
        return
    END FUNCTION ppiclf_iglmin

    real*8 FUNCTION ppiclf_vlmin(vec,n)

        ! 
        ! Vars:
        ! 
        ! TLJ changed dimension of VEC
        REAL*8 VEC(n),tmin
        integer*4 i, n
        !
        TMIN = 99.0E20
        DO I=1,N
            TMIN = MIN(TMIN,VEC(I))
        END DO
        ppiclf_VLMIN = TMIN
        return
    END FUNCTION ppiclf_vlmin

    real*8 FUNCTION ppiclf_vlmax(vec,n)
        ! 
        ! Vars:
        ! 
        ! TLJ changed dimension of VEC
        REAL*8 VEC(n),tmax
        integer*4 i, n
        !
        TMAX =-99.0E20
        do i=1,n
            TMAX = MAX(TMAX,VEC(I))
        enddo
        ppiclf_VLMAX = TMAX
        return
    END FUNCTION ppiclf_vlmax

    SUBROUTINE ppiclf_copy(a,b,n)
        ! 
        ! Vars:
        ! 
        ! TLJ changed dimension of A and B
        real*8 a(n),b(n)
        integer*4 i,n
        !

        do i=1,n
            a(i)=b(i)
        enddo

        return
    END SUBROUTINE ppiclf_copy

    SUBROUTINE ppiclf_icopy(a,b,n)
        ! 
        ! Vars:
        ! 
        ! TLJ changed dimension of A and B
        INTEGER*4 A(n), B(n)
        integer*4 i,n
        !
        DO I = 1, N
            A(I) = B(I)
        END DO
        return
    END SUBROUTINE ppiclf_icopy

    SUBROUTINE ppiclf_chcopy(a,b,n)
        ! 
        ! Vars:
        ! 
        ! TLJ changed A and B dimenions
        CHARACTER*1 A(n), B(n)
        integer*4 i,n
        !
        DO I = 1, N
            A(I) = B(I)
        END DO
        return
    END SUBROUTINE ppiclf_chcopy
    
    SUBROUTINE ppiclf_exittr(stringi,rdata,idata)
        ! 
        ! Vars:
        ! 
        character(len=*) stringi
        character(len=fixed_str_len) stringo
        character*25 s25
        integer*4 ilen, ierr, k, idata
        real*8 rdata
        !
        call ppiclf_blank(stringo,fixed_str_len)
        call ppiclf_chcopy(stringo,stringi,len(stringi))
        ilen = len(stringi) !ppiclf_indx1(stringo,'$')
        write(s25,25) rdata,idata
        25 format(1x,1p1e14.6,i10)
        call ppiclf_chcopy(stringo(ilen + 1:),s25,25)

        ! changed to always print the error message, not just on rank 0
        write(6,11) ppiclf_nid, stringo(1:ilen+24)
        11 format('PPICLF: ERROR on rank ',i2,":",a)

        !     call mpi_finalize (ierr)
        call mpi_abort(ppiclf_comm, 1, ierr)

        return
    END SUBROUTINE ppiclf_exittr

    SUBROUTINE ppiclf_printsri(stringi,rdata,idata)
        !
        ! Vars:
        !
        character(len=*) stringi
        character(len=fixed_str_len) stringo
        character*25 s25
        integer*4 ilen, idata, k, ierr
        real*8 rdata
        !
#ifdef TEST
        return
#endif
        call ppiclf_blank(stringo,fixed_str_len)
        call ppiclf_chcopy(stringo,stringi,len(stringi))
        ilen = ppiclf_indx1(stringo,'$')
        write(s25,25) rdata,idata
        25 format(1x,1p1e14.6,i10)
        call ppiclf_chcopy(stringo(ilen:),s25,25)

        call mpi_barrier(ppiclf_comm,ierr)

        if (ppiclf_nid.eq.0) write(6,12) stringo(1:ilen+24)
        12 format('PPICLF: ',a)

        call mpi_barrier(ppiclf_comm,ierr)

        return
    END SUBROUTINE ppiclf_printsri

    SUBROUTINE ppiclf_printsi(stringi,idata)
        !
        ! Vars:
        !
        character(len=*) stringi
        character(len=fixed_str_len) stringo
        character*10 s10
        integer*4 ilen, idata, k, ierr
        !
#ifdef TEST
        return
#endif
        call ppiclf_blank(stringo,fixed_str_len)
        call ppiclf_chcopy(stringo,stringi,len(stringi))
        ilen = ppiclf_indx1(stringo,'$')
        write(s10,10) idata
        10 format(1x,i9)
        call ppiclf_chcopy(stringo(ilen:),s10,10)

        call mpi_barrier(ppiclf_comm,ierr)

        if (ppiclf_nid.eq.0) write(6,13) stringo(1:ilen+9)
        13 format('PPICLF: ',a)

        call mpi_barrier(ppiclf_comm,ierr)

        return
    END SUBROUTINE ppiclf_printsi

    SUBROUTINE ppiclf_printsr(stringi,rdata)
        !
        ! Vars:
        !
        character(len=*) stringi
        character(len=fixed_str_len) stringo
        character*15 s15
        integer*4 ilen, k, ierr
        real*8 rdata
        !
#ifdef TEST
        return
#endif
        call ppiclf_blank(stringo,fixed_str_len)
        call ppiclf_chcopy(stringo,stringi,len(stringi))
        ilen = ppiclf_indx1(stringo,'$')
        write(s15,15) rdata
        15 format(1x,1p1e14.6)
        call ppiclf_chcopy(stringo(ilen:),s15,15)

        call mpi_barrier(ppiclf_comm,ierr)

        if (ppiclf_nid.eq.0) write(6,14) stringo(1:ilen+14)
        14 format('PPICLF: ',a)

        call mpi_barrier(ppiclf_comm,ierr)

        return
    END SUBROUTINE ppiclf_printsr

    SUBROUTINE ppiclf_prints(stringi)
        !
        ! Vars:
        !
        character(len=*) stringi
        character(len=fixed_str_len) stringo
        integer*4 ilen, k, ierr
        !
#ifdef TEST
        return
#endif
        call ppiclf_blank(stringo,fixed_str_len)
        call ppiclf_chcopy(stringo,stringi,len(stringi))
        ilen = ppiclf_indx1(stringo,'$')

        call mpi_barrier(ppiclf_comm,ierr)

        if (ppiclf_nid.eq.0) write(6,21) stringo(1:ilen-1)
        21 format('PPICLF: ',a)

        call mpi_barrier(ppiclf_comm,ierr)

        return
    END SUBROUTINE ppiclf_prints

    SUBROUTINE PPICLF_BLANK(A,N)
        ! 
        ! Vars:
        !
        ! TLJ changed dimension of A
        CHARACTER(len=*), intent(out) :: A
        CHARACTER*1 BLNK
        SAVE        BLNK
        DATA        BLNK /' '/
        integer*4 i,n
        !
        !
        DO I=1,N
            A(I:I)=BLNK
        END DO
        RETURN
    END SUBROUTINE PPICLF_BLANK

    ! finds the first index of a character(S2) in a string(S1)
    INTEGER*4 FUNCTION PPICLF_INDX1(S1,S2)
        ! 
        ! Vars:
        !
        CHARACTER(len=*) S1
        character(len=1) S2
        integer*4 n1, i
        !
        N1 = LEN(S1)
        PPICLF_INDX1=0
        ! IF (N1.LT.1) return
        !
        DO I=1,N1
            IF (S1(I:I).EQ.S2) THEN
                PPICLF_INDX1=I
                return
            END IF
        END DO
        !
        return
    END FUNCTION PPICLF_INDX1

    character*132 FUNCTION PPICLF_CHSTR(S1,indx1)
        ! 
        ! Vars:
        !
        ! TLJ modified, but not sure why I had to
        CHARACTER S1
        INTEGER indx1
        !
        PPICLF_CHSTR = S1(1:indx1)
        return
    END FUNCTION PPICLF_CHSTR

    SUBROUTINE ppiclf_byte_open_mpi(fnamei,mpi_fh,ifro,ierr)
        !
        ! Vars:
        !
        character fnamei*(*)
        logical ifro
        CHARACTER*1 BLNK
        DATA BLNK/' '/
        character*132 fname
        character*1   fname1(132)
        equivalence  (fname1,fname)
        integer*4 imode, ierr, mpi_fh
        !
        imode = MPI_MODE_WRONLY+MPI_MODE_CREATE
        if(ifro) then
            imode = MPI_MODE_RDONLY 
        endif

        call MPI_file_open(ppiclf_comm,fnamei,imode, MPI_INFO_NULL,mpi_fh,ierr)

        return
    END SUBROUTINE ppiclf_byte_open_mpi

    SUBROUTINE ppiclf_byte_read_mpi(buf,icount,mpi_fh,ierr)
        !
        ! Vars:
        !
        real*4 buf(1)          ! buffer
        integer*4 iout, icount, mpi_fh, ierr
        !

        iout = icount ! icount is in 4-byte words
        call MPI_file_read_all(mpi_fh,buf,iout,MPI_REAL, MPI_STATUS_IGNORE,ierr)

        return
    END SUBROUTINE ppiclf_byte_read_mpi

    SUBROUTINE ppiclf_byte_write_mpi(buf,icount,iorank,mpi_fh,ierr)
        !
        ! Vars:
        !
        real*4 buf(1)          ! buffer
        integer*4 icount, iorank, mpi_fh, ierr, iout
        !

        iout = icount ! icount is in 4-byte words
        if(iorank.ge.0 .and. ppiclf_nid.ne.iorank) iout = 0
        call MPI_file_write_all(mpi_fh,buf,iout,MPI_REAL, MPI_STATUS_IGNORE,ierr)

        return
    END SUBROUTINE ppiclf_byte_write_mpi

    SUBROUTINE ppiclf_byte_close_mpi(mpi_fh,ierr)
        !
        ! Vars:
        !
        integer*4 mpi_fh, ierr
        !

        call MPI_file_close(mpi_fh,ierr)

        return
    END SUBROUTINE ppiclf_byte_close_mpi

    SUBROUTINE ppiclf_byte_set_view(ioff_in,mpi_fh)
        !
        ! Vars:
        !
        integer*8 ioff_in
        integer*4 mpi_fh, ierr
        !
        call MPI_file_set_view(mpi_fh,ioff_in,MPI_BYTE,MPI_BYTE,'native',MPI_INFO_NULL,ierr)

        return
    END SUBROUTINE ppiclf_byte_set_view

    SUBROUTINE ppiclf_bcast(buf,len)
        !
        ! Vars:
        !
        !real*4 buf(1)
        type(*) buf(*)
        integer*4 len, ierr
        !

        call mpi_bcast (buf,len,mpi_byte,0,ppiclf_comm,ierr)

        return
    END SUBROUTINE ppiclf_bcast
end module ppiclf_op
