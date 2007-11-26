!----------------------------------------------------------------------
!
!  module: tfsf
!
!  field sources using total-field scattered-field.
!
!----------------------------------------------------------------------


module tfsf

  use constant
  use grid  
  use fdtd

  implicit none
  save

  ! --- Constants

  character(len=STRLNG) :: pfxtfsf = "tfsf"

  ! --- Variables

  integer :: k_inc
  integer :: j, ioffset, joffset, ilength, jlength
  real(8) :: om, t0, sig, latcon

  ! --- Fields

  real(8), allocatable :: Hy_inc(:,:), Ex_inc(:,:), Hx_inc(:,:), Ey_inc(:,:)

contains

  subroutine InitializeTFSF()

    implicit none

    integer :: i,j, offi, offj
    real(8) :: Hzinch, Eyinch

    integer :: ios
    character(len=255) :: file, str, tfsffn, label

    ! Read input parameters and file

    file = cat2(pfxtfsf,sfxin)
    
    open(UNITTMP,FILE=file,STATUS='unknown')
    do
       read(UNITTMP,IOSTAT=ios,FMT=*) str
       if(ios .ne. 0) exit
       if(str(1:5).eq. '#TFSF') then
          read(UNITTMP,*) label
          read(UNITTMP,*) tfsffn
          read(UNITTMP,*) ioffset, joffset, k_inc
          read(UNITTMP,*) ilength, jlength
          read(UNITTMP,*) latcon
          read(UNITTMP,*) om, t0, sig
          exit
       endif
    enddo
    close(UNITTMP)

    if ( k_inc .gt. KEND .or. k_inc .lt. KBEG) return   ! MPI CHECK: are we in or what?

    om=om*2.0*PI/latcon

    t0=t0*2*PI/om
    sig=sig*2*PI/om

!    write(6,*) "TFSF-Test: ", ioffset, ilength, joffset, jlength, k_inc, latcon, om, t0, sig
!    write(6,*) "TFSF-Test: ", tfsffn, label

    allocate(Ex_inc(ioffset:ioffset+ilength,joffset:joffset+jlength))
    allocate(Ey_inc(ioffset:ioffset+ilength,joffset:joffset+jlength))
    allocate(Hx_inc(ioffset:ioffset+ilength,joffset:joffset+jlength))
    allocate(Hy_inc(ioffset:ioffset+ilength,joffset:joffset+jlength))
 
    i = ioffset
    j = joffset
    
    Hy_inc(i,j) = 1.0
    Hx_inc(i,j) = 1.0
    Ex_inc(i,j) = 1.0
    Ey_inc(i,j) = 1.0

    open(UNITTMP,FILE=tfsffn,STATUS='unknown',IOSTAT=ios) 
    
    if ( ios .ne. 0 ) goto 10


! changed loop order of i and j
! Changed input/output order from Hy, Hx, Ex, Ey -> Ex, Ey, Hx, Hy

    do j=joffset, joffset+jlength-1
       do i=ioffset, ioffset+ilength-1
          read(UNITTMP,*,END=10,ERR=10) Ex_inc(i,j), Ey_inc(i,j), Hx_inc(i,j), Hy_inc(i,j) 
       end do
    end do
        
10  continue

    close(UNITTMP)

    offi = max(i,ioffset+1)
    offj = max(j,joffset+1)

    do j=offj, joffset+jlength-1
       do i=offi, ioffset+ilength-1
          Ex_inc(i,j) = Ex_inc(offi-1,offj)
          Ey_inc(i,j) = Ey_inc(offi-1,offj)
          Hx_inc(i,j) = Hx_inc(offi-1,offj)
          Hy_inc(i,j) = Hy_inc(offi-1,offj)
       end do
    end do

!DEBUG >>>

!    do j=joffset, joffset+jlength-1
!       do i=ioffset, ioffset+ilength-1
!          Ex_inc(i,j) = 1.0
!          Ey_inc(i,j) = 1.0
!          Hx_inc(i,j) = 1.0
!          Hy_inc(i,j) = 1.0 
!       end do
!    end do

!<<< DEBUG
    ! Test output 

    open(UNITTMP,FILE=cat2(label,'.dat'),STATUS='unknown')
    
    do j=0, joffset-1
       do i=0, IMAX
          write(UNITTMP,*) 0.0, 0.0, 0.0, 0.0, 1.0/epsinv(i,j,k_inc)
       end do
       write(UNITTMP,*)
    end do
                                                                                                   
    do j=joffset, joffset+jlength-1
       do i=0, ioffset-1
          write(UNITTMP,*) 0.0, 0.0, 0.0, 0.0, 1.0/epsinv(i,j,k_inc)
       end do
       do i=ioffset, ioffset+ilength-1
          write(UNITTMP,*) Ex_inc(i,j), Ey_inc(i,j), Hx_inc(i,j), Hy_inc(i,j), 1.0/epsinv(i,j,k_inc)
       end do
       do i=ioffset+ilength, IMAX
          write(UNITTMP,*) 0.0, 0.0, 0.0, 0.0, 1.0/epsinv(i,j,k_inc)
       end do
       write(UNITTMP,*)
    end do
                                                                                                   
    do j=joffset+jlength, JMAX
       do i=0, IMAX
          write(UNITTMP,*) 0.0, 0.0, 0.0, 0.0, 1.0/epsinv(i,j,k_inc)
       end do
       write(UNITTMP,*)
    end do

    close(UNITTMP)

  end subroutine InitializeTFSF


  subroutine NewTFSF_E()
 
   implicit none

    real(8) :: H_time
    integer :: i,j

    if ( k_inc .gt. KEND .or. k_inc .lt. KBEG) return   ! MPI CHECK: are we in or what?

    H_time=sin(om*(GT+0.5*DT-t0))*exp(-0.5*(GT+0.5*DT-t0)**2/sig**2)

    do j=joffset, joffset+jlength-1
       do i=ioffset, ioffset+ilength-1
          Ex(i,j,k_inc)=Ex(i,j,k_inc)+DT/SZ*H_time*(Hy_inc(i,j)+Hy_inc(i+1,j))*epsinv(i,j,k_inc)
          Ey(i,j,k_inc)=Ey(i,j,k_inc)+DT/SZ*H_time*(Hx_inc(i,j)+Hx_inc(i,j+1))*epsinv(i,j,k_inc)
       end do
    end do
    
  end subroutine NewTFSF_E

  subroutine NewTFSF_H()

    implicit none

    real(8) :: E_time
    integer :: i,j
    
    if ( k_inc .gt. KEND .or. k_inc .lt. KBEG) return   ! MPI CHECK: are we in or what?

    E_time=sin(om*(GT-t0))*exp(-0.5*(GT-t0)**2/sig**2)

    do j=joffset, joffset+jlength-1
       do i=ioffset, ioffset+ilength-1
          Hy(i,j,k_inc-1)=Hy(i,j,k_inc-1)+DT/SZ*E_time*(Ex_inc(i,j)+Ex_inc(i+1,j))
          Hx(i,j,k_inc-1)=Hx(i,j,k_inc-1)+DT/SZ*E_time*(Ey_inc(i,j)+Ey_inc(i,j+1))
       end do
    end do
    
  end subroutine NewTFSF_H


  subroutine FinalizeTFSF

    implicit none

    if ( k_inc .gt. KEND .or. k_inc .lt. KBEG) return   ! MPI CHECK: are we in or what?

    deallocate(Ex_inc)
    deallocate(Hy_inc)
    deallocate(Ey_inc)
    deallocate(Hx_inc)

  end subroutine FinalizeTFSF


end module tfsf
