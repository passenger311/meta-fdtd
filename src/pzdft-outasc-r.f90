!----------------------------------------------------------------------
!
!  module: pzdft-output-r
!
!  dft output module.
!
!  subs:
!
!  
!
!----------------------------------------------------------------------


module pzdft_output

  use constant
  use strings
  use mpiworld
  use output
  use grid  
  use fdtd

  implicit none
  save

  ! --- Constants

  character(len=255), parameter :: pfxoutput = 'output'
  integer, parameter :: PARTSMAXGPL = 100


  ! --- Variables

  type(T_OutBas) :: gpl(PARTSMAXGPL)
  integer  partsgpl
  real(8), allocatable :: DataGpl(:)
  integer DataIndxGpl(PARTSMAXGPL+1)

  integer :: xoff, yoff 
  integer :: kt
  integer :: ndf  
  real(8) :: fl      
  real(8) :: fu      
  real(8) :: latcon  

  character(STRLNG) :: pztfn

  ! --- Fields

  real(8), allocatable :: Eyt_r(:,:,:)
  real(8), allocatable :: Eyt_i(:,:,:)
  real(8), allocatable :: Hxt_r(:,:,:)
  real(8), allocatable :: Hxt_i(:,:,:)

  real(8), allocatable :: Ext_r(:,:,:)
  real(8), allocatable :: Ext_i(:,:,:)
  real(8), allocatable :: Hyt_r(:,:,:)
  real(8), allocatable :: Hyt_i(:,:,:)


contains


  subroutine InitializePzDftOutput

    implicit none

    integer ::  ios,n,i,err
    character(len=STRLNG) :: file, str

    call ReadConfig
    call Initialize

  contains
  
    subroutine ReadConfig

      implicit none
      
      integer ios
      character(len=STRLNG) :: file, str
      
      ndf=-1
      
      file = cat2(pfxoutput,sfxin)
      
      ! Read Data
      open(UNITTMP,FILE=file,STATUS='unknown')
      do
         read(UNITTMP,IOSTAT=ios,FMT=*) str
         if(ios .ne. 0) exit
         if(str(1:4).eq. '#PzDft') then
            read(UNITTMP,*) pztfn
            read(UNITTMP,*) fl, fu
            read(UNITTMP,*) ndf
            read(UNITTMP,*) kt, xoff, yoff
            read(UNITTMP,*) latcon
            exit
         endif
      enddo
      close(UNITTMP)
      
    end subroutine ReadConfig

    subroutine Initialize

      implicit none


    end subroutine Initialize


  end subroutine InitializePzDftOutput


  subroutine FinalizePzDftOutput
    implicit none


  end subroutine FinalizePzDftOutput


  subroutine WritePzDftOutput(ncyc)


    implicit none
    
    integer :: ncyc
    integer :: n
    logical :: ret

    
    do n=1, PARTSGPL
       
       call WriteOutput(n, ncyc, ret)
       
       if ( ret ) then

          ! output
          select case (gpl(n)%Mode(1:2))
          case('Px')
             call LoadPx(n)
             call WriteData(n,UNITTMP) 
          case('Py')
             call LoadPy(n)
             call WriteData(n,UNITTMP) 
          case('Pz')
             call LoadPz(n)
             call WriteData(n,UNITTMP) 
          case('En')
             call LoadEn(n)
             call WriteData(n,UNITTMP) 
          case('EH')
             call WriteEH(UNITTMP)
          case('Ex')
             call WriteComp(Ex,UNITTMP)
          case('Ey')
             call WriteComp(Ey,UNITTMP)
          case('Ez')
             call WriteComp(Ez,UNITTMP)
          case('Hx')
             call WriteComp(Hx,UNITTMP)
          case('Hy')
             call WriteComp(Hy,UNITTMP)
          case('Hz')
             call WriteComp(Hz,UNITTMP)
          case('Di')         
             call WriteDi(UNITTMP)
          end select
          ! skip all the others
          call CloseOutput
       endif
    enddo
       
  contains

  end subroutine WritePzDftOutput


  subroutine InitPzDFT()

    use constant
    use grid
    implicit none

    integer ios
    character(len=STRLNG) :: file, str

    ndf=-1

    file = cat2(pfxoutput,sfxin)

    ! Read Data
    open(UNITTMP,FILE=file,STATUS='unknown')
    do
       read(UNITTMP,IOSTAT=ios,FMT=*) str
       if(ios .ne. 0) exit
       if(str(1:4).eq. '#PZT') then
          read(UNITTMP,*) pztfn
          read(UNITTMP,*) fl, fu
          read(UNITTMP,*) ndf
          read(UNITTMP,*) kt, xoff, yoff
          read(UNITTMP,*) latcon
          exit
       endif
    enddo
    close(UNITTMP)

    allocate(Ext_r(xoff:IMAX-xoff,yoff:JMAX-yoff,0:ndf))
    allocate(Ext_i(xoff:IMAX-xoff,yoff:JMAX-yoff,0:ndf))
    allocate(Hyt_r(xoff:IMAX-xoff,yoff:JMAX-yoff,0:ndf))
    allocate(Hyt_i(xoff:IMAX-xoff,yoff:JMAX-yoff,0:ndf))

    allocate(Eyt_r(xoff:IMAX-xoff,yoff:JMAX-yoff,0:ndf))
    allocate(Eyt_i(xoff:IMAX-xoff,yoff:JMAX-yoff,0:ndf))
    allocate(Hxt_r(xoff:IMAX-xoff,yoff:JMAX-yoff,0:ndf))
    allocate(Hxt_i(xoff:IMAX-xoff,yoff:JMAX-yoff,0:ndf))

    Ext_r=0.0
    Ext_i=0.0
    Hyt_r=0.0
    Hyt_i=0.0

    Eyt_r=0.0
    Eyt_i=0.0
    Hxt_r=0.0
    Hxt_i=0.0

  end subroutine InitPzDFT

  subroutine PzDFT(ncyc)

    use constant
    use grid
    implicit none
    
    integer :: ncyc, i, j, k, m
    real(8) :: cosn, sinn, hilf, fac

    if (ndf .ne. -1) then
       do m=0, ndf

          cosn=cos(((1.0*m)/(1.0*ndf)*(fu-fl)+fl)*6.2831853/latcon*DT*ncyc)*DT 
          sinn=sin(((1.0*m)/(1.0*ndf)*(fu-fl)+fl)*6.2831853/latcon*DT*ncyc)*DT 
       
          do j=yoff, JMAX-yoff
             do i=xoff, IMAX-xoff

                hilf=0.5*(Ex(i,j,kt)+Ex(i-1,j,kt))

                Ext_r(i,j,m)=Ext_r(i,j,m)+ hilf*cosn
                Ext_i(i,j,m)=Ext_i(i,j,m)+ hilf*sinn

                hilf=0.25*(Hy(i,j,kt)+Hy(i-1,j,kt)+Hz(i,j-1,kt)+Hz(i-1,j-1,kt))
             
                Hyt_r(i,j,m)=Hyt_r(i,j,m)+ hilf*cosn
                Hyt_i(i,j,m)=Hyt_i(i,j,m)+ hilf*sinn

                hilf=0.5*(Ey(i,j,kt)+Ey(i,j-1,kt))
                
                Eyt_r(i,j,m)=Eyt_r(i,j,m)+ hilf*cosn
                Eyt_i(i,j,m)=Eyt_i(i,j,m)+ hilf*sinn
                
                hilf=0.25*(Hx(i,j,kt)+Hx(i,j-1,kt)+Hx(i,j,kt-1)+Hx(i,j-1,kt-1))
                
                Hxt_r(i,j,m)=Hxt_r(i,j,m)+ hilf*cosn
                Hxt_i(i,j,m)=Hxt_i(i,j,m)+ hilf*sinn

             end do
          end do
       end do
    end if

  end subroutine PzDFT

  subroutine OutPzDFT()

    use constant
    use grid
    implicit none
    
    integer :: i,j,k,m
    real(8) :: pzt(0:ndf)

    
    if (ndf .ne. -1) then

       pzt=0.0

       do m=0, ndf
          do j=yoff, JMAX-yoff
             do i=xoff, IMAX-xoff

                pzt(m)=pzt(m)+Ext_r(i,j,m)*Hyt_r(i,j,m)+ & 
                     Ext_i(i,j,m)*Hyt_i(i,j,m)-Eyt_r(i,j,m)*Hxt_r(i,j,m)- &
                     Eyt_i(i,j,m)*Hxt_i(i,j,m)

             end do
          end do
       end do
    
       open(UNITTMP,FILE=pztfn)

       do m=0, ndf

          write(UNITTMP,*) (1.0*m)/(1.0*ndf)*(fu-fl)+fl, pzt(m)

       end do

       close(UNITTMP)

    end if

  end subroutine OutPzDFT

end module output
