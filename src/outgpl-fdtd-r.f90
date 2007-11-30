!----------------------------------------------------------------------
!
!  module: fdtd-outgpl-r
!
!  ascii outgpl module.
!
!  subs:
!
!  
!
!----------------------------------------------------------------------

! ---------------------------------------------------------------------
!     Supported Data: Ex, Ey, Ez, Hx, Hy, Hz, Di, Px, Py, Pz, En
! ---------------------------------------------------------------------
!     Contained Subroutines:
!
!     InitOutgplParameters
!     WriteOutgplHeader
!     InitOutgpl             used in max3d.f90
!     Outgpl(ncyc)           used in max3d.f90
!            WriteEH
!            WriteComp
!            WriteDi
!            WriteData
!     DataPrepOutgpl(ncyc)   used in max3d.f90 (between StepH and StepE)
!     LoadPx
!     LoadPy
!     LoadPz

! ---------------------------------------------------------------------

!  outgpl mode 'abcd':
!  ab = component(s) 'Ex', 'Ey', 'Ez', 'Hx', 'Hy','Hz','EH','Di',
!                     'En', 'Px', 'Py', oder 'Pz'
!  c = 'E':  one file (for all time steps)
!  c = 'M':  multiple files
!  d = 'R':  spatially resolved outgpl
!  d = 'S':  spatially integrated outgpl

!  Spatial and temporal localization of the components in (i,j,k,ncyc):
!  Ex-Ez: (i+1/2,j,k,GT) - (i,j,k+1/2,GT)  
!  Hx-Hz: (i,j+1/2,k+1/2,GT-0.5*DT) - (i+1/2,j+1/2,k,GT-0.5*DT)
!  EH = (Ex,Ey,Ez,Hx,Hy,Hz): as above
!  Di = diectric constant: (i,j,k)
!  Energy En:  electric part EnE: (i,j,k,GT)
!              magnetic part EnM: (i,j,k,GT-0.5*DT)
!              En = EnE + EnM
!  Px-Pz:  (i+1/2,j,k,GT-0.5*DT) - (i,j,k+1/2,GT-0.5*DT)  
! ---------------------------------------------------------------------


module outgpl_fdtd

  use constant
  use strings
  use mpiworld
  use grid  
  use outgpl
  use fdtd

  implicit none
  save


  ! --- Variables

  real(8), allocatable :: DataGpl(:)
  integer DataIndxGpl(MAXOBJGPL+1)

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


  subroutine InitializeOutgplFdtd

    implicit none

    call Initialize    

  contains
  

    subroutine Initialize

      implicit none

      integer :: i, n, err

      i = 1
      DataIndxGpl(1)=1
      do n=1, numobjgpl
         if( objgpl(n)%Mode(1:2) .eq. 'En' .or. &
              objgpl(n)%Mode(1:2) .eq. 'Px' .or. &
              objgpl(n)%Mode(1:2) .eq. 'Py' .or. &
              objgpl(n)%Mode(1:2) .eq. 'Pz' ) then
            i = i + objgpl(n)%NumNodes               
         endif
         DataIndxGpl(n+1)=i
      enddo
      allocate(DataGpl(1:i),STAT=err)
      if(err .ne. 0) then
         write(6,*) 'Error 1 in InitOutGpl()'
         stop
         
      endif

    end subroutine Initialize


  end subroutine InitializeOutgplFdtd


  subroutine FinalizeOutgplFdtd
    implicit none

    deallocate(DataGpl)

  end subroutine FinalizeOutgplFdtd


  subroutine WritOutgplFdtd(ncyc)

    implicit none
    
    integer :: ncyc
    type (T_OUTBAS) :: gpl
    integer :: n
    logical :: ret

    
    do n=1, numobjgpl
       
       gpl = objgpl(n)

       call OpenObjOutgpl(gpl, ncyc, ret)
       
       if ( ret ) then

          ! outgpl
          select case (gpl%Mode(1:2))
          case('Px')
             call LoadPx(gpl)
             call WriteData(gpl) 
          case('Py')
             call LoadPy(gpl)
             call WriteData(gpl) 
          case('Pz')
             call LoadPz(gpl)
             call WriteData(gpl) 
          case('En')
             call LoadEn(gpl)
             call WriteData(gpl) 
          case('EH')
             call WriteEH(gpl)
          case('Ex')
             call WriteComp(gpl,Ex)
          case('Ey')
             call WriteComp(gpl,Ey)
          case('Ez')
             call WriteComp(gpl,Ez)
          case('Hx')
             call WriteComp(gpl,Hx)
          case('Hy')
             call WriteComp(gpl,Hy)
          case('Hz')
             call WriteComp(gpl,Hz)
          case('Di')         
             call WriteDi(gpl)
          end select
          ! skip all the others

          call CloseObjOutgpl

       endif

    enddo
       
  contains

    ! **************************************************************** !   

    subroutine WriteEH(gpl)

      implicit none
      type (T_OUTBAS) :: gpl
      integer i,j,k

      do k=gpl%ks, gpl%ke, gpl%dk  
         do j=gpl%js, gpl%je, gpl%dj      
            do i=gpl%is, gpl%ie, gpl%di 
               write(UNITTMP,'(6E15.6E3))') Ex(i,j,k),  Ey(i,j,k), &
                    Ez(i,j,k), Hx(i,j,k), Hy(i,j,k), Hz(i,j,k)
            enddo
            if(gpl%is .ne. gpl%ie) write(UNITTMP,*)
         enddo
         if(gpl%js .ne. gpl%je) write(UNITTMP,*)
      enddo
    end subroutine WriteEH


    ! **************************************************************** !

    subroutine WriteComp(gpl,Comp)

      implicit none

      type (T_OUTBAS) :: gpl
      real(8), dimension(IMIN:KMAX,JMIN:KMAX,KMIN:KMAX) :: Comp
      integer i,j,k

      do k=gpl%ks, gpl%ke, gpl%dk  
         do j=gpl%js, gpl%je, gpl%dj      
            do i=gpl%is, gpl%ie, gpl%di 
               write(UNITTMP,'(E15.6))') Comp(i,j,k)
            enddo
            if(gpl%is .ne. gpl%ie) write(UNITTMP,*) 
         enddo
         if(gpl%js .ne. gpl%je) write(UNITTMP,*)
      enddo

    end subroutine WriteComp


    ! **************************************************************** !

    subroutine WriteDi(gpl)

      implicit none

      type (T_OUTBAS) :: gpl
      integer i,j,k

      write(6,*) gpl%is, gpl%ie
      write(6,*) gpl%js, gpl%je
      write(6,*) gpl%ks, gpl%ke 
      do k=gpl%ks, gpl%ke, gpl%dk
         do j=gpl%js, gpl%je, gpl%dj      
            do i=gpl%is, gpl%ie, gpl%di 
               write(UNITTMP,'(E15.6E3))') 1.0/epsinv(i,j,k)
            enddo
            if(gpl%is .ne. gpl%ie) write(UNITTMP,*)
         enddo
         if(gpl%js .ne. gpl%je) write(UNITTMP,*)
      enddo
      
    end subroutine WriteDi 


    ! **************************************************************** !

    subroutine WriteData(gpl)

      implicit none
  
      type (T_OUTBAS) :: gpl
      integer :: i,j,k,n,m
      real(8) :: sum
  
      sum = 0.0
      m = DataIndxGpl(gpl%idx)
      if( gpl%Mode(4:4) .eq. 'S' ) then      ! integration 	
         do k=gpl%ks, gpl%ke, gpl%dk  
            do j=gpl%js, gpl%je, gpl%dj      
               do i=gpl%is, gpl%ie, gpl%di 
                  sum = sum+DataGpl(m)
                  m = m+1
               enddo
            enddo
         enddo
         write(UNITTMP,*) sum
      else                                     ! direct outgpl
         do k=gpl%ks, gpl%ke, gpl%dk  
            do j=gpl%js, gpl%je, gpl%dj      
               do i=gpl%is, gpl%ie, gpl%di 
                  write(UNITTMP,'(E15.6E3))') DataGpl(m)
                  m = m+1
               enddo
               if(gpl%is .ne. gpl%ie) write(UNITTMP,*)
            enddo
            if(gpl%js .ne. gpl%je) write(UNITTMP,*)
         enddo
      endif
    end subroutine WriteData

  end subroutine WritOutgplFdtd


  ! ***************************************************************** !

  subroutine PrepareOutgplFdtd(ncyc)

    implicit none

    integer ncyc, n, i
    type (T_OUTBAS) :: gpl

    ! Loop over all outgpl units
    do n=1, numobjgpl

       gpl = objgpl(n)

       ! outgpl ?
       if(mod(ncyc, gpl%dn) .eq. 0 .and. &
            ncyc .ge. gpl%ns       .and. &
            ncyc .le. gpl%ne ) then

          ! First part of Px,y,z calculation
          select case (gpl%Mode(1:2))
          case('Px')
             DataGpl(DataIndxGpl(gpl%idx):DataIndxGpl(gpl%idx+1))=0.0
             call LoadPx(gpl)
          case('Py')
             DataGpl(DataIndxGpl(gpl%idx):DataIndxGpl(gpl%idx+1))=0.0
             call LoadPy(gpl)
          case('Pz')
             DataGpl(DataIndxGpl(gpl%idx):DataIndxGpl(gpl%idx+1))=0.0
             call LoadPz(gpl)
          end select
          
       endif
    enddo
  end subroutine PrepareOutgplFdtd


  subroutine LoadPx(gpl)
    ! Calculates x-component of the Poynting vector P
    ! Stores result in DataGpl
    ! Localization: Px[i,j,k] = Px(i+1/2,j,k)
    implicit none
    type (T_OUTBAS) :: gpl
    integer :: i,j,k,n,m
    real(8) :: val
    ! Code
    m = DataIndxGpl(gpl%idx)
    do k=gpl%ks, gpl%ke, gpl%dk 
       do j=gpl%js, gpl%je, gpl%dj 
          do i=gpl%is, gpl%ie, gpl%di
             val = 0.125*((Ey(i,j,k)+Ey(i+1,j,k))*Hz(i,j,k) &
                  + (Ey(i,j-1,k)+Ey(i+1,j-1,k))* Hz(i,j-1,k) &
                  - (Ez(i,j,k)+Ez(i+1,j,k))*Hy(i,j,k) &
                  - (Ez(i,j,k-1)+Ez(i+1,j,k-1))* Hy(i,j,k-1)) 
             DataGpl(m)=DataGpl(m)+val/(4.0*PI) 
             m = m+1
          enddo
       enddo
    enddo
  end subroutine LoadPx


  subroutine LoadPy(gpl)
    ! Calculates y-component of the Poynting vector P
    ! Stores result in DataGpl
    ! Localization: Py[i,j,k] = Py(i,j+1/2,k)
    implicit none
    type (T_OUTBAS) :: gpl
    integer :: i,j,k,n,m
    real(8) :: val
    ! Code
    m = DataIndxGpl(gpl%idx)
    do k=gpl%ks, gpl%ke, gpl%dk
       do j=gpl%js, gpl%je, gpl%dj
          do i=gpl%is, gpl%ie, gpl%di
               val = 0.125*( (Ez(i,j,k)+Ez(i,j+1,k))*Hx(i,j,k) &
                    + (Ez(i,j,k-1)+Ez(i,j+1,k-1))*Hx(i,j,k-1) &
                    - (Ex(i,j,k)+Ex(i,j+1,k))*Hz(i,j,k) &
                    - (Ex(i-1,j,k)+Ex(i-1,j+1,k))*Hz(i-1,j,k))
             DataGpl(m)=DataGpl(m)+val/(4.0*PI)
             m = m+1
          enddo
      enddo
    enddo
  end subroutine LoadPy


  subroutine LoadPz(gpl)
    ! Calculates z-component of the Poynting vector P
    ! Stores result in DataGpl
    ! Localization: Pz[i,j,k] = Pz(i,j,k+1/2)
    implicit none
    type (T_OUTBAS) :: gpl
    integer :: i,j,k,n, m
    real(8) :: val
    ! Code
    m = DataIndxGpl(gpl%idx)
    do k=gpl%ks, gpl%ke, gpl%dk
       do j=gpl%js, gpl%je, gpl%dj
          do i=gpl%is, gpl%ie, gpl%di
             val = 0.125*( (Ex(i,j,k)+Ex(i,j,k+1))*Hy(i,j,k) &
                  +(Ex(i-1,j,k)+Ex(i-1,j,k+1))*Hy(i-1,j,k) &
                  - (Ey(i,j,k)+Ey(i,j,k+1))*Hx(i,j,k) &
                  - (Ey(i,j-1,k)+Ey(i,j-1,k+1))* Hx(i,j-1,k))
             DataGpl(m)=DataGpl(m)+val/(4.0*PI)
             m = m+1
          enddo
      enddo
    enddo
  end subroutine LoadPz


  subroutine LoadEn(gpl)
    ! Calculates energy density En
    ! Stores result in DataGpl
    ! Localization: En[i,j,k] = En(i,j,k)
    implicit none
    type (T_OUTBAS) :: gpl
    integer :: i,j,k,n,m
    real(8) :: EEx,EEy,EEz,EHx,EHy,EHz,eps
    ! Start Code
    m = DataIndxGpl(gpl%idx)
    do k=gpl%ks, gpl%ke, gpl%dk
       do j=gpl%js, gpl%je, gpl%dj
          do i=gpl%is, gpl%ie, gpl%di
             eps = 1.0/epsinv(i,j,k)
             EEx = 0.5*eps*(Ex(i,j,k)**2+Ex(i-1,j,k)**2)
             EEy = 0.5*eps*(Ey(i,j,k)**2+Ey(i,j-1,k)**2)
             EEz = 0.5*eps*(Ez(i,j,k)**2+Ez(i,j,k-1)**2)
             EHx = 0.25*(Hx(i,j,k)**2+Hx(i,j,k-1)**2 + &
                         Hx(i,j-1,k)**2+Hx(i,j-1,k-1)**2)
             EHy = 0.25*(Hy(i,j,k)**2+Hy(i-1,j,k)**2 + &
                         Hy(i,j,k-1)**2+Hy(i-1,j,k-1)**2)
             EHz = 0.25*(Hz(i,j,k)**2+Hz(i-1,j,k)**2 + &
                         Hz(i,j-1,k)**2+Hz(i-1,j-1,k)**2)
             DataGpl(m)=(EEx+EEy+EEz+EHx+EHy+EHz)/(4.0*PI)
             m = m+1
          enddo
       enddo
    enddo
  end subroutine LoadEn

end module outgpl_fdtd
