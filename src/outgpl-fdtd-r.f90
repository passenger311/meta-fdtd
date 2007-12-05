!-*- F90 -*------------------------------------------------------------
!
!  module: outgpl-fdtd / max3d
!
!  this module handles GPL output of data related to the fdtd module.
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
  use outlist
  use reglist
  use buflist
  use mpiworld
  use grid 
  use fdtd

  implicit none
  save

contains

!----------------------------------------------------------------------

  subroutine InitializeOutgplFdtd
    
  end subroutine InitializeOutgplFdtd

!----------------------------------------------------------------------

  subroutine FinalizeOutgplFdtd

  end subroutine FinalizeOutgplFdtd

!----------------------------------------------------------------------


  subroutine WriteDataOutgplFdtdObj(out, ncyc)

    type (T_OUT) :: out
    integer :: ncyc

    logical :: ret

    select case (out%fn)
    case('Px')
       call LoadPx(out)
       call WriteData(out) 
    case('Py')
       call LoadPy(out)
       call WriteData(out) 
    case('Pz')
       call LoadPz(out)
       call WriteData(out) 
    case('En')
       call LoadEn(out)
       call WriteData(out) 
    case('EH')
       call WriteEH(out)
    case('Ex')
       call WriteComp(out,Ex)
    case('Ey')
       call WriteComp(out,Ey)
    case('Ez')
       call WriteComp(out,Ez)
    case('Hx')
       call WriteComp(out,Hx)
    case('Hy')
       call WriteComp(out,Hy)
    case('Hz')
       call WriteComp(out,Hz)
    case('Di')         
       call WriteDi(out)
    end select
    
  contains

    ! **************************************************************** !   

    subroutine WriteEH(out)

      implicit none
      type (T_OUT) :: out
      integer i,j,k

      do k=out%ks, out%ke, out%dk  
         do j=out%js, out%je, out%dj      
            do i=out%is, out%ie, out%di 
               write(UNITTMP,'(6E15.6E3))') Ex(i,j,k),  Ey(i,j,k), &
                    Ez(i,j,k), Hx(i,j,k), Hy(i,j,k), Hz(i,j,k)
            enddo
            if(out%is .ne. out%ie) write(UNITTMP,*)
         enddo
         if(out%js .ne. out%je) write(UNITTMP,*)
      enddo
    end subroutine WriteEH


    ! **************************************************************** !

    subroutine WriteComp(out,Comp)

      implicit none

      type (T_OUT) :: out
      real(8), dimension(IMIN:KMAX,JMIN:KMAX,KMIN:KMAX) :: Comp
      integer i,j,k

      do k=out%ks, out%ke, out%dk  
         do j=out%js, out%je, out%dj      
            do i=out%is, out%ie, out%di 
               write(UNITTMP,'(E15.6))') Comp(i,j,k)
            enddo
            if(out%is .ne. out%ie) write(UNITTMP,*) 
         enddo
         if(out%js .ne. out%je) write(UNITTMP,*)
      enddo

    end subroutine WriteComp


    ! **************************************************************** !

    subroutine WriteDi(out)

      implicit none

      type (T_OUT) :: out
      integer i,j,k

      write(6,*) out%is, out%ie
      write(6,*) out%js, out%je
      write(6,*) out%ks, out%ke 
      do k=out%ks, out%ke, out%dk
         do j=out%js, out%je, out%dj      
            do i=out%is, out%ie, out%di 
               write(UNITTMP,'(E15.6E3))') 1.0/epsinv(i,j,k)
            enddo
            if(out%is .ne. out%ie) write(UNITTMP,*)
         enddo
         if(out%js .ne. out%je) write(UNITTMP,*)
      enddo
      
    end subroutine WriteDi 


    ! **************************************************************** !

    subroutine WriteData(out)

      implicit none
  
      type (T_OUT) :: out
      integer :: i,j,k,n,m
      real(8) :: sum
  
      sum = 0.0
      m = DataIndxOut(out%idx)
      if( out%Mode(4:4) .eq. 'S' ) then      ! integration 	
         do k=out%ks, out%ke, out%dk  
            do j=out%js, out%je, out%dj      
               do i=out%is, out%ie, out%di 
                  sum = sum+DataOut(m)
                  m = m+1
               enddo
            enddo
         enddo
         write(UNITTMP,*) sum
      else                                     ! direct outgpl
         do k=out%ks, out%ke, out%dk  
            do j=out%js, out%je, out%dj      
               do i=out%is, out%ie, out%di 
                  write(UNITTMP,'(E15.6E3))') DataOut(m)
                  m = m+1
               enddo
               if(out%is .ne. out%ie) write(UNITTMP,*)
            enddo
            if(out%js .ne. out%je) write(UNITTMP,*)
         enddo
      endif
    end subroutine WriteData

  end subroutine WritOutgplFdtd


  ! ***************************************************************** !

  subroutine PrepareOutgplFdtd(ncyc)

    implicit none

    integer ncyc, n, i
    type (T_OUT) :: out

    ! Loop over all outgpl units
    do n=1, numoutlist

       out = outlist(n)

       ! outgpl ?
       if(mod(ncyc, out%dn) .eq. 0 .and. &
            ncyc .ge. out%ns       .and. &
            ncyc .le. out%ne ) then

          ! First part of Px,y,z calculation
          select case (out%Mode(1:2))
          case('Px')
             DataOut(DataIndxOut(out%idx):DataIndxOut(out%idx+1))=0.0
             call LoadPx(out)
          case('Py')
             DataOut(DataIndxOut(out%idx):DataIndxOut(out%idx+1))=0.0
             call LoadPy(out)
          case('Pz')
             DataOut(DataIndxOut(out%idx):DataIndxOut(out%idx+1))=0.0
             call LoadPz(out)
          end select
          
       endif
    enddo
  end subroutine PrepareOutgplFdtd


  subroutine LoadPx(out)
    ! Calculates x-component of the Poynting vector P
    ! Stores result in DataOut
    ! Localization: Px[i,j,k] = Px(i+1/2,j,k)
    implicit none
    type (T_OUT) :: out
    integer :: i,j,k,n,m
    real(8) :: val
    ! Code
    m = DataIndxOut(out%idx)
    do k=out%ks, out%ke, out%dk 
       do j=out%js, out%je, out%dj 
          do i=out%is, out%ie, out%di
             val = 0.125*((Ey(i,j,k)+Ey(i+1,j,k))*Hz(i,j,k) &
                  + (Ey(i,j-1,k)+Ey(i+1,j-1,k))* Hz(i,j-1,k) &
                  - (Ez(i,j,k)+Ez(i+1,j,k))*Hy(i,j,k) &
                  - (Ez(i,j,k-1)+Ez(i+1,j,k-1))* Hy(i,j,k-1)) 
             DataOut(m)=DataOut(m)+val/(4.0*PI) 
             m = m+1
          enddo
       enddo
    enddo
  end subroutine LoadPx


  subroutine LoadPy(out)
    ! Calculates y-component of the Poynting vector P
    ! Stores result in DataOut
    ! Localization: Py[i,j,k] = Py(i,j+1/2,k)
    implicit none
    type (T_OUT) :: out
    integer :: i,j,k,n,m
    real(8) :: val
    ! Code
    m = DataIndxOut(out%idx)
    do k=out%ks, out%ke, out%dk
       do j=out%js, out%je, out%dj
          do i=out%is, out%ie, out%di
               val = 0.125*( (Ez(i,j,k)+Ez(i,j+1,k))*Hx(i,j,k) &
                    + (Ez(i,j,k-1)+Ez(i,j+1,k-1))*Hx(i,j,k-1) &
                    - (Ex(i,j,k)+Ex(i,j+1,k))*Hz(i,j,k) &
                    - (Ex(i-1,j,k)+Ex(i-1,j+1,k))*Hz(i-1,j,k))
             DataOut(m)=DataOut(m)+val/(4.0*PI)
             m = m+1
          enddo
      enddo
    enddo
  end subroutine LoadPy


  subroutine LoadPz(out)
    ! Calculates z-component of the Poynting vector P
    ! Stores result in DataOut
    ! Localization: Pz[i,j,k] = Pz(i,j,k+1/2)
    implicit none
    type (T_OUT) :: out
    integer :: i,j,k,n, m
    real(8) :: val
    ! Code
    m = DataIndxOut(out%idx)
    do k=out%ks, out%ke, out%dk
       do j=out%js, out%je, out%dj
          do i=out%is, out%ie, out%di
             val = 0.125*( (Ex(i,j,k)+Ex(i,j,k+1))*Hy(i,j,k) &
                  +(Ex(i-1,j,k)+Ex(i-1,j,k+1))*Hy(i-1,j,k) &
                  - (Ey(i,j,k)+Ey(i,j,k+1))*Hx(i,j,k) &
                  - (Ey(i,j-1,k)+Ey(i,j-1,k+1))* Hx(i,j-1,k))
             DataOut(m)=DataOut(m)+val/(4.0*PI)
             m = m+1
          enddo
      enddo
    enddo
  end subroutine LoadPz


  subroutine LoadEn(out)
    ! Calculates energy density En
    ! Stores result in DataOut
    ! Localization: En[i,j,k] = En(i,j,k)
    implicit none
    type (T_OUT) :: out
    integer :: i,j,k,n,m
    real(8) :: EEx,EEy,EEz,EHx,EHy,EHz,eps
    ! Start Code
    m = DataIndxOut(out%idx)
    do k=out%ks, out%ke, out%dk
       do j=out%js, out%je, out%dj
          do i=out%is, out%ie, out%di
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
             DataOut(m)=(EEx+EEy+EEz+EHx+EHy+EHz)/(4.0*PI)
             m = m+1
          enddo
       enddo
    enddo
  end subroutine LoadEn

end module outgpl_fdtd

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
