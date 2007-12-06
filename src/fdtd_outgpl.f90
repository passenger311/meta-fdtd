!-*- F90 -*------------------------------------------------------------
!
!  module: fdtd_outgpl / meta3
!
!  this module handles GPL output of data related to the fdtd module.
!
!  subs:
!
!  InitializeFdtdOutgplObj
!  FinalizeFdtdOutgplObj
!  LoadDataFdtdOutgplObj
!  WriteDataFdtdOutgplObj
!
!----------------------------------------------------------------------


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

  subroutine InitializeFdtdOutgplObj(out)

    type (T_OUT) :: out
    type (T_BUF) :: buf

    if ( out%fn .eq. 'Px' .or. out%fn .eq. 'Py' .or. out%fn .eq. 'Pz' .or. &
         out%fn .eq. 'En' ) then

       ! allocate additional buffer for Load functions here

       buf = CreateBufObj(objreg(out%regidx))
       out%bufidx = buf%idx

    endif

  end subroutine InitializeFdtdOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeFdtdOutgplObj

    type (T_OUT) :: out

    if ( out%fn .eq. 'Px' .or. out%fn .eq. 'Py' .or. out%fn .eq. 'Pz' .or. &
         out%fn .eq. 'En' ) then

       ! destroy additional buffer for Load functions here

       call DestroyBufObj(objbuf(out%bufidx))
       
    endif

  end subroutine FinalizeFdtdOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataFdtdOutgplObj(out, ncyc)

    type (T_OUT) :: out
    integer :: ncyc

    select case (out%fn)
    case('Px')
       call LoadBufPx(out)
       call WriteBufData(out) 
    case('Py')
       call LoadPy(out)
       call WriteBufData(out) 
    case('Pz')
       call LoadBufPz(out)
       call WriteBufData(out) 
    case('En')
       call LoadBufEn(out)
       call WriteBufData(out) 
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

    subroutine WriteEH(out)

      type (T_OUT) :: out

      M4_REGLOOP_DECL(reg,p,i,j,k,w)  
      reg = regobj(out%regidx)

      M4_REGLOOP_WRITE(reg,p,i,j,k,w, 

         UNITTMP,6E15.6E3, `Ex(i,j,k),Ey(i,j,k),Ez(i,j,k),Hx(i,j,k),Hy(i,j,k),Hz(i,j,k)',
         
         `')

    end subroutine WriteEH


    ! **************************************************************** !

    subroutine WriteComp(out,Comp)

      type (T_OUT) :: out
      M4_FTYPE, dimension(IMIN:KMAX,JMIN:KMAX,KMIN:KMAX) :: Comp
 
      M4_REGLOOP_DECL(reg,p,i,j,k,w)  
      reg = regobj(out%regidx)

      M4_REGLOOP_WRITE(reg,p,i,j,k,w, 

         UNITTMP,E15.6E3, `Comp(i,j,k)',
         
         `')

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

  end subroutine WritFdtdOutgpl


  ! ***************************************************************** !

  subroutine PrepareFdtdOutgpl(ncyc)

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
  end subroutine PrepareFdtdOutgpl


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
