!----------------------------------------------------------------------
!
!  module: fdtd-out-r
!
!  ascii output module.
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
!     InitOutputParameters
!     WriteOutputHeader
!     InitOutput             used in max3d.f90
!     Output(ncyc)           used in max3d.f90
!            WriteEH
!            WriteComp
!            WriteDi
!            WriteData
!     DataPrepOutput(ncyc)   used in max3d.f90 (between StepH and StepE)
!     LoadPx
!     LoadPy
!     LoadPz

! ---------------------------------------------------------------------

!  output mode 'abcd':
!  ab = component(s) 'Ex', 'Ey', 'Ez', 'Hx', 'Hy','Hz','EH','Di',
!                     'En', 'Px', 'Py', oder 'Pz'
!  c = 'E':  one file (for all time steps)
!  c = 'M':  multiple files
!  d = 'R':  spatially resolved output
!  d = 'S':  spatially integrated output

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


module fdtd_out

  use constant
  use strings
  use mpiworld
  use grid  
  use fdtd

  implicit none
  save

  ! --- Constants

  character(len=255), parameter :: pfxoutput = 'output'
  integer, parameter :: PARTSMAXGPL = 100

  ! --- Types

  type T_OUTBAS 

     ! Output files:
     character(STRLNG) fn

     ! Output mode:
     character(len=10) :: Mode

     ! Spatial area
     integer ns, ne, dn
     integer is, ie, di
     integer js, je, dj
     integer ks, ke, dk

     ! Other
     integer NumNodes

  end type T_OUTBAS

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


  subroutine CreateOutFDTD

    use constant
    implicit none

    integer ::  ios,n,i,err
    character(len=STRLNG) :: file, str

    call ReadConfig
    call Initialize
    call WriteHeader
    

  contains
  
    subroutine ReadConfig
      
      implicit none

      character(len=STRLNG) :: file

      file=cat2(pfxoutput,sfxin)
      
      open(UNITTMP,FILE=file,STATUS='unknown')
      n = 0
      do
         read(UNITTMP,IOSTAT=ios,FMT=*) str
         if(ios .ne. 0) exit
         if(str(1:4).eq. '#GPL') then
            n = n+1
            call ReadConfigObject(gpl(n),UNITTMP)
            if(n .ge. PARTSMAXGPL) exit
         endif
      enddo
      close(UNITTMP)
      PARTSGPL=n  
      
    end subroutine ReadConfig


    subroutine WriteHeader

      implicit none

      do n=1, PARTSGPL
         open(UNITTMP,FILE=gpl(n)%fn,STATUS='unknown')      
         call WriteHeaderObject(gpl(n),UNITTMP,'Produced by DataOutGpl()')
         close(UNITTMP)           
      enddo

    end subroutine WriteHeader

    subroutine Initialize

      implicit none

      i = 1
      DataIndxGpl(1)=1
      do n=1, PARTSGPL
         if( gpl(n)%Mode(1:2) .eq. 'En' .or. &
              gpl(n)%Mode(1:2) .eq. 'Px' .or. &
              gpl(n)%Mode(1:2) .eq. 'Py' .or. &
              gpl(n)%Mode(1:2) .eq. 'Pz' ) then
            i = i + gpl(n)%NumNodes               
         endif
         DataIndxGpl(n+1)=i
      enddo
      allocate(DataGpl(1:i),STAT=err)
      if(err .ne. 0) then
         write(6,*) 'Error 1 in InitOutGpl()'
         stop
         
      endif

    end subroutine Initialize


    subroutine ReadConfigObject(this, funit)

      ! Initialisiert Modul OutBas	
      implicit none
      
      type (T_OUTBAS) :: this
      integer funit, isteps, jsteps, ksteps

      ! Read Modul (Type) Data
      read(funit,*) this%fn, this%Mode
      read(funit,*) this%ns, this%ne, this%dn
      read(funit,*) this%is, this%ie, this%di
      read(funit,*) this%js, this%je, this%dj
      read(funit,*) this%ks, this%ke, this%dk
      
      this%fn = cat2(this%fn, mpi_sfxout)
      
      ! check ranges	
      this%is = max(IBEG, this%is); 
      this%ie = min(IEND, this%ie);
      
      this%js = max(JBEG, this%js); 
      this%je = min(JEND, this%je);
      
      this%ks = max(KBEG, this%ks); 
      !    this%ks = max(KMIN, this%ks); ! >>>DEBUG 
      this%ke = min(KEND, this%ke);
      !    this%ke = min(Kmax, this%ke); ! >>>DEBUG
      
      isteps = max(int((this%ie-this%is+this%di)/this%di),0)
      jsteps = max(int((this%je-this%js+this%dj)/this%dj),0)
      jsteps = max(int((this%ke-this%ks+this%dk)/this%dk),0)
      this%NumNodes = isteps*jsteps*ksteps
      
    end subroutine ReadConfigObject

    subroutine WriteHeaderObject(this,funit, TitleStr)
      
      implicit none
      
      type (T_OUTBAS) :: this
      integer :: funit
      character(len=*) :: TitleStr
      real(8) :: ts,te,xs,xe,ys,ye,zs,ze
      
      ! Do some Preparations
      ts = this%ns*DT+GT
      te = this%ne*DT+GT
      xs = this%is*SX
      xe = this%ie*SX
      ys = this%js*SY
      ye = this%je*SY
      zs = this%ks*SZ
      ze = this%ke*SZ
      
      ! Write Data in File Header
      write(funit,'(A1,A30)') '#',TitleStr
      write(funit,'(A1,A30)') '#','label1'
      write(funit,'(A1,A30)') '#','label2'
      write(funit,'(A1,A30)') '#','label3'
      write(funit,'(A1,A30)') '#','label4'
      write(funit,'(A1,A30)') '#','label5'
      write(funit,'(A1,A30)') '#','label6'
      write(funit,'(A1,A30,A4)') '#',this%fn(1:27),this%Mode
      write(funit,'(A1,3I8)') '#',this%ns,this%ne,this%dn
      write(funit,'(A1,2E15.6E3)') '#',ts,te
      write(funit,'(A1,3I8)') '#',this%is,this%ie,this%di
      write(funit,'(A1,2E15.6E3)') '#',xs,xe
      write(funit,'(A1,3I8)') '#',this%js,this%je,this%dj
      write(funit,'(A1,2E15.6E3)') '#',ys,ye
      write(funit,'(A1,3I8)') '#',this%ks,this%ke,this%dk
      write(funit,'(A1,2E15.6E3)') '#',zs,ze

    end subroutine WriteHeaderObject


  end subroutine CreateOutFDTD


  subroutine DestroyOutFDTD
    implicit none

    deallocate(DataGpl)

  end subroutine DestroyOutFDTD


  subroutine WriteOutFDTD(ncyc)

    implicit none
    
    integer ncyc, n, ios
    character(STRLNG) fn, fnh
    
    do n=1, PARTSGPL
       
       ! output ?
       if(mod(ncyc, gpl(n)%dn) .eq. 0 .and. &
            ncyc .ge. gpl(n)%ns       .and. &
            ncyc .le. gpl(n)%ne ) then
          
          ! open file
          if( gpl(n)%Mode(3:3) .eq. 'M') then
             fn = gpl(n)%fn
             write(fnh,*) ncyc     
             fnh=adjustl(fnh)
             fn((len_trim(fn)+1):STRLNG) = fnh         
             open(UNITTMP,IOSTAT=ios, FILE=fn) 
          else
             fn = gpl(n)%fn
             open(UNITTMP,IOSTAT=ios,POSITION= "APPEND", FILE=fn) 
          endif

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

          close(UNITTMP)
       endif
    enddo
       
  contains

    ! **************************************************************** !   
    subroutine WriteEH(funit)
      ! output of all field commponents
      implicit none
      integer funit,i,j,k
      ! Start Code
      do k=gpl(n)%ks, gpl(n)%ke, gpl(n)%dk  
         do j=gpl(n)%js, gpl(n)%je, gpl(n)%dj      
            do i=gpl(n)%is, gpl(n)%ie, gpl(n)%di 
               write(funit,'(6E15.6E3))') Ex(i,j,k),  Ey(i,j,k), &
                    Ez(i,j,k), Hx(i,j,k), Hy(i,j,k), Hz(i,j,k)
            enddo
            if(gpl(n)%is .ne. gpl(n)%ie) write(funit,*)
         enddo
         if(gpl(n)%js .ne. gpl(n)%je) write(funit,*)
      enddo
    end subroutine WriteEH


    ! **************************************************************** !
    subroutine WriteComp(Comp,funit)
      ! output of a single field component
      implicit none
      real(8), dimension(IMIN:KMAX,JMIN:KMAX,KMIN:KMAX) :: Comp
      integer funit,i,j,k
      ! Start Code
      do k=gpl(n)%ks, gpl(n)%ke, gpl(n)%dk  
         do j=gpl(n)%js, gpl(n)%je, gpl(n)%dj      
            do i=gpl(n)%is, gpl(n)%ie, gpl(n)%di 
               write(funit,'(E15.6))') Comp(i,j,k)
            enddo
            if(gpl(n)%is .ne. gpl(n)%ie) write(funit,*) 
         enddo
         if(gpl(n)%js .ne. gpl(n)%je) write(funit,*)
      enddo

    end subroutine WriteComp


    ! **************************************************************** !
    subroutine WriteDi(funit)
      ! output of the dielectric constant
      implicit none
      integer funit,i,j,k
      ! Start Code
      write(6,*) gpl(n)%is, gpl(n)%ie
      write(6,*) gpl(n)%js, gpl(n)%je
      write(6,*) gpl(n)%ks, gpl(n)%ke 
      do k=gpl(n)%ks, gpl(n)%ke, gpl(n)%dk
         do j=gpl(n)%js, gpl(n)%je, gpl(n)%dj      
            do i=gpl(n)%is, gpl(n)%ie, gpl(n)%di 
               write(funit,'(E15.6E3))') 1.0/epsinv(i,j,k)
            enddo
            if(gpl(n)%is .ne. gpl(n)%ie) write(funit,*)
         enddo
         if(gpl(n)%js .ne. gpl(n)%je) write(funit,*)
      enddo
      
    end subroutine WriteDi 


    ! **************************************************************** !
    subroutine WriteData(n,funit)
      ! output of data from field DataGpl 
      implicit none
      integer :: funit,i,j,k,n,m
      real(8) :: sum
      ! Start Code 
      sum = 0.0
      m = DataIndxGpl(n)
      if( gpl(n)%Mode(4:4) .eq. 'S' ) then      ! integration 	
         do k=gpl(n)%ks, gpl(n)%ke, gpl(n)%dk  
            do j=gpl(n)%js, gpl(n)%je, gpl(n)%dj      
               do i=gpl(n)%is, gpl(n)%ie, gpl(n)%di 
                  sum = sum+DataGpl(m)
                  m = m+1
               enddo
            enddo
         enddo
         write(funit,*) sum
      else                                     ! direct output
         do k=gpl(n)%ks, gpl(n)%ke, gpl(n)%dk  
            do j=gpl(n)%js, gpl(n)%je, gpl(n)%dj      
               do i=gpl(n)%is, gpl(n)%ie, gpl(n)%di 
                  write(funit,'(E15.6E3))') DataGpl(m)
                  m = m+1
               enddo
               if(gpl(n)%is .ne. gpl(n)%ie) write(funit,*)
            enddo
            if(gpl(n)%js .ne. gpl(n)%je) write(funit,*)
         enddo
      endif
    end subroutine WriteData

  end subroutine WriteOutFDTD


  ! ***************************************************************** !

  subroutine PrepareOutFDTD(ncyc)

    implicit none

    ! Variables
    integer ncyc, n, i

    ! Loop over all output units
    do n=1, PARTSGPL

       ! output ?
       if(mod(ncyc, gpl(n)%dn) .eq. 0 .and. &
            ncyc .ge. gpl(n)%ns       .and. &
            ncyc .le. gpl(n)%ne ) then

          ! First part of Px,y,z calculation
          select case (gpl(n)%Mode(1:2))
          case('Px')
             DataGpl(DataIndxGpl(n):DataIndxGpl(n+1))=0.0
             call LoadPx(n)
          case('Py')
             DataGpl(DataIndxGpl(n):DataIndxGpl(n+1))=0.0
             call LoadPy(n)
          case('Pz')
             DataGpl(DataIndxGpl(n):DataIndxGpl(n+1))=0.0
             call LoadPz(n)
          end select
          
       endif
    enddo
  end subroutine PrepareOutFDTD


  subroutine LoadPx(n)
    ! Calculates x-component of the Poynting vector P
    ! Stores result in DataGpl
    ! Localization: Px[i,j,k] = Px(i+1/2,j,k)
    use constant
    use grid
    implicit none
    integer :: i,j,k,n,m
    real(8) :: val
    ! Code
    m = DataIndxGpl(n)
    do k=gpl(n)%ks, gpl(n)%ke, gpl(n)%dk 
       do j=gpl(n)%js, gpl(n)%je, gpl(n)%dj 
          do i=gpl(n)%is, gpl(n)%ie, gpl(n)%di
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


  subroutine LoadPy(n)
    ! Calculates y-component of the Poynting vector P
    ! Stores result in DataGpl
    ! Localization: Py[i,j,k] = Py(i,j+1/2,k)
    use constant
    use grid
    implicit none
    integer :: i,j,k,n,m
    real(8) :: val
    ! Code
    m = DataIndxGpl(n)
    do k=gpl(n)%ks, gpl(n)%ke, gpl(n)%dk
       do j=gpl(n)%js, gpl(n)%je, gpl(n)%dj
          do i=gpl(n)%is, gpl(n)%ie, gpl(n)%di
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


  subroutine LoadPz(n)
    ! Calculates z-component of the Poynting vector P
    ! Stores result in DataGpl
    ! Localization: Pz[i,j,k] = Pz(i,j,k+1/2)
    use constant
    use grid
    implicit none
    integer :: i,j,k,n, m
    real(8) :: val
    ! Code
    m = DataIndxGpl(n)
    do k=gpl(n)%ks, gpl(n)%ke, gpl(n)%dk
       do j=gpl(n)%js, gpl(n)%je, gpl(n)%dj
          do i=gpl(n)%is, gpl(n)%ie, gpl(n)%di
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


  subroutine LoadEn(n)
    ! Calculates energy density En
    ! Stores result in DataGpl
    ! Localization: En[i,j,k] = En(i,j,k)
    use constant
    use grid
    implicit none
    integer :: i,j,k,n,m
    real(8) :: EEx,EEy,EEz,EHx,EHy,EHz,eps
    ! Start Code
    m = DataIndxGpl(n)
    do k=gpl(n)%ks, gpl(n)%ke, gpl(n)%dk
       do j=gpl(n)%js, gpl(n)%je, gpl(n)%dj
          do i=gpl(n)%is, gpl(n)%ie, gpl(n)%di
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

end module fdtd_out
