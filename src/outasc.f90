!----------------------------------------------------------------------
!
!  module: outasc-r
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


module outasc

  use constant
  use strings
  use mpiworld
  use grid  

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
  integer :: partsgpl


contains


  subroutine InitializeOutAsc

    use constant
    implicit none

    integer ::  ios,n,i,err
    character(len=STRLNG) :: file, str

    call ReadConfig
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
         if(str(1:4).eq. '#ASC') then
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
         call WriteHeaderObject(gpl(n),UNITTMP,'Ascii Output')
         close(UNITTMP)           
      enddo

    end subroutine WriteHeader


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
      this%ke = min(KEND, this%ke);
      
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
!      write(funit,'(A1,A30)') '#','label1'
!      write(funit,'(A1,A30)') '#','label2'
!      write(funit,'(A1,A30)') '#','label3'
!      write(funit,'(A1,A30)') '#','label4'
!      write(funit,'(A1,A30)') '#','label5'
!      write(funit,'(A1,A30)') '#','label6'
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


  end subroutine InitializeOutAsc


  subroutine FinalizeOutAsc
    implicit none


  end subroutine FinalizeOutAsc

  subroutine WriteOutAsc(n, ncyc, ret)

    integer :: n, ncyc
    logical :: ret

    integer :: ios
    character(len=STRLNG) :: fn, fnh

    ret = .False.
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
       ret = .True.
    endif
    
  end subroutine WriteOutAsc

  subroutine CloseOutAsc

    close(UNITTMP)

  end subroutine CloseOutAsc


end module outasc
