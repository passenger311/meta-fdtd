!----------------------------------------------------------------------
!
!  module: outgpl-r
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
!     InitOutgpl             used in meta3.f90
!     Outgpl(ncyc)           used in meta3.f90
!            WriteEH
!            WriteComp
!            WriteDi
!            WriteData
!     DataPrepOutgpl(ncyc)   used in meta3.f90 (between StepH and StepE)
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


module outgpl

  use constant
  use strings
  use mpiworld
  use grid  

  implicit none
  save

  ! --- Constants

  character(len=255), parameter :: pfxoutgpl = 'outgpl'
  integer, parameter :: MAXOBJGPL = 100

  ! --- Types

  type T_OUTBAS 

     ! Outgpl files:
     character(STRLNG) fn

     ! Outgpl mode:
     character(len=10) :: Mode

     ! Spatial area
     integer ns, ne, dn
     integer is, ie, di
     integer js, je, dj
     integer ks, ke, dk

     ! Other
     integer idx
     integer NumNodes

  end type T_OUTBAS

  ! --- Variables

  type(T_OutBas) :: objgpl(MAXOBJGPL)
  integer :: numobjgpl


contains


  subroutine InitializeOutgpl

    integer ::  ios,n,i,err
    character(len=STRLNG) :: file, str

    call ReadOutgpl
    call WriteHeaderOutgpl

  end subroutine InitializeOutgpl

  
  subroutine ReadOutgpl
    
    character(len=STRLNG) :: file, str
    integer :: n, ios
    
    file=cat2(pfxoutgpl,sfxin)
    
    open(UNITTMP,FILE=file,STATUS='unknown')
    n = 0
    do
       read(UNITTMP,IOSTAT=ios,FMT=*) str
       if(ios .ne. 0) exit
       if(str(1:4).eq. '(OUTGPL') then
          n = n+1
          call ReadOutObjgpl(objgpl(n))
          objgpl(n)%idx = n  
          if(n .ge. MAXOBJGPL) exit
       endif
    enddo
    close(UNITTMP)
    numobjgpl=n  
    
  end subroutine ReadOutgpl
  

  subroutine WriteHeaderOutgpl
    
    type(T_OutBas) :: objgpl(MAXOBJGPL)
    integer :: n

    do n=1, numobjgpl
       open(UNITTMP,FILE=objgpl(n)%fn,STATUS='unknown')      
       call WriteHeaderOutObjgpl(objgpl(n))
       close(UNITTMP)           
    enddo
    
  end subroutine WriteHeaderOutgpl
  

  subroutine ReadOutObjgpl(gpl)
    
    type (T_OUTBAS) :: gpl
    integer UNITTMP, isteps, jsteps, ksteps
    
    ! Read Modul (Type) Data
    read(UNITTMP,*) gpl%fn, gpl%Mode
    read(UNITTMP,*) gpl%ns, gpl%ne, gpl%dn
    read(UNITTMP,*) gpl%is, gpl%ie, gpl%di
    read(UNITTMP,*) gpl%js, gpl%je, gpl%dj
    read(UNITTMP,*) gpl%ks, gpl%ke, gpl%dk
    
    gpl%fn = cat2(gpl%fn, mpi_sfxout)
    
    ! check ranges	
    gpl%is = max(IBEG, gpl%is); 
    gpl%ie = min(IEND, gpl%ie);
    
    gpl%js = max(JBEG, gpl%js); 
    gpl%je = min(JEND, gpl%je);
    
    gpl%ks = max(KBEG, gpl%ks); 
    gpl%ke = min(KEND, gpl%ke);
    
    isteps = max(int((gpl%ie-gpl%is+gpl%di)/gpl%di),0)
    jsteps = max(int((gpl%je-gpl%js+gpl%dj)/gpl%dj),0)
    jsteps = max(int((gpl%ke-gpl%ks+gpl%dk)/gpl%dk),0)
    gpl%NumNodes = isteps*jsteps*ksteps
    
  end subroutine ReadOutObjgpl
  
  subroutine WriteHeaderOutObjgpl(gpl)
    
    type (T_OUTBAS) :: gpl
    real(8) :: ts,te,xs,xe,ys,ye,zs,ze
    
    ts = gpl%ns*DT+GT
    te = gpl%ne*DT+GT
    xs = gpl%is*SX
    xe = gpl%ie*SX
    ys = gpl%js*SY
    ye = gpl%je*SY
    zs = gpl%ks*SZ
    ze = gpl%ke*SZ
    
    ! Write Data in File Header
    write(UNITTMP,'(A1,A30)') '# Generated by Outgpl module'
    !      write(UNITTMP,'(A1,A30)') '#','label1'
    !      write(UNITTMP,'(A1,A30)') '#','label2'
    !      write(UNITTMP,'(A1,A30)') '#','label3'
    !      write(UNITTMP,'(A1,A30)') '#','label4'
    !      write(UNITTMP,'(A1,A30)') '#','label5'
    !      write(UNITTMP,'(A1,A30)') '#','label6'
    write(UNITTMP,'(A1,A30,A4)') '#',gpl%fn(1:27),gpl%Mode
    write(UNITTMP,'(A1,3I8)') '#',gpl%ns,gpl%ne,gpl%dn
    write(UNITTMP,'(A1,2E15.6E3)') '#',ts,te
    write(UNITTMP,'(A1,3I8)') '#',gpl%is,gpl%ie,gpl%di
    write(UNITTMP,'(A1,2E15.6E3)') '#',xs,xe
    write(UNITTMP,'(A1,3I8)') '#',gpl%js,gpl%je,gpl%dj
    write(UNITTMP,'(A1,2E15.6E3)') '#',ys,ye
    write(UNITTMP,'(A1,3I8)') '#',gpl%ks,gpl%ke,gpl%dk
    write(UNITTMP,'(A1,2E15.6E3)') '#',zs,ze
    
  end subroutine WriteHeaderOutObjgpl
  

  subroutine FinalizeOutgpl
   
  end subroutine FinalizeOutgpl

  subroutine OpenOutObjgpl(gpl, ncyc, ret)
        
    type (T_OUTBAS) :: gpl
    integer :: ncyc
    logical :: ret

    integer :: ios
    character(len=STRLNG) :: fn, fnh

    ret = .False.
    ! outgpl ?
    if(mod(ncyc, gpl%dn) .eq. 0 .and. &
         ncyc .ge. gpl%ns       .and. &
         ncyc .le. gpl%ne ) then
       
       ! open file
       if( gpl%Mode(3:3) .eq. 'M') then
          fn = gpl%fn
          write(fnh,*) ncyc     
          fnh=adjustl(fnh)
          fn((len_trim(fn)+1):STRLNG) = fnh         
          open(UNITTMP,IOSTAT=ios, FILE=fn) 
       else
          fn = gpl%fn
          open(UNITTMP,IOSTAT=ios,POSITION= "APPEND", FILE=fn) 
       endif
       ret = .True.
    endif
    
  end subroutine OpenOutObjgpl

  subroutine CloseOutObjgpl

    close(UNITTMP)

  end subroutine CloseOutObjgpl


end module outgpl