!----------------------------------------------------------------------
!
!  module: outset
!
!  ascii (gnuplot) output module
!
!  subs:
!
!  InitializeOutsetObj
!  FinalizeOutsetObj
!  OpenOutsetObj
!  CloseOutsetObj
!  WriteHeaderOutsetObj
!  WriteDataOutsetObj
!  
!----------------------------------------------------------------------


!======================================================================
!
!

module outset

  use constant
  use strings
  use mpiworld
  use reglist
  use outlist
  use grid

! outset modules 
  use fdtd_outset

  implicit none
  private
  save
  
  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'OUTSET'

  ! --- Public Methods

  public :: InitializeOutsetObj
  public :: FinalizeOutsetObj
  public :: OpenOutsetObj
  public :: CloseOutsetObj
  public :: WriteHeaderOutsetObj
  public :: WriteDataOutsetObj

  ! --- Public Data

  ! --- Constants

  character(len=STRLNG), parameter :: outsetsfx = '.set'

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializeOutsetObj(out)

    type(T_OUT) :: out


! override out%snap as vtk does not support multiple frames in one file

    out%snap = .true.

! call various finalization methods
    select case ( out%modl ) 

    case ("FDTD")
       call  InitializeFdtdOutsetObj(out)
    end select

    if ( .not. out%snap ) then
       call OpenOutsetObj(out, 0, .true.)
       call CloseOutsetObj(out)
    endif


! if not in snapshot mode: write a header once at initialization 

  end subroutine InitializeOutsetObj


!----------------------------------------------------------------------

  subroutine FinalizeOutsetObj(out)

    type(T_OUT) :: out

    select case ( out%modl ) 

    case ("FDTD")
       call FinalizeFdtdOutsetObj(out)
    end select

  end subroutine FinalizeOutsetObj


!----------------------------------------------------------------------

  subroutine OpenOutsetObj(out, ncyc, writehdr)

    type(T_OUT) :: out
    integer :: ncyc
    logical :: writehdr
    character(len=STRLNG) :: name, stepstr
 
    if ( out%snap ) then
       stepstr = TRIM(i2str(ncyc))
       name = cat5(out%filename,"_",TRIM(i2str(ncyc)),mpi_sfx,outsetsfx)
    else
       name = cat3(out%filename,mpi_sfx,outsetsfx)
    endif
   
    if ( writehdr ) then 
       
       M4_WRITE_DBG({"opening ",TRIM(out%filename),", as ",TRIM(name),", unit ", &
            TRIM(i2str(out%funit)), " (new)"})
       open(out%funit,FILE=name,STATUS="UNKNOWN",POSITION="REWIND")
       call WriteHeaderOutsetObj(out,ncyc)

    else

       M4_WRITE_DBG({"opening ",TRIM(out%filename),", as ",TRIM(name),", unit ", &
            TRIM(i2str(out%funit)), " (append)"})
       open(out%funit,FILE=name,STATUS="UNKNOWN",POSITION="APPEND")

    endif

  end subroutine OpenOutsetObj

!----------------------------------------------------------------------

  subroutine CloseOutsetObj(out)

    type(T_OUT) :: out

    M4_WRITE_DBG({"closing ",TRIM(out%filename),", unit ", TRIM(i2str(out%funit))})

    close(out%funit)

  end subroutine CloseOutsetObj

!----------------------------------------------------------------------

  subroutine WriteHeaderOutsetObj(out,ncyc)

    type(T_OUT) :: out
    integer :: ncyc
    character(len=STRLNG) :: dataset
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    reg = regobj(out%regidx)
    
  
    M4_WRITE_DBG({"write header ",TRIM(out%filename)})
    
    write(out%funit,*) '! (SET-) FIELD DATA FILE '
    write(out%funit,*) '! META: M4_VERSION(), M4_FLAVOUR()'
    write(out%funit,*) '! ISBOX: ',reg%isbox           
    write(out%funit,*) '! NUMNODES: ',reg%numnodes 
    write(out%funit,*) '! IRANGE: ',reg%is,reg%ie,reg%di 
    write(out%funit,*) '! JRANGE: ',reg%js,reg%je,reg%dj
    write(out%funit,*) '! KRANGE: ',reg%ks,reg%ke,reg%dk


  end subroutine WriteHeaderOutsetObj


!----------------------------------------------------------------------

  subroutine WriteDataOutsetObj(out, ncyc, mode)

    type(T_OUT) :: out
    integer :: ncyc
    logical :: mode  
    
    M4_WRITE_DBG({"write data ",TRIM(out%filename),", numnodes ", TRIM(i2str(out%numnodes))})

    select case ( out%modl ) 

! call various output methods

    case ("FDTD")

       if ( mode ) then 
          call OpenOutsetObj(out, ncyc, out%snap)
!          if ( out%numnodes .gt. 0 ) write(out%funit,*)
       end if
       call WriteDataFdtdOutsetObj(out,mode)
       if ( mode ) call CloseOutsetObj(out)
    
    case default
       if ( mode ) then 
          call OpenOutsetObj(out, ncyc, out%snap)
          write(out%funit,*) "OUTPUT MODULE NOT IMPLEMENTED" 
       end if
       if ( mode ) call CloseOutsetObj(out)

    end select

  end subroutine WriteDataOutsetObj
  
!----------------------------------------------------------------------

end module outset

!
! Authors:  J.Hamm 
! Modified: 15/02/2007
!
!======================================================================
