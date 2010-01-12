!----------------------------------------------------------------------
!
!  module: outR
!
!  ascii (R) output module
!
!  subs:
!
!  InitializeOutRObj
!  FinalizeOutRObj
!  OpenOutRObj
!  CloseOutRObj
!  WriteHeaderOutRObj
!  WriteDataOutRObj
!  
!----------------------------------------------------------------------


!======================================================================
!
!

module outR

  use constant
  use strings
  use mpiworld
  use reglist
  use outlist
  use out_calc

! outr modules
   use fdtd_outr
M4_FOREACH_OUTR({use }, {_OUTR
  })

  implicit none
  private
  save
  
  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'OUTR'

  ! --- Public Methods

  public :: InitializeOutRObj
  public :: FinalizeOutRObj
  public :: OpenOutRObj
  public :: CloseOutRObj
  public :: WriteHeaderOutRObj
  public :: WriteDataOutRObj

  ! --- Public Data

  ! --- Constants

  character(len=STRLNG), parameter :: outRsfx = '.r'
  character(len=STRLNG), parameter :: outinfosfx = '.info'

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializeOutRObj(out)

    type(T_OUT) :: out

! call various finalization methods
    select case ( out%modl ) 
    
    case ("FDTD")
       call  InitializeFdtdOutRObj(out)

    M4_FOREACH_OUTR2({case ("},{")
       call Initialize},{OutRObj(out)
             })
    end select

    if ( .not. out%snap ) then
       call OpenOutRObj(out, 0, .true.)
       call CloseOutRObj(out)
    endif


! if not in snapshot mode: write a header once at initialization 

  end subroutine InitializeOutRObj


!----------------------------------------------------------------------

  subroutine FinalizeOutRObj(out)

    type(T_OUT) :: out

    select case ( out%modl ) 

! call various finalization methods   
    case ("FDTD")
       call FinalizeFdtdOutRObj(out)
       M4_FOREACH_OUTR2({case ("},{")
             call Finalize},{OutRObj(out)
             })
    end select

  end subroutine FinalizeOutRObj


!----------------------------------------------------------------------

  subroutine OpenOutRObj(out, ncyc, writehdr)

    type(T_OUT) :: out
    integer :: ncyc
    logical :: writehdr
    character(len=STRLNG) :: namer, nameinfo, stepstr
 
    if ( out%snap ) then
       stepstr = TRIM(i2str(ncyc))
       namer = cat5(out%filename,"_",TRIM(i2str(ncyc)),mpi_sfx,outRsfx)
       nameinfo = cat5(out%filename,"_",TRIM(i2str(ncyc)),mpi_sfx,outinfosfx)
    else
       namer = cat3(out%filename,mpi_sfx,outRsfx)
       nameinfo = cat3(out%filename,mpi_sfx,outinfosfx)
    endif
   
    if ( writehdr ) then 
       
       M4_WRITE_DBG({"opening ",TRIM(out%filename),", as ",TRIM(nameinfo),", unit ", &
            TRIM(i2str(out%funit)), " (new)"})
       open(out%funit,FILE=nameinfo,STATUS="UNKNOWN",POSITION="REWIND")
       call WriteHeaderOutRObj(out,ncyc)
       close(out%funit)
       open(out%funit,FILE=namer,STATUS="UNKNOWN",POSITION="REWIND")

    else

       M4_WRITE_DBG({"opening ",TRIM(out%filename),", as ",TRIM(namer),", unit ", &
            TRIM(i2str(out%funit)), " (append)"})
       open(out%funit,FILE=namer,STATUS="UNKNOWN",POSITION="APPEND")

    endif

  end subroutine OpenOutRObj

!----------------------------------------------------------------------

  subroutine CloseOutRObj(out)

    type(T_OUT) :: out

    M4_WRITE_DBG({"closing ",TRIM(out%filename),", unit ", TRIM(i2str(out%funit))})

    close(out%funit)

  end subroutine CloseOutRObj

!----------------------------------------------------------------------

  subroutine WriteHeaderOutRObj(out,ncyc)

    type(T_OUT) :: out
    integer :: ncyc
    type(T_REG) :: reg

    reg = regobj(out%regidx)
    
    M4_WRITE_DBG({"write header ",TRIM(out%filename)})

    write(out%funit,*) 'snapshot mode tnum tmin tmax'
    if(out%snap) then
       write(out%funit,*) out%snap, out%mode, 1, ncyc, ncyc
    else
       write(out%funit,*) out%snap, out%mode, out%numsteps, out%ns, out%ne
    endif

  end subroutine WriteHeaderOutRObj


!----------------------------------------------------------------------

  subroutine WriteDataOutRObj(out, ncyc, mode)

    type(T_OUT) :: out
    integer :: ncyc
    logical :: mode  
    
    M4_WRITE_DBG({"write data ",TRIM(out%filename),", numnodes ", TRIM(i2str(out%numnodes))})


    select case ( out%modl ) 

! call various output methods

   case ("FDTD")

       if ( mode ) then 
          call OpenOutRObj(out, ncyc, out%snap)
          if ( out%numnodes .gt. 1 .and. out%mode .ne. 'S' ) write(out%funit,*)
       end if
       call WriteDataFdtdOutRObj(out,mode)
       if ( mode ) call CloseOutRObj(out)
    
    M4_FOREACH_OUTR2({case ("},{")
       if ( mode ) then 
          call OpenOutRObj(out, ncyc, out%snap)
          if ( out%numnodes .gt. 1 ) write(out%funit,*) 
       endif
       call WriteData},{OutRObj(out,mode)
       if ( mode ) call CloseOutRObj(out)
       })

    case default
       if ( mode ) then 
          call OpenOutRObj(out, ncyc, out%snap)
          write(out%funit,*) "OUTPUT MODULE NOT IMPLEMENTED" 
       end if
       if ( mode ) call CloseOutRObj(out)

    end select

  end subroutine WriteDataOutRObj


!----------------------------------------------------------------------


end module outR

!
! Authors:  Andreas Pusch 
! Modified: 24/09/2009
!
!======================================================================
