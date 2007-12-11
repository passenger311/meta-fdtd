!----------------------------------------------------------------------
!
!  module: outgpl
!
!  ascii (gnuplot) output module
!
!  subs:
!
!  InitializeOutgplObj
!  FinalizeOutgplObj
!  OpenOutgplObj
!  CloseOutgplObj
!  WriteHeaderOutgplObj
!  WriteDataOutgplObj
!  
!----------------------------------------------------------------------


!======================================================================
!
!

module outgpl

  use constant
  use strings
  use mpiworld
  use reglist
  use outlist

! outgpl modules 
  use fdtd_outgpl
M4_FOREACH_OUTGPL({use }, {
  })

  implicit none
  private
  save
  
  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'OUTGPL'

  ! --- Public Methods

  public :: InitializeOutgplObj
  public :: FinalizeOutgplObj
  public :: OpenOutgplObj
  public :: CloseOutgplObj
  public :: WriteHeaderOutgplObj
  public :: WriteDataOutgplObj

  ! --- Public Data

  ! --- Constants

  character(len=STRLNG), parameter :: outgplsfx = '.gpl'

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializeOutgplObj(out)

    type(T_OUT) :: out

! call various finalization methods
    select case ( out%modl ) 

    case ("FDTD")
       call  InitializeFdtdOutgplObj(out)
       call OpenOutgplObj(out, 0, .not. out%snap)
       call CloseOutgplObj(out)

    M4_FOREACH_OUTGPL2({case ("},{")
       call Initialize},{OutgplObj(UNITTMP)
       call OpenOutgplObj(out, 0, .not. out%snap)
       call CloseOutgplObj(out)
             })
    end select

! if not in snapshot mode: write a header once at initialization 

  end subroutine InitializeOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeOutgplObj(out)

    type(T_OUT) :: out

    select case ( out%modl ) 

! call various finalization methods
    case ("FDTD")
       call FinalizeFdtdOutgplObj(out)
       M4_FOREACH_OUTGPL2({case ("},{")
             call Finalize},{OutgplObj
             })
    end select

  end subroutine FinalizeOutgplObj


!----------------------------------------------------------------------

  subroutine OpenOutgplObj(out, ncyc, writehdr)

    type(T_OUT) :: out
    integer :: ncyc
    logical :: writehdr
    character(len=STRLNG) :: name, stepstr
 
    if ( out%snap ) then
       stepstr = TRIM(i2str(ncyc))
       name = cat5(out%filename,"_",TRIM(i2str(ncyc)),mpi_sfx,outgplsfx)
    else
       name = cat3(out%filename,mpi_sfx,outgplsfx)
    endif
   
    if ( writehdr ) then 
       
       M4_WRITE_DBG({"opening ",TRIM(out%filename),", as ",TRIM(name),", unit ", &
            TRIM(i2str(out%funit)), " (new)"})
       open(out%funit,FILE=name,STATUS="UNKNOWN",POSITION="REWIND")
       call WriteHeaderOutgplObj(out,ncyc)

    else

       M4_WRITE_DBG({"opening ",TRIM(out%filename),", as ",TRIM(name),", unit ", &
            TRIM(i2str(out%funit)), " (append)"})
       open(out%funit,FILE=name,STATUS="UNKNOWN",POSITION="APPEND")

    endif

  end subroutine OpenOutgplObj

!----------------------------------------------------------------------

  subroutine CloseOutgplObj(out)

    type(T_OUT) :: out

    M4_WRITE_DBG({"closing ",TRIM(out%filename),", unit ", TRIM(i2str(out%funit))})

    close(out%funit)

  end subroutine CloseOutgplObj

!----------------------------------------------------------------------

  subroutine WriteHeaderOutgplObj(out,ncyc)

    type(T_OUT) :: out
    integer :: ncyc
    type(T_REG) :: reg
    integer :: numn 

    numn = (out%ne - out%ns + 1) / out%dn + 1
    
    reg = regobj(out%regidx)
    
    M4_WRITE_DBG({"write header ",TRIM(out%filename)})

    write(out%funit,*) '# ',out%fmt                 ! format
    write(out%funit,*) '# ',out%snap                ! snapshot mode
    write(out%funit,*) '# ',out%modl                ! module
    write(out%funit,*) '# ',out%fn                  ! function
    write(out%funit,*) '# ',out%mode                ! mode
    if ( out%snap ) then
       write(out%funit,*) '# ',1                    ! number of time points
       write(out%funit,*) '# ',ncyc,ncyc,1          ! time frame
    else
       write(out%funit,*) '# ',numn                 ! number of time points
       write(out%funit,*) '# ',out%ns,out%ne,out%dn ! time frame
    endif
    write(out%funit,*) '# ',reg%isbox               ! mode
    write(out%funit,*) '# ',reg%numnodes            ! number of points
    write(out%funit,*) '# ',reg%is,reg%ie,reg%di    ! space box
    write(out%funit,*) '# ',reg%js,reg%je,reg%dj
    write(out%funit,*) '# ',reg%ks,reg%ke,reg%dk

  end subroutine WriteHeaderOutgplObj


!----------------------------------------------------------------------

  subroutine WriteDataOutgplObj(out, ncyc, mode)

    type(T_OUT) :: out
    integer :: ncyc
    logical :: mode  
    
    M4_WRITE_DBG({"write data ",TRIM(out%filename),", numnodes ", TRIM(i2str(out%numnodes))})

    
    select case ( out%modl ) 

! call various output methods

    case ("FDTD")

       if ( mode ) call OpenOutgplObj(out, ncyc, out%snap)
       call WriteDataFdtdOutgplObj(out,ncyc,mode)
       if ( mode ) call CloseOutgplObj(out)
    
    M4_FOREACH_OUTGPL2({case ("},{")
       if ( mode ) call OpenOutgplObj(out, ncyc, out%snap)
       call WriteData},{OutgplObj(out,ncyc,mode)
       if ( mode ) call CloseOutgplObj(out)
       })

    end select



  end subroutine WriteDataOutgplObj
  
!----------------------------------------------------------------------

end module outgpl

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
