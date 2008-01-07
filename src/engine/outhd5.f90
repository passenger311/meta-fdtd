!----------------------------------------------------------------------
!
!  module: outhd5
!
!  ascii (gnuplot) output module
!
!  subs:
!
!  InitializeOuthd5Obj
!  FinalizeOuthd5Obj
!  OpenOuthd5Obj
!  CloseOuthd5Obj
!  WriteHeaderOuthd5Obj
!  WriteDataOuthd5Obj
!  
!----------------------------------------------------------------------


!======================================================================
!
!

module outhd5

  use constant
  use strings
  use mpiworld
  use reglist
  use outlist

! outhd5 modules 
  use fdtd_outhd5
M4_FOREACH_OUTHD5({use }, {
  })

  implicit none
  private
  save
  
  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'OUTHD5'

  ! --- Public Methods

  public :: InitializeOuthd5Obj
  public :: FinalizeOuthd5Obj
  public :: OpenOuthd5Obj
  public :: CloseOuthd5Obj
  public :: WriteHeaderOuthd5Obj
  public :: WriteDataOuthd5Obj

  ! --- Public Data

  ! --- Constants

  character(len=STRLNG), parameter :: outhd5sfx = '.hd5'

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializeOuthd5Obj(out)

    type(T_OUT) :: out

! call various finalization methods
    select case ( out%modl ) 

    case ("FDTD")
       call  InitializeFdtdOuthd5Obj(out)

    M4_FOREACH_OUTHD52({case ("},{")
       call Initialize},{Outhd5Obj(UNITTMP)
             })
    end select

    if ( .not. out%snap ) then
       call OpenOuthd5Obj(out, 0, .true.)
       call CloseOuthd5Obj(out)
    endif


! if not in snapshot mode: write a header once at initialization 

  end subroutine InitializeOuthd5Obj


!----------------------------------------------------------------------

  subroutine FinalizeOuthd5Obj(out)

    type(T_OUT) :: out

    select case ( out%modl ) 

! call various finalization methods
    case ("FDTD")
       call FinalizeFdtdOuthd5Obj(out)
       M4_FOREACH_OUTHD52({case ("},{")
             call Finalize},{Outhd5Obj
             })
    end select

  end subroutine FinalizeOuthd5Obj


!----------------------------------------------------------------------

  subroutine OpenOuthd5Obj(out, ncyc, writehdr)

    type(T_OUT) :: out
    integer :: ncyc
    logical :: writehdr
    character(len=STRLNG) :: name, stepstr
 
    if ( out%snap ) then
       stepstr = TRIM(i2str(ncyc))
       name = cat5(out%filename,"_",TRIM(i2str(ncyc)),mpi_sfx,outhd5sfx)
    else
       name = cat3(out%filename,mpi_sfx,outhd5sfx)
    endif
   
    if ( writehdr ) then 
       
       M4_WRITE_DBG({"opening ",TRIM(out%filename),", as ",TRIM(name),", unit ", &
            TRIM(i2str(out%funit)), " (new)"})
       open(out%funit,FILE=name,STATUS="UNKNOWN",POSITION="REWIND")
       call WriteHeaderOuthd5Obj(out,ncyc)

    else

       M4_WRITE_DBG({"opening ",TRIM(out%filename),", as ",TRIM(name),", unit ", &
            TRIM(i2str(out%funit)), " (append)"})
       open(out%funit,FILE=name,STATUS="UNKNOWN",POSITION="APPEND")

    endif

  end subroutine OpenOuthd5Obj

!----------------------------------------------------------------------

  subroutine CloseOuthd5Obj(out)

    type(T_OUT) :: out

    M4_WRITE_DBG({"closing ",TRIM(out%filename),", unit ", TRIM(i2str(out%funit))})

    close(out%funit)

  end subroutine CloseOuthd5Obj

!----------------------------------------------------------------------

  subroutine WriteHeaderOuthd5Obj(out,ncyc)

    type(T_OUT) :: out
    integer :: ncyc
    type(T_REG) :: reg

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
       write(out%funit,*) '# ',out%numsteps         ! number of time points
       write(out%funit,*) '# ',out%ns,out%ne,out%dn ! time frame
    endif
    write(out%funit,*) '# ',reg%isbox               ! mode
    write(out%funit,*) '# ',reg%numnodes            ! number of points
    write(out%funit,*) '# ',reg%is,reg%ie,reg%di    ! space box
    write(out%funit,*) '# ',reg%js,reg%je,reg%dj
    write(out%funit,*) '# ',reg%ks,reg%ke,reg%dk

  end subroutine WriteHeaderOuthd5Obj


!----------------------------------------------------------------------

  subroutine WriteDataOuthd5Obj(out, ncyc, mode)

    type(T_OUT) :: out
    integer :: ncyc
    logical :: mode  
    
    M4_WRITE_DBG({"write data ",TRIM(out%filename),", numnodes ", TRIM(i2str(out%numnodes))})


    select case ( out%modl ) 

! call various output methods

    case ("FDTD")

       if ( mode ) then 
          call OpenOuthd5Obj(out, ncyc, out%snap)
          if ( out%numnodes .gt. 0 ) write(out%funit,*)
       end if
       call WriteDataFdtdOuthd5Obj(out,mode)
       if ( mode ) call CloseOuthd5Obj(out)
    
    M4_FOREACH_OUTHD52({case ("},{")
       if ( mode ) then 
          call OpenOuthd5Obj(out, ncyc, out%snap)
          if ( out%numnodes .gt. 0 ) write(out%funit,*) 
       endif
       call WriteData},{Outhd5Obj(out,mode)
       if ( mode ) call CloseOuthd5Obj(out)
       })

    case default
       if ( mode ) then 
          call OpenOuthd5Obj(out, ncyc, out%snap)
          write(out%funit,*) "OUTPUT MODULE NOT IMPLEMENTED" 
       end if
       if ( mode ) call CloseOuthd5Obj(out)

    end select



  end subroutine WriteDataOuthd5Obj
  
!----------------------------------------------------------------------

end module outhd5

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
