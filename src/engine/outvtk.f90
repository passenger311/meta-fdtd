!----------------------------------------------------------------------
!
!  module: outvtk
!
!  ascii (gnuplot) output module
!
!  subs:
!
!  InitializeOutvtkObj
!  FinalizeOutvtkObj
!  OpenOutvtkObj
!  CloseOutvtkObj
!  WriteHeaderOutvtkObj
!  WriteDataOutvtkObj
!  
!----------------------------------------------------------------------


!======================================================================
!
!

module outvtk

  use constant
  use strings
  use mpiworld
  use reglist
  use outlist
  use grid

! outvtk modules 
  use fdtd_outvtk
M4_FOREACH_OUTVTK({use }, {
  })

  implicit none
  private
  save
  
  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'OUTVTK'

  ! --- Public Methods

  public :: InitializeOutvtkObj
  public :: FinalizeOutvtkObj
  public :: OpenOutvtkObj
  public :: CloseOutvtkObj
  public :: WriteHeaderOutvtkObj
  public :: WriteDataOutvtkObj

  ! --- Public Data

  ! --- Constants

  character(len=STRLNG), parameter :: vtusfx = '.vtu' ! unstructured grid
  character(len=STRLNG), parameter :: vtrsfx = '.vtr' ! rectilinear grid

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializeOutvtkObj(out)

    type(T_OUT) :: out


! override out%snap as vtk does not support multiple frames in one file

    out%snap = .true.

! call various finalization methods
    select case ( out%modl ) 

    case ("FDTD")
       call  InitializeFdtdOutvtkObj(out)

    M4_FOREACH_OUTVTK2({case ("},{")
       call Initialize},{OutvtkObj(UNITTMP)
             })
    end select

    if ( .not. out%snap ) then
       call OpenOutvtkObj(out, 0, .true.)
       call CloseOutvtkObj(out)
    endif


! if not in snapshot mode: write a header once at initialization 

  end subroutine InitializeOutvtkObj


!----------------------------------------------------------------------

  subroutine FinalizeOutvtkObj(out)

    type(T_OUT) :: out

    select case ( out%modl ) 

! call various finalization methods
    case ("FDTD")
       call FinalizeFdtdOutvtkObj(out)
       M4_FOREACH_OUTVTK2({case ("},{")
             call Finalize},{OutvtkObj
             })
    end select

  end subroutine FinalizeOutvtkObj


!----------------------------------------------------------------------

  subroutine OpenOutvtkObj(out, ncyc, writehdr)

    type(T_OUT) :: out
    integer :: ncyc
    logical :: writehdr
    character(len=STRLNG) :: name, stepstr, outvtksfx
 
    if ( regobj(out%regidx)%isbox ) then
       outvtksfx = ".vtr"
    else
       outvtksfx = ".vtu"
    endif

    if ( out%snap ) then
       stepstr = TRIM(i2str(ncyc))
       name = cat5(out%filename,"_",TRIM(i2str(ncyc)),mpi_sfx,outvtksfx)
    else
       name = cat3(out%filename,mpi_sfx,outvtksfx)
    endif
   
    if ( writehdr ) then 
       
       M4_WRITE_DBG({"opening ",TRIM(out%filename),", as ",TRIM(name),", unit ", &
            TRIM(i2str(out%funit)), " (new)"})
       open(out%funit,FILE=name,STATUS="UNKNOWN",POSITION="REWIND")
       call WriteHeaderOutvtkObj(out,ncyc)

    else

       M4_WRITE_DBG({"opening ",TRIM(out%filename),", as ",TRIM(name),", unit ", &
            TRIM(i2str(out%funit)), " (append)"})
       open(out%funit,FILE=name,STATUS="UNKNOWN",POSITION="APPEND")

    endif

  end subroutine OpenOutvtkObj

!----------------------------------------------------------------------

  subroutine CloseOutvtkObj(out)

    type(T_OUT) :: out

    M4_WRITE_DBG({"closing ",TRIM(out%filename),", unit ", TRIM(i2str(out%funit))})

    close(out%funit)

  end subroutine CloseOutvtkObj

!----------------------------------------------------------------------

  subroutine WriteHeaderOutvtkObj(out,ncyc)

    type(T_OUT) :: out
    integer :: ncyc
    character(len=STRLNG) :: dataset
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    reg = regobj(out%regidx)
    
    if ( reg%isbox ) then
       dataset = "STRUCTURED_GRID"
    else
       dataset = "UNSTRUCTURED_GRID"
    endif

    M4_WRITE_DBG({"write header ",TRIM(out%filename)})
    
    write(out%funit,*) '# vtk DataFile Version 2.0'
    write(out%funit,*) "META: M4_VERSION(), M4_FLAVOUR()"
    write(out%funit,*) 'ASCII'
    if ( out%mode .eq. 'S' ) then
       write(out%funit,*) 'FIELD Sum 1'
       return
    endif
    if ( reg%isbox ) then
       write(out%funit,*) 'DATASET RECTILINEAR_GRID'
       write(out%funit,*) 'DIMENSIONS ', (reg%ie-reg%is)/reg%di+1,(reg%je-reg%js)/reg%dj+1,(reg%ke-reg%ks)/reg%dk+1
       write(out%funit,*) 'X_COORDINATES float'
       write(out%funit,"(E15.6E3)") (SX*i, i=IBEG,IEND,1)
       write(out%funit,*) 'Y_COORDINATES float'
       write(out%funit,"(E15.6E3)") (SY*i, i=JBEG,JEND,1)
       write(out%funit,*) 'Z_COORDINATES float'
       write(out%funit,"(E15.6E3)") (SZ*i, i=KBEG,KEND,1)
    else
       write(out%funit,*) 'DATASET UNSTRUCTURED_GRID'
       write(out%funit,*) 'POINTS ',reg%numnodes,' float'
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       write(out%funit,"(M4_SDIM({I5}))") M4_DIM123({i},{i,j},{i,j,k})
       })
       write(out%funit,*) 'CELLS ',reg%numnodes,2*reg%numnodes
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       write(out%funit,*) "1 ",p
       })
       write(out%funit,*) 'CELLS_TYPES ',reg%numnodes
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       write(out%funit,*) "1"
       })
       write(out%funit,*) 'POINT_DATA ',reg%numnodes
    endif

  end subroutine WriteHeaderOutvtkObj


!----------------------------------------------------------------------

  subroutine WriteDataOutvtkObj(out, ncyc, mode)

    type(T_OUT) :: out
    integer :: ncyc
    logical :: mode  
    
    M4_WRITE_DBG({"write data ",TRIM(out%filename),", numnodes ", TRIM(i2str(out%numnodes))})


    select case ( out%modl ) 

! call various output methods

    case ("FDTD")

       if ( mode ) then 
          call OpenOutvtkObj(out, ncyc, out%snap)
          if ( out%numnodes .gt. 0 ) write(out%funit,*)
       end if
       call WriteDataFdtdOutvtkObj(out,mode)
       if ( mode ) call CloseOutvtkObj(out)
    
    M4_FOREACH_OUTVTK2({case ("},{")
       if ( mode ) then 
          call OpenOutvtkObj(out, ncyc, out%snap)
          if ( out%numnodes .gt. 0 ) write(out%funit,*) 
       endif
       call WriteData},{OutvtkObj(out,mode)
       if ( mode ) call CloseOutvtkObj(out)
       })

    case default
       if ( mode ) then 
          call OpenOutvtkObj(out, ncyc, out%snap)
          write(out%funit,*) "OUTPUT MODULE NOT IMPLEMENTED" 
       end if
       if ( mode ) call CloseOutvtkObj(out)

    end select



  end subroutine WriteDataOutvtkObj
  
!----------------------------------------------------------------------

end module outvtk

!
! Authors:  J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
