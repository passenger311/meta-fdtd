!-*- F90 -*------------------------------------------------------------
!
!  module: diagebal_outgpl / meta
!
!  this module handles GPL output of data related to the diagebal module.
!
!----------------------------------------------------------------------


module diagebal_outgpl

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use fdtd_calc
  use diagebal
  use out_calc

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'DIAGEBAL_OUTGPL'

 ! --- Public Methods

  public :: InitializeDiagebalOutgplObj
  public :: FinalizeDiagebalOutgplObj
  public :: WriteDataDiagebalOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeDiagebalOutgplObj(out)

    type (T_OUT) :: out
    
    out%numnodes = 1

  end subroutine InitializeDiagebalOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeDiagebalOutgplObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeDiagebalOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataDiagebalOutgplObj(out, mode)

    type (T_OUT) :: out
    logical :: mode

    if ( .not. mode ) return

    M4_WRITE_DBG({"write data ",TRIM(out%filename), " ",TRIM(out%fn)})

    select case (out%fn)
    case('En') ! differential
       call WriteValues(out, 1 )
    case('EnI') ! time integrated
       call WriteValues(out, 2 )
    case('DS') ! differential
       call WriteValues(out, 3 )
    case('DSI') ! time integrated
       call WriteValues(out, 4 )
    case('UR') ! position weighted energy desity 
       call WriteValues(out, 5 )
    case('URR') ! position weighted energy desity 
       call WriteValues(out, 6 )
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED" 
    end select
    
  contains

    ! **************************************************************** !

    subroutine WriteValues(out, mode)

      type (T_OUT) :: out
      integer :: mode
      type (T_DIAGEBAL) :: diag

      diag = diagebalobj(out%objidx)
    
      M4_WRITE_DBG({"WriteValues!"})
      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})


      select case ( mode ) 

      case ( 1 )
         write(out%funit,"(5E20.11E3)") diag%dudt, diag%ds, diag%jeabs, diag%khabs, diag%res
      case ( 2 )
         write(out%funit,"(5E20.11E3)") diag%sumdudt, diag%sumds, diag%sumje, diag%sumkh, diag%sumres
      case ( 3 ) 
         write(out%funit,"(4E20.11E3)") diag%ds, diag%dsx, diag%dsy, diag%dsz
      case ( 4 ) 
         write(out%funit,"(4E20.11E3)") diag%sumds, diag%sumdsx, diag%sumdsy, diag%sumdsz
      case ( 5 ) 
         write(out%funit,"(4E22.13E3)") diag%sumdudtr(1)/diag%sumdudt, &
              diag%sumdudtr(2)/diag%sumdudt, diag%sumdudtr(3)/diag%sumdudt
      case ( 6 ) 
         write(out%funit,"(4E20.11E3)") diag%sumdudtrr(1)/diag%sumdudt, &
              diag%sumdudtrr(2)/diag%sumdudt, diag%sumdudtrr(3)/diag%sumdudt
      end select

    end subroutine WriteValues

  end subroutine WriteDataDiagebalOutgplObj


end module diagebal_outgpl

!
! Authors:  J.Hamm
! Modified: 1/02/2008
!
!======================================================================
