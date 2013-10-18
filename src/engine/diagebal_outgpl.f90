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
         write(out%funit,"(5E15.6E3)") dble(diag%dudt), dble(diag%ds), dble(diag%jeabs), &
                dble(diag%khabs), dble(diag%res)
      case ( 2 )
         write(out%funit,"(5E15.6E3)") dble(diag%sumdudt), dble(diag%sumds), &
                dble(diag%sumje), dble(diag%sumkh), dble(diag%sumres)
      case ( 3 ) 
         write(out%funit,"(4E15.6E3)") dble(diag%ds), dble(diag%dsx), dble(diag%dsy), &
                dble(diag%dsz)
      case ( 4 ) 
         write(out%funit,"(4E15.6E3)") dble(diag%sumds), dble(diag%sumdsx), &
                dble(diag%sumdsy), dble(diag%sumdsz)
      end select

    end subroutine WriteValues

  end subroutine WriteDataDiagebalOutgplObj


end module diagebal_outgpl

!
! Authors:  J.Hamm
! Modified: 1/02/2008
!
!======================================================================
