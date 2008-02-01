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
    case('D') ! differential
       call WriteValues(out, .false.)
    case('I') ! time integrated
       call WriteValues(out,.true.)
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED" 
    end select
    
  contains

    ! **************************************************************** !

    subroutine WriteValues(out, sum)

      type (T_OUT) :: out
      logical :: sum
      type (T_DIAGEBAL) :: diag

      diag = diagebalobj(out%objidx)
    
      M4_WRITE_DBG({"WriteValues!"})
      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      if ( sum ) then

         write(out%funit,"(4E15.6E3)") diag%sumdudt, diag%sumdivs, diag%sumjekh, diag%sumres
         
      else

         write(out%funit,"(4E15.6E3)") diag%dudt, diag%divs, diag%jekh, diag%res

      end if

    end subroutine WriteValues

  end subroutine WriteDataDiagebalOutgplObj


end module diagebal_outgpl

!
! Authors:  J.Hamm
! Modified: 1/02/2008
!
!======================================================================
