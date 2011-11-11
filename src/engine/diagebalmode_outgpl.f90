!-*- F90 -*------------------------------------------------------------
!
!  module: diagebalmode_outgpl / meta
!
!  this module handles GPL output of data related to the diagebalmode
!  module.
!
!----------------------------------------------------------------------


module diagebalmode_outgpl

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use fdtd_calc
  use diagebalmode
  use out_calc
  use matlorentz

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'DIAGEBALMODE_OUTGPL'

 ! --- Public Methods

  public :: InitializeDiagebalModeOutgplObj
  public :: FinalizeDiagebalModeOutgplObj
  public :: WriteDataDiagebalModeOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeDiagebalModeOutgplObj(out)

    type (T_OUT) :: out
    
    out%numnodes = 1    

  end subroutine InitializeDiagebalModeOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeDiagebalModeOutgplObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeDiagebalModeOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataDiagebalModeOutgplObj(out, mode)

    type (T_OUT) :: out
    logical :: mode

    if ( .not. mode ) return

    select case (out%fn)
    case('En') ! differential
       call WriteValuesDiagebalMode(out, 1 )
    case('EnI') ! time integrated
       call WriteValuesDiagebalMode(out, 2 )
    case('DS') ! differential
       call WriteValuesDiagebalMode(out, 3 )
    case('DSI') ! time integrated
       call WriteValuesDiagebalMode(out, 4 )
    case('E')
       call WriteSumDiagebalMode(out, 1 )   
    case('H')
       call WriteSumDiagebalMode(out, 2 )   
    case('J')
       call WriteSumDiagebalMode(out, 3 )
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED" 
    end select
    
  end subroutine WriteDataDiagebalModeOutgplObj

  subroutine WriteSumDiagebalMode(out, mode)
    
    type (T_OUT) :: out
    integer :: mode, c
    type (T_DIAGEBALMODE) :: diag
    real(kind=8) fld_sum(0:2)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    fld_sum = 0.
    
    diag = diagebalmodeobj(out%objidx)
    M4_MODOBJ_GETREG(out,reg)

    M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
    do c=0,2
       select case ( mode ) 
       case (1)
          fld_sum(c) = fld_sum(c) + abs(diag%buf_E(diag%h_pos_E,i,j,k,c))
       case (2)
          fld_sum(c) = fld_sum(c) + abs(diag%buf_H(diag%h_pos_H,i,j,k,c))
       case (3)
          fld_sum(c) = fld_sum(c) + abs(diag%buf_J(diag%h_pos_J,i,j,k,c,0))
       end select
    end do
    })

    write(out%funit,"(5E15.6E3)") fld_sum(0), fld_sum(1), fld_sum(2)
    
  end subroutine WriteSumDiagebalMode

  subroutine WriteValuesDiagebalMode(out, mode)
    
    type (T_OUT) :: out
    integer :: mode
    type (T_DIAGEBALMODE) :: diag
    
    diag = diagebalmodeobj(out%objidx)
    
    M4_WRITE_DBG({"WriteValues!"})
    M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})
    
    select case ( mode ) 
       
    case ( 1 )
       write(out%funit,"(5E15.6E3)") diag%dudt, diag%ds, diag%je, diag%res
    case ( 2 )
       write(out%funit,"(5E15.6E3)") diag%sumdudt, diag%sumds,  diag%sumje, diag%sumres
    case ( 3 ) 
       write(out%funit,"(4E15.6E3)") diag%ds, diag%dsx, diag%dsy, diag%dsz
    case ( 4 ) 
       write(out%funit,"(4E15.6E3)") diag%sumds, diag%sumdsx, diag%sumdsy, diag%sumdsz
    end select

  end subroutine WriteValuesDiagebalMode
        
end module diagebalmode_outgpl

!
! Authors:  F.Renn, J.Hamm
! Modified: 1/08/2011
!
!======================================================================
