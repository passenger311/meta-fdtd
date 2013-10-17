!-*- F90 -*------------------------------------------------------------
!
!  module: diagebalmode_outvtk / meta
!
!  this module handles VTK output of data related to the diagebalmode module.
!
!----------------------------------------------------------------------


module diagebalmode_outvtk

  use constant
  use strings
  use reglist
  use fdtd
  use outlist
  use buflist
  use mpiworld
  use grid 
  use diagebalmode
  use out_calc

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'DIAGEBALMODE_OUTVTK'

 ! --- Public Methods

  public :: InitializeDiagebalmodeOutvtkObj
  public :: FinalizeDiagebalmodeOutvtkObj
  public :: WriteDataDiagebalmodeOutvtkObj

contains

!----------------------------------------------------------------------

  subroutine InitializeDiagebalmodeOutvtkObj(out)

    type (T_OUT) :: out

  end subroutine InitializeDiagebalmodeOutvtkObj


!----------------------------------------------------------------------

  subroutine FinalizeDiagebalmodeOutvtkObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeDiagebalmodeOutvtkObj


!----------------------------------------------------------------------

  subroutine WriteDataDiagebalmodeOutvtkObj(out, mode)

    type (T_OUT) :: out
    type (T_BUF) :: buf
    logical :: mode

    if ( .not. mode ) return

    M4_WRITE_DBG({"write data ",TRIM(out%filename), " ",TRIM(out%fn)})

    select case (out%fn)
    case('En')
       call WriteValuesDiagebalmode(out, 1)
    case('E')
       call WriteValuesDiagebalmode(out, 2)
    case('H')
       call WriteValuesDiagebalmode(out, 3)
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED"
    end select

  contains

    ! **************************************************************** !

    subroutine WriteValuesDiagebalmode(out, mode)

      type (T_OUT) :: out
      integer :: mode
      real(kind=8) :: en
      M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
      type(T_DIAGEBALMODE) :: diag

      M4_WRITE_DBG({"WriteScalars!"})

      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      diag = diagebalmodeobj(out%objidx)

      select case (mode)
      case(1)
         write(out%funit,"(A,A)") "SCALARS scalars float 1"
         write(out%funit,*) "LOOKUP_TABLE default"
      case (2)
         write(out%funit,"(A)") "VECTORS vectors float"
      case (3)
         write(out%funit,"(A)") "VECTORS vectors float"
      end select

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
      select case (mode)
      case(1)
         en = ( M4_VOLEX(i,j,k) / epsinvx(i,j,k) * &
              (dble(diag%buf_E(diag%h_pos_E,i,j,k,0))**2+&
              dimag(diag%buf_E(diag%h_pos_E,i,j,k,0))**2) + &
              M4_VOLEY(i,j,k) / epsinvy(i,j,k) * &
              (dble(diag%buf_E(diag%h_pos_E,i,j,k,1))**2+&
              dimag(diag%buf_E(diag%h_pos_E,i,j,k,1))**2) + &
              M4_VOLEZ(i,j,k) / epsinvz(i,j,k) * &
              (dble(diag%buf_E(diag%h_pos_E,i,j,k,2))**2+&
              dimag(diag%buf_E(diag%h_pos_E,i,j,k,2))**2) + &
              M4_VOLHX(i,j,k) / M4_MUINVX(i,j,k) * &
              (dble(diag%buf_H(diag%h_pos_H,i,j,k,0))**2+&
              dimag(diag%buf_H(diag%h_pos_H,i,j,k,0))**2) + &
              M4_VOLHY(i,j,k) / M4_MUINVY(i,j,k) * &
              (dble(diag%buf_H(diag%h_pos_H,i,j,k,1))**2+&
              dimag(diag%buf_H(diag%h_pos_H,i,j,k,1))**2) + &
              M4_VOLHZ(i,j,k) / M4_MUINVZ(i,j,k) * &
              (dble(diag%buf_H(diag%h_pos_H,i,j,k,2))**2+&
              dimag(diag%buf_H(diag%h_pos_H,i,j,k,2))**2) )
         write(out%funit,"(4E15.6E3)") real(en)
      case (2)
         write(out%funit,"(3E14.6E2)") abs(diag%buf_E(diag%h_pos_E,i,j,k,0)), &
              abs(diag%buf_E(diag%h_pos_E,i,j,k,1)), &
              abs(diag%buf_E(diag%h_pos_E,i,j,k,2))
      case (3)
         write(out%funit,"(3E14.6E2)") abs(diag%buf_H(diag%h_pos_H,i,j,k,0)), &
              abs(diag%buf_H(diag%h_pos_H,i,j,k,1)), &
              abs(diag%buf_H(diag%h_pos_H,i,j,k,2))
      end select
      } )

    end subroutine WriteValuesDiagebalmode

  end subroutine WriteDataDiagebalmodeOutvtkObj


end module diagebalmode_outvtk

!
! Authors:  J.Hamm, F.Renn, A. Pusch

! Modified: 08/08/2011
!
!======================================================================
