!-*- F90 -*------------------------------------------------------------
!
!  module: fdtd_outset / meta
!
!  this module handles SET output of data related to the fdtd module.
!
!----------------------------------------------------------------------


module fdtd_outset

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use fdtd_calc
  use out_calc

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'FDTD_OUTSET'

 ! --- Public Methods

  public :: InitializeFdtdOutsetObj
  public :: FinalizeFdtdOutsetObj
  public :: WriteDataFdtdOutsetObj

contains

!----------------------------------------------------------------------

  subroutine InitializeFdtdOutsetObj(out)

    type (T_OUT) :: out
    type (T_BUF) :: buf

  end subroutine InitializeFdtdOutsetObj


!----------------------------------------------------------------------

  subroutine FinalizeFdtdOutsetObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeFdtdOutsetObj


!----------------------------------------------------------------------


  subroutine WriteDataFdtdOutsetObj(out, mode)

    type (T_OUT) :: out
    type (T_BUF) :: buf
    logical :: mode

    M4_WRITE_DBG({"write data ",TRIM(out%filename), " ",TRIM(out%fn)})

    select case (out%fn)
    case('EH')
       call WriteEH(out, mode)
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED" 
    end select
    
  contains

    ! **************************************************************** !

   subroutine WriteEH(out, mode)

      type (T_OUT) :: out
      logical :: mode
      M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  

      if ( .not. mode ) return

      M4_WRITE_DBG({"WriteEH!"})
      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      sx = 0.0
      sy = 0.0
      sz = 0.0

      reg = regobj(out%regidx)
      
      write(out%funit,"(A)") "(SET"

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      write(out%funit,"(6E15.6E3)") Ex(i,j,k),Ey(i,j,k),Ez(i,j,k),Hx(i,j,k),Hy(i,j,k),Hz(i,j,k)
      
      },{}, {} )

     write(out%funit,"(A)") ")SET"
   
   end subroutine WriteEH


    ! *************************************************************** !

  end subroutine WriteDataFdtdOutsetObj


end module fdtd_outset

!
! Authors:  J.Hamm

! Modified: 4/12/2007
!
!======================================================================
