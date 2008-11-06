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
    logical :: mode

    select case (out%fn)
    case('EH')
       call WriteEH(out,mode)
    case('NRef')
       call WriteN(out,mode)
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED" 
    end select

  end subroutine WriteDataFdtdOutsetObj

  subroutine WriteEH(out,mode)
 
    type (T_OUT) :: out
    logical :: mode

    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  

    M4_WRITE_DBG({"write data ",TRIM(out%filename), " ",TRIM(out%fn)})

    if ( .not. mode ) return
    
    M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})
    
    reg = regobj(out%regidx)
    
     M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
    
    write(out%funit,"(6E15.6E3)") Ex(i,j,k),Ey(i,j,k),Ez(i,j,k),Hx(i,j,k),Hy(i,j,k),Hz(i,j,k)
    
    },{}, {} )
    
    write(out%funit,"(A)") ")SET"

  end subroutine WriteEH

  subroutine WriteN(out,mode)
 
    type (T_OUT) :: out
    logical :: mode
    real(kind=8) :: eps

    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))  

    M4_WRITE_DBG({"write data ",TRIM(out%filename), " ",TRIM(out%fn)})

    if ( .not. mode ) return
    
    M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})
    
    reg = regobj(out%regidx)
    
     M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
    
     eps = 1./6. * ( 1./epsinvx(i,j,k) + 1./epsinvx(i-1,j,k) + &
          1./epsinvy(i,j,k) + 1./epsinvy(i,j-1,k) + &
          1./epsinvz(i,j,k) + 1./epsinvz(i,j,k-1) )

    write(out%funit,"(E15.6E3)") sqrt(eps)
    
    },{}, {} )
    
    write(out%funit,"(A)") ")SET"

  end subroutine WriteN


end module fdtd_outset

!
! Authors:  J.Hamm

! Modified: 4/12/2007
!
!======================================================================
