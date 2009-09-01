!-*- F90 -*------------------------------------------------------------
!
!  module: mat3lvl_outgpl / meta
!
!  this module handles GPL output of data related to the mat3lvl module.
!
!----------------------------------------------------------------------


module mat3lvl_outgpl

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use out_calc
  use mat3lvl

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MAT3LVL_OUTGPL'

 ! --- Public Methods

  public :: InitializeMat3lvlOutgplObj
  public :: FinalizeMat3lvlOutgplObj
  public :: WriteDataMat3lvlOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeMat3lvlOutgplObj(out)

    type (T_OUT) :: out

  end subroutine InitializeMat3lvlOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeMat3lvlOutgplObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeMat3lvlOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataMat3lvlOutgplObj(out, mode)

    type (T_OUT) :: out
    type (T_BUF) :: buf
    logical :: mode

    M4_WRITE_DBG({"write data ",TRIM(out%filename), " ",TRIM(out%fn)})

    select case (out%fn)
    case('N')
       call WriteValues(out, 1)
    case('P')
       call WriteValues(out, 2)
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED" 
    end select
    
  contains

    ! **************************************************************** !

    subroutine WriteValues(out, mode)

      type (T_OUT) :: out
      integer :: mode

      M4_REGLOOP_DECL(reg,p,i,j,k,w(3))  
      real(kind=8) :: val, sum
      type(T_MAT3LVL) :: mat
      integer :: m

      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      mat = mat3lvlobj(out%objidx)

      sum = 0.

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      select case ( mode ) 
      case ( 1 )

         val = mat%N(p)

         if ( out%mode .ne. 'S' ) then
 
            if ( reg%isbox ) then
               write(out%funit,"(E15.6E3)") real(val,8)
            else
               write(out%funit,"(M4_SDIM({I5}),(E15.6E3))") M4_DIM123({i},{i,j},{i,j,k}),real(val,8)
            endif
         
         else
            sum = sum + val
         endif

      case ( 2 )

         val = mat%P(p,mat%cyc)

         if ( out%mode .ne. 'S' ) then

            if ( reg%isbox ) then
               write(out%funit,"(3E15.6E3)") real(val,8)
            else
               write(out%funit,"(M4_SDIM({I5}),(3E15.6E3))") M4_DIM123({i},{i,j},{i,j,k}),real(val,8)
            endif

         else
            sum = sum + val
         end if

      end select      

      },{if ( reg%is .ne. reg%ie ) write(out%funit,*)}, {if ( reg%js .ne. reg%je ) write(out%funit,*)} )
   
      if ( out%mode .eq. 'S' ) then
         write(out%funit,"(E15.6E3)") real(sum,8)
      endif

    end subroutine WriteValues

  end subroutine WriteDataMat3lvlOutgplObj


end module mat3lvl_outgpl


!
! Authors:  J.Hamm
! Modified: 28/05/2008
!
!======================================================================
