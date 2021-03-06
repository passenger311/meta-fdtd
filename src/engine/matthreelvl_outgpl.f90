!-*- F90 -*------------------------------------------------------------
!
!  module: matthreelvl_outgpl / meta
!
!  this module handles GPL output of data related to the matthreelvl module.
!
!----------------------------------------------------------------------


module matthreelvl_outgpl

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use out_calc
  use matthreelvl

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MATTHREELVL_OUTGPL'

 ! --- Public Methods

  public :: InitializeMatthreelvlOutgplObj
  public :: FinalizeMatthreelvlOutgplObj 
  public :: WriteDataMatthreelvlOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeMatthreelvlOutgplObj(out)

    type (T_OUT) :: out

  end subroutine InitializeMatthreelvlOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeMatthreelvlOutgplObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeMatthreelvlOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataMatthreelvlOutgplObj(out, mode)

    type (T_OUT) :: out
    type (T_BUF) :: buf
    logical :: mode

    if ( .not. mode ) return

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
      real(kind=8) :: val1, val2, val3, sum1, sum2, sum3
      type(T_MATTHREELVL) :: mat
      integer :: m

      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      mat = matthreelvlobj(out%objidx)

      sum1 = 0.
      sum2 = 0.
      sum3 = 0.

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      select case ( mode ) 
      case ( 1 )

         val1 = mat%rho11(p)
         val2 = mat%rho22(p)
         val3 = mat%rho33(p)

         if ( out%mode .ne. 'S' ) then
 
            if ( reg%isbox ) then
               write(out%funit,"(3E15.6E3)") dble(val1), dble(val2), dble(val3)
            else
               write(out%funit,"(M4_SDIM({I5}),(3E15.6E3))") &
                    M4_DIM123({i},{i,j},{i,j,k}),dble(val1), dble(val2), dble(val3)
            endif
         
         else
            sum1 = sum1 + val1
            sum2 = sum2 + val2
            sum3 = sum3 + val3
         endif

      case ( 2 )

         val1 = mat%rho12(p)
         val2 = mat%rho13(p)
         val3 = mat%rho23(p)

         if ( out%mode .ne. 'S' ) then

            if ( reg%isbox ) then
               write(out%funit,"(3E15.6E3)") dble(val1), dble(val2), dble(val3)
            else
               write(out%funit,"(M4_SDIM({I5}),(3E15.6E3))") &
                    M4_DIM123({i},{i,j},{i,j,k}),dble(val1), dble(val2), dble(val3)
            endif

         else
            sum1 = sum1 + val1
            sum2 = sum2 + val2
            sum3 = sum3 + val3
         end if

      end select      

      },{
      
      if ( reg%is .ne. reg%ie ) write(out%funit,*)}, {if ( reg%js .ne. reg%je ) write(out%funit,*)

      } )
   
      if ( out%mode .eq. 'S' ) then
         write(out%funit,"(3E15.6E3)") dble(sum1), dble(sum2), dble(sum3)
      endif

    end subroutine WriteValues

  end subroutine WriteDataMatthreelvlOutgplObj


end module matthreelvl_outgpl


!
! Authors:  J.Hamm
! Modified: 28/05/2008
!
!======================================================================
