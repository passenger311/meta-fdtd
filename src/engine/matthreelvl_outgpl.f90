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
      integer :: m, q
      type(T_REG) :: reg2


      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      mat = matthreelvlobj(out%objidx)
      M4_MODOBJ_GETREG(mat,reg2)

      sum1 = 0.
      sum2 = 0.
      sum3 = 0.

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      if ( reg2%mask(i,j,k) .ne. 0 ) then 

     q = reg2%mask(i,j,k)

      select case ( mode ) 
      case ( 1 )

         val1 = mat%rho11(q)
         val2 = mat%rho22(q)
         val3 = mat%rho33(q)

         if ( out%mode .ne. 'S' ) then
 
            if ( reg%isbox ) then
               write(out%funit,"(3E15.6E3)") real(val1,8), real(val2,8), real(val3,8)
            else
               write(out%funit,"(M4_SDIM({I5}),(3E15.6E3))") &
                    M4_DIM123({i},{i,j},{i,j,k}),real(val1,8),real(val2,8),real(val3,8)
            endif
         
         else

            sum1 = sum1 + val1
            sum2 = sum2 + val2
            sum3 = sum3 + val3

         endif

      case ( 2 )

         val1 = mat%rho12(q)
         val2 = mat%rho13(q)
         val3 = mat%rho23(q)

         if ( out%mode .ne. 'S' ) then

            if ( reg%isbox ) then
               write(out%funit,"(3E15.6E3)") real(val1,8), real(val2,8), real(val3,8)
            else
               write(out%funit,"(M4_SDIM({I5}),(3E15.6E3))") &
                    M4_DIM123({i},{i,j},{i,j,k}),real(val1,8),real(val2,8),real(val3,8)
            endif

         else

            sum1 = sum1 + val1
            sum2 = sum2 + val2
            sum3 = sum3 + val3

         end if

      end select      

      endif

      },{
      
      if ( reg%is .ne. reg%ie ) write(out%funit,*)}, {if ( reg%js .ne. reg%je ) write(out%funit,*)

      } )
   
      if ( out%mode .eq. 'S' ) then
         write(out%funit,"(3E15.6E3)") real(sum1,8), real(sum2,8), real(sum3,8)
      endif

    end subroutine WriteValues

  end subroutine WriteDataMatthreelvlOutgplObj


end module matthreelvl_outgpl


!
! Authors:  J.Hamm
! Modified: 28/05/2008
!
!======================================================================
