!-*- F90 -*------------------------------------------------------------
!
!  module: matthreelvl_outR / meta
!
!  this module handles R output of data related to the matthreelvl module.
!
!----------------------------------------------------------------------


module matthreelvl_outR

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

  public :: InitializeMatthreelvlOutRObj
  public :: FinalizeMatthreelvlOutRObj 
  public :: WriteDataMatthreelvlOutRObj

contains

!----------------------------------------------------------------------

  subroutine InitializeMatthreelvlOutRObj(out)

    type (T_OUT) :: out

  end subroutine InitializeMatthreelvlOutRObj


!----------------------------------------------------------------------

  subroutine FinalizeMatthreelvlOutRObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeMatthreelvlOutRObj


!----------------------------------------------------------------------


  subroutine WriteDataMatthreelvlOutRObj(out, mode)

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
      real(kind=8) :: val1, val2, val3, sum1, sum2, sum3, val4, val5, val6, sum4, sum5, sum6
      type(T_MATTHREELVL) :: mat
      integer :: m

      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      mat = matthreelvlobj(out%objidx)

      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      sum4 = 0.
      sum5 = 0.
      sum6 = 0.

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      select case ( mode ) 
      case ( 1 )

         val1 = mat%rho11(p)
         val2 = mat%rho22(p)
         val3 = mat%rho33(p)

         if ( out%mode .ne. 'S' ) then
               write(out%funit,"(M4_SDIM({I5}),(3E15.6E3))") &
                    M4_DIM123({i},{i,j},{i,j,k}),dble(val1),dble(val2),dble(val3)
         
         else
            sum1 = sum1 + val1
            sum2 = sum2 + val2
            sum3 = sum3 + val3
         endif

      case ( 2 )

         val1 = dble(mat%rho12(p))
         val2 = dble(mat%rho13(p))
         val3 = dble(mat%rho23(p))
         val4 = dimag(mat%rho12(p))
         val5 = dimag(mat%rho13(p))
         val6 = dimag(mat%rho23(p))

         if ( out%mode .ne. 'S' ) then
                       
            write(out%funit,"(M4_SDIM({I5}),(6E15.6E3))") &
                    M4_DIM123({i},{i,j},{i,j,k}),dble(val1),dble(val2),dble(val3)& 
                    ,dble(val4),dble(val5),dble(val6)

         else
            sum1 = sum1 + val1
            sum2 = sum2 + val2
            sum3 = sum3 + val3
            sum4 = sum4 + val4
            sum5 = sum5 + val5
            sum6 = sum6 + val6
         end if

      end select      

      },{
      
      if ( reg%is .ne. reg%ie ) write(out%funit,*)}, {if ( reg%js .ne. reg%je ) write(out%funit,*)

      } )
   
      if ( out%mode .eq. 'S' ) then
         write(out%funit,"(3E15.6E3)") dble(sum1), dble(sum2), dble(sum3), &
              dble(sum4), dble(sum5), dble(sum6)
      endif

    end subroutine WriteValues

  end subroutine WriteDataMatthreelvlOutRObj


end module matthreelvl_outR


!
! Authors:  J.Hamm, Andreas Pusch
! Modified: 24/09/2009
!
!======================================================================
