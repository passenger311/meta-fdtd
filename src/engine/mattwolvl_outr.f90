!-*- F90 -*------------------------------------------------------------
!
!  module: mattwolvl_outR / meta
!
!  this module handles R output of data related to the mattwolvl module.
!
!----------------------------------------------------------------------


module mattwolvl_outR

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use out_calc
  use mattwolvl

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MATTWOLVL_OUTGPL'

 ! --- Public Methods

  public :: InitializeMattwolvlOutRObj
  public :: FinalizeMattwolvlOutRObj 
  public :: WriteDataMattwolvlOutRObj

contains

!----------------------------------------------------------------------

  subroutine InitializeMattwolvlOutRObj(out)

    type (T_OUT) :: out

  end subroutine InitializeMattwolvlOutRObj


!----------------------------------------------------------------------

  subroutine FinalizeMattwolvlOutRObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeMattwolvlOutRObj


!----------------------------------------------------------------------


  subroutine WriteDataMattwolvlOutRObj(out, mode)

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
      type(T_MATTWOLVL) :: mat
      integer :: m

      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      mat = mattwolvlobj(out%objidx)

      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      sum4 = 0.
      sum5 = 0.
      sum6 = 0.

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      select case ( mode ) 
      case ( 1 )

         val1 = mat%inversion(p)

         if ( out%mode .ne. 'S' ) then
               write(out%funit,"(M4_SDIM({I5}),(1E15.6E3))") &
                    M4_DIM123({i},{i,j},{i,j,k}),dble(val1)
         
         else
            sum1 = sum1 + val1
         endif

      case ( 2 )

         val1 = dble(mat%rho12(p))

         if ( out%mode .ne. 'S' ) then
                       
            write(out%funit,"(M4_SDIM({I5}),(2E15.6E3))") &
                    M4_DIM123({i},{i,j},{i,j,k}),dble(val1) 

         else
            sum1 = sum1 * val1
         end if

      end select      

      },{
      
      if ( reg%is .ne. reg%ie ) write(out%funit,*)}, {if ( reg%js .ne. reg%je ) write(out%funit,*)

      } )
   
      if ( out%mode .eq. 'S' ) then
         write(out%funit,"(2E15.6E3)") dble(sum1)
      endif

    end subroutine WriteValues

  end subroutine WriteDataMattwolvlOutRObj


end module mattwolvl_outR


!
! Authors:  J.Hamm, Andreas Pusch
! Modified: 13/08/2010
!
!======================================================================
