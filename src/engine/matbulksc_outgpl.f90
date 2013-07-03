!-*- F90 -*------------------------------------------------------------
!
!  module: matbulksc_outgpl / meta
!
!  this module handles GPL output of data related to the matbulksc module.
!
!----------------------------------------------------------------------


module matbulksc_outgpl

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use out_calc
  use matbulksc

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MATBULKSC_OUTGPL'

 ! --- Public Methods

  public :: InitializeMatBulkscOutgplObj
  public :: FinalizeMatBulkscOutgplObj
  public :: WriteDataMatBulkscOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeMatBulkscOutgplObj(out)

    type (T_OUT) :: out

  end subroutine InitializeMatBulkscOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeMatBulkscOutgplObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeMatBulkscOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataMatBulkscOutgplObj(out, mode)

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
      real(kind=8) :: val1, val2, val3, sum(3)
      type(T_MATBULKSC) :: mat
      integer :: m,c

      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      mat = matbulkscobj(out%objidx)

      c = 0
      sum = 0.

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      select case ( mode ) 
      case ( 1 )

         val1 = mat%N(p)

         if ( out%mode .ne. 'S' ) then
 
            if ( reg%isbox ) then
               write(out%funit,"(1E15.6E3)") real(val1,8)
            else
               write(out%funit,"(M4_SDIM({I5}),(1E15.6E3))") M4_DIM123({i},{i,j},{i,j,k}), real(val1,8)
            endif
         
         else
            sum(1) = sum(1) + val1
         endif

      case ( 2 )

         val1 = mat%Psum(1,mat%cyc,p)
         val2 = mat%Psum(2,mat%cyc,p)
         val3 = mat%Psum(3,mat%cyc,p)

         if ( out%mode .ne. 'S' ) then

            if ( reg%isbox ) then
               write(out%funit,"(3E15.6E3)") real(val1,8), real(val2,8), real(val3,8)
           else
               write(out%funit,"(M4_SDIM({I5}),(3E15.6E3))") M4_DIM123({i},{i,j},{i,j,k}),real(val1,8), &
                    real(val2,8), real(val3,8)
            endif

         else
            sum(1) = sum(1) + val1
            sum(2) = sum(2) + val2
            sum(3) = sum(3) + val3
         end if

      end select      
      
      c = c + 1

      },{if ( reg%is .ne. reg%ie ) write(out%funit,*)},{if ( reg%js .ne. reg%je ) write(out%funit,*)} )
   
      if ( out%mode .eq. 'S' ) then
         if ( mode == 1 ) then
            write(out%funit,"(1E15.6E3)") real(sum(1),8)/c
         else
            write(out%funit,"(3E15.6E3)") real(sum,8)/c  
         endif
      endif

    end subroutine WriteValues

  end subroutine WriteDataMatBulkscOutgplObj


end module matbulksc_outgpl


!
! Authors:  S.Wuestner J.Wood
! Modified: 03/07/2013
!
!======================================================================
