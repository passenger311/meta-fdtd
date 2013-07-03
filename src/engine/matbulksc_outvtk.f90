!-*- F90 -*------------------------------------------------------------
!
!  module: matbulksc_outvtk / meta
!
!  this module handles VTK output of data related to the matbulksc module.
!
!----------------------------------------------------------------------


module matbulksc_outvtk

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use matbulksc
  use out_calc

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MATBULKSC_OUTVTK'

 ! --- Public Methods

  public :: InitializeMatbulkscOutvtkObj
  public :: FinalizeMatbulkscOutvtkObj
  public :: WriteDataMatbulkscOutvtkObj

contains

!----------------------------------------------------------------------

  subroutine InitializeMatbulkscOutvtkObj(out)

    type (T_OUT) :: out

  end subroutine InitializeMatbulkscOutvtkObj


!----------------------------------------------------------------------

  subroutine FinalizeMatbulkscOutvtkObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeMatbulkscOutvtkObj


!----------------------------------------------------------------------

  subroutine WriteDataMatbulkscOutvtkObj(out, mode)

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
      real(kind=8) :: val1, val2, val3, sum
      type(T_MATBULKSC) :: mat

      M4_WRITE_DBG({"WriteScalars!"})

      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      mat = matbulkscobj(out%objidx)

      sum = 0.

      select case ( mode )
      case ( 1 )
         write(out%funit,"(A,A)") "SCALARS scalars float 1"
         write(out%funit,*) "LOOKUP_TABLE default"
      case ( 2 )
         write(out%funit,"(A,A)") "SCALARS scalars float 3"
         write(out%funit,*) "LOOKUP_TABLE default"
      end select

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      select case ( mode )
      case ( 1 )

         val1 = mat%N(p)

         if ( out%mode .ne. 'S' ) then

           write(out%funit,"(1E15.6E3)") real(val1,8)

         else

            sum = sum + val1

         endif

      case ( 2 )

         val1 = mat%Psum(1,mat%cyc,p)
         val2 = mat%Psum(2,mat%cyc,p)
         val3 = mat%Psum(3,mat%cyc,p)

         if ( out%mode .ne. 'S' ) then

            write(out%funit,"(3E15.6E3)") real(val1,8), real(val2,8), real(val3,8)

         else

            sum = sum + val1

         end if

      end select

      },{

      if ( reg%is .ne. reg%ie ) write(out%funit,*)}, {if ( reg%js .ne. reg%je ) write(out%funit,*)

      } )

      if ( out%mode .eq. 'S' ) then
         write(out%funit,"(E15.6E3)") real(sum,8)
      endif

    end subroutine WriteValues

  end subroutine WriteDataMatbulkscOutvtkObj


end module matbulksc_outvtk

!
! Authors:  S.Wuestner, J.Wood

! Modified: 03/07/2013
!
!======================================================================
