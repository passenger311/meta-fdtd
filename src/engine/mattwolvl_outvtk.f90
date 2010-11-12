!-*- F90 -*------------------------------------------------------------
!
!  module: mattwolvl_outvtk / meta
!
!  this module handles VTK output of data related to the mattwolvl module.
!
!----------------------------------------------------------------------


module mattwolvl_outvtk

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use mattwolvl
  use out_calc

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MATTWOLVL_OUTVTK'

 ! --- Public Methods

  public :: InitializeMattwolvlOutvtkObj
  public :: FinalizeMattwolvlOutvtkObj
  public :: WriteDataMattwolvlOutvtkObj

contains

!----------------------------------------------------------------------

  subroutine InitializeMattwolvlOutvtkObj(out)

    type (T_OUT) :: out

  end subroutine InitializeMattwolvlOutvtkObj


!----------------------------------------------------------------------

  subroutine FinalizeMattwolvlOutvtkObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeMattwolvlOutvtkObj


!----------------------------------------------------------------------

  subroutine WriteDataMattwolvlOutvtkObj(out, mode)

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
      type(T_MATTWOLVL) :: mat

      M4_WRITE_DBG({"WriteScalars!"})

      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      mat = mattwolvlobj(out%objidx)

      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      write(out%funit,"(A,A)") "SCALARS scalars float 3"
      write(out%funit,*) "LOOKUP_TABLE default"

!      write(out%funit,"(A)") "VECTORS vectors float"

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      select case ( mode )
      case ( 1 )

         val1 = mat%inversion(p)

         if ( out%mode .ne. 'S' ) then

           write(out%funit,"(1E15.6E3)") real(val1,8)

         else

            sum1 = sum1 + val1

         endif

      case ( 2 )

         val1 = real(mat%rho12(p))

         if ( out%mode .ne. 'S' ) then

            write(out%funit,"(1E15.6E3)") real(val1,8)

         else

            sum1 = sum1 + val1

         end if

      end select

      },{

      if ( reg%is .ne. reg%ie ) write(out%funit,*)}, {if ( reg%js .ne. reg%je ) write(out%funit,*)

      } )

      if ( out%mode .eq. 'S' ) then
         write(out%funit,"(1E15.6E3)") real(sum1,8)
      endif

    end subroutine WriteValues

  end subroutine WriteDataMattwolvlOutvtkObj


end module mattwolvl_outvtk

!
! Authors:  J.Hamm, S.Wuestner, A.Pusch

! Modified: 13/08/2010
!
!======================================================================
