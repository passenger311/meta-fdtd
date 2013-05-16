!-*- F90 -*------------------------------------------------------------
!
!  module: matlvfourlvlep_outvtk / meta
!
!  this module handles VTK output of data related to the matlvfourlvlep module.
!
!----------------------------------------------------------------------


module matlvfourlvlep_outvtk

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use matlvfourlvlep
  use out_calc

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MATLVFOURLVLEP_OUTVTK'

 ! --- Public Methods

  public :: InitializeMatlvfourlvlepOutvtkObj
  public :: FinalizeMatlvfourlvlepOutvtkObj
  public :: WriteDataMatlvfourlvlepOutvtkObj

contains

!----------------------------------------------------------------------

  subroutine InitializeMatlvfourlvlepOutvtkObj(out)

    type (T_OUT) :: out

  end subroutine InitializeMatlvfourlvlepOutvtkObj


!----------------------------------------------------------------------

  subroutine FinalizeMatlvfourlvlepOutvtkObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeMatlvfourlvlepOutvtkObj


!----------------------------------------------------------------------

  subroutine WriteDataMatlvfourlvlepOutvtkObj(out, mode)

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
      real(kind=8) :: val1, val2, val3, val4, val5, val6, sum
      type(T_MATLVFOURLVLEP) :: mat

      M4_WRITE_DBG({"WriteScalars!"})

      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      mat = matlvfourlvlepobj(out%objidx)

      sum = 0.

      select case ( mode )
      case ( 1 )
         write(out%funit,"(A,A)") "SCALARS scalars float 4"
         write(out%funit,*) "LOOKUP_TABLE default"
      case ( 2 )
         write(out%funit,"(A,A)") "SCALARS scalars float 6"
         write(out%funit,*) "LOOKUP_TABLE default"
      end select

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      select case ( mode )
      case ( 1 )

         val1 = mat%N0(p)
         val2 = mat%N1(p)
         val3 = mat%N2(p)
         val4 = mat%N3(p)

         if ( out%mode .ne. 'S' ) then

           write(out%funit,"(4E15.6E3)") real(val1,8), real(val2,8), real(val3,8), real(val4,8)

         else

            sum = sum + val4

         endif

      case ( 2 )

         val1 = mat%Pax(p,mat%cyc)
         val2 = mat%Pay(p,mat%cyc)
         val3 = mat%Paz(p,mat%cyc)
  !       val4 = mat%Pbx(p,mat%cyc)
!         val5 = mat%Pby(p,mat%cyc)
!         val6 = mat%Pbz(p,mat%cyc)

         if ( out%mode .ne. 'S' ) then

            write(out%funit,"(6E15.6E3)") real(val1,8), real(val2,8), real(val3,8)!, &
!            real(val4,8), real(val5,8), real(val6,8)

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

  end subroutine WriteDataMatlvfourlvlepOutvtkObj


end module matlvfourlvlep_outvtk

!
! Authors:  J.Hamm, S.Wuestner, A. Pusch

! Modified: 08/08/2011
!
!======================================================================
