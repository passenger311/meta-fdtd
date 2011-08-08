!-*- F90 -*------------------------------------------------------------
!
!  module: matfourlvl_outgpl / meta
!
!  this module handles GPL output of data related to the matfourlvl module.
!
!----------------------------------------------------------------------


module matfourlvl_outgpl

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use out_calc
  use matfourlvl

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MATFOURLVL_OUTGPL'

 ! --- Public Methods

  public :: InitializeMatFourlvlOutgplObj
  public :: FinalizeMatFourlvlOutgplObj
  public :: WriteDataMatFourlvlOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeMatFourlvlOutgplObj(out)

    type (T_OUT) :: out

  end subroutine InitializeMatFourlvlOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeMatFourlvlOutgplObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeMatFourlvlOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataMatFourlvlOutgplObj(out, mode)

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
      real(kind=8) :: val1, val2, val3, val4, val5, val6, sum(6)
      type(T_MATFOURLVL) :: mat
      integer :: m,c

      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      mat = matfourlvlobj(out%objidx)

      c = 0
      sum = 0.

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      select case ( mode ) 
      case ( 1 )

         val1 = mat%N0(p)
         val2 = mat%N1(p)
         val3 = mat%N2(p)
         val4 = mat%N3(p)

         if ( out%mode .ne. 'S' ) then
 
            if ( reg%isbox ) then
               write(out%funit,"(4E15.6E3)") real(val1,8), real(val2,8), real(val3,8), real(val4,8)
            else
               write(out%funit,"(M4_SDIM({I5}),(4E15.6E3))") M4_DIM123({i},{i,j},{i,j,k}), real(val1,8), &
                    real(val2,8),real(val3,8),real(val4,8)
            endif
         
         else
            sum(1) = sum(1) + val1
            sum(2) = sum(2) + val2
            sum(3) = sum(3) + val3
            sum(4) = sum(4) + val4
         endif

      case ( 2 )

         val1 = mat%Pax(p,mat%cyc)
         val2 = mat%Pay(p,mat%cyc)
         val3 = mat%Paz(p,mat%cyc)
         val4 = mat%Pbx(p,mat%cyc)
         val5 = mat%Pby(p,mat%cyc)
         val6 = mat%Pbz(p,mat%cyc)

         if ( out%mode .ne. 'S' ) then

            if ( reg%isbox ) then
               write(out%funit,"(6E15.6E3)") real(val1,8), real(val2,8), real(val3,8), real(val4,8), &
                    real(val5,8),real(val6,8)
           else
               write(out%funit,"(M4_SDIM({I5}),(6E15.6E3))") M4_DIM123({i},{i,j},{i,j,k}),real(val1,8), &
                    real(val2,8), real(val3,8), real(val4,8), real(val5,8), real(val6,8)
            endif

         else
            sum(1) = sum(1) + val1
            sum(2) = sum(2) + val2
            sum(3) = sum(3) + val3
            sum(4) = sum(4) + val4
            sum(5) = sum(5) + val5
            sum(6) = sum(6) + val6
         end if

      end select      
      
      c = c + 1

      },{if ( reg%is .ne. reg%ie ) write(out%funit,*)},{if ( reg%js .ne. reg%je ) write(out%funit,*)} )
   
      if ( out%mode .eq. 'S' ) then
         if ( mode == 1 ) then
            write(out%funit,"(4E15.6E3)") real(sum(1:4),8)/c
         else
            write(out%funit,"(6E15.6E3)") real(sum,8)/c  
         endif
      endif

    end subroutine WriteValues

  end subroutine WriteDataMatFourlvlOutgplObj


end module matfourlvl_outgpl


!
! Authors:  J.Hamm, A. Pusch
! Modified: 11/01/2010
!
!======================================================================
