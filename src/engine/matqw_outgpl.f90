!-*- F90 -*------------------------------------------------------------
!
!  module: matqw_outgpl / meta
!
!  this module handles GPL output of data related to the matfourlvl module.
!
!----------------------------------------------------------------------


module matqw_outgpl

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use out_calc
  use matqw

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MATQW_OUTGPL'

 ! --- Public Methods

  public :: InitializeMatQWOutgplObj
  public :: FinalizeMatQWOutgplObj
  public :: WriteDataMatQWOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeMatQWOutgplObj(out)

    type (T_OUT) :: out

  end subroutine InitializeMatQWOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeMatQWOutgplObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeMatQWOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataMatQWOutgplObj(out, mode)

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
    case('M')
       call WriteValues(out, 3)
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED" 
    end select
    
  contains

    ! **************************************************************** !

    subroutine WriteValues(out, mode)

      type (T_OUT) :: out
      integer :: mode,pk

      M4_REGLOOP_DECL(reg,p,i,j,k,w(3))  
      real(kind=8) :: val1, val2, val3, val4, val5, val6, sum, val7
      type(T_MATQW) :: mat
      integer :: m

      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      mat = matqwobj(out%objidx)

      sum = 0.

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      select case ( mode ) 
      case ( 1 )
         ! microscopic densities
         do pk=1, mat%kCount
            val1 = pk
            val2 = mat%fe(p,pk)
            val3 = mat%fh(p,pk)
            val4 = mat%fe_eq(p,pk)
            val5 = mat%fh_eq(p,pk)
            if ( out%mode .ne. 'S' ) then
 
               if ( reg%isbox ) then
                  write(out%funit,"(5E15.6E3)") dble(val1), dble(val2), dble(val3), dble(val4), dble(val5)
               else
                  write(out%funit,"(M4_SDIM({I5}),(5E15.6E3))") M4_DIM123({i},{i,j},{i,j,k}), dble(val1), &
                    dble(val2), dble(val3), dble(val4), dble(val5)
               endif
         
            else
               sum = sum + val2 * mat%pW(pk)
            endif
         end do

      case ( 2 )
         do pk=1, mat%kCount
            ! microscopic polarisations
            val1 = pk
            val2 = mat%Px(p,pk,mat%cyc)
            val3 = mat%Py(p,pk,mat%cyc)
            val4 = mat%Pz(p,pk,mat%cyc)
            val5 = mat%gnx(p,pk)
            val6 = mat%gny(p,pk)
            val7 = mat%gnz(p,pk)
            if ( out%mode .ne. 'S' ) then

               if ( reg%isbox ) then
                  write(out%funit,"(7E15.6E3)") dble(val1), dble(val2), dble(val3), dble(val4), &
                       dble(val5), dble(val6), dble(val7)
               else
                  write(out%funit,"(M4_SDIM({I5}),(7E15.6E3))") M4_DIM123({i},{i,j},{i,j,k}),dble(val1), &
                       dble(val2), dble(val3), dble(val4), dble(val5), dble(val6), dble(val7)
               endif

            else
               sum = sum + val1
            end if
         end do
      case ( 3)
         !macroscopic values
         val1 = mat%n(p)
         val2 = mat%G(p)
         val3 = mat%P_ma(p)
         if ( out%mode .ne. 'S' ) then

            if ( reg%isbox ) then
               write(out%funit,"(3E15.6E3)") dble(val1), dble(val2), dble(val3)
            else
               write(out%funit,"(M4_SDIM({I5}),(3E15.6E3))") M4_DIM123({i},{i,j},{i,j,k}),dble(val1), &
                    dble(val2), dble(val3)
            endif

         else
            sum = sum + val1
         end if
      end select

      },{if ( reg%is .ne. reg%ie ) write(out%funit,*)}, {if ( reg%js .ne. reg%je ) write(out%funit,*)} )
   
      if ( out%mode .eq. 'S' ) then
         write(out%funit,"(E15.6E3)") dble(sum)
      endif

    end subroutine WriteValues

  end subroutine WriteDataMatQWOutgplObj


end module matqw_outgpl


!
! Authors:  J.Hamm, A. Pusch
! Modified: 19/04/2010
!
!======================================================================
