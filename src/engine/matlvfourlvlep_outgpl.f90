!-*- F90 -*------------------------------------------------------------
!
!  module: MatLVFourlvlep_outgpl / meta
!
!  this module handles GPL output of data related to the MatLVFourlvlep module.
!
!----------------------------------------------------------------------


module MatLVFourlvlep_outgpl

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use out_calc
  use MatLVFourlvlep

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MatLVFourlvlep_OUTGPL'

 ! --- Public Methods

  public :: InitializeMatLVFourlvlepOutgplObj
  public :: FinalizeMatLVFourlvlepOutgplObj
  public :: WriteDataMatLVFourlvlepOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeMatLVFourlvlepOutgplObj(out)

    type (T_OUT) :: out

  end subroutine InitializeMatLVFourlvlepOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeMatLVFourlvlepOutgplObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeMatLVFourlvlepOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataMatLVFourlvlepOutgplObj(out, mode)

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
    case('R')
       call WriteValues(out, 3)
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED" 
    end select
    
  contains

    ! **************************************************************** !

    subroutine WriteValues(out, mode)

      type (T_OUT) :: out
      integer :: mode

      M4_REGLOOP_DECL(reg,p,i,j,k,w(3))  
      real(kind=8) :: val1, val2, val3, val4, val5, val6, sum(6),lEx,lEy,lEz,pema,pena
      type(T_MatLVFourlvlep) :: mat
      integer :: m,n,c

      M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

      reg = regobj(out%regidx)
      mat = MatLVFourlvlepobj(out%objidx)

      c = 0
      sum = 0.
      m = mat%cyc
      n = mod(m,2)+1

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      select case ( mode ) 
      case ( 1 )

         val1 = mat%N0(p)
         val2 = mat%N1(p)
         val3 = mat%N2(p)
         val4 = mat%N3(p)

         if ( out%mode .ne. 'S' ) then
 
            if ( reg%isbox ) then
               write(out%funit,"(4E15.6E3)") dble(val1), dble(val2), dble(val3), dble(val4)
            else
               write(out%funit,"(M4_SDIM({I5}),(4E15.6E3))") M4_DIM123({i},{i,j},{i,j,k}), dble(val1), &
                    dble(val2), dble(val3), dble(val4)
            endif
         
         else
            sum(1) = sum(1) + val1
            sum(2) = sum(2) + val2
            sum(3) = sum(3) + val3
            sum(4) = sum(4) + val4
         endif

      case ( 2 )

         val1 = mat%Pax(p,m)
         val2 = mat%Pay(p,m)
         val3 = mat%Paz(p,m)
         val4 = 0!mat%Pbx(p,m)
         val5 = 0!mat%Pby(p,m)
         val6 = 0!mat%Pbz(p,m)

         if ( out%mode .ne. 'S' ) then

            if ( reg%isbox ) then
               write(out%funit,"(3E15.6E3)") dble(val1), dble(val2), dble(val3) !, dble(val4), &
!                    dble(val5), dble(val6)
           else
               write(out%funit,"(M4_SDIM({I5}),(3E15.6E3))") M4_DIM123({i},{i,j},{i,j,k}),dble(val1), &
                    dble(val2), dble(val3) !, dble(val4), dble(val5), dble(val6)
            endif

         else
            sum(1) = sum(1) + val1
            sum(2) = sum(2) + val2
            sum(3) = sum(3) + val3
            sum(4) = sum(4) + val4
            sum(5) = sum(5) + val5
            sum(6) = sum(6) + val6
         end if

      case ( 3 )

         lEx = Ex(i,j,k) * ( mat%epsLFE * ( 2. + 1./epsinvx(i,j,k) ) / 3. + ( 1. - mat%epsLFE ) )
         lEy = Ey(i,j,k) * ( mat%epsLFE * ( 2. + 1./epsinvy(i,j,k) ) / 3. + ( 1. - mat%epsLFE ) )
         lEz = Ez(i,j,k) * ( mat%epsLFE * ( 2. + 1./epsinvz(i,j,k) ) / 3. + ( 1. - mat%epsLFE ) )

         pema = 2 * mat%Pax(p,m) * lEx + mat%Pay(p,m) * lEy + mat%Paz(p,m) * lEz
         pena = 2 * mat%Pax(p,n) * lEx + mat%Pay(p,n) * lEy + mat%Paz(p,n) * lEz

         val1 = mat%rp * mat%N0(p)
         val2 = ( mat%gamma30 + mat%gamma32 ) * mat%N3(p)
         val3 = - mat%x2fac1 * ( pema - pena ) - mat%x2fac2 * ( pena + pema ) + mat%gamma21 * mat%N2(p)
         val4 = mat%gamma10 * mat%N1(p)

         if ( out%mode .ne. 'S' ) then

            if ( reg%isbox ) then
               write(out%funit,"(2E15.6E3)") dble(val1), dble(val2), dble(val3), dble(val4)
           else
               write(out%funit,"(M4_SDIM({I5}),(2E15.6E3))") M4_DIM123({i},{i,j},{i,j,k}),dble(val1), &
                    dble(val2), dble(val3), dble(val4)
            endif

         else
            sum(1) = sum(1) + val1
            sum(2) = sum(2) + val2
            sum(3) = sum(3) + val3
            sum(4) = sum(4) + val4
         end if

      end select      
      
      c = c + 1

      },{if ( reg%is .ne. reg%ie ) write(out%funit,*)},{if ( reg%js .ne. reg%je ) write(out%funit,*)} )
   
      if ( out%mode .eq. 'S' ) then
         select case ( mode )
         case ( 1 )
            write(out%funit,"(4E15.6E3)") dble(sum(1:4))/c
         case ( 2 )
            write(out%funit,"(3E15.6E3)") dble(sum(1:3))/c  
         case ( 3 )
            write(out%funit,"(4E15.6E3)") dble(sum(1:4))/c
         end select
      endif

    end subroutine WriteValues

  end subroutine WriteDataMatLVFourlvlepOutgplObj


end module MatLVFourlvlep_outgpl


!
! Authors:  J.Hamm, A. Pusch
! Modified: 07/08/2011
!
!======================================================================
