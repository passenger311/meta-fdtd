


























 














































































!-*- F90 -*------------------------------------------------------------
!
!  module: mattwolvl_outgpl / meta
!
!  this module handles GPL output of data related to the mattwolvl module.
!
!----------------------------------------------------------------------


module mattwolvl_outgpl

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

  public :: InitializeMattwolvlOutgplObj
  public :: FinalizeMattwolvlOutgplObj 
  public :: WriteDataMattwolvlOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeMattwolvlOutgplObj(out)

    type (T_OUT) :: out

  end subroutine InitializeMattwolvlOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeMattwolvlOutgplObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeMattwolvlOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataMattwolvlOutgplObj(out, mode)

    type (T_OUT) :: out
    type (T_BUF) :: buf
    logical :: mode

    

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

      
 type(T_REG) :: reg
 integer :: p,i,j,k
 real(kind=8) :: w(3)
  
      real(kind=8) :: val1, val2, val3, sum1, sum2, sum3
      type(T_MATTWOLVL) :: mat
      integer :: m

      

      reg = regobj(out%regidx)
      mat = mattwolvlobj(out%objidx)

      sum1 = 0.
      sum2 = 0.
      sum3 = 0.

      
! --- START M4: REGLOOP_EXPR
! reg,p,i,j,k,w,expr,spac1, spac2

if ( reg%numnodes .gt. 0 ) then

if ( reg%isbox ) then
p = 0
do k = reg%ks, reg%ke, reg%dk
do j = reg%js, reg%je, reg%dj
do i = reg%is, reg%ie, reg%di
p = p + 1	     

if ( .not. reg%isref ) then
if ( reg%compressval ) then 
w(1:reg%numval) = reg%val(1:reg%numval,reg%valptr(p))
else
w(1:reg%numval) = reg%val(1:reg%numval,p)
endif
endif 



      select case ( mode ) 
      case ( 1 )

         val1 = mat%inversion(p)

         if ( out%mode .ne. 'S' ) then
 
            if ( reg%isbox ) then
               write(out%funit,"(1E15.6E3)") real(val1,8)
            else
               write(out%funit,"(1I5,(1E15.6E3))") &
                    i,real(val1,8)
            endif
         
         else
            sum1 = sum1 + val1
         endif

      case ( 2 )

         val1 = mat%rho12(p)

         if ( out%mode .ne. 'S' ) then

            if ( reg%isbox ) then
               write(out%funit,"(1E15.6E3)") real(val1,8)
            else
               write(out%funit,"(1I5,(1E15.6E3))") &
                    i,real(val1,8)
            endif

         else
            sum1 = sum1 + val1
         end if

      end select      

      
	   
enddo !i

      
      if ( reg%is .ne. reg%ie ) write(out%funit,*)
enddo !j

enddo !k

else !isbox


do k = reg%ks, reg%ke, reg%dk ! -> k

do j = reg%js, reg%je, reg%dj ! -> j
!$OMP PARALLEL DO PRIVATE(p,w)
do i = reg%is, reg%ie, reg%di ! -> i
p = reg%mask(i,j,k) ! -> p
if ( p .gt. 0 ) then

if ( .not. reg%isref ) then
if ( reg%compressval ) then 
w(1:reg%numval) = reg%val(1:reg%numval,reg%valptr(p))
else
w(1:reg%numval) = reg%val(1:reg%numval,p)
endif
endif



      select case ( mode ) 
      case ( 1 )

         val1 = mat%inversion(p)

         if ( out%mode .ne. 'S' ) then
 
            if ( reg%isbox ) then
               write(out%funit,"(1E15.6E3)") real(val1,8)
            else
               write(out%funit,"(1I5,(1E15.6E3))") &
                    i,real(val1,8)
            endif
         
         else
            sum1 = sum1 + val1
         endif

      case ( 2 )

         val1 = mat%rho12(p)

         if ( out%mode .ne. 'S' ) then

            if ( reg%isbox ) then
               write(out%funit,"(1E15.6E3)") real(val1,8)
            else
               write(out%funit,"(1I5,(1E15.6E3))") &
                    i,real(val1,8)
            endif

         else
            sum1 = sum1 + val1
         end if

      end select      

      
	   
endif 
enddo ! i
enddo ! j
enddo ! k

endif !isbox

endif !numnodes > 0

! --- END M4: REGLOOP_EXPR

   
      if ( out%mode .eq. 'S' ) then
         write(out%funit,"(1E15.6E3)") real(sum1,8)
      endif

    end subroutine WriteValues

  end subroutine WriteDataMattwolvlOutgplObj


end module mattwolvl_outgpl


!
! Authors:  A.Pusch, J.Hamm
! Modified: 13/08/2010
!
!======================================================================
