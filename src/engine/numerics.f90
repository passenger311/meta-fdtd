
!-*- F90 -*------------------------------------------------------------
!
!  module: numerics
!
!  some generic numerical algorithms
!
!----------------------------------------------------------------------


! =====================================================================
!

module numerics

  use constant
 
  private
  save


  character(len=20), private, parameter :: modname = 'NUMERICS'

  public :: root_find

contains

! ********** SUB : ROOT_FIND

! newton raphson root finding algorithm a la numerical recipies

! find zero <xr> of function < valfunc(x) - val > inside the interval [<xl>,<xr>] 
! <xr> : initial guess on entry, root position on exit. 
! <iter> : if -1 on exit root_find failed, else it is the number of iterations needed
! <xacc> : absolute x accuracy to reach 


  recursive subroutine root_find(xr, val, val_func, xl, xh, xacc, iter)

    real(8) :: xr, val, xl, xh, xacc
    real(8), external :: val_func
    integer :: iter
    
!    integer :: i
    real(8) :: x1, x2, dx, dx_old
    real(8) :: yr, y1, y2, dy, yr_dx

    integer, parameter :: max_it = 200 

    iter = -1

    y1 = val - val_func(xl)
    y2 = val - val_func(xh)
    
    if ( y1 * y2 .ge. 0 ) then 

       M4_WRITE_WARN({"root_find: no zero or more than one!"})
       return 

    end if

    if ( y1 .lt. 0 ) then

       x1 = xl    ! orient search so that f(x1) < 0 !
       x2 = xh
       
    else

       x1 = xh
       x2 = xl

    endif

!    xr = 0.5*( x1 + x2 )   ! initial guess
    dx_old = abs(x2 - x1)
    dx = dx_old

   yr = val - val_func(xr)
   yr_dx = val - val_func(xr+xacc)
   dy = (yr_dx - yr)/(xacc)


!    dy = (y1 - y2)/(xl - xh)

    do iter = 1, max_it
       
!       write(6,*) x1, x2

       if ( (( xr - x2 )*dy - yr)*(( xr - x1 )*dy - yr) .ge. 0 .or. &
            abs( 2.*yr ) .gt. abs(dx_old * dy) ) then

! switch to bisect if newton is not decreasing fast enough or
! leaves the intervall       

          dx_old = dx
          dx = 0.5 * (x2 - x1)
          xr = x1 + dx

          if ( x1 .eq. xr ) return 

      
       else

! newton step
          
          dx_old = dx
          dx = yr/dy

          if ( xr - dx .eq. xr ) return

          xr = xr - dx
 
       end if

       if ( abs(dx) .lt. xacc ) return 

       yr = val - val_func(xr)
       yr_dx = val - val_func(xr+xacc)
       dy = (yr_dx - yr)/(xacc)

!       yr = val - val_func(xr)
!       y2 = val - val_func(x2)
!       y1 = val - val_func(x1)
!       dy = (y2 - y1)/(x2 - x1)

       if ( yr .lt. 0 ) then
          
          x1 = xr
          
       else
          
          x2 = xr
          
       end if

    end do

    iter = -1
    M4_WRITE_WARN({"root_find: maximum iterations exceeded!"})
    return

  end subroutine root_find


end module numerics

!
! Authors:  J.Hamm
! Modified: 25/02/2007
!
!======================================================================


