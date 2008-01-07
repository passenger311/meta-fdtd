! -*-F90-*-------------------------------------------------------------
!
!  module: utils / meta
!
!  utility functions.
!
! ---------------------------------------------------------------------
 

!======================================================================
!
! 
!

module utils

  implicit none
  save

contains

  integer function checktick()
                                                                               
    !     .. Parameters ..
    integer, parameter :: n = 20

    !     .. Local Scalars ..
    real(kind=8) :: t1, t2
    integer :: i,j,jmin

    !     .. Local Arrays ..
    real(kind=8) :: timesfound(n)

    !     .. External Functions ..
    real(kind=8) :: mysecond
    external mysecond

    !     .. Intrinsic Functions ..
    intrinsic :: max,min,nint

    do i = 0, n-1
       t2 = mysecond()
       
       do while (t2.ne.t1)
          t2 = mysecond()
       end do
       
       t1 = t2
       timesfound(i) = t1
       
    end do
    
    jmin = 1000000
    do i = 2,n
       j = nint((timesfound(i)-timesfound(i-1))*1d6)
       jmin = min(jmin,max(j,0))
    end do
       
    if (jmin.GT.0) then
       checktick = jmin
    else
       write(6,*) 'Your clock granularity appears to be less ', &
            'than one microsecond'
       checktick = 1
    end if
    return
    
  end function checktick

end module utils

!
! Authors:  J.Hamm
! Modified: 7/1/2008
!
!======================================================================


