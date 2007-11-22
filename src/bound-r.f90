!----------------------------------------------------------------------
!
!  module: bound(-r) / max3d
!
!  set metallic boundary conditions.
!
!----------------------------------------------------------------------
 
module bound
 
  use grid
  use constant
  use mpistart
  implicit none
  save
  
contains

  subroutine PEC
    
    implicit none
 
    integer :: i,j,k


! top, bottom

    do k=KBEG, KMAX
       do i=IBEG, IMAX
          Ez(i,0,k)=0.0
          Ex(i,0,k)=0.0
          Ez(i,JMAX,k)=0.0
          Ex(i,JMAX,k)=0.0
       enddo
       do j=JBEG, JMAX
          Ez(0,j,k)=0.0
          Ey(0,j,k)=0.0
          Ez(IMAX,j,k)=0.0
          Ey(IMAX,j,k)=0.0
       enddo
    enddo


! left

    if ( myrank .eq. 0 ) then 

       do j=JBEG, JMAX
          do i=IBEG, IMAX
             Ey(i,j,KBEG)=0.0
             Ex(i,j,KBEG)=0.0
          enddo
       enddo
    
    endif

! right

    if ( myrank .eq. numproc-1 ) then 

       do j=JBEG, JMAX
          do i=IBEG, IMAX
             Ey(i,j,KMAX)=0.0
             Ex(i,j,KMAX)=0.0
          enddo
       enddo
    
    endif

  end subroutine PEC

end module bound
