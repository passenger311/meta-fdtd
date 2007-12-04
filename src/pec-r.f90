!-*- F90 -*------------------------------------------------------------
!
!  module: pec / max3d
!
!  set metallic boundary conditions.
!
!  subs:
!
!  InitializePEC
!  FinalizePEC
!  SetPEC
!
!----------------------------------------------------------------------
 
!======================================================================
!
!

module pec
 
  use constant
  use mpiworld
  use grid
  use fdtd

  implicit none
  save
  

contains

!----------------------------------------------------------------------

  subroutine InitializePEC

  end subroutine InitializePEC

!----------------------------------------------------------------------

  subroutine FinalizePEC

  end subroutine FinalizePEC

!----------------------------------------------------------------------

  subroutine SetPEC
    
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

  end subroutine SetPEC

!----------------------------------------------------------------------

end module pec

!
! Authors:  S.Scholz, A.Klaedtke, J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
