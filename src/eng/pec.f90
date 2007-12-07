!-*- F90 -*------------------------------------------------------------
!
!  module: pec / meta3
!
!  set metallic boundary conditions.
!
!  subs:
!
!  InitializePec
!  FinalizePec
!  SetPec
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
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'pec'

  ! --- Public Methods

  public :: InitializePec
  public :: FinalizePec
  public :: SetPec

  ! --- Public Data

  ! --- Constants

  ! --- Data

contains

!----------------------------------------------------------------------

  subroutine InitializePec

  end subroutine InitializePec

!----------------------------------------------------------------------

  subroutine FinalizePec

  end subroutine FinalizePec

!----------------------------------------------------------------------

  subroutine SetPec
    
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

  end subroutine SetPec

!----------------------------------------------------------------------

end module pec

!
! Authors:  S.Scholz, A.Klaedtke, J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
