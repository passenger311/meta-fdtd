
!-*- F90 -*------------------------------------------------------------
!
!  module: timer
!
!  semi-portable timer module.
!
!----------------------------------------------------------------------


! =====================================================================
!

module timer

  use constant
 
  private
  save

  character(len=20), private, parameter :: modname = 'TIMER'

  public :: CheckTimer, ResetTimer, DisplayTotalTimer, DisplayMcpsTimer
  public :: StartTimer, StopTimer


  integer :: clockgranularity
  real(kind=8) :: timervalue(MAXTIMER,2)

contains

!----------------------------------------------------------------------

  integer function checktick()

    integer, parameter :: n = 20

    real(kind=8) :: t1, t2
    integer :: i,j,jmin

    real(kind=8) :: timesfound(n)

    real(kind=8),external :: mysecond

      i = 0

   10 t2 = mysecond()
      IF (t2.EQ.t1) GO TO 10

      t1 = t2
      i = i + 1
      timesfound(i) = t1
      IF (i.LT.n) GO TO 10

      jmin = 1000000
      DO 20 i = 2,n
          j = nint((timesfound(i)-timesfound(i-1))*1d6)
          jmin = min(jmin,max(j,0))
   20 CONTINUE

      IF (jmin.GT.0) THEN
          checktick = jmin
      ELSE
          checktick = 1
      END IF

      RETURN
    
  end function checktick


!----------------------------------------------------------------------

  
  subroutine CheckTimer()

    clockgranularity = checktick()

    if ( clockgranularity  .ne. 1 ) then
       M4_FATAL_ERROR({"Clock granularity is bigger than 1 microsecond!"})
    end if

  end subroutine CheckTimer


!----------------------------------------------------------------------

  subroutine ResetTimer(i)
    
    integer :: i
    real(kind=8),external :: mysecond

    if ( i .lt. 1 .or. i .gt. 10 ) return

    timervalue(i,1) = 0
    timervalue(i,2) = mysecond()

  end subroutine ResetTimer
  
!----------------------------------------------------------------------

  subroutine StartTimer(i)
    
    integer :: i
    real(kind=8),external :: mysecond

    if ( i .lt. 1 .or. i .gt. 10 ) return

    timervalue(i,2) = mysecond()

  end subroutine StartTimer
  

!----------------------------------------------------------------------

  subroutine StopTimer(i)
    
    integer :: i
    real(kind=8) :: now, laptime
    real(kind=8),external :: mysecond

    if ( i .lt. 1 .or. i .gt. 10 ) return

    now = mysecond()

    laptime = now - timervalue(i,2)

    timervalue(i,1) = timervalue(i,1) + laptime
    timervalue(i,2) = now

  end subroutine StopTimer

!----------------------------------------------------------------------

  subroutine DisplayTotalTimer(str,i)

    character(len=*) :: str
    integer :: i

    M4_WRITE_FMT_INFO({A,F10.6,A},{str,timervalue(i,1)," secs"}) 

  end subroutine DisplayTotalTimer


!----------------------------------------------------------------------

  subroutine DisplayMcpsTimer(str,i,cells)

    character(len=*) :: str
    integer :: i,cells

    M4_WRITE_FMT_INFO({A,F10.6,A},{str,cells/(timervalue(i,1)*1.e6)," mcps"}) 

  end subroutine DisplayMcpsTimer


!----------------------------------------------------------------------

end module timer

!
! Authors:  J.Hamm
! Modified: 30/04/2008
!
!======================================================================


