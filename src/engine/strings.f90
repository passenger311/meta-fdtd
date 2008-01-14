! -*-F90-*-------------------------------------------------------------
!
!  module: strings / meta
!
!  fortran smarter string handling.
!
! ---------------------------------------------------------------------
 

!======================================================================
!
!
!

module strings
 
  implicit none
  save

contains

!
! ---- join <n> strings <strarr> to outgpl string
!

  character(len=255) function join(strarr, n)

    integer :: n
    character(len=*), dimension(n) :: strarr

    integer :: l, i, j, p
    character(len=1024) :: str

    p = 1

    do i = 1, n

       l = strlen(strarr(i))

       do j = 1, l
          
          str(p:p) = strarr(i)(j:j)
          p = p + 1

       end do

    end do

    join = str

  end function join


!
!
! ---- split string <str> to <elem> substrings
!

  subroutine charsplit(sc, str, strarr, n, elem)

    character(len=*) :: str
    character(len=*) :: sc
    integer :: n, elem
    character(len=*), dimension(n) :: strarr

    logical :: splitit, first
    character(len=1) :: c
    integer ::  i, j, s

    j = 0
    s = 0

    first = .true.
 
    do i = 1, len(str)

       c = str(i:i) 
       if ( c .ne. sc  ) then
          if ( splitit .or. first ) then
             j = 1
             s = s + 1
             strarr(s) = '' 
             first = .false.
          endif
          strarr(s)(j:j) = c
          j = j + 1
          splitit = .false.
       else
          splitit = .true.
       endif

    enddo
    elem = s

  end subroutine charsplit


!
! ---- split string <str> to <elem> substrings at whitespace characters
!

  subroutine wssplit(str, strarr, n, elem)

    character(len=*) :: str
    integer :: n, elem
    character(len=*), dimension(n) :: strarr

    logical :: splitit, first
    character(len=1) :: c
    integer :: i, j, s

    j = 0
    s = 0

    first = .true.

    do i = 1, len(str)

       c = str(i:i) 
       if ( c .ne. ' ' .and. iachar(c) .ne. 9 ) then
          if ( splitit .or. first ) then
             j = 1
             s = s + 1
             strarr(s) = '' 
             first = .false.
          endif
          strarr(s)(j:j) = c
          j = j + 1
          splitit = .false.
       else
          splitit = .true.
       endif

    enddo
    elem = s

  end subroutine wssplit


!
! ---- stripticks,  remove ticks at begin and end of string
!


  subroutine stripticks(str)

    character(len=*) :: str

    integer :: l ,i 

    l = strlen(str)

    if ( ( str(1:1) .eq. '"' .and. str(l:l) .eq. '"' ) .or. &
         ( str(1:1) .eq. "'" .and. str(l:l) .eq. "'" ) ) then

       str(l:l) = '' 
       
       do i = 1, l-2

          str(i:i) = str(i+1:i+1)
          
       end do

       str(l-1:l-1) = '' 

    end if

  end subroutine stripticks


!
!       ---- cat<n>, cat up to 9 strings
!

  character(len=255) function cat9(s1,s2,s3,s4,s5,s6,s7,s8,s9)

    implicit none

    character(len=*) :: s1,s2,s3,s4,s5,s6,s7,s8,s9
        
    cat9 = cat2(cat2(cat2(cat2(cat2(cat2(cat2( &
         cat2(s1,s2),s3),s4),s5),s6),s7),s8),s9)

    return
  end function cat9
      
  character(len=255) function cat8(s1,s2,s3,s4,s5,s6,s7,s8)
    
    implicit none

    character(len=*) :: s1,s2,s3,s4,s5,s6,s7,s8
    
    cat8 = cat2(cat2(cat2(cat2(cat2(cat2( &
         cat2(s1,s2),s3),s4),s5),s6),s7),s8)

    return
  end function cat8

  character(len=255) function cat7(s1,s2,s3,s4,s5,s6,s7)

    implicit none
    
    character(len=*) :: s1,s2,s3,s4,s5,s6,s7

    cat7 = cat2(cat2(cat2(cat2(cat2(cat2(s1,s2),s3),s4),s5),s6),s7)

    return
  end function cat7

  character(len=255) function cat6(s1,s2,s3,s4,s5,s6)
    
    implicit none
    
    character(len=*) :: s1,s2,s3,s4,s5,s6

    cat6 = cat2(cat2(cat2(cat2(cat2(s1,s2),s3),s4),s5),s6)

    return
  end function cat6

  character(len=255) function cat5(s1,s2,s3,s4,s5)

    implicit none

    character(len=*) :: s1,s2,s3,s4,s5
    
    cat5 = cat2(cat2(cat2(cat2(s1,s2),s3),s4),s5)

    return
  end function cat5

  character(len=255) function cat4(s1,s2,s3,s4)

    implicit none

    character(len=*) :: s1,s2,s3,s4

    cat4 = cat2(cat2(cat2(s1,s2),s3),s4)

    return
  end function cat4

  character(len=255) function cat3(s1,s2,s3)

    implicit none

    character(len=*) :: s1,s2,s3

    cat3 = cat2(cat2(s1,s2),s3)

    return
  end function cat3


  character(len=255) function cat2(s1,s2)

    implicit none

    character(len=*) :: s1,s2

    cat2 = TRIM(ADJUSTL(s1)) // ADJUSTL(s2)

    return
  end function cat2

  integer function ilen(i)

    implicit none

    integer :: i
    
    character(len=255) :: str
    
    str = i2str(i)
    ilen = strlen(str)
    
    return
    
  end function ilen


  ! ---- i2str, convert integer to str

  character(len=255) function i2str(i)

    implicit none

    integer :: i

    character(len=1024) :: str

    write(str,*) i
!    i2str = wipe(str)
    i2str = TRIM(ADJUSTL(str))
    return

  end function i2str

  ! ---- f2str, convert float to str

  character(len=255) function f2str(f)

    implicit none

    real(kind=8) :: f

    character(len=1024) :: str

    write(str,"(E8.3E2)") f
!    i2str = wipe(str)
    f2str = TRIM(ADJUSTL(str))
    return

  end function f2str



  ! ---- ltrim, removes leading spaces

  character(len=255) function ltrim(str)

    implicit none

    character(len=*) :: str

    ltrim = TRIM(ADJUSTL(str))

    return
  end function ltrim


! ---- effective string length

  integer function strlen(str)

    implicit none

    character(len=*) :: str

    character(len=1) :: ch
    integer :: i,l

    if ( iachar(str(1:1)) .eq. 0 ) then
       strlen = 0
       return
    endif
    
    l = 0
    do i=1,len(str)    
       ch = str(i:i)
       if (  (ch .ne. ' ' ) .and. &
            ( iachar(ch) .ne. 9 ) ) l = i
       if ( ch .eq. '~') then
          l = i-1
          exit
       endif
    enddo

    strlen = l
    return

  end function strlen

end module strings

!
! Authors:  J.Hamm
! Modified: 4/12/2007
!
!======================================================================


