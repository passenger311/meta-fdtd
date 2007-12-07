! -*-F90-*-------------------------------------------------------------
!
!  module: strings / meta3
!
!  fortran smart string handling
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
        
    cat9 = cat(cats(cats(cats(cats(cats(cats( &
         cats(s1,s2),s3),s4),s5),s6),s7),s8),s9)

    return
  end function cat9
      
  character(len=255) function cat8(s1,s2,s3,s4,s5,s6,s7,s8)
    
    implicit none

    character(len=*) :: s1,s2,s3,s4,s5,s6,s7,s8
    
    cat8 = cat(cats(cats(cats(cats(cats( &
         cats(s1,s2),s3),s4),s5),s6),s7),s8)

    return
  end function cat8

  character(len=255) function cat7(s1,s2,s3,s4,s5,s6,s7)

    implicit none
    
    character(len=*) :: s1,s2,s3,s4,s5,s6,s7

    cat7 = cat(cats(cats(cats(cats(cat(s1,s2),s3),s4),s5),s6),s7)

    return
  end function cat7

  character(len=255) function cat6(s1,s2,s3,s4,s5,s6)
    
    implicit none
    
    character(len=*) :: s1,s2,s3,s4,s5,s6

    cat6 = cat(cats(cats(cats(cats(s1,s2),s3),s4),s5),s6)

    return
  end function cat6

  character(len=255) function cat5(s1,s2,s3,s4,s5)

    implicit none

    character(len=*) :: s1,s2,s3,s4,s5
    
    cat5 = cat(cats(cats(cats(s1,s2),s3),s4),s5)

    return
  end function cat5

  character(len=255) function cat4(s1,s2,s3,s4)

    implicit none

    character(len=*) :: s1,s2,s3,s4

    cat4 = cat(cats(cats(s1,s2),s3),s4)

    return
  end function cat4

  character(len=255) function cat3(s1,s2,s3)

    implicit none

    character(len=*) :: s1,s2,s3

    cat3 = cat(cats(s1,s2),s3)

    return
  end function cat3


  character(len=255) function cat2(s1,s2)

    implicit none

    character(len=*) :: s1,s2

    cat2 = cat(s1,s2)

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


  ! ---- wipe, removes *all* spaces

  character(len=255) function wipe(str)

    implicit none

    character(len=*) :: str

    character(len=1) :: c
    character(len=1024) :: nstr
    integer :: i,j,l

    l = strlen(str)

    j = 1
    do i=1,l
       c = str(i:i) 
       if ( c .ne. ' ' .and. iachar(c) .ne. 9 ) then
          nstr(j:j) = c
          j = j + 1
       endif
    enddo
    do i=j,len(nstr)
       nstr(i:i) = ' '
    enddo
    
    wipe = nstr

    return

  end function wipe



  ! ---- ltrim, removes leading spaces

  character(len=255) function ltrim(str)

    implicit none

    character(len=*) :: str

    character(len=1) :: c
    character(len=1024) :: nstr 
    integer :: i, j
    logical :: flag

    j = 1
    
    flag = .false.
    nstr = ''
    
    do i = 1, len(str)
       
       c = str(i:i)
       
       if ( ( c .ne. ' ' .and. iachar(c) .ne. 9 ) &
            .or. flag ) then
          
          nstr(j:j) = str(i:i)  
          j = j + 1
          flag = .true.
       endif

    enddo

    ltrim = nstr

    return
  end function ltrim


  ! ---- cat 2 strings, add string terminator '~'

  character(len=255) function cats(str1,str2)

    implicit none

    character(len=*) :: str1,str2

    character(len=1024) :: str
    integer :: l1,l2

    l1 = strlen(str1)
    l2 = strlen(str2)
    str = str1(1:l1)
    str(l1+1:l1+l2) = str2(1:l2)
    str(l1+l2+1:l1+l2+1) = '~'
    
    cats = str

    return
  end function cats


! ---- cat 2 strings

  character(len=255) function cat(str1,str2)

    implicit none

    character(len=*) :: str1,str2

    character(len=1024) :: str
    integer :: l1,l2

    l1 = strlen(str1)
    l2 = strlen(str2)
    str = str1(1:l1)
    str(l1+1:l1+l2) = str2(1:l2)

    cat = str
    
    return
  end function cat


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


