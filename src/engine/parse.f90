! -*-F90-*-------------------------------------------------------------
!
!  module: parse / meta
!
!  subroutines to parse input files
!
! ---------------------------------------------------------------------
 

!======================================================================
!
!
!

module parse
 
  use constant
  use strings

  implicit none
  save

  character(len=20), private, parameter :: modname = 'PARSE'


contains


! ---------------------------------------------------------------------

  subroutine readline(unit, lcount, eof, line)
    
    integer :: unit, lcount
    logical :: eof
    character(len=*) :: line

    character(len=LINELNG) :: rline, skiptill = ""
    integer :: ios, i
    logical :: wipe, skip = .false.

    do

       line = ""

       read(unit,"(A159)",iostat=ios) line
       
       if ( ios .ne. 0 ) then
          eof = .true.
          exit
       else
          eof = .false. 
          lcount = lcount + 1
       end if
       
       rline = ltrim(ADJUSTL(line))

       if ( line(1:1) .eq. "!" ) cycle

       if ( rline(1:2) .eq. "(!" ) then
          skiptill = cat2(")",rline(3:))
          skip = .true.
          cycle
       end if
       
       if ( skip ) then 
          if ( rline .eq. skiptill ) then 
             skiptill = ""  
             skip = .false.
          end if
          cycle
       else
          ! cut comment
          wipe = .false.
          do i = 1, LINELNG
             if ( wipe ) then
                rline(i:i) = ' ' 
             else
                if ( rline(i:i) .eq. '!' ) then 
                   rline(i:i) = ' ' 
                   wipe = .true.
                end if
             endif
          end do
          
       endif
       
       if ( .not. skip .and. rline .eq. "" ) cycle

       exit

    end do

    line = rline

  end subroutine readline


! ---------------------------------------------------------------------

 subroutine getlogical(line,val, err)

   character(len=*) :: line
   logical :: err
   logical :: val

   integer :: i ,s, e, ios

   s = -1
   e = -1

   if ( err ) return

   do i = 1, LINELNG

      if ( line(i:i) .ne. ' ' .and. s .eq. -1 ) s = i      
      if ( line(i:i) .eq. ' ' .and. s .ne. -1 ) then 
         e = i-1
         read( line(s:e),*, iostat=ios ) val
         if ( ios .ne. 0 ) then 
            err =.true.
         else
            line(s:e) = "" ! wipe
         end if
         return
      end if
   end do
   
   val = .false.
   err = .true. 

 end subroutine getlogical


! ---------------------------------------------------------------------


  subroutine getint(line, val, err)

   character(len=*) :: line
   logical :: err
   integer :: val

   integer :: i ,s, e, ios

   s = -1
   e = -1

   if ( err ) return

   do i = 1, LINELNG

      if ( line(i:i) .ne. ' ' .and. s .eq. -1 ) s = i      
      if ( line(i:i) .eq. ' ' .and. s .ne. -1 ) then 
         e = i-1
         read( line(s:e),*, iostat=ios ) val
         if ( ios .ne. 0 ) then 
            err =.true.
         else
            line(s:e) = "" ! wipe
         end if
         return
      end if
   end do
   
   val = 0
   err = .true. 

 end subroutine getint


! ---------------------------------------------------------------------


 subroutine getfloat(line, val, err)

   character(len=*) :: line
   logical :: err
   real(kind=8) :: val

   integer :: i ,s, e, ios

   s = -1
   e = -1

   if ( err ) return

   do i = 1, LINELNG

      if ( line(i:i) .ne. ' ' .and. s .eq. -1 ) s = i      
      if ( line(i:i) .eq. ' ' .and. s .ne. -1 ) then 
         e = i-1
         read( line(s:e),*, iostat=ios ) val
         if ( ios .ne. 0 ) then 
            err =.true.
         else
            line(s:e) = "" ! wipe
         end if
         return
      end if
   end do
   
   val = .0
   err = .true. 

 end subroutine getfloat


! ---------------------------------------------------------------------


  subroutine getcomplex(line, val, err)

   character(len=*) :: line
   logical :: err
   complex(kind=8) :: val

   integer :: i ,s, e, ios

   s = -1
   e = -1

   if ( err ) return
 

   do i = 1, LINELNG

      if ( line(i:i) .ne. ' ' .and. s .eq. -1 ) s = i      
      if ( line(i:i) .eq. ' ' .and. s .ne. -1 ) then 
         e = i-1
         read( line(s:e),*, iostat=ios ) val
         if ( ios .ne. 0 ) then 
            err =.true.
         else
            line(s:e) = "" ! wipe
         end if
         return
      end if
   end do
   
   val = .0
   err = .true. 

 end subroutine getcomplex


! ---------------------------------------------------------------------


 subroutine getstring(line, val, err)

   character(len=*) :: line
   logical :: err
   character(len=*) :: val
 
   integer :: i ,s, e, ios

   s = -1
   e = -1

   if ( err ) return

   do i = 1, LINELNG

      if ( line(i:i) .ne. ' ' .and. s .eq. -1 ) s = i      
      if ( line(i:i) .eq. ' ' .and. s .ne. -1 ) then 
         e = i-1
         read( line(s:e),*, iostat=ios ) val
         if ( ios .ne. 0 ) then 
            err =.true.
         else
            line(s:e) = "" ! wipe
         end if
         return
      end if
   end do
   
   val = ""
   err = .true. 

 end subroutine getstring

! ---------------------------------------------------------------------

 subroutine getints(line, val, num, err)

   character(len=*) :: line
   integer :: num, val(num)
   logical :: err

   integer :: i ,s, e, c, ios

   s = -1
   e = -1
   c = 0

   if ( err ) return

   do i = 1, LINELNG
      
      if ( line(i:i) .ne. ' ' .and. s .eq. -1 ) s = i      
      if ( line(i:i) .eq. ' ' .and. s .ne. -1 ) then 
         e = i-1
         c = c + 1
         read( line(s:e),*, iostat=ios ) val(c)
         if ( ios .ne. 0 ) then 
            err =.true.
         else
            line(s:e) = "" ! wipe
            s = -1
            e = -1
         end if
         if ( c .eq. num ) return
      end if
   end do
   
   err = .true. 

 end subroutine getints

! ---------------------------------------------------------------------

 subroutine getintvec(line, val, num, delim, err)

   character(len=*) :: line
   integer :: num, val(num)
   character :: delim
   logical :: err
   integer :: l,i ,s, e, c, ios
   
   s = -1
   e = -1
   c = 0

   if ( err ) return

   do i = 1, LINELNG

      if ( line(i:i) .eq. delim ) then 
         line(i:i) = " "
         return
      end if
      
      if ( line(i:i) .ne. ' ' .and. s .eq. -1 ) s = i      
      if ( line(i:i) .eq. ' ' .and. s .ne. -1 ) then 
         e = i-1
         c = c + 1
         if ( c .gt. num ) then
            err = .true. 
            return
         end if
         read( line(s:e),*, iostat=ios ) val(c)
         if ( ios .ne. 0 ) then 
            err =.true.
            return
         else
            line(s:e) = "" ! wipe
            s = -1
            e = -1
         end if
      end if
   end do
   
 end subroutine getintvec

! ---------------------------------------------------------------------

 subroutine getfloatvec(line, val, num, delim, err)

   character(len=*) :: line
   integer :: num
   real(kind=8) :: val(num)
   character :: delim
   logical :: err
   integer :: i ,s, e, c, ios
   
   s = -1
   e = -1
   c = 0

   if ( err ) return

   do i = 1, LINELNG

      if ( line(i:i) .eq. delim ) then 
         line(i:i) = " "
         return
      end if
      
      if ( line(i:i) .ne. ' ' .and. s .eq. -1 ) s = i      
      if ( line(i:i) .eq. ' ' .and. s .ne. -1 ) then 
         e = i-1
         c = c + 1
         if ( c .gt. num ) then
            err = .true. 
            return
         end if
         read( line(s:e),*, iostat=ios ) val(c)
         if ( ios .ne. 0 ) then 
            err =.true.
            return
         else
            line(s:e) = "" ! wipe
            s = -1
            e = -1
         end if
      end if
   end do
   
 end subroutine getfloatvec


! ---------------------------------------------------------------------

 subroutine getfloats(line, val, num, err)

   character(len=*) :: line
   integer :: num
   real(kind=8) :: val(num)
   logical :: err

   integer :: i ,s, e, c, ios

   s = -1
   e = -1
   c = 0

   if ( err ) return

   do i = 1, LINELNG
      
      if ( line(i:i) .ne. ' ' .and. s .eq. -1 ) s = i      
      if ( line(i:i) .eq. ' ' .and. s .ne. -1 ) then 
         e = i-1
         c = c + 1
         read( line(s:e),*, iostat=ios ) val(c)
         if ( ios .ne. 0 ) then 
            err =.true.
         else
            line(s:e) = "" ! wipe
            s = -1
            e = -1
         end if
         if ( c .eq. num ) return
      end if
   end do
   
   err = .true. 

 end subroutine getfloats

! ---------------------------------------------------------------------

 subroutine gettoken(line, token, err)

   character(len=*) :: token
   logical :: err
   character(len=*) :: line
   character(len=80) :: val
 
   integer :: i ,s, e, ios

   s = -1
   e = -1

   if ( err ) return

   do i = 1, LINELNG

      if ( line(i:i) .ne. ' ' .and. s .eq. -1 ) s = i      
      if ( line(i:i) .eq. ' ' .and. s .ne. -1 ) then 
         e = i-1
         read( line(s:e),*, iostat=ios ) val
         if ( ios .ne. 0 .or. val .ne. token ) then 
            err =.true.
         else
            line(s:e) = "" ! wipe
         end if
         return
      end if
   end do
   
   val = ""
   err = .true. 

 end subroutine gettoken

! ---------------------------------------------------------------------

  subroutine readtoken(unit, lcount, token)
    
    integer :: unit, lcount
    character(len=*) :: token
    logical :: err = .false.

    character(len=LINELNG) :: line
    logical :: eof

    call readline(unit, lcount, eof, line)
    if ( eof ) err = .true.
    call gettoken(line, token, err)
    if ( line .ne. "" ) err = .true.
    M4_SYNTAX_ERROR(err,lcount,token)

  end subroutine readtoken

! ---------------------------------------------------------------------

  subroutine readint(unit, lcount, val)
    
    integer :: unit, lcount
    integer :: val
    logical :: err = .false.

    character(len=LINELNG) :: line
    logical :: eof

    call readline(unit, lcount, eof, line)
    if ( eof ) err = .true.
    call getint(line, val, err)
    if ( line .ne. "" ) err = .true.
    M4_SYNTAX_ERROR(err,lcount,{"1 INTEGER"})

  end subroutine readint

! ---------------------------------------------------------------------

  subroutine readlogical(unit, lcount, val)
    
    integer :: unit, lcount
    logical :: val
    logical :: err = .false.

    character(len=LINELNG) :: line
    logical :: eof

    call readline(unit, lcount, eof, line)
    if ( eof ) err = .true.
    call getlogical(line, val, err)
    if ( line .ne. "" ) err = .true.
    M4_SYNTAX_ERROR(err,lcount,{"1 INTEGER"})

  end subroutine readlogical

! ---------------------------------------------------------------------

  subroutine readfloat(unit, lcount, val)
    
    integer :: unit, lcount
    real(kind=8) :: val
    logical :: err = .false.

    character(len=LINELNG) :: line
    logical :: eof

    call readline(unit, lcount, eof, line)
    if ( eof ) err = .true.
    call getfloat(line, val, err)
    if ( line .ne. "" ) err = .true.
    M4_SYNTAX_ERROR(err,lcount,{"1 FLOAT"})

  end subroutine readfloat

! ---------------------------------------------------------------------

  subroutine readints(unit, lcount, val, num)
    
    integer :: unit, lcount
    integer :: num
    integer :: val(num)
    logical :: err = .false.

    character(len=LINELNG) :: line
    logical :: eof

    call readline(unit, lcount, eof, line)
    if ( eof ) err = .true.
    call getints(line, val, num, err)
    if ( line .ne. "" ) err = .true.
    M4_SYNTAX_ERROR(err,lcount,{TRIM(i2str(num)), " INTEGERS"})

  end subroutine readints


! ---------------------------------------------------------------------

  subroutine readintvec(unit, lcount, val, num)
    
    integer :: unit, lcount
    integer :: num
    integer :: val(num)
    logical :: err = .false.

    character(len=LINELNG) :: line
    logical :: eof

    call readline(unit, lcount, eof, line)
    if ( eof ) err = .true.
    call getintvec(line, val, num, ":", err)
    if ( line .ne. "" ) err = .true.
    M4_SYNTAX_ERROR(err,lcount,{TRIM(i2str(num)), " INTEGERS"})

  end subroutine readintvec


! ---------------------------------------------------------------------

  subroutine readfloats(unit, lcount, val, num)
    
    integer :: unit, lcount
    integer :: num
    real(kind=8) :: val(num)
    logical :: err = .false.

    character(len=LINELNG) :: line
    logical :: eof

    call readline(unit, lcount, eof, line)
    if ( eof ) err = .true.
    call getfloats(line, val, num, err)
    if ( line .ne. "" ) err = .true.
    M4_SYNTAX_ERROR(err,lcount,{TRIM(i2str(num)), " FLOATS"})

  end subroutine readfloats


! ---------------------------------------------------------------------

  
  character(len=LINELNG) function ltrim(str) 

    character(len=*) :: str

    character :: c 
    integer :: i, n 
    logical :: copy = .false.
    
    ltrim = ""

    n = 1
    do i = 1, LINELNG
       c = str(i:i)
       if ( iachar(c) .eq. 9 ) c = ' '
       if ( c .ne. ' ' .or. copy ) then 
          ltrim(n:n) = c  
          n = n + 1
          copy = .true.
       end if
    end do

  end function ltrim

! ---------------------------------------------------------------------


end module parse

!
! Authors:  J.Hamm
! Modified: 13/02/2007
!
!======================================================================


