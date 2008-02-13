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
 
  use strings

  implicit none
  save

contains


! ---------------------------------------------------------------------


  subroutine readline(unit, rline, lcount, eof, skip, skiptill)
    
    integer :: unit, lcount
    logical :: skip, eof
    character(len=255) :: rline, skiptill

    integer :: ios, i
    logical :: wipe
    character(len=255) :: line

    line = ""

    if ( skiptill .eq. "" ) then
       skip = .false.
    end if

    read(unit,"(A255)",iostat=ios) line

    if ( ios .ne. 0 ) then
       eof = .true.
    else
       eof = .false. 
       lcount = lcount + 1
    end if

    rline =  TRIM(ADJUSTL(line))

    if ( rline(1:2) .eq. "(!" ) then
       skiptill = cat2(")",rline(3:))
       skip = .true.
    end if
    
    if ( skip ) then 
       if ( rline .eq. skiptill ) skiptill = ""  
    else
       ! cut comment
       wipe = .false.
       do i = 1, len(rline)
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

    if ( .not. skip .and. rline .eq. "" ) then
       skip = .true.
       skiptill = ""
    end if


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

   do i = 1, len(line)

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

   do i = 1, len(line)

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

   do i = 1, len(line)

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
 

   do i = 1, len(line)

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

   do i = 1, len(line)

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

   do i = 1, len(line)
      
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
   integer :: i ,s, e, c, ios
   
   val = 0

   s = -1
   e = -1
   c = 0

   if ( err ) return

   do i = 1, len(line)

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
         else
            line(s:e) = "" ! wipe
            s = -1
            e = -1
         end if
      end if
   end do
   
   err = .true. 

 end subroutine getintvec

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

   do i = 1, len(line)
      
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

 subroutine expectstring(line, str, err)

  character(len=*) :: line
  logical :: err
  character(len=*) :: str
  character(len=80) :: val
 
   integer :: i ,s, e, ios

   s = -1
   e = -1

   if ( err ) return

   do i = 1, len(line)

      if ( line(i:i) .ne. ' ' .and. s .eq. -1 ) s = i      
      if ( line(i:i) .eq. ' ' .and. s .ne. -1 ) then 
         e = i-1
         read( line(s:e),*, iostat=ios ) val
         if ( ios .ne. 0 .or. val .ne. str ) then 
            err =.true.
         else
            line(s:e) = "" ! wipe
         end if
         return
      end if
   end do
   
   val = ""
   err = .true. 

 end subroutine expectstring


! ---------------------------------------------------------------------




end module parse

!
! Authors:  J.Hamm
! Modified: 13/02/2007
!
!======================================================================


