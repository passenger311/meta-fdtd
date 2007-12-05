define(`M4_REGLOOP_DECL',`
 type(T_REG) :: $1
 integer :: $2,$3,$4,$5
 real(kind=8) :: $6
')

define(`M4_REGLOOP_EXPR',`
! --- START M4: REG_LOOP
! 1 = reg ( T_REG )
! 2 = p ( integer )
! 3 = i ( integer )
! 4 = j ( integer ) 
! 5 = k ( integer )
! 6 = w ( real*8 )
! 7 = expr

if ( $1%islist ) then
 
 do $2 = 1, $1%pe

  $3 = $1%i(p)
  $4 = $1%j(p)
  $5 = $1%k(p)
  $6 = $1%val(p)

  $7

 enddo

else 
 
 $2 = 0
 if ( $1%isbox ) then
  $6 = 1.0
  do $5 = $1%ks, $1%ke, $1%dk
   do $4 = $1%js, $1%je, $1%dj
    do $3 = $1%is, $1%ie, $1%di
     $2 = $2 + 1	     

     $7
	   
    enddo
   enddo
  enddo

 else 

  do $5 = $1%ks, $1%ke, $1%dk
   do $4 = $1%js, $1%je, $1%dj
    do $3 = $1%is, $1%ie, $1%di
     if ( $1%mask(i,j,k) .gt. 0 ) then
      $6 = $1%val($1%mask(i,j,k)) 
      $2 = $2 + 1	     

      $7
	   
     endif
    enddo
   enddo
  enddo

 endif

endif

! --- END M4: REG_LOOP

')