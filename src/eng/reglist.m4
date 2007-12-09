
define({M4_GET_REG_AND_TERMINATOR},{
! read regobj information

    read(funit,*) string
    if ( string .eq. "(REG" ) then
       M4_WRITE_DBG({"got (REG -> ReadRegObj"})
       call ReadRegObj(reg, funit)
       $1%regidx = reg%idx
    else
       M4_FATAL_ERROR({"NO REGION DEFINED: ",$2})
    end if

! get terminator

    read(funit,*) string
    M4_WRITE_DBG({"read terminator: ", TRIM(string)})

    if ( string .ne. ")" ) then
       M4_FATAL_ERROR({"BAD TERMINATOR: ",$2})
    end if
})

define({M4_REGLOOP_DECL},{
 type(T_REG) :: $1
 integer :: $2,$3,$4,$5
 real(kind=8) :: $6
})

define({M4_REGLOOP_EXPR},{
! --- START M4: REGLOOP_EXPR
! reg,p,i,j,k,w,expr,spac1, spac2

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
    $8
   enddo
   $9
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

! --- END M4: REGLOOP_EXPR

})

define({M4_REGLOOP_WRITE},{
M4_REGLOOP_EXPR($1,$2,$3,$4,$5,$6,{
if ( $1%isbox ) then
write($7,"($8)") $9
else
write($7,"(I5,I5,I5,$8)") $3,$4,$5,$9
endif
},{if ( $1%is .ne. $1%ie ) write($7,*)}, {if ( $1%js .ne. $1%je ) write($7,*)} )
})