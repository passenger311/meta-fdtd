
define({M4_REGLOOP_DECL},{
 type(T_REG) :: $1
 integer :: $2,$3,$4,$5
 real(kind=8) :: $6
})

define({M4_REGLOOP_EXPR},{
! --- START M4: REGLOOP_EXPR
! reg,p,i,j,k,w,expr,spac1, spac2

if ( $1%numnodes .gt. 0 ) then

if ( $1%isbox ) then
$2 = 0
do $5 = $1%ks, $1%ke, $1%dk
do $4 = $1%js, $1%je, $1%dj
do $3 = $1%is, $1%ie, $1%di
$2 = $2 + 1	     

if ( $1%compressval ) then 
$6(:) = $1%val(:,$1%valptr($2))
else
$6(:) = $1%val(:,$2)
endif

$7
	   
enddo !i
$8
enddo !j
M4_IFELSE_3D({$9})
enddo !k

else !isbox

if ( $1%islist ) then ! list loop mode
 
do $2 = 1, $1%pe ! -> p

$3 = $1%i($2) ! -> i
$4 = $1%j($2) ! -> j
$5 = $1%k($2) ! -> k

if ( $1%compressval ) then 
$6(:) = $1%val(:,$1%valptr($2))
else
$6(:) = $1%val(:,$2)
endif

$7

enddo ! p

else ! mask loop mode

do $5 = $1%ks, $1%ke, $1%dk ! -> k
do $4 = $1%js, $1%je, $1%dj ! -> j
do $3 = $1%is, $1%ie, $1%di ! -> i
$2 = $1%mask($3,$4,$5) ! -> p
if ( $2 .gt. 0 ) then

if ( $1%compressval ) then 
$6(:) = $1%val(:,$1%valptr($2))
else
$6(:) = $1%val(:,$2)
endif

$7
	   
endif 
enddo ! i
enddo ! j
enddo ! k

endif !islist

endif !isbox

endif !numnodes > 0

! --- END M4: REGLOOP_EXPR
})

define({M4_REGLOOP_WRITE},{
M4_REGLOOP_EXPR($1,$2,$3,$4,$5,$6,{
if ( $1%isbox ) then
write($7,"($8)") $9
else
write($7,"(M4_SDIM({I5}),$8)") M4_DIM123({$3},{$3,$4},{$3,$4,$5}),$9
endif
},{if ( $1%is .ne. $1%ie ) write($7,*)}, {if ( $1%js .ne. $1%je ) write($7,*)} )
})


define({M4_DOMREGLOOP_EXPR},{

! --- START M4: DOMREGLOOP_EXPR
! dom,reg,p,i,j,k,w,expr

if ( $2%numnodes .gt. 0 ) then

if ( $2%isbox ) then
$3 = 0
do $6 = max($1%ks,$2%ks), min($1%ke,$2%ke), $2%dk
do $5 = max($1%js,$2%js), min($1%je,$2%je), $2%dj
do $4 = max($1%is,$2%is), min($1%ie,$2%ie), $2%di

if ( $1%mask($4,$5,$6) .gt. 0 ) then

$3 = $3 + 1	     

if ( $2%compressval ) then 
$7(:) = $2%val(:,$2%valptr($3))
else
$7(:) = $2%val(:,$3)
endif

$8

endif
	   
enddo !i
enddo !j
enddo !k

else !isbox

do $6 = max($1%ks,$2%ks), min($1%ke,$2%ke), $2%dk ! -> k
do $5 = max($1%ks,$2%js), min($1%je,$2%je), $2%dj ! -> j
do $4 = max($1%ks,$2%is), min($1%ie,$2%ie), $2%di ! -> i
$3 = $2%mask($4,$5,$6) ! -> p
if ( $3 .gt. 0 .and. $1%mask($4,$5,$6) .gt. 0 ) then

if ( $2%compressval ) then 
$7(:) = $2%val(:,$2%valptr($3))
else
$7(:) = $2%val(:,$3)
endif

$8
	   
endif 
enddo ! i
enddo ! j
enddo ! k

endif !isbox

endif !numnodes > 0

! --- END M4: DOMREGLOOP_EXPR
})