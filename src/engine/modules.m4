define({M4_MODHEAD_DECL},{

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = '$1'

  ! --- Public Methods

  public :: Initialize$1
  public :: Finalize$1
  public :: StepE$1
  public :: StepH$1
  public :: Read$1Obj
  public :: Echo$1Obj

  ! --- Public Data

  public :: T_$1
  public :: $1obj
  public :: num$1obj

  ! --- Constants

  integer, parameter :: MAX$1OBJ = $2

  ! --- Types

  type T_$1

     character(len=STRLNG) :: type = "$1"
     integer :: regidx, idx       		! regobj index
     $3	

  end type T_$1

  ! --- Fields

  type(T_$1) :: $1obj(MAX$1OBJ) 
  integer :: num$1obj

})


define({M4_MATHEAD_DECL},{
public :: SumJE$1
public :: SumKH$1
M4_MODHEAD_DECL({$1},{$2},{$3})})

define({M4_SRCHEAD_DECL},{
public :: SumJE$1
public :: SumKH$1
M4_MODHEAD_DECL({$1},{$2},{$3})})

define({M4_DIAGHEAD_DECL},{
M4_MODHEAD_DECL({$1},{$2},{$3})})

define({M4_LUMPEDHEAD_DECL},{
M4_MODHEAD_DECL({$1},{$2},{$3})})

define({M4_MODLOOP_DECL},{
    integer :: m4_n
    type(T_$1) :: $2   
})

define({M4_MODLOOP_EXPR},{
  do m4_n = 1, num$1obj

       $2 = $1obj(m4_n)
       
       $3

       $1obj(m4_n) = $2

  end do
})

define({M4_MODOBJ_GETREG},{
  $2 = regobj($1%regidx)
})


define({M4_MODREAD_DECL},{
    integer:: $2,$3 ! funit, lcount	
    type(T_$1) :: $4 
    type(T_REG) :: $5 ! reg
    type(T_OUT) :: $6 ! out
    character(len=LINELNG) :: m4_string
    logical :: m4_eof

})



define({M4_MODREAD_EXPR},{

 num$1obj = num$1obj + 1
 $4 = $1obj(num$1obj)
 $4%idx = num$1obj	

 $8

! read regobj information

    call readline($2,$3, m4_eof, m4_string)
    M4_EOF_ERROR({m4_eof},$3)
    M4_WRITE_DBG({"got token: ",TRIM(m4_string)})
    M4_SYNTAX_ERROR({m4_string .ne. "(REG"},$3,{"(REG"})
    M4_WRITE_DBG({"got token (REG -> ReadRegObj"})
    call ReadRegObj($5, fdtdreg, $2, $3, $6)
    $4%regidx = reg%idx

! get terminator

    do	
      call readline($2,$3, m4_eof, m4_string)
      M4_EOF_ERROR({m4_eof},$3)
      select case (m4_string)
      case("(OUT") 
	M4_WRITE_DBG({"got token (OUT -> ReadOutObj"})
       call ReadOutObj($7, $5, $2, $3, modname)
       outobj($7%idx)%objidx = num$1obj
       $4%regidx = $5%idx
      case(")$1")
	exit
      case default
	M4_SYNTAX_ERROR(.true.,lcount,{")$1"})
      end select
    enddo


! write back
    $1obj(num$1obj) = $4		
})
