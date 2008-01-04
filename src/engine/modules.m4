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
    integer:: $2 ! funit	
    type(T_$1) :: $3 
    type(T_REG) :: $4 ! reg
    type(T_OUT) :: $5 ! out
    character(len=STRLNG) :: m4_string, m4_linestr, m4_skiptill

})



define({M4_MODREAD_EXPR},{

 num$1obj = num$1obj + 1
 $3 = $1obj(num$1obj)
 $3%idx = num$1obj	

 $7

! read regobj information

    read($2,*) m4_string
    if ( m4_string .eq. "(REG" ) then
       M4_WRITE_DBG({"got token (REG -> ReadRegObj"})
       call ReadRegObj($4, fdtdreg, $2, $5)
       $3%regidx = reg%idx
    else
       M4_FATAL_ERROR({"NO REGION DEFINED: Read$1Obj"})
    end if

! get terminator

    m4_skiptill = ""
    do	
      read($2,*) m4_linestr
      m4_string = TRIM(ADJUSTL(m4_linestr))

      if ( m4_skiptill .ne. "" ) then 
         M4_WRITE_DBG({"skipping line ",TRIM(m4_string)})
         if ( m4_string .eq. m4_skiptill ) m4_skiptill = ""  
         cycle              
      endif

      select case (m4_string)
      case("(OUT") 
	M4_WRITE_DBG({"got token (OUT -> ReadOutObj"})
       call ReadOutObj($6, $4, modname, $2)
       $3%regidx = $4%idx
      case default	
	if ( m4_string(1:2) .eq. "(!" ) then
           m4_skiptill = cat2(")",m4_string(3:))
           M4_WRITE_DBG({"got token (! -> skiptill = ", TRIM(m4_skiptill)})  
           cycle
        end if

        M4_WRITE_DBG({"read terminator: ", TRIM(m4_string)})
        if ( m4_string(1:1) .ne. ")" ) then
          M4_FATAL_ERROR({"BAD TERMINATOR: Read$1Obj/$1"})
        end if
        exit
      end select
    enddo


! write back
    $1obj(num$1obj) = $3		
})
