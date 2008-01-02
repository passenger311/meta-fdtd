!-*- F90 -*------------------------------------------------------------
!
!  module: mat2lvl / meta3
!
!  Dummy material module.
!
!  subs:
!
!    InitializeMat2lvl
!    FinalizeMat2lvl
!    ReadMat2lvlObj
!    StepEMat2lvl
!    StepHMat2lvl
!
!----------------------------------------------------------------------


! =====================================================================
!
! The Mat2lvl module allows to add sources to the electromagnetic 
! field equations


module mat2lvl

  use constant
  use mpiworld
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MODHEAD_DECL({MAT2LVL},100,{

  ! enter mat data structure here:
  integer :: something
  M4_FTYPE, pointer, dimension(:) :: n

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMat2lvlObj(funit)

    M4_MODREAD_DECL({MAT2LVL}, funit,mat,reg,out)

    M4_WRITE_DBG({". enter ReadMat2lvlObj"})
    
    M4_MODREAD_EXPR({MAT2LVL}, funit,mat,reg,out,{ 

    ! read mat parameters here, as defined in mat data structure
    read(funit,*) mat%something
    M4_WRITE_DBG({"somthing = ",mat%something})

    })

    M4_WRITE_DBG({". exit ReadMat2lvlObj"})

  end subroutine ReadMat2lvlObj

!----------------------------------------------------------------------

  subroutine InitializeMat2lvl

    integer :: err
    type(T_REG) :: reg
    M4_MODLOOP_DECL({MAT2LVL},mat) 

    M4_WRITE_DBG(". enter InitializeMat2lvl")
    M4_MODLOOP_EXPR({MAT2LVL},mat,{
    
       ! initialize mat object here
       M4_IFELSE_DBG({call EchoMat2lvlObj(mat)})

       M4_MODOBJ_GETREG(mat,reg)
       allocate(mat%n(1:reg%numnodes),stat = err)
       M4_ALLOC_ERROR(err,"ReadMat2lvlObj")

       mat%n = mat%something
 
    })
    M4_WRITE_DBG(". exit InitializeMat2lvl")

  end subroutine InitializeMat2lvl

!----------------------------------------------------------------------

  subroutine FinalizeMat2lvl

    M4_MODLOOP_DECL({MAT2LVL},mat)
    M4_WRITE_DBG(". enter FinalizeMat2lvl")
    M4_MODLOOP_EXPR({MAT2LVL},mat,{

       ! finalize mat object here
    deallocate(mat%n)



    })
    M4_WRITE_DBG(". exit FinalizeMat2lvl")

  end subroutine FinalizeMat2lvl

!----------------------------------------------------------------------

  subroutine StepHMat2lvl(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MAT2LVL},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w)

    M4_MODLOOP_EXPR({MAT2LVL},mat,{

       ! this loops over all mat structures, setting mat

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
           
       ! this loops over all points of the region 
       mat%n(p) = w* Ex(i,j,k)

       })      
    })
  
  end subroutine StepHMat2lvl


!----------------------------------------------------------------------


  subroutine StepEMat2lvl(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MAT2LVL},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w)

    M4_MODLOOP_EXPR({MAT2LVL},mat,{

       ! this loops over all mat structures, setting mat

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! this loops over all points of the region 

       })      
    })

  end subroutine StepEMat2lvl


!----------------------------------------------------------------------

   subroutine EchoMat2lvlObj(mat)

    type(T_MAT2LVL) :: mat

    M4_WRITE_INFO({"--- mat2lvl # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"something = ",mat%something })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMat2lvlObj

  
!----------------------------------------------------------------------

end module mat2lvl

! =====================================================================


