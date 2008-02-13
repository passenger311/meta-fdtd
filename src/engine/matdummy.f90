!-*- F90 -*------------------------------------------------------------
!
!  module: matdummy / meta
!
!  Dummy material module.
!
!  subs:
!
!    InitializeMatDummy
!    FinalizeMatDummy
!    ReadMatDummyObj
!    StepEMatDummy
!    StepHMatDummy
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatDummy module allows to add sources to the electromagnetic 
! field equations


module matdummy

  use constant
  use mpiworld
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MODHEAD_DECL({MATDUMMY},100,{

  ! enter mat data structure here:
  integer :: something

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatDummyObj(funit,lcount)

    M4_MODREAD_DECL({MATDUMMY}, funit,lcount,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatDummyObj")
    
    M4_MODREAD_EXPR({MATDUMMY},funit,lcount,mat,reg,0,out,{ 

    ! read mat parameters here, as defined in mat data structure
    read(funit,*) mat%something
 
    })

    M4_WRITE_DBG(". exit ReadMatDummyObj")

  end subroutine ReadMatDummyObj

!----------------------------------------------------------------------

  subroutine InitializeMatDummy

    M4_MODLOOP_DECL({MATDUMMY},mat) 
    M4_WRITE_DBG(". enter InitializeMatDummy")
    M4_MODLOOP_EXPR({MATDUMMY},mat,{
    
       ! initialize mat object here
       M4_IFELSE_DBG({call EchoMatDummyObj(mat)})
 
    })
    M4_WRITE_DBG(". exit InitializeMatDummy")

  end subroutine InitializeMatDummy

!----------------------------------------------------------------------

  subroutine FinalizeMatDummy

    M4_MODLOOP_DECL({MATDUMMY},mat)
    M4_WRITE_DBG(". enter FinalizeMatDummy")
    M4_MODLOOP_EXPR({MATDUMMY},mat,{

       ! finalize mat object here

    })
    M4_WRITE_DBG(". exit FinalizeMatDummy")

  end subroutine FinalizeMatDummy

!----------------------------------------------------------------------

  subroutine StepHMatDummy(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATDUMMY},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    M4_MODLOOP_EXPR({MATDUMMY},mat,{

       ! this loops over all mat structures, setting mat

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! this loops over all points of the region 

       })      
    })
  
  end subroutine StepHMatDummy


!----------------------------------------------------------------------


  subroutine StepEMatDummy(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATDUMMY},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    M4_MODLOOP_EXPR({MATDUMMY},mat,{

       ! this loops over all mat structures, setting mat

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! this loops over all points of the region 

       })      
    })

  end subroutine StepEMatDummy


!----------------------------------------------------------------------

   subroutine EchoMatDummyObj(mat)

    type(T_MATDUMMY) :: mat

    M4_WRITE_INFO({"--- matdummy # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"something = ",mat%something })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatDummyObj

  
!----------------------------------------------------------------------

end module matdummy

! =====================================================================


