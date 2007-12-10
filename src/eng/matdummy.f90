!-*- F90 -*------------------------------------------------------------
!
!  module: matdummy / meta3
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

  subroutine ReadMatDummyObj(funit)

    M4_MODREAD_DECL({MATDUMMY}, funit,mat,reg,out)
    
    ! read mat parameters here, as defined in mat data structure
    read(funit,*) mat%something
 
    M4_MODREAD_END({MATDUMMY}, funit,mat,reg,out)

  end subroutine ReadMatDummyObj

!----------------------------------------------------------------------

  subroutine InitializeMatDummy

    M4_MODLOOP_DECL({MATDUMMY},mat)
    M4_MODLOOP_EXPR({MATDUMMY},mat,{
    
       ! initialize mat object here

    })

  end subroutine InitializeMatDummy

!----------------------------------------------------------------------

  subroutine FinalizeMatDummy

    M4_MODLOOP_DECL({MATDUMMY},mat)
    M4_MODLOOP_EXPR({MATDUMMY},mat,{

       ! finalize mat object here

    })

  end subroutine FinalizeMatDummy

!----------------------------------------------------------------------

  subroutine StepHMatDummy(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATDUMMY},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w)

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
    M4_REGLOOP_DECL(reg,p,i,j,k,w)

    M4_MODLOOP_EXPR({MATDUMMY},mat,{

       ! this loops over all mat structures, setting mat

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! this loops over all points of the region 

       })      
    })

  end subroutine StepEMatDummy

  
  
!----------------------------------------------------------------------

end module matdummy

! =====================================================================


