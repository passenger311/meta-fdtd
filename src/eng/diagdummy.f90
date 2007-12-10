!-*- F90 -*------------------------------------------------------------
!
!  module: diagdummy / meta3
!
!  Dummy diagnostics module.
!
!  subs:
!
!    InitializeDiagDummy
!    FinalizeDiagDummy
!    ReadDiagDummyObj
!    StepEDiagDummy
!    StepHDiagDummy
!
!----------------------------------------------------------------------


! =====================================================================
!
! The DiagDummy module allows to add sources to the electromagnetic 
! field equations


module diagdummy

  use constant
  use mpiworld
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MODHEAD_DECL({DIAGDUMMY},100,{

  ! enter diag data structure here:
  integer :: something

  })

contains

!----------------------------------------------------------------------

  subroutine ReadDiagDummyObj(funit)

    M4_MODREAD_DECL({DIAGDUMMY}, funit,diag,reg,out)
    
    ! read diag parameters here, as defined in diag data structure
    read(funit,*) diag%something
 
    M4_MODREAD_END({DIAGDUMMY}, funit,diag,reg,out)

  end subroutine ReadDiagDummyObj

!----------------------------------------------------------------------

  subroutine InitializeDiagDummy

    M4_MODLOOP_DECL({DIAGDUMMY},diag)
    M4_MODLOOP_EXPR({DIAGDUMMY},diag,{
    
       ! initialize diag object here

    })

  end subroutine InitializeDiagDummy

!----------------------------------------------------------------------

  subroutine FinalizeDiagDummy

    M4_MODLOOP_DECL({DIAGDUMMY},diag)
    M4_MODLOOP_EXPR({DIAGDUMMY},diag,{

       ! finalize diag object here

    })

  end subroutine FinalizeDiagDummy

!----------------------------------------------------------------------

  subroutine StepHDiagDummy(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({DIAGDUMMY},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w)

    M4_MODLOOP_EXPR({DIAGDUMMY},diag,{

       ! this loops over all diag structures, setting diag

       M4_MODOBJ_GETREG(diag,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! this loops over all points of the region 

       })      
    })
  
  end subroutine StepHDiagDummy


!----------------------------------------------------------------------


  subroutine StepEDiagDummy(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({DIAGDUMMY},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w)

    M4_MODLOOP_EXPR({DIAGDUMMY},diag,{

       ! this loops over all diag structures, setting diag

       M4_MODOBJ_GETREG(diag,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! this loops over all points of the region 

       })      
    })

  end subroutine StepEDiagDummy

  
  
!----------------------------------------------------------------------

end module diagdummy

! =====================================================================


