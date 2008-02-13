!-*- F90 -*------------------------------------------------------------
!
!  module: diagdummy / meta
!
!  Dummy diagnostics module.
!
!  CF,2D
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

  subroutine ReadDiagDummyObj(funit,lcount)

    M4_MODREAD_DECL({DIAGDUMMY}, funit,lcount,diag,reg,out)

    M4_WRITE_DBG(". enter ReadMatDummyObj")
    
    M4_MODREAD_EXPR({DIAGDUMMY}, funit,diag,reg,0,out, {

    ! read diag parameters here, as defined in diag data structure
    read(funit,*) diag%something
 
    })

    M4_WRITE_DBG(". exit ReadMatDummyObj")

  end subroutine ReadDiagDummyObj

!----------------------------------------------------------------------

  subroutine InitializeDiagDummy

    M4_MODLOOP_DECL({DIAGDUMMY},diag)
    M4_WRITE_DBG(". enter InitializeMatDummy")
    M4_MODLOOP_EXPR({DIAGDUMMY},diag,{
    
       ! initialize diag object here

       M4_IFELSE_DBG({call EchoDiagDummyObj(diag)})

    })
    M4_WRITE_DBG(". exit InitializeMatDummy")

  end subroutine InitializeDiagDummy

!----------------------------------------------------------------------

  subroutine FinalizeDiagDummy

    M4_MODLOOP_DECL({DIAGDUMMY},diag)
    M4_WRITE_DBG(". enter FinalizeMatDummy")
    M4_MODLOOP_EXPR({DIAGDUMMY},diag,{

       ! finalize diag object here

    })
    M4_WRITE_DBG(". exit FinalizeMatDummy")

  end subroutine FinalizeDiagDummy

!----------------------------------------------------------------------

  subroutine StepHDiagDummy(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({DIAGDUMMY},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

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
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    M4_MODLOOP_EXPR({DIAGDUMMY},diag,{

       ! this loops over all diag structures, setting diag

       M4_MODOBJ_GETREG(diag,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! this loops over all points of the region 

       })      
    })

  end subroutine StepEDiagDummy

!----------------------------------------------------------------------

   subroutine EchoDiagDummyObj(diag)

    type(T_DIAGDUMMY) :: diag
 
    M4_WRITE_INFO({"--- diagdummy # ",&
         TRIM(i2str(diag%idx))," ", TRIM(diag%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"something = ",diag%something })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(diag%regidx))
    

  end subroutine EchoDiagDummyObj
  
!----------------------------------------------------------------------

end module diagdummy

! =====================================================================


