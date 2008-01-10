!-*- F90 -*------------------------------------------------------------
!
!  module: matjdrude / meta
!
!  JDrude material module.
!
!  subs:
!
!    InitializeMatJDrude
!    FinalizeMatJDrude
!    ReadMatJDrudeObj
!    StepEMatJDrude
!    StepHMatJDrude
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatJDrude module allows to add sources to the electromagnetic 
! field equations
!
! StepHMatJDrude: update eq. J(n+1/2) = c1 * J(n-1/2) + c2 * E(n)
! StepEMatJDrude: update eq. E(n+1)* = E(n+1) + epsinv * DT * J(n+1/2)


module matjdrude

  use constant
  use mpiworld
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MODHEAD_DECL({MATJDRUDE},100,{

  ! input parameters
  real(kind=8) :: lambdapl, gamma

  real(kind=8) :: omegapl
  ! coefficients
  real(kind=8) :: c1, c2
  

  ! current field: J
  M4_FTYPE, dimension(:), pointer :: Jx, Jy, Jz

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatJDrudeObj(funit)

    M4_MODREAD_DECL({MATJDRUDE}, funit,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatJDrudeObj")
    
    M4_MODREAD_EXPR({MATJDRUDE},funit,mat,reg,3,out,{ 

    ! read mat parameters here, as defined in mat data structure
    read(funit,*) mat%lambdapl
    read(funit,*) mat%gamma
    
    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatJDrudeObj")

  end subroutine ReadMatJDrudeObj

!----------------------------------------------------------------------

  subroutine InitializeMatJDrude

    type (T_REG) :: reg
    integer :: err
    M4_MODLOOP_DECL({MATJDRUDE},mat) 
    M4_WRITE_DBG(". enter InitializeMatJDrude")
    M4_MODLOOP_EXPR({MATJDRUDE},mat,{
    
       ! initialize mat object here
       M4_IFELSE_DBG({call EchoMatJDrudeObj(mat)})

       mat%omegapl = 2. * PI * 1. / ( mat%lambdapl * DT )

       mat%c1 = ( 2. - DT * mat%gamma ) / ( 2. + DT * mat%gamma )
       mat%c2 = ( 2. * DT ) / ( 2. + DT * mat%gamma ) * mat%omegapl**2
       
       reg = regobj(mat%regidx)

       allocate(mat%Jx(reg%numnodes),mat%Jy(reg%numnodes),mat%Jz(reg%numnodes), stat = err)
       M4_ALLOC_ERROR(err,"InitializeMatJDrude")

       mat%Jx = 0.
       mat%Jy = 0.
       mat%Jz = 0.
 
    })
    M4_WRITE_DBG(". exit InitializeMatJDrude")

  end subroutine InitializeMatJDrude

!----------------------------------------------------------------------

  subroutine FinalizeMatJDrude

    M4_MODLOOP_DECL({MATJDRUDE},mat)
    M4_WRITE_DBG(". enter FinalizeMatJDrude")
    M4_MODLOOP_EXPR({MATJDRUDE},mat,{

    ! finalize mat object here
    deallocate(mat%Jx,mat%Jy,mat%Jz)

    })
    M4_WRITE_DBG(". exit FinalizeMatJDrude")

  end subroutine FinalizeMatJDrude

!----------------------------------------------------------------------

  subroutine StepHMatJDrude(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATJDRUDE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATJDRUDE},mat,{

       ! this loops over all mat structures, setting mat

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

       mat%Jx(p) = mat%c1 * mat%Jx(p) + mat%c2 * w(1) * Ex(i,j,k)
       mat%Jy(p) = mat%c1 * mat%Jy(p) + mat%c2 * w(2) * Ey(i,j,k)
       mat%Jz(p) = mat%c1 * mat%Jz(p) + mat%c2 * w(3) * Ez(i,j,k)
       
       ! this loops over all points of the region 

       })      
    })
  
  end subroutine StepHMatJDrude


!----------------------------------------------------------------------


  subroutine StepEMatJDrude(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATJDRUDE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATJDRUDE},mat,{

       ! this loops over all mat structures, setting mat

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! this loops over all points of the region 

       Ex(i,j,k) = Ex(i,j,k) - epsinvx(i,j,k) * DT * mat%Jx(p)
       Ey(i,j,k) = Ey(i,j,k) - epsinvy(i,j,k) * DT * mat%Jy(p)
       Ez(i,j,k) = Ez(i,j,k) - epsinvz(i,j,k) * DT * mat%Jz(p)

       })      
    })

  end subroutine StepEMatJDrude


!----------------------------------------------------------------------

   subroutine EchoMatJDrudeObj(mat)

    type(T_MATJDRUDE) :: mat

    M4_WRITE_INFO({"--- matjdrude # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"omegapl = ",mat%omegapl })
    M4_WRITE_INFO({"gamma = ",mat%gamma })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatJDrudeObj

  
!----------------------------------------------------------------------

end module matjdrude

! =====================================================================


