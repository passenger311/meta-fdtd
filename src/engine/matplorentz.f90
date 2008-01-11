!-*- F90 -*------------------------------------------------------------
!
!  module: matplorentz / meta
!
!  Plorentz material module.
!
!  subs:
!
!    InitializeMatPlorentz
!    FinalizeMatPlorentz
!    ReadMatPlorentzObj
!    StepEMatPlorentz
!    StepHMatPlorentz
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatPlorentz module allows to add sources to the electromagnetic 
! field equations
!
! StepHMatPlorentz: update eq. P(n+1) = c1 * P(n) + c2 * Pold(n) + c3 * E(n)
! StepEMatPlorentz: update eq. E(n+1)* = E(n+1) - epsinv * (P(n+1) - P(n))


module matplorentz

  use constant
  use mpiworld
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MODHEAD_DECL({MATPLORENTZ},100,{

  ! input parameters
  real(kind=8) :: omegal, gamma, deltaepsl

  ! coefficients
  real(kind=8) :: c1, c2, c3
  real(kind+8) :: Pxbuffer, Pybuffer, Pzbuffer
  

  ! current field: P
  M4_FTYPE, dimension(:), pointer :: Px, Py, Pz
  M4_FTYPE, dimension(:), pointer :: Pxold, Pyold, Pzold

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatPlorentzObj(funit)

    M4_MODREAD_DECL({MATPLORENTZ},funit,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatPlorentzObj")
    
    M4_MODREAD_EXPR({MATPLORENTZ},funit,mat,reg,3,out,{ 

    ! read mat parameters here, as defined in mat data structure
    read(funit,*) mat%omegal
    read(funit,*) mat%gamma
    read(funit,*) mat%deltaepsl
   
    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatPlorentzObj")

  end subroutine ReadMatPlorentzObj

!----------------------------------------------------------------------

  subroutine InitializeMatPlorentz

    type (T_REG) :: reg
    integer :: err
    M4_MODLOOP_DECL({MATPLORENTZ},mat) 
    M4_WRITE_DBG(". enter InitializeMatPlorentz")
    M4_MODLOOP_EXPR({MATPLORENTZ},mat,{
    
       ! initialize mat object here
       M4_IFELSE_DBG({call EchoMatPlorentzObj(mat)})

       mat%c1 = ( 2. - omegal**2 * DT**2 ) / ( 1. + DT * mat%gamma )
       mat%c2 = ( -1. + DT * mat%gamma ) / ( 1. + DT * mat%gamma )
       mat%c3 = DT**2 * mat%omegapl**2 * mat%deltaepsl / ( 1. + DT * mat%gamma )

       mat%Pxbuffer = 0.
       mat%Pybuffer = 0.
       mat%Pzbuffer = 0.
       
       reg = regobj(mat%regidx)

       allocate(mat%Px(reg%numnodes),mat%Py(reg%numnodes),mat%Pz(reg%numnodes), &
            mat%Pxold(reg%numnodes),mat%Pyold(reg%numnodes),mat%Pzold(reg%numnodes), stat = err)
       M4_ALLOC_ERROR(err,"InitializeMatPlorentz")

       mat%Px = 0.
       mat%Py = 0.
       mat%Pz = 0.
       mat%Pxold = 0.
       mat%Pyold = 0.
       mat%Pzold = 0.
 
    })
    M4_WRITE_DBG(". exit InitializeMatPlorentz")

  end subroutine InitializeMatPlorentz

!----------------------------------------------------------------------

  subroutine FinalizeMatPlorentz

    M4_MODLOOP_DECL({MATPLORENTZ},mat)
    M4_WRITE_DBG(". enter FinalizeMatPlorentz")
    M4_MODLOOP_EXPR({MATPLORENTZ},mat,{

    ! finalize mat object here
    deallocate(mat%Px,mat%Py,mat%Pz,mat%Pxold,mat%Pyold,mat%Pzold)

    })
    M4_WRITE_DBG(". exit FinalizeMatPlorentz")

  end subroutine FinalizeMatPlorentz

!----------------------------------------------------------------------

  subroutine StepHMatPlorentz(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATPLORENTZ},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATPLORENTZ},mat,{

       ! this loops over all mat structures, setting mat

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

       Pxbuffer = mat%Px(p)
       Pybuffer = mat%Py(p)
       Pzbuffer = mat%Pz(p)

       mat%Px(p) = mat%c1 * mat%Px(p) + mat%c2 * mat%Pxold(p) + mat%c3 * Ex(i,j,k)
       mat%Py(p) = mat%c1 * mat%Py(p) + mat%c2 * mat%Pyold(p) + mat%c3 * Ey(i,j,k)
       mat%Pz(p) = mat%c1 * mat%Pz(p) + mat%c2 * mat%Pzold(p) + mat%c3 * Ez(i,j,k)
       
       mat%Pxold(p) = Pxbuffer
       mat%Pyold(p) = Pybuffer
       mat%Pzold(p) = Pzbuffer

       ! this loops over all points of the region 

       })      
    })
  
  end subroutine StepHMatPlorentz


!----------------------------------------------------------------------


  subroutine StepEMatPlorentz(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATPLORENTZ},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATPLORENTZ},mat,{

       ! this loops over all mat structures, setting mat

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! this loops over all points of the region 

       Ex(i,j,k) = Ex(i,j,k) - epsinvx(i,j,k) * w(1) * ( mat%Px(p) - mat%Pxold(p) )
       Ey(i,j,k) = Ey(i,j,k) - epsinvy(i,j,k) * w(2) * ( mat%Py(p) - mat%Pyold(p) )
       Ez(i,j,k) = Ez(i,j,k) - epsinvz(i,j,k) * w(3) * ( mat%Pz(p) - mat%Pzold(p) )

       })      
    })

  end subroutine StepEMatPlorentz


!----------------------------------------------------------------------

   subroutine EchoMatPlorentzObj(mat)

    type(T_MATPLORENTZ) :: mat

    M4_WRITE_INFO({"--- matplorentz # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"omegapl = ",mat%omegapl })
    M4_WRITE_INFO({"gamma = ",mat%gamma })
    M4_WRITE_INFO({"deltaepsl = ",mat%deltaepsl })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatPlorentzObj

  
!----------------------------------------------------------------------

end module matplorentz

! =====================================================================


