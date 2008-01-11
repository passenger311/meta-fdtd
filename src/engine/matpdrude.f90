!-*- F90 -*------------------------------------------------------------
!
!  module: matpdrude / meta
!
!  Pdrude material module.
!
!  subs:
!
!    InitializeMatPdrude
!    FinalizeMatPdrude
!    ReadMatPdrudeObj
!    StepEMatPdrude
!    StepHMatPdrude
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatPdrude module allows to add sources to the electromagnetic 
! field equations
!
! StepHMatPdrude: update eq. P(n+1) = c1 * P(n) + c2 * Pold(n) + c3 * E(n)
! StepEMatPdrude: update eq. E(n+1)* = E(n+1) - epsinv * (P(n+1) + P(n))


module matpdrude

  use constant
  use mpiworld
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MODHEAD_DECL({MATPDRUDE},100,{

  ! input parameters
  real(kind=8) :: omegapl, gamma

  ! coefficients
  real(kind=8) :: c1, c2, c3
  real(kind+8) :: Pxbuffer, Pybuffer, Pzbuffer
  

  ! current field: P
  M4_FTYPE, dimension(:), pointer :: Px, Py, Pz
  M4_FTYPE, dimension(:), pointer :: Pxold, Pyold, Pzold

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatPdrudeObj(funit)

    M4_MODREAD_DECL({MATPDRUDE}, funit,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatPdrudeObj")
    
    M4_MODREAD_EXPR({MATPDRUDE},funit,mat,reg,3,out,{ 

    ! read mat parameters here, as defined in mat data structure
    read(funit,*) mat%omegapl
    read(funit,*) mat%gamma
    
    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatPdrudeObj")

  end subroutine ReadMatPdrudeObj

!----------------------------------------------------------------------

  subroutine InitializeMatPdrude

    type (T_REG) :: reg
    integer :: err
    M4_MODLOOP_DECL({MATPDRUDE},mat) 
    M4_WRITE_DBG(". enter InitializeMatPdrude")
    M4_MODLOOP_EXPR({MATPDRUDE},mat,{
    
       ! initialize mat object here
       M4_IFELSE_DBG({call EchoMatPdrudeObj(mat)})

       mat%c1 = 4. / ( 2. + DT * mat%gamma )
       mat%c2 = ( -2. + DT * mat%gamma ) / ( 2. + DT * mat%gamma )
       mat%c3 = ( 2. * DT**2 * mat%omegapl**2) / ( 2. + DT * mat%gamma )

       mat%Pxbuffer = 0.
       mat%Pybuffer = 0.
       mat%Pzbuffer = 0.
       
       reg = regobj(mat%regidx)

       allocate(mat%Px(reg%numnodes),mat%Py(reg%numnodes),mat%Pz(reg%numnodes), &
            mat%Pxold(reg%numnodes),mat%Pyold(reg%numnodes),mat%Pzold(reg%numnodes), stat = err)
       M4_ALLOC_ERROR(err,"InitializeMatPdrude")

       mat%Px = 0.
       mat%Py = 0.
       mat%Pz = 0.
       mat%Pxold = 0.
       mat%Pyold = 0.
       mat%Pzold = 0.
 
    })
    M4_WRITE_DBG(". exit InitializeMatPdrude")

  end subroutine InitializeMatPdrude

!----------------------------------------------------------------------

  subroutine FinalizeMatPdrude

    M4_MODLOOP_DECL({MATPDRUDE},mat)
    M4_WRITE_DBG(". enter FinalizeMatPdrude")
    M4_MODLOOP_EXPR({MATPDRUDE},mat,{

    ! finalize mat object here
    deallocate(mat%Px,mat%Py,mat%Pz,mat%Pxold,mat%Pyold,mat%Pzold)

    })
    M4_WRITE_DBG(". exit FinalizeMatPdrude")

  end subroutine FinalizeMatPdrude

!----------------------------------------------------------------------

  subroutine StepHMatPdrude(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATPDRUDE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATPDRUDE},mat,{

       ! this loops over all mat structures, setting mat

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

       Pxbuffer = mat%Px(p)
       Pybuffer = mat%Py(p)
       Pzbuffer = mat%Pz(p)

       mat%Px(p) = mat%c1 * mat%Px(p) + mat%c2 * mat%Pxold(p) + mat%c3 * Ex(i,j,k)
       mat%Py(p) = mat%c1 * mat%Py(p) + mat%c2 * mat%Pyold(p) + mat%c3 * Ey(i,j,k)
       mat%Pz(p) = mat%c1 * mat%Pz(p) + mat%c2 * mat%Pyold(p) + mat%c3 * Ez(i,j,k)
       
       mat%Pxold(p) = Pxbuffer
       mat%Pyold(p) = Pybuffer
       mat%Pzold(p) = Pzbuffer

       ! this loops over all points of the region 

       })      
    })
  
  end subroutine StepHMatPdrude


!----------------------------------------------------------------------


  subroutine StepEMatPdrude(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATPDRUDE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATPDRUDE},mat,{

       ! this loops over all mat structures, setting mat

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! this loops over all points of the region 

       Ex(i,j,k) = Ex(i,j,k) - epsinvx(i,j,k) * ( mat%Px(p) - mat%Pxold(p) )
       Ey(i,j,k) = Ey(i,j,k) - epsinvy(i,j,k) * ( mat%Py(p) - mat%Pyold(p) )
       Ez(i,j,k) = Ez(i,j,k) - epsinvz(i,j,k) * ( mat%Pz(p) - mat%Pzold(p) )

       })      
    })

  end subroutine StepEMatPdrude


!----------------------------------------------------------------------

   subroutine EchoMatPdrudeObj(mat)

    type(T_MATPDRUDE) :: mat

    M4_WRITE_INFO({"--- matpdrude # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"omegapl = ",mat%omegapl })
    M4_WRITE_INFO({"gamma = ",mat%gamma })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatPdrudeObj

  
!----------------------------------------------------------------------

end module matpdrude

! =====================================================================


