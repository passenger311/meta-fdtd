!-*- F90 -*------------------------------------------------------------
!
!  module: matlorentz / meta
!
!  JLorentz material module.
!
!  subs:
!
!    InitializeMatlorentz
!    FinalizeMatlorentz
!    ReadMatlorentzObj
!    StepEMatlorentz
!    StepHMatlorentz
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatLorentz module calculates the reponse of a Lorentz pole
!
! d/dt d/dt P + 2 * gammal * d/dt P + omegal**2 P = deltaepsl * omegal**2 * E
! d/dt E = d/dt E* - epsinv * d/dt P 
!
! where E* is the electric field as calculated without the sources.  
!
! 2. order method
!
! StepHMatPlorentz: update eq. P(n+1) = c1 * P(n) + c2 * P(n-1) + c3 * E(n)
! StepEMatPlorentz: update eq. E(n+1) = E(n+1)* - epsinv * (P(n+1) - P(n))
!

module matlorentz

  use constant
  use mpiworld
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MODHEAD_DECL({MATLORENTZ},100,{

  ! input parameters
  real(kind=8) :: lambdal    ! vac. plasma wavelength and abs. length
  real(kind=8) :: gammal     ! resonance width
  real(kind=8) :: deltaepsl  ! delta epsilon

  real(kind=8) :: omegal

  ! coefficients
  real(kind=8) :: c1, c2, c3
  
  ! polarisation field 
  M4_FTYPE, dimension(:,:), pointer :: Px, Py, Pz

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatLorentzObj(funit)

    M4_MODREAD_DECL({MATLORENTZ}, funit,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatLorentzObj")
    
    M4_MODREAD_EXPR({MATLORENTZ},funit,mat,reg,3,out,{ 

    ! read mat parameters here, as defined in mat data structure
    read(funit,*) mat%lambdal
    read(funit,*) mat%gammal
    read(funit,*) mat%deltaepsl

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatLorentzObj")

  end subroutine ReadMatLorentzObj

!----------------------------------------------------------------------

  subroutine InitializeMatLorentz

    type (T_REG) :: reg
    integer :: err
    M4_MODLOOP_DECL({MATLORENTZ},mat) 
    M4_WRITE_DBG(". enter InitializeMatLorentz")
    M4_MODLOOP_EXPR({MATLORENTZ},mat,{
    
       ! initialize mat object here

       mat%omegal = 2. * PI * 1. / ( mat%lambdal * DT )
!       mat%gammal = 2. / ( mat%abslenl * DT )

       reg = regobj(mat%regidx)

       allocate(mat%Px(reg%numnodes,2),mat%Py(reg%numnodes,2),mat%Pz(reg%numnodes,2), stat = err)
       M4_ALLOC_ERROR(err,"InitializeMatLorentz")

       mat%Px = 0.
       mat%Py = 0.
       mat%Pz = 0.

       mat%c1 = ( 2. - mat%omegal**2 * DT**2 ) / ( 1. + DT * mat%gammal )
       mat%c2 = ( -1. + DT * mat%gammal ) / ( 1. + DT * mat%gammal )
       mat%c3 = DT**2 * mat%omegal**2 * mat%deltaepsl / ( 1. + DT * mat%gammal )

       M4_IFELSE_DBG({call EchoMatLorentzObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatLorentz")

  end subroutine InitializeMatLorentz

!----------------------------------------------------------------------

  subroutine FinalizeMatLorentz

    M4_MODLOOP_DECL({MATLORENTZ},mat)
    M4_WRITE_DBG(". enter FinalizeMatLorentz")
    M4_MODLOOP_EXPR({MATLORENTZ},mat,{

    ! finalize mat object here
    deallocate(mat%Px,mat%Py,mat%Pz)

    })
    M4_WRITE_DBG(". exit FinalizeMatLorentz")

  end subroutine FinalizeMatLorentz

!----------------------------------------------------------------------

  subroutine StepHMatLorentz(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATLORENTZ},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATLORENTZ},mat,{

    ! this loops over all mat structures, setting mat

      M4_MODOBJ_GETREG(mat,reg)

      n = mod(ncyc-1,2) + 1
      m = mod(ncyc,2) + 1
    
      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      ! calculate P(n+1) from P(n),P(n-1) and E(n)

      ! before: J(*,m) is P(n-1), J(*,n) is P(n)

      mat%Px(p,m) = mat%c1 * mat%Px(p,n) + mat%c2 * mat%Px(p,m) + mat%c3 * Ex(i,j,k)
      mat%Py(p,m) = mat%c1 * mat%Py(p,n) + mat%c2 * mat%Py(p,m) + mat%c3 * Ey(i,j,k)
      mat%Pz(p,m) = mat%c1 * mat%Pz(p,n) + mat%c2 * mat%Pz(p,m) + mat%c3 * Ez(i,j,k)
      
      ! after: J(*,m) is now P(n+1)
       
      ! m and n will be flipped in the next timestep!

      })      

    })
  
  end subroutine StepHMatLorentz


!----------------------------------------------------------------------


  subroutine StepEMatLorentz(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATLORENTZ},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATLORENTZ},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1,2) + 1
       m = mod(ncyc,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

       Ex(i,j,k) = Ex(i,j,k) - w(1) * epsinvx(i,j,k) * ( mat%Px(p,m) - mat%Px(p,n) )
       Ey(i,j,k) = Ey(i,j,k) - w(2) * epsinvy(i,j,k) * ( mat%Py(p,m) - mat%Py(p,n) )
       Ez(i,j,k) = Ez(i,j,k) - w(3) * epsinvz(i,j,k) * ( mat%Pz(p,m) - mat%Pz(p,n) )

       })      

    })

  end subroutine StepEMatLorentz


!----------------------------------------------------------------------

   subroutine EchoMatLorentzObj(mat)

    type(T_MATLORENTZ) :: mat

    M4_WRITE_INFO({"--- matlorentz # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"lambdal = ",mat%lambdal })
    M4_WRITE_INFO({"omegal = ",mat%omegal })
    M4_WRITE_INFO({"deltaepsl = ",mat%deltaepsl })
    M4_WRITE_INFO({"gammal = ",mat%gammal })
    M4_WRITE_INFO({"c1 = ",mat%c1 })
    M4_WRITE_INFO({"c2 = ",mat%c2 })
    M4_WRITE_INFO({"c3 = ",mat%c3 })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatLorentzObj

  
!----------------------------------------------------------------------

end module matlorentz

! Authors:  K.Boehringer, J.Hamm 
! Modified: 14/1/2008
!
! =====================================================================


