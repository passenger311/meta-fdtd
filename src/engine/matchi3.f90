!-*- F90 -*------------------------------------------------------------
!
!  module: matchi3 / meta
!
!  Chi3 material module.
!
!  subs:
!
!    InitializeMatChi3
!    FinalizeMatChi3
!    ReadMatChi3Obj
!    StepEMatChi3
!    StepHMatChi3
!    SumJEKHMatChi3
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatChi3 module calculates the reponse of a Chi3 pole consisting of 
! a combined Raman/Kerr nonlinearity:
!
! d/dt d/dt S + 2 * gammar * d/dt S + omegar**2 S = omegar**2 * chi3r * E^2
!
! P = chi3k * E^3 + S * E
!
! define:
!
! E' = E + epsinv * P 
!
! and consider StepE and all previous StepEMat to solve for E'. E' is 
! E'(n) when StepEMatChi3 is entered and must be reverted to E(n) before
! entering StepH. StepHMatChi3 in turn must change E(n) back to E'(n) 
! for consistency.
!
! Discretization 
!
! StepEMatChi3 (at n)
! 
! A) solve iteratively:
! 
! E(n) = E'(n) / ( 1 + epsinv * S(n) + epsinv * chi3k * E(n)^2 )
!
! StepHMatChi3 (at n+1/2)
!
! A) similar to matlorentz do:
!
! S(n+1) = c1 * S(n) + c2 * S(n-1) + c3 * E(n)^2
!
! B) revert E(n) to E'(n):
!
! E'(n) = E(n) + epsinv * P(n)
!
! with P(n) = chi3k * E(n)^3 + S(n) * E(n)
!
!
! Note on Scaling:
!
! invlambdar scales like an inverse length, gammar as an inverse time
! chi3r + chi3k scale as follows:
!
! chi3 = ( dx[SI]^3 * c[SI]^2 / eps0[SI] ) * chi3(SI)
!


module matchi3

  use constant
  use parse
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MATHEAD_DECL({MATCHI3},MAXMATOBJ,{

  ! input parameters
  real(kind=8) :: lambdarinv    ! inv. vac. plasma wavelength and abs. length
  real(kind=8) :: gammar        ! resonance width
  real(kind=8) :: chi3r         ! chi3 for raman
  real(kind=8) :: chi3k         ! chi3 for kerr

  real(kind=8) :: omegar

  integer :: maxit

  ! coefficients
  real(kind=8) :: c1, c2, c3
  
  ! raman response field 
  M4_FTYPE, dimension(:,:), pointer :: Sx, Sy, Sz

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatChi3Obj(funit,lcount)

    M4_MODREAD_DECL({MATCHI3}, funit,lcount,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatChi3Obj")
    
    M4_MODREAD_EXPR({MATCHI3},funit,lcount,mat,reg,3,out,{ 

    ! read mat parameters here, as defined in mat data structure
    call readfloat(funit,lcount, mat%lambdarinv)
    call readfloat(funit,lcount, mat%gammar)
    call readfloat(funit,lcount, mat%chi3r)
    call readfloat(funit,lcount, mat%chi3k)
    call readint(funit,lcount, mat%maxit)

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatChi3Obj")

  end subroutine ReadMatChi3Obj

!----------------------------------------------------------------------

  subroutine InitializeMatChi3

    type (T_REG) :: reg
    integer :: err
    M4_MODLOOP_DECL({MATCHI3},mat) 
    M4_WRITE_DBG(". enter InitializeMatChi3")
    M4_MODLOOP_EXPR({MATCHI3},mat,{
    
       ! initialize mat object here

       mat%omegar = 2. * PI * mat%lambdarinv
!       mat%gammal = 2. / ( mat%abslenl * DT )

       reg = regobj(mat%regidx)

       allocate(mat%Sx(reg%numnodes,2),mat%Sy(reg%numnodes,2),mat%Sz(reg%numnodes,2), stat = err)
       M4_ALLOC_ERROR(err,"InitializeMatChi3")

       mat%Sx = 0.
       mat%Sy = 0.
       mat%Sz = 0.
       
       mat%c1 = ( 2. - mat%omegar**2 * DT**2 ) / ( 1. + DT * mat%gammar )
       mat%c2 = ( -1. + DT * mat%gammar ) / ( 1. + DT * mat%gammar )
       mat%c3 = DT**2 * mat%omegar**2 * mat%chi3r / ( 1. + DT * mat%gammar )

       M4_IFELSE_DBG({call EchoMatChi3Obj(mat)},{call DisplayMatChi3Obj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatChi3")

  end subroutine InitializeMatChi3

!----------------------------------------------------------------------

  subroutine FinalizeMatChi3

    M4_MODLOOP_DECL({MATCHI3},mat)
    M4_WRITE_DBG(". enter FinalizeMatChi3")
    M4_MODLOOP_EXPR({MATCHI3},mat,{

    ! finalize mat object here
    deallocate(mat%Sx,mat%Sy,mat%Sz)

    })
    M4_WRITE_DBG(". exit FinalizeMatChi3")

  end subroutine FinalizeMatChi3

!----------------------------------------------------------------------

  subroutine StepHMatChi3(ncyc)

    integer :: ncyc, m, n
    M4_FTYPE :: exh, eyh, ezh
    M4_MODLOOP_DECL({MATCHI3},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATCHI3},mat,{

    ! this loops over all mat structures, setting mat

      M4_MODOBJ_GETREG(mat,reg)

      n = mod(ncyc-1+2,2) + 1
      m = mod(ncyc+2,2) + 1
    
      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      exh = Ex(i,j,k)
      eyh = Ey(i,j,k)
      ezh = Ez(i,j,k)

      ! A) calculate S(n+1) from S(n),S(n-1) and E(n)

M4_IFELSE_TM({
      mat%Sx(p,m) = mat%c1 * mat%Sx(p,n) + mat%c2 * mat%Sx(p,m) + mat%c3 * exh**2
      mat%Sy(p,m) = mat%c1 * mat%Sy(p,n) + mat%c2 * mat%Sy(p,m) + mat%c3 * eyh**2
})
M4_IFELSE_TE({
      mat%Sz(p,m) = mat%c1 * mat%Sz(p,n) + mat%c2 * mat%Sz(p,m) + mat%c3 * ezh**2
})      

      ! NOTE: m and n will be flipped in the next timestep!

      ! B) revert E to E'

      ! E'(n) = E(n) + epsinv * P(n)
      ! with P(n) = chi3k * E(n)^3 + S(n) * E(n)

M4_IFELSE_TM({
      Ex(i,j,k) = exh + epsinvx(i,j,k) * ( mat%chi3k * exh**3 + mat%Sx(p,n) * exh )
      Ey(i,j,k) = eyh + epsinvy(i,j,k) * ( mat%chi3k * eyh**3 + mat%Sy(p,n) * eyh )
})
M4_IFELSE_TE({
      Ez(i,j,k) = ezh + epsinvz(i,j,k) * ( mat%chi3k * ezh**3 + mat%Sz(p,n) * ezh )
})      


      })      

    })
  
  end subroutine StepHMatChi3


!----------------------------------------------------------------------

  subroutine StepEMatChi3(ncyc)

    integer :: ncyc, m, n, it
    M4_FTYPE :: exh, eyh, ezh, exo, eyo, ezo, sxh, syh, szh
    real(kind=8) :: eix, eiy, eiz
    M4_MODLOOP_DECL({MATCHI3},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATCHI3},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
! A) solve iteratively:
! 
! E(n) = E'(n) / ( 1 + epsinv * ( S(n) + chi3k * E(n)^2 ) )

       exh = Ex(i,j,k)
       eyh = Ey(i,j,k)
       ezh = Ez(i,j,k)

       exo = exh
       eyo = eyh
       ezo = ezh

       sxh = mat%Sx(p,m)
       syh = mat%Sy(p,m)
       szh = mat%Sz(p,m)

       eix = epsinvx(i,j,k)
       eiy = epsinvy(i,j,k)
       eiz = epsinvz(i,j,k)

       do it = 1, mat%maxit

M4_IFELSE_TM({
         exh = exo / ( 1 + eix * ( sxh + mat%chi3k * exh**2 ) )
         eyh = eyo / ( 1 + eiy * ( syh + mat%chi3k * eyh**2 ) )
})
M4_IFELSE_TE({
         ezh = ezo / ( 1 + eiz * ( szh + mat%chi3k * ezh**2 ) )
})


       end do
       
       Ex(i,j,k) = exh
       Ey(i,j,k) = eyh
       Ez(i,j,k) = ezh

       })      
       


    })

  end subroutine StepEMatChi3

!----------------------------------------------------------------------

  real(kind=8) function SumJEMatChi3(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc
   
    SumJEMatChi3 = 0.
    
  end function SumJEMatChi3

!----------------------------------------------------------------------

  real(kind=8) function SumKHMatChi3(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc

    SumKHMatChi3 = 0.

  end function SumKHMatChi3
 
!----------------------------------------------------------------------

  subroutine DisplayMatChi3Obj(mat)

    type(T_MATCHI3) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" lambdarinv=",TRIM(f2str(mat%lambdarinv,5)),&
    	" gammar=",TRIM(f2str(mat%gammar,5)),&
    	" chi3r=",TRIM(f2str(mat%chi3r,5)),&
    	" chi3k=",TRIM(f2str(mat%chi3r,5))
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatChi3Obj

!----------------------------------------------------------------------

   subroutine EchoMatChi3Obj(mat)

    type(T_MATCHI3) :: mat

    M4_WRITE_INFO({"--- matchi3 # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"lambdarinv = ",mat%lambdarinv })
    M4_WRITE_INFO({"omegar = ",mat%omegar })
    M4_WRITE_INFO({"gammar = ",mat%gammar })
    M4_WRITE_INFO({"chi3r = ",mat%chi3r })
    M4_WRITE_INFO({"chi3k = ",mat%chi3k })
    M4_WRITE_INFO({"c1 = ",mat%c1 })
    M4_WRITE_INFO({"c2 = ",mat%c2 })
    M4_WRITE_INFO({"c3 = ",mat%c3 })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatChi3Obj

  
!----------------------------------------------------------------------

end module matchi3

! Authors:  J.Hamm 
! Modified: 1/10/2009
!
! =====================================================================


