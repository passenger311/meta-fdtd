!-*- F90 -*------------------------------------------------------------
!
!  module: matkerr / meta
!
!  Kerr material module.
!
!  subs:
!
!    InitializeMatKerr
!    FinalizeMatKerr
!    ReadMatKerrObj
!    StepEMatKerr
!    StepHMatKerr
!    SumJEMatKerr
!    SumKHMatKerr
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatKerr module calculates the reponse of a Kerr nonlinearity resulting
! in an isotropic permittivity
!
! P = chi3k * |E|^2 E
!
! define:
!
! E' = E + epsinv * P 
!
! and consider StepE and all previous StepEMat to solve for E'. E' is 
! E'(n) when StepEMatKerr is entered and must be reverted to E(n) before
! entering StepH. StepHMatKerr in turn must change E(n) back to E'(n) 
! for consistency.
!
! Discretization 
!
! StepEMatKerr (at n)
! 
! A) solve iteratively:
! 
! E(n) = E'(n) / ( 1 + epsinv * chi3k * E(n)^2 )
!
! StepHMatKerr (at n+1/2)
!
! A) revert E(n) to E'(n):
!
! E'(n) = E(n) + epsinv * P(n)
!
! with P(n) = chi3k * |E(n)|^2 E(n)
!
!
! Note on Scaling:
!
! invlambdar scales like an inverse length, gammar as an inverse time
! chi3k scale as follows:
!
! chi3k = c[SI]^2 / ( dx[SI]^3 eps0[SI] ) * chi3(SI)
!


module matkerr

  use constant
  use parse
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MATHEAD_DECL({MATKERR},MAXMATOBJ,{

  ! input parameters
  real(kind=8) :: chi3k         ! chi3 for kerr

  integer :: maxit
  
  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatKerrObj(funit,lcount)

    M4_MODREAD_DECL({MATKERR}, funit,lcount,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatKerrObj")
    
    M4_MODREAD_EXPR({MATKERR},funit,lcount,mat,reg,3,out,{ 

    ! read mat parameters here, as defined in mat data structure
    call readfloat(funit,lcount, mat%chi3k)
    call readint(funit,lcount, mat%maxit)
    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatKerrObj")

  end subroutine ReadMatKerrObj

!----------------------------------------------------------------------

  subroutine InitializeMatKerr

    type (T_REG) :: reg
    integer :: err 
    M4_WRITE_DBG(". enter InitializeMatKerr")

! No initialization necessary since there is no memory

    M4_WRITE_DBG(". exit InitializeMatKerr")

  end subroutine InitializeMatKerr

!----------------------------------------------------------------------

  subroutine FinalizeMatKerr

    M4_WRITE_DBG(". enter FinalizeMatKerr")
    M4_WRITE_DBG(". exit FinalizeMatKerr")

  end subroutine FinalizeMatKerr

!----------------------------------------------------------------------

  subroutine StepHMatKerr(ncyc)

    integer :: ncyc, m, n
    M4_FTYPE :: exh, eyh, ezh, esq
    M4_MODLOOP_DECL({MATKERR},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATKERR},mat,{

    ! this loops over all mat structures, setting mat

      M4_MODOBJ_GETREG(mat,reg)

      n = mod(ncyc-1+2,2) + 1
      m = mod(ncyc+2,2) + 1
    
      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      exh = Ex(i,j,k)
      eyh = Ey(i,j,k)
      ezh = Ez(i,j,k)

      esq = 0

M4_IFELSE_TM({
      esq = esq + exh*exh + eyh*eyh
})
M4_IFELSE_TE({
      esq = esq + ezh*ezh
})
 
M4_IFELSE_TM({
      Ex(i,j,k) = exh + epsinvx(i,j,k) * ( mat%chi3k * esq*exh )!+ mat%S(p,n) * exh )
      Ey(i,j,k) = eyh + epsinvy(i,j,k) * ( mat%chi3k * esq*eyh )!+ mat%S(p,n) * eyh )
})
M4_IFELSE_TE({
      Ez(i,j,k) = ezh + epsinvz(i,j,k) * ( mat%chi3k * esq*ezh )! + mat%S(p,n) * ezh )
})      
     
      })      

    })
  
  end subroutine StepHMatKerr


!----------------------------------------------------------------------

  subroutine StepEMatKerr(ncyc)

    integer :: ncyc, m, n, it
    M4_FTYPE :: exh, eyh, ezh, exo, eyo, ezo, esq
    real(kind=8) :: eix, eiy, eiz
    M4_MODLOOP_DECL({MATKERR},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATKERR},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
! A) solve iteratively:
! 
! E(n) = E'(n) / ( 1 + epsinv * chi3k * E(n)^2 )
       esq = 0
M4_IFELSE_TM({
       exo = Ex(i,j,k)
       eyo = Ey(i,j,k)
       eix = epsinvx(i,j,k)
       eiy = epsinvy(i,j,k)
       exh = exo
       eyh = eyo
       esq = exh*exh + eyh*eyh + esq
})
M4_IFELSE_TE({
       ezo = Ez(i,j,k)
       eiz = epsinvz(i,j,k)
       ezh = ezo
       esq = esq + ezh*ezh
})
       do it = 1, mat%maxit
         
M4_IFELSE_TM({
         exh = exo / ( 1 + eix * ( mat%chi3k * esq ) )
         eyh = eyo / ( 1 + eiy * ( mat%chi3k * esq ) )
})
M4_IFELSE_TE({
         ezh = ezo / ( 1 + eiz * ( mat%chi3k * esq ) )
}) 
         esq = 0  
M4_IFELSE_TM({
         esq = esq + exh*exh + eyh*eyh
})
M4_IFELSE_TE({
         esq = esq + ezh*ezh
})
       end do

M4_IFELSE_TM({
       Ex(i,j,k) = exh
       Ey(i,j,k) = eyh
})
M4_IFELSE_TE({
       Ez(i,j,k) = ezh

       })      
})
       
    })

  end subroutine StepEMatKerr

!----------------------------------------------------------------------

  subroutine SumJEMatKerr(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    logical :: mode
    integer :: ncyc, idx
    real(kind=8) :: sum(MAXEBALCH)
   
  end subroutine SumJEMatKerr

!----------------------------------------------------------------------

  subroutine SumKHMatKerr(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc, idx
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH)

  end subroutine SumKHMatKerr
 
!----------------------------------------------------------------------

  subroutine DisplayMatKerrObj(mat)

    type(T_MATKERR) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" chi3k=",TRIM(f2str(mat%chi3k,5))
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatKerrObj

!----------------------------------------------------------------------

   subroutine EchoMatKerrObj(mat)

    type(T_MATKERR) :: mat

    M4_WRITE_INFO({"--- matkerr # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"chi3k = ",mat%chi3k })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatKerrObj

  
!----------------------------------------------------------------------

end module matkerr

! Authors:  A. Pusch 
! Modified: 06/09/2012
! =====================================================================


