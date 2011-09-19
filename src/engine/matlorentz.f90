
!-*- F90 -*------------------------------------------------------------
!
!  module: matlorentz / meta
!
!  Lorentz material module.
!
!  subs:
!
!    InitializeMatLorentz
!    FinalizeMatLorentz
!    ReadMatLorentzObj
!    StepEMatLorentz
!    StepHMatLorentz
!    SumJEMatLorentz
!    SumKHMatLorentz
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatLorentz module calculates the reponse of a Lorentz pole
!
! a * d/dt d/dt P + b * d/dt P + c * P = d * E
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
  use checkpoint
  use parse
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MATHEAD_DECL({MATLORENTZ},MAXMATOBJ,{

  ! input parameters
!  real(kind=8) :: lambdalinv    ! inv. vac. plasma wavelength and abs. length
!  real(kind=8) :: gammal        ! resonance width
!  real(kind=8) :: deltaepsl     ! delta epsilon

  real(kind=8) :: a, b, c, d	 ! parameters in extended Lorentzian


  real(kind=8) :: omegal

  ! coefficients
  real(kind=8) :: c1, c2, c3
  ! coefficients energy balance
  real(kind=8) :: d1, d2, d3, d4
  
  ! polarisation field 
  M4_FTYPE, dimension(:,:), pointer :: Px, Py, Pz

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatLorentzObj(funit,lcount)

    M4_MODREAD_DECL({MATLORENTZ}, funit,lcount,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatLorentzObj")
    
    M4_MODREAD_EXPR({MATLORENTZ},funit,lcount,mat,reg,3,out,{ 

    ! read mat parameters here, as defined in mat data structure
!    call readfloat(funit,lcount, mat%lambdalinv)
!    call readfloat(funit,lcount, mat%gammal)
!    call readfloat(funit,lcount, mat%deltaepsl)
    call readfloat(funit,lcount, mat%a)
    call readfloat(funit,lcount, mat%b)
    call readfloat(funit,lcount, mat%c)
    call readfloat(funit,lcount, mat%d)

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatLorentzObj")

  end subroutine ReadMatLorentzObj

!----------------------------------------------------------------------

  subroutine InitializeMatLorentz

    integer :: err
    M4_MODLOOP_DECL({MATLORENTZ},mat) 
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    M4_WRITE_DBG(". enter InitializeMatLorentz")
    M4_MODLOOP_EXPR({MATLORENTZ},mat,{
    
       ! initialize mat object here

!       mat%omegal = 2. * PI * mat%lambdalinv
       mat%c = mat%c * ( 2. * PI ) ** 2
       mat%d = mat%d * ( 2. * PI ) ** 2

!!       mat%gammal = 2. / ( mat%abslenl * DT )

       reg = regobj(mat%regidx)

       allocate(mat%Px(reg%numnodes,2),mat%Py(reg%numnodes,2),mat%Pz(reg%numnodes,2), stat = err)
       M4_ALLOC_ERROR(err,"InitializeMatLorentz")

       mat%Px = 0.
       mat%Py = 0.
       mat%Pz = 0.
       
       mat%c1 = ( 4. * mat%a - 2. * mat%c * DT**2 ) / ( 2. * mat%a + DT * mat%b )
       mat%c2 = ( -2. * mat%a + DT * mat%b ) / ( 2. * mat%a + DT * mat%b )
       mat%c3 = 2. * mat%d * DT**2 / ( 2. * mat%a + DT * mat%b )

       mat%d1 = 0.5 * mat%b / ( 2. * mat%d * DT ) * (  - 16. * mat%a**2  + &
                2. * mat%c * DT**2 * ( 2. * mat%a - mat%b * DT ) ) / &
                ( 4. * mat%a**2 - mat%b**2 * DT**2 )
       mat%d2 = - 0.5 * mat%b / ( 2. * mat%d * DT ) * (  - 16. * mat%a**2  + &
                2. * mat%c * DT**2 * ( 2. * mat%a + mat%b * DT ) ) / &
                ( 4. * mat%a**2 - mat%b**2 * DT**2 )
       mat%d3 = 2. * mat%a / ( 2. * mat%a + mat%b * DT )
       mat%d4 = 2. * mat%a / ( 2. * mat%a - mat%b * DT )

! load from checkpoint file

       if ( load_state .and. detail_level .ge. 2 ) then

          read(UNITCHK) mat%Px, mat%Py, mat%Pz

       end if

       M4_IFELSE_DBG({call EchoMatLorentzObj(mat)},{call DisplayMatLorentzObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatLorentz")

  end subroutine InitializeMatLorentz

!----------------------------------------------------------------------

  subroutine FinalizeMatLorentz

    M4_MODLOOP_DECL({MATLORENTZ},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    M4_WRITE_DBG(". enter FinalizeMatLorentz")
    M4_MODLOOP_EXPR({MATLORENTZ},mat,{

       M4_MODOBJ_GETREG(mat,reg)

! save to checkpoint file

       if ( save_state .and. detail_level .ge. 2 ) then

          write(UNITCHK) mat%Px, mat%Py, mat%Pz

       end if

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

      n = mod(ncyc-1+2,2) + 1
      m = mod(ncyc+2,2) + 1
    
      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      ! calculate P(n+1) from P(n),P(n-1) and E(n)

      ! before: J(*,m) is P(n-1), J(*,n) is P(n)

M4_IFELSE_TM({
      mat%Px(p,m) = mat%c1 * mat%Px(p,n) + mat%c2 * mat%Px(p,m) + mat%c3 * Ex(i,j,k)
      mat%Py(p,m) = mat%c1 * mat%Py(p,n) + mat%c2 * mat%Py(p,m) + mat%c3 * Ey(i,j,k)
})
M4_IFELSE_TE({
      mat%Pz(p,m) = mat%c1 * mat%Pz(p,n) + mat%c2 * mat%Pz(p,m) + mat%c3 * Ez(i,j,k)
})      

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

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

M4_IFELSE_TM({
       Ex(i,j,k) = Ex(i,j,k) - w(1) * epsinvx(i,j,k) * ( mat%Px(p,m) - mat%Px(p,n) )
       Ey(i,j,k) = Ey(i,j,k) - w(2) * epsinvy(i,j,k) * ( mat%Py(p,m) - mat%Py(p,n) )
})
M4_IFELSE_TE({
       Ez(i,j,k) = Ez(i,j,k) - w(3) * epsinvz(i,j,k) * ( mat%Pz(p,m) - mat%Pz(p,n) )
})
       })      

    })

  end subroutine StepEMatLorentz

!----------------------------------------------------------------------

  subroutine SumJEMatLorentz(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    logical :: mode
    real(kind=8) :: d34, sum(MAXEBALCH), sum1, sum2
    integer :: ncyc, m, n, idx
   
    M4_MODLOOP_DECL({MATLORENTZ},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATLORENTZ},mat,{

    sum1 = 0
    sum2 = 0

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       if ( mode ) then
          d34 = mat%d4
       else
          d34 = mat%d3
       endif

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

       if ( mask(i,j,k) ) then

          sum1 = sum1 + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * w(1) * ( mat%d1 * mat%Px(p,m) + mat%d2 * mat%Px(p,n) + d34 * Ex(i,j,k)) * &
               ( mat%Px(p,m) - mat%Px(p,n) ) / DT + &
               M4_VOLEY(i,j,k) * w(2) * ( mat%d1 * mat%Py(p,m) + mat%d2 * mat%Py(p,n) + d34 * Ey(i,j,k)) * &
               ( mat%Py(p,m) - mat%Py(p,n) ) / DT +},{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * w(3) * ( mat%d1 * mat%Pz(p,m) + mat%d2 * mat%Pz(p,n) + d34 * Ez(i,j,k)) * &
               ( mat%Pz(p,m) - mat%Pz(p,n) ) / DT  },{0.  }) &
               )

          sum2 = sum2 + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * w(1) * Ex(i,j,k) * ( mat%Px(p,m) - mat%Px(p,n) ) / DT + &
               M4_VOLEY(i,j,k) * w(2) * Ey(i,j,k) * ( mat%Py(p,m) - mat%Py(p,n) ) / DT +},{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * w(3) * Ez(i,j,k) * ( mat%Pz(p,m) - mat%Pz(p,n) ) / DT  },{0.  }) &
               )

       endif

       })      

    sum(idx) = sum(idx) + sum1
    sum(idx+1) = sum(idx+1) + sum2 - sum1
    idx = idx + NUMEBALCH

    })

    
  end subroutine SumJEMatLorentz

!----------------------------------------------------------------------

  subroutine SumKHMatLorentz(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc, idx
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH)

    M4_MODLOOP_DECL({MATLorentz},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATLorentz},mat,{

    idx = idx + NUMEBALCH

    })

  end subroutine SumKHMatLorentz
 
!----------------------------------------------------------------------

  subroutine DisplayMatLorentzObj(mat)

    type(T_MATLORENTZ) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" a=",TRIM(f2str(mat%a,5)),&
    	" b=",TRIM(f2str(mat%b,5)),&
    	" c=",TRIM(f2str(mat%c,5)),&
    	" d=",TRIM(f2str(mat%d,5))
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatLorentzObj

!----------------------------------------------------------------------

   subroutine EchoMatLorentzObj(mat)

    type(T_MATLORENTZ) :: mat

    M4_WRITE_INFO({"--- matlorentz # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"a = ",mat%a })
    M4_WRITE_INFO({"b = ",mat%b })
    M4_WRITE_INFO({"c = ",mat%c })
    M4_WRITE_INFO({"d = ",mat%d })
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
! Changed : 7/07/2011 S.Wuestner
!
! =====================================================================


