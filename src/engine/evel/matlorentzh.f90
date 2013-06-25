
!-*- F90 -*------------------------------------------------------------
!
!  module: matlorentzh / meta
!
!  Lorentzh material module.
!
!  subs:
!
!    InitializeMatLorentzh
!    FinalizeMatLorentzh
!    ReadMatLorentzhObj
!    StepEMatLorentzh
!    StepHMatLorentzh
!    SumJEMatLorentzh
!    SumKHMatLorentzh
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatLorentzh module calculates the reponse of a Lorentzh pole
!
! a * d/dt d/dt P + b * d/dt P + c * P = d * E
! d/dt E = d/dt E* - epsinv * d/dt P 
!
! where E* is the electric field as calculated without the sources.  
!
! 2. order method
!
! StepHMatPlorentzh: update eq. P(n+1) = c1 * P(n) + c2 * P(n-1) + c3 * E(n)
! StepEMatPlorentzh: update eq. E(n+1) = E(n+1)* - epsinv * (P(n+1) - P(n))
!

module matlorentzh

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

  M4_MATHEAD_DECL({MATLORENTZH},MAXMATOBJ,{

  ! input parameters
!  real(kind=8) :: lambdalinv    ! inv. vac. plasma wavelength and abs. length
!  real(kind=8) :: gammal        ! resonance width
!  real(kind=8) :: deltaepsl     ! delta epsilon

  real(kind=8) :: a, b, c, d	 ! parameters in extended Lorentzhian


  real(kind=8) :: omegal

  ! coefficients
  real(kind=8) :: c1, c2, c3
  ! coefficients energy balance
  real(kind=8) :: d1, d2, d3, d4
  
  ! polarisation field 
  M4_FTYPE, dimension(:,:), pointer :: Mx, My, Mz

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatLorentzhObj(funit,lcount)

    M4_MODREAD_DECL({MATLORENTZH}, funit,lcount,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatLorentzhObj")
    
    M4_MODREAD_EXPR({MATLORENTZH},funit,lcount,mat,reg,3,out,{ 

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

    M4_WRITE_DBG(". exit ReadMatLorentzhObj")

  end subroutine ReadMatLorentzhObj

!----------------------------------------------------------------------

  subroutine InitializeMatLorentzh

    integer :: err
    M4_MODLOOP_DECL({MATLORENTZH},mat) 
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    M4_WRITE_DBG(". enter InitializeMatLorentzh")
    M4_MODLOOP_EXPR({MATLORENTZH},mat,{
    
       ! initialize mat object here

!       mat%omegal = 2. * PI * mat%lambdalinv
       mat%c = mat%c * ( 2. * PI ) ** 2
       mat%d = mat%d * ( 2. * PI ) ** 2

!!       mat%gammal = 2. / ( mat%abslenl * DT )

       reg = regobj(mat%regidx)

       allocate(mat%Mx(reg%numnodes,2),mat%My(reg%numnodes,2),mat%Mz(reg%numnodes,2), stat = err)
       M4_ALLOC_ERROR(err,"InitializeMatLorentzh")

       mat%Mx = 0.
       mat%My = 0.
       mat%Mz = 0.
       
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

          read(UNITCHK) mat%Mx, mat%My, mat%Mz

       end if

       M4_IFELSE_DBG({call EchoMatLorentzhObj(mat)},{call DisplayMatLorentzhObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatLorentzh")

  end subroutine InitializeMatLorentzh

!----------------------------------------------------------------------

  subroutine FinalizeMatLorentzh

    M4_MODLOOP_DECL({MATLORENTZH},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    M4_WRITE_DBG(". enter FinalizeMatLorentzh")
    M4_MODLOOP_EXPR({MATLORENTZH},mat,{

       M4_MODOBJ_GETREG(mat,reg)

! save to checkpoint file

       if ( save_state .and. detail_level .ge. 2 ) then

          write(UNITCHK) mat%Mx, mat%My, mat%Mz

       end if

    ! finalize mat object here
    deallocate(mat%Mx,mat%My,mat%Mz)

    })
    M4_WRITE_DBG(". exit FinalizeMatLorentzh")

  end subroutine FinalizeMatLorentzh

!----------------------------------------------------------------------

  subroutine StepHMatLorentzh(ncyc)


    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATLORENTZH},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATLORENTZH},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

M4_IFELSE_TM({
       Hx(i,j,k) = Hx(i,j,k) - w(1) * M4_MUINVX(i,j,k) * ( mat%Mx(p,m) - mat%Mx(p,n) )
       Hy(i,j,k) = Hy(i,j,k) - w(2) * M4_MUINVY(i,j,k) * ( mat%My(p,m) - mat%My(p,n) )
})
M4_IFELSE_TE({
       Hz(i,j,k) = Hz(i,j,k) - w(3) * M4_MUINVZ(i,j,k) * ( mat%Mz(p,m) - mat%Mz(p,n) )
})
       })      

    })

  end subroutine StepHMatLorentzh


!----------------------------------------------------------------------

  subroutine StepEMatLorentzh(ncyc)


    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATLORENTZH},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATLORENTZH},mat,{

    ! this loops over all mat structures, setting mat

      M4_MODOBJ_GETREG(mat,reg)

      n = mod(ncyc-1+2,2) + 1
      m = mod(ncyc+2,2) + 1
    
      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      ! calculate P(n+1) from P(n),P(n-1) and E(n)

      ! before: J(*,m) is P(n-1), J(*,n) is P(n)

M4_IFELSE_TM({
      mat%Mx(p,m) = mat%c1 * mat%Mx(p,n) + mat%c2 * mat%Mx(p,m) + mat%c3 * Hx(i,j,k)
      mat%My(p,m) = mat%c1 * mat%My(p,n) + mat%c2 * mat%My(p,m) + mat%c3 * Hy(i,j,k)
})
M4_IFELSE_TE({
      mat%Mz(p,m) = mat%c1 * mat%Mz(p,n) + mat%c2 * mat%Mz(p,m) + mat%c3 * Hz(i,j,k)
})      

      ! after: J(*,m) is now P(n+1)
       
      ! m and n will be flipped in the next timestep!

      })      

    })
  
  end subroutine StepEMatLorentzh

!----------------------------------------------------------------------

  subroutine SumJEMatLorentzh(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc, idx
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH)

    M4_MODLOOP_DECL({MATLorentzh},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATLorentzh},mat,{

    idx = idx + NUMEBALCH

    })
    
  end subroutine SumJEMatLorentzh

!----------------------------------------------------------------------

  subroutine SumKHMatLorentzh(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    logical :: mode
    real(kind=8) :: d34, sum(MAXEBALCH), sum1, sum2, val1, val2
    real(kind=8) :: sum1x, sum1y, sum1z, sum1xx, sum1yy, sum1zz
    integer :: ncyc, m, n, idx
   
    M4_MODLOOP_DECL({MATLORENTZH},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATLORENTZH},mat,{

    sum1 = 0.
    sum2 = 0.
    sum1x = 0.
    sum1y = 0.
    sum1z = 0.
    sum1xx = 0.
    sum1yy = 0.
    sum1zz = 0.

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

          val1 = ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * w(1) * ( mat%d1 * mat%Mx(p,m) + mat%d2 * mat%Mx(p,n) + d34 * Hx(i,j,k)) * &
               ( mat%Mx(p,m) - mat%Mx(p,n) ) / DT + &
               M4_VOLEY(i,j,k) * w(2) * ( mat%d1 * mat%My(p,m) + mat%d2 * mat%My(p,n) + d34 * Hy(i,j,k)) * &
               ( mat%My(p,m) - mat%My(p,n) ) / DT +},{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * w(3) * ( mat%d1 * mat%Mz(p,m) + mat%d2 * mat%Mz(p,n) + d34 * Hz(i,j,k)) * &
               ( mat%Mz(p,m) - mat%Mz(p,n) ) / DT  },{0.  }) &
               )

          sum1 = sum1 + val1
          sum1x = sum1x + i*val1
          sum1y = sum1y + j*val1
          sum1z = sum1z + k*val1
          sum1xx = sum1xx + i*i*val1 ! Corrected from sum1xx = sum1x + i*i*val1
          sum1yy = sum1yy + j*j*val1
          sum1zz = sum1zz + k*k*val1


          val2 = ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * w(1) * Hx(i,j,k) * ( mat%Mx(p,m) - mat%Mx(p,n) ) / DT + &
               M4_VOLEY(i,j,k) * w(2) * Hy(i,j,k) * ( mat%My(p,m) - mat%My(p,n) ) / DT +},{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * w(3) * Hz(i,j,k) * ( mat%Mz(p,m) - mat%Mz(p,n) ) / DT  },{0.  }) &
               )

          sum2 = sum2 + val2

       endif

       })      

    sum(idx) = sum(idx) + sum1
    sum(idx+1) = sum(idx+1) + sum2 - sum1

    sum(idx+3) = sum(idx+3) + sum1x
    sum(idx+4) = sum(idx+4) + sum1y
    sum(idx+5) = sum(idx+5) + sum1z
    sum(idx+6) = sum(idx+6) + sum1xx
    sum(idx+7) = sum(idx+7) + sum1yy
    sum(idx+8) = sum(idx+8) + sum1zz

    idx = idx + NUMEBALCH

    })

  end subroutine SumKHMatLorentzh
 
!----------------------------------------------------------------------

  subroutine DisplayMatLorentzhObj(mat)

    type(T_MATLORENTZH) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" a=",TRIM(f2str(mat%a,5)),&
    	" b=",TRIM(f2str(mat%b,5)),&
    	" c=",TRIM(f2str(mat%c,5)),&
    	" d=",TRIM(f2str(mat%d,5))
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatLorentzhObj

!----------------------------------------------------------------------

   subroutine EchoMatLorentzhObj(mat)

    type(T_MATLORENTZH) :: mat

    M4_WRITE_INFO({"--- matlorentzh # ",&
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
    

  end subroutine EchoMatLorentzhObj

  
!----------------------------------------------------------------------

end module matlorentzh

! Authors:  K.Boehringer, J.Hamm 
! Modified: 14/1/2008
! Changed : 7/07/2011 S.Wuestner
!
! =====================================================================


