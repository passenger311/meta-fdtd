!-*- F90 -*------------------------------------------------------------
!
!  module: matdrude / meta
!
!  Drude material module.
!
!  subs:
!
!    InitializeMatDrude
!    FinalizeMatDrude
!    ReadMatDrudeObj
!    StepEMatDrude
!    StepHMatDrude
!    SumJEMatDrude
!    SumKHMatDrude
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatDrude calculates the material response of a Drude pole. 
!
! d/dt J + gammapl * J = omegapl**2 * E
! E = E* + J
!
! where E* is the electric field as calculated without the sources.  
!
! A first order (calculating J) or second order approach (calculating
! P) can be chosen.
!
! 1. order method
!
! StepHMatDrude: update eq. J(n+1/2) = c1 * J(n-1/2) + c2 * E(n)
! StepEMatDrude: update eq. E(n+1)* = E(n+1) - epsinv * DT * J(n+1/2)
!
! 2. order method (see matlorentz)
!
! StepHMatPdrude: update eq. P(n+1) = c1 * P(n) + c2 * P(n-1) + c3 * E(n)
! StepEMatPdrude: update eq. E(n+1)* = E(n+1) - epsinv * (P(n+1) - P(n))
!

module matdrude

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

  M4_MATHEAD_DECL({MATDRUDE},MAXMATOBJ,{

  ! input parameters
  real(kind=8) :: lambdaplinv ! inv. vac. plasma wavelength [2 pi c]
  real(kind=8) :: gammapl     ! current damping [1/dt]
  integer :: order            ! use 1. or 2. order solver?

  real(kind=8) :: omegapl

  ! coefficients
  real(kind=8) :: c1, c2, c3
  

  ! current field: J (or Polarisation P) 
  M4_FTYPE, dimension(:,:), pointer :: Jx, Jy, Jz

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatDrudeObj(funit,lcount)

    M4_MODREAD_DECL({MATDRUDE}, funit,lcount,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatDrudeObj")
    
    M4_MODREAD_EXPR({MATDRUDE},funit,lcount,mat,reg,3,out,{ 

    ! read mat parameters here, as defined in mat data structure

    call readfloat(funit,lcount,mat%lambdaplinv)
    call readfloat(funit,lcount,mat%gammapl)
    call readint(funit,lcount,mat%order)

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatDrudeObj")

  end subroutine ReadMatDrudeObj

!----------------------------------------------------------------------

  subroutine InitializeMatDrude

    integer :: err
    M4_MODLOOP_DECL({MATDRUDE},mat) 
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    M4_WRITE_DBG(". enter InitializeMatDrude")
    M4_MODLOOP_EXPR({MATDRUDE},mat,{
    
       ! initialize mat object here

       mat%omegapl = 2. * PI * mat%lambdaplinv
!       mat%gammapl = 2. / ( mat%abslenpl * DT )

       reg = regobj(mat%regidx)

       if ( mat%order .lt. 1 .or. mat%order .gt. 2 ) then
          M4_FATAL_ERROR({"ORDER PARAMETER MUST BE 1 OR 2"})
       endif

       allocate(mat%Jx(reg%numnodes,mat%order),mat%Jy(reg%numnodes,mat%order),mat%Jz(reg%numnodes,mat%order), stat = err)
       M4_ALLOC_ERROR(err,"InitializeMatDrude")

       mat%Jx = 0.
       mat%Jy = 0.
       mat%Jz = 0.


       if ( mat%order .eq. 1 ) then

          mat%c1 = ( 2. - DT * mat%gammapl ) / ( 2. + DT * mat%gammapl )
          mat%c2 = ( 2. * DT ) / ( 2. + DT * mat%gammapl ) * mat%omegapl**2
          mat%c3 = 0.
          
       else

          mat%c1 = 4. / ( 2. + DT * mat%gammapl )
          mat%c2 = ( -2. + DT * mat%gammapl ) / ( 2. + DT * mat%gammapl )
          mat%c3 = ( 2. * DT**2 * mat%omegapl**2) / ( 2. + DT * mat%gammapl )

       endif

! load from checkpoint file

       if ( load_state .and. detail_level .ge. 2 ) then

          read(UNITCHK) mat%Jx, mat%Jy, mat%Jz

       end if

       M4_IFELSE_DBG({call EchoMatDrudeObj(mat)},{call DisplayMatDrudeObj(mat)})


    })
    M4_WRITE_DBG(". exit InitializeMatDrude")

  end subroutine InitializeMatDrude

!----------------------------------------------------------------------

  subroutine FinalizeMatDrude

    M4_MODLOOP_DECL({MATDRUDE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    M4_WRITE_DBG(". enter FinalizeMatDrude")
    M4_MODLOOP_EXPR({MATDRUDE},mat,{

       M4_MODOBJ_GETREG(mat,reg)

! save to checkpoint file

       if ( save_state .and. detail_level .ge. 2 ) then

          write(UNITCHK) mat%Jx, mat%Jy, mat%Jz

       end if

    ! finalize mat object here
    deallocate(mat%Jx,mat%Jy,mat%Jz)

    })
    M4_WRITE_DBG(". exit FinalizeMatDrude")

  end subroutine FinalizeMatDrude

!----------------------------------------------------------------------

  subroutine StepHMatDrude(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATDRUDE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATDRUDE},mat,{

    ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

    if ( mat%order .eq. 1 ) then ! 1. order equation

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

       ! calculate J(n+1/2) from J(n-1/2) and E(n)

M4_IFELSE_TM({
       mat%Jx(p,1) = mat%c1 * mat%Jx(p,1) + mat%c2 * Ex(i,j,k)
       mat%Jy(p,1) = mat%c1 * mat%Jy(p,1) + mat%c2 * Ey(i,j,k)
})
M4_IFELSE_TE({
       mat%Jz(p,1) = mat%c1 * mat%Jz(p,1) + mat%c2 * Ez(i,j,k)
})

       })      

    else ! 2. order equation

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

       ! calculate P(n+1) from P(n),P(n-1) and E(n)

       ! before: J(*,m) is P(n-1), J(*,n) is P(n)

M4_IFELSE_TM({
       mat%Jx(p,m) = mat%c1 * mat%Jx(p,n) + mat%c2 * mat%Jx(p,m) + mat%c3 * Ex(i,j,k)
       mat%Jy(p,m) = mat%c1 * mat%Jy(p,n) + mat%c2 * mat%Jy(p,m) + mat%c3 * Ey(i,j,k)
})
M4_IFELSE_TE({
       mat%Jz(p,m) = mat%c1 * mat%Jz(p,n) + mat%c2 * mat%Jz(p,m) + mat%c3 * Ez(i,j,k)
})

       ! after: J(*,m) is now P(n+1)
       
       ! m and n will be flipped in the next timestep!

       })      

    endif
    })
  
  end subroutine StepHMatDrude


!----------------------------------------------------------------------


  subroutine StepEMatDrude(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATDRUDE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATDRUDE},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

    if ( mat%order .eq. 1 ) then ! 1. order equation
          
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and J(n+1/2)

M4_IFELSE_TM({
       Ex(i,j,k) = Ex(i,j,k) -  w(1) * epsinvx(i,j,k) * DT * mat%Jx(p,1)
       Ey(i,j,k) = Ey(i,j,k) -  w(2) * epsinvy(i,j,k) * DT * mat%Jy(p,1)
})
M4_IFELSE_TE({
       Ez(i,j,k) = Ez(i,j,k) -  w(3) * epsinvz(i,j,k) * DT * mat%Jz(p,1)
})

       })      

    else ! 2. order equation

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

M4_IFELSE_TM({
       Ex(i,j,k) = Ex(i,j,k) - w(1) * epsinvx(i,j,k) * ( mat%Jx(p,m) - mat%Jx(p,n) )
       Ey(i,j,k) = Ey(i,j,k) - w(2) * epsinvy(i,j,k) * ( mat%Jy(p,m) - mat%Jy(p,n) )
})
M4_IFELSE_TE({
       Ez(i,j,k) = Ez(i,j,k) - w(3) * epsinvz(i,j,k) * ( mat%Jz(p,m) - mat%Jz(p,n) )
})

       })      


    end if

    })

  end subroutine StepEMatDrude

!----------------------------------------------------------------------

  subroutine SumJEMatDrude(mask, ncyc, sum, idx, mode)


    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH)
    integer :: ncyc, m, n, idx
   
    M4_MODLOOP_DECL({MATDRUDE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATDRUDE},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

       if ( mask(i,j,k) ) then


          if ( mat%order .eq. 1 ) then ! 1. order equation

             sum(idx+1) = sum(idx+1) + ( &
M4_IFELSE_TM({    M4_VOLEX(i,j,k) * w(1) * dble(Ex(i,j,k)) * dble(mat%Jx(p,1)) + },{0. +}) &
M4_IFELSE_TM({    M4_VOLEY(i,j,k) * w(2) * dble(Ey(i,j,k)) * dble(mat%Jy(p,1)) + },{0. +}) &
M4_IFELSE_TE({    M4_VOLEZ(i,j,k) * w(3) * dble(Ez(i,j,k)) * dble(mat%Jz(p,1))   },{0.  }) &
                  )
             
          else

             sum(idx+1) = sum(idx+1) + M4_VOL(i,j,k) * ( &
M4_IFELSE_TM({    M4_VOLEX(i,j,k) * w(1) * dble(Ex(i,j,k)) * dble( mat%Jx(p,m) - mat%Jx(p,n) ) / DT + },{0. + }) &
M4_IFELSE_TM({    M4_VOLEY(i,j,k) * w(2) * dble(Ey(i,j,k)) * dble( mat%Jy(p,m) - mat%Jy(p,n) ) / DT + },{0. + }) &
M4_IFELSE_TE({    M4_VOLEZ(i,j,k) * w(3) * dble(Ez(i,j,k)) * dble( mat%Jz(p,m) - mat%Jz(p,n) ) / DT   },{0.   }) &
                  )

          end if

       endif

       })      

    idx = idx + NUMEBALCH

    })
    
  end subroutine SumJEMatDrude

!----------------------------------------------------------------------

  subroutine SumKHMatDrude(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc, idx
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH)

    M4_MODLOOP_DECL({MATDRUDE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATDRUDE},mat,{

    idx = idx + NUMEBALCH

    })

  end subroutine SumKHMatDrude
 
!----------------------------------------------------------------------

  subroutine DisplayMatDrudeObj(mat)

    type(T_MATDRUDE) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" lambdapl=",TRIM(f2str(mat%lambdaplinv,5)),&
    	" gammapl=",TRIM(f2str(mat%gammapl,5))
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatDrudeObj
  
!----------------------------------------------------------------------

   subroutine EchoMatDrudeObj(mat)

    type(T_MATDRUDE) :: mat

    M4_WRITE_INFO({"--- matdrude # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"lambdapl = ",mat%lambdaplinv })
    M4_WRITE_INFO({"omegapl = ",mat%omegapl })
    M4_WRITE_INFO({"gammapl = ",mat%gammapl })
    M4_WRITE_INFO({"c1 = ",mat%c1 })
    M4_WRITE_INFO({"c2 = ",mat%c2 })
    M4_WRITE_INFO({"c3 = ",mat%c3 })
    M4_WRITE_INFO({"order = ",mat%order })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatDrudeObj

  
!----------------------------------------------------------------------

end module matdrude

! Authors:  K.Boehringer, J.Hamm 
! Modified: 14/1/2008
! Changed: 7/07/2011 S.Wuestner
!
! =====================================================================


