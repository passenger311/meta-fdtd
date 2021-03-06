!-*- F90 -*------------------------------------------------------------
!
!  module: matdebye / meta
!
!  Debye material module.
!
!  subs:
!
!    InitializeMatDebye
!    FinalizeMatDebye
!    ReadMatDebyeObj
!    StepEMatDebye
!    StepHMatDebye
!    SumJEMatDebye
!    SumKHMatDebye
!
!----------------------------------------------------------------------
! =====================================================================
!
! The MatDebye module calculates the reponse of a Debye pole
!
! taud d/dt P + P = deltaepsd * E
!
! where E* is the electric field as calculated without the sources.  
!
! 1. order method
!
! StepHMatPdebye: update eq. P(n+1/2) = c1 * P(n-1/2) + c2 * E(n)
! StepEMatPdebye: update eq. E(n+1) = a1 * E(n+1)* - a2 * P(n+1/2) - a3 * P(n-1/2)
!

module matdebye

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

  M4_MATHEAD_DECL({MATDEBYE},MAXMATOBJ,{

  ! input parameters
  real(kind=8) :: taud       ! relaxation time
  real(kind=8) :: deltaepsd  ! delta epsilon

  real(kind=8) :: omegad

  ! coefficients
  real(kind=8) :: c1, c2
  
  ! polarisation field 
  M4_FTYPE, dimension(:,:), pointer :: Px, Py, Pz

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatDebyeObj(funit,lcount)

    M4_MODREAD_DECL({MATDEBYE}, funit,lcount,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatDebyeObj")
    
    M4_MODREAD_EXPR({MATDEBYE},funit,lcount,mat,reg,3,out,{ 

    ! read mat parameters here, as defined in mat data structure
    call readfloat(funit,lcount, mat%taud)
    call readfloat(funit,lcount, mat%deltaepsd)

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatDebyeObj")

  end subroutine ReadMatDebyeObj

!----------------------------------------------------------------------

  subroutine InitializeMatDebye

    type (T_REG) :: reg
    integer :: err
    M4_MODLOOP_DECL({MATDEBYE},mat) 
    M4_WRITE_DBG(". enter InitializeMatDebye")
    M4_MODLOOP_EXPR({MATDEBYE},mat,{
    
       ! initialize mat object here

       reg = regobj(mat%regidx)

       allocate(mat%Px(reg%numnodes,2),mat%Py(reg%numnodes,2),mat%Pz(reg%numnodes,2), stat = err)
       M4_ALLOC_ERROR(err,"InitializeMatDebye")

       mat%Px = 0.
       mat%Py = 0.
       mat%Pz = 0.

       mat%c1 = ( 2. * mat%taud - DT ) / ( 2. * mat%taud + DT )
       mat%c2 = 2. * DT * mat%deltaepsd / ( 2. * mat%taud + DT )

! load from checkpoint file

       if ( load_state .and. detail_level .ge. 2 ) then

          write(UNITCHK) mat%Px, mat%Py, mat%Pz

       end if

       M4_IFELSE_DBG({call EchoMatDebyeObj(mat)},{call DisplayMatDebyeObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatDebye")

  end subroutine InitializeMatDebye

!----------------------------------------------------------------------

  subroutine FinalizeMatDebye

    M4_MODLOOP_DECL({MATDEBYE},mat)
    M4_WRITE_DBG(". enter FinalizeMatDebye")
    M4_MODLOOP_EXPR({MATDEBYE},mat,{

! save to checkpoint file

       if ( save_state .and. detail_level .ge. 2 ) then

          write(UNITCHK) mat%Px, mat%Py, mat%Pz

       end if

    ! finalize mat object here
    deallocate(mat%Px,mat%Py,mat%Pz)

    })
    M4_WRITE_DBG(". exit FinalizeMatDebye")

  end subroutine FinalizeMatDebye

!----------------------------------------------------------------------

  subroutine StepHMatDebye(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATDEBYE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATDEBYE},mat,{

    ! this loops over all mat structures, setting mat

      M4_MODOBJ_GETREG(mat,reg)

      m = mod(ncyc+2,2) + 1
      n = mod(ncyc-1+2,2) + 1

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

      ! calculate P(n+1/2) from P(n-1/2) and E(n)

M4_IFELSE_TM({
      mat%Px(p,m) = mat%c1 * mat%Px(p,n) + mat%c2 * Ex(i,j,k)
      mat%Py(p,m) = mat%c1 * mat%Py(p,n) + mat%c2 * Ey(i,j,k)
})
M4_IFELSE_TE({
      mat%Pz(p,m) = mat%c1 * mat%Pz(p,n) + mat%c2 * Ez(i,j,k)
})      

      })      

    })
  
  end subroutine StepHMatDebye


!----------------------------------------------------------------------


  subroutine StepEMatDebye(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATDEBYE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8) :: dx,dy,dz

    M4_MODLOOP_EXPR({MATDEBYE},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

    m = mod(ncyc+2,2) + 1
    n = mod(ncyc-1+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1/2),P(n-1/2)

       dx = 2. + w(1) * epsinvx(i,j,k) * mat%c2
       dy = 2. + w(2) * epsinvy(i,j,k) * mat%c2
       dz = 2. + w(3) * epsinvz(i,j,k) * mat%c2

M4_IFELSE_TM({
       Ex(i,j,k) = 1./dx * ( 2. * Ex(i,j,k) - &
            w(1) * epsinvx(i,j,k) * ( mat%c1 * mat%Px(p,m) - mat%Px(p,n) ) )
       Ey(i,j,k) = 1./dy * ( 2. * Ey(i,j,k) - &
            w(2) * epsinvy(i,j,k) * ( mat%c1 * mat%Py(p,m) - mat%Py(p,n) ) )
})
M4_IFELSE_TE({
       Ez(i,j,k) = 1./dz * ( 2. * Ez(i,j,k) - &
            w(3) * epsinvz(i,j,k) * ( mat%c1 * mat%Pz(p,m) + mat%Pz(p,n) ) )
})

       })      

    })

  end subroutine StepEMatDebye

!----------------------------------------------------------------------

  subroutine SumJEMatDebye(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH)
    integer :: ncyc, m, n, idx
   
    M4_MODLOOP_DECL({MATDEBYE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATDEBYE},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

       if ( mask(i,j,k) ) then

          sum(idx+1) = sum(idx+1) + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * w(1) * dble(Ex(i,j,k)) * dble( mat%Px(p,m) - mat%Px(p,n) ) / DT +},{0. +}) &
M4_IFELSE_TM({ M4_VOLEY(i,j,k) * w(2) * dble(Ey(i,j,k)) * dble( mat%Py(p,m) - mat%Py(p,n) ) / DT +},{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * w(3) * dble(Ez(i,j,k)) * dble( mat%Pz(p,m) - mat%Pz(p,n) ) / DT  },{0.  }) &
               )

       endif

       })      

    idx = idx + NUMEBALCH

    })
    
  end subroutine SumJEMatDebye

!----------------------------------------------------------------------

  subroutine SumKHMatDebye(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc, idx
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH)

    M4_MODLOOP_DECL({MATDEBYE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATDEBYE},mat,{

    idx = idx + NUMEBALCH

    })

  end subroutine SumKHMatDebye
 
!----------------------------------------------------------------------

  subroutine DisplayMatDebyeObj(mat)

    type(T_MATDEBYE) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" taud=",TRIM(f2str(mat%taud,5)),&
    	" deltaepsd=",TRIM(f2str(mat%deltaepsd,5))
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatDebyeObj
  
!----------------------------------------------------------------------

   subroutine EchoMatDebyeObj(mat)

    type(T_MATDEBYE) :: mat

    M4_WRITE_INFO({"--- matdebye # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"deltaepsd = ",mat%deltaepsd })
    M4_WRITE_INFO({"taud = ",mat%taud })
    M4_WRITE_INFO({"c1 = ",mat%c1 })
    M4_WRITE_INFO({"c2 = ",mat%c2 })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatDebyeObj

  
!----------------------------------------------------------------------

end module matdebye

! Authors:  J.Hamm 
! Modified: 14/1/2008
! Changed: 7/07/2011 S.Wuestner
!
! =====================================================================


