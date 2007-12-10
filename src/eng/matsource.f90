!-*- F90 -*------------------------------------------------------------
!
!  module: matsource / meta3
!
!  feed electromagnetic field with source 
!
!  subs:
!
!    InitializeMatSource
!    FinalizeMatSource
!    ReadMatSourceObj
!    StepEMatSource
!    StepHMatSource
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatSource module allows to add sources to the electromagnetic 
! field equations


module matsource

  use constant
  use mpiworld
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save


  M4_MODHEAD_DECL({MATSOURCE},100,{

     real(kind=8) :: lambda0       ! vacuum wavelength in units of [dx]
     real(kind=8) :: dlambda0      ! spectral width of vac.wave in units of [dx]
     real(kind=8) :: a0            ! gaussian start value as fraction of peak
     real(kind=8) :: ampl          ! amplitude = 1.
     logical      :: cw            ! go over to cw after peak?
     logical      :: esource       ! e (or h) field source
     real(kind=8) :: vec(3)        ! vector

     real(kind=8) :: gamma         ! some internal values ...
     real(kind=8) :: omega0, omega1 
     real(kind=8) :: npeak, nmax
     real(kind=8) :: es

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatSourceObj(funit)


    M4_MODREAD_DECL({MATSOURCE}, funit,mat,reg,out)
    
    M4_WRITE_DBG(". enter ReadMatSourceObj/matsource")

    ! read parameters here, as defined in mat data structure

    read(funit,*) mat%lambda0     ! vacuum wavelength in units of [dx]
    M4_WRITE_DBG({"read lambda0: ", mat%lambda0})
    read(funit,*) mat%dlambda0    ! spectral width of vac.wave in units of [dx]
    M4_WRITE_DBG({"read dlambda0: ", mat%dlambda0})
    read(funit,*) mat%a0          ! gaussian start value as fraction of peak
    M4_WRITE_DBG({"read a0: ", mat%a0})
    read(funit,*) mat%ampl        ! amplitude = 1.
    M4_WRITE_DBG({"read ampl: ", mat%ampl})
    read(funit,*) mat%cw          ! go over to cw after peak
    M4_WRITE_DBG({"read cw: ", mat%cw})
    read(funit,*) mat%esource     ! electric/magnetic field source
    M4_WRITE_DBG({"read esource: ", mat%esource})
    read(funit,*) mat%vec(1),mat%vec(2), mat%vec(3) ! vector components
    M4_WRITE_DBG({"read vec(1:3): ", mat%vec(1),mat%vec(2), mat%vec(3) })
   
    ! read regions and output structures

    M4_MODREAD_END({MATSOURCE},funit,mat,reg,out)

  end subroutine ReadMatSourceObj

!----------------------------------------------------------------------

  subroutine InitializeMatSource

    M4_MODLOOP_DECL({MATSOURCE},mat)
    M4_MODLOOP_EXPR({MATSOURCE},mat,{
    
       mat%omega0 = 2. * PI * 1. / ( mat%lambda0 * DT )
       mat%omega1 = 2. * PI * 1. / ( ( mat%lambda0 + mat%dlambda0 ) * DT )
       mat%gamma = (mat%omega0 - mat%omega1 ) / log(2.0)
       mat%npeak =  sqrt ( - log(mat%a0) / mat%gamma**2 )
       mat%nmax = mat%npeak

    })

  end subroutine InitializeMatSource

!----------------------------------------------------------------------

  subroutine FinalizeMatSource

    M4_MODLOOP_DECL({MATSOURCE},mat)
    M4_MODLOOP_EXPR({MATSOURCE},mat,{

       ! finalize mat object here

    })

  end subroutine FinalizeMatSource

!----------------------------------------------------------------------

  subroutine StepEMatSource(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATSOURCE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w)
    real(kind=8) :: es

    M4_MODLOOP_EXPR({MATSOURCE},mat,{

       if ( .not. mat%esource ) return
       es = 1.0
       if ( mat%cw .and. ncyc .lt. mat%nmax ) then
          es =  exp ( - mat%gamma**2 * ( 1.0 * ncyc - mat%npeak )**2 )
       endif

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
          Ex(i,j,k) = Ex(i,j,k) + mat%vec(1) * es  * mat%ampl * cos(mat%omega0*ncyc) * DT
          Ey(i,j,k) = Ey(i,j,k) + mat%vec(2) * es  * mat%ampl * cos(mat%omega0*ncyc) * DT
          Ez(i,j,k) = Ez(i,j,k) + mat%vec(3) * es  * mat%ampl * cos(mat%omega0*ncyc) * DT

       })      
    })

  end subroutine StepEMatSource

!----------------------------------------------------------------------


  subroutine StepHMatSource(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATSOURCE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w)
    real(kind=8) :: es

    M4_MODLOOP_EXPR({MATSOURCE},mat,{

       if ( mat%esource ) return
       es = 1.0
       if ( mat%cw .and. ncyc .lt. mat%nmax ) then
          es =  exp ( - mat%gamma**2 * ( 1.0 * ncyc - mat%npeak )**2 )
       endif

       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
    
          Hx(i,j,k) = Hx(i,j,k) + mat%vec(1) * es  * mat%ampl * cos(mat%omega0*ncyc) * DT
          Hy(i,j,k) = Hy(i,j,k) + mat%vec(2) * es  * mat%ampl * cos(mat%omega0*ncyc) * DT
          Hz(i,j,k) = Hz(i,j,k) + mat%vec(3) * es  * mat%ampl * cos(mat%omega0*ncyc) * DT
       
       })

    })
    
  end subroutine StepHMatSource

 !----------------------------------------------------------------------

end module matsource

! =====================================================================


