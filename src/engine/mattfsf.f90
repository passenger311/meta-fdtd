!-*- F90 -*------------------------------------------------------------
!
!  module: mattfsf / meta
!
!  total-field-scattered-field source module. 
!
!----------------------------------------------------------------------


! =====================================================================
!
! The mattfsf module is similar to to the matsource module but feeds
! an incident EH field instead of a driving current. 
!
! Each Yee cell is considered a TFSF source of its own. The incident
! field componenents decide over the angle into which the plane wave
! is directed [see tavlov p 212 for calculation].
!

module mattfsf

  use constant
  use mpiworld
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save


  M4_MATHEAD_DECL({MATTFSF},MAXMATOBJ,{

     real(kind=8) :: lambdainv0    ! inv vacuum wavelength in units of [dx]

     logical      :: tmode         ! time or wavelength mode
! read these is frequency mode     
     real(kind=8) :: dlambdainv    ! spectral width of vac.wave in units of [dx]
     real(kind=8) :: a0            ! gaussian start value as fraction of peak [dt]
! read these is time mode     
     real(kind=8) :: dn            ! time width of gaussian [dt]
     real(kind=8) :: n0            ! time offset of peak [dt]

     real(kind=8) :: noffs,ncw     ! offset in time domain and duration of cw [dt]    

     real(kind=8) :: n1,nend       ! some values used internally ...
     real(kind=8) :: gamma
     real(kind=8) :: omega0, omega1, domega 
     real(kind=8) :: amp, wavefct

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatTfsfObj(funit,lcount)


    M4_MODREAD_DECL({MATTFSF}, funit,lcount,mat,reg,out)
    integer :: v(2)

    M4_WRITE_DBG(". enter ReadMatTfsfObj")

    M4_MODREAD_EXPR({MATTFSF}, funit,lcount,mat,reg, 6, out,{ 

    ! read parameters here, as defined in mat data structure


    call readfloat(funit, lcount, mat%lambdainv0)  ! inv. vacuum wavelength in units of [2 pi c]
    call readints(funit,lcount,v,2)                ! time offset and cw activity [dt] 
    mat%noffs = v(1)
    mat%ncw = v(2)   
    call readlogical(funit,lcount,mat%tmode)       ! time mode? (or wavelength mode)
    if ( mat%tmode ) then
       
       call readfloat(funit,lcount,mat%n0)           ! peak of gaussian in time domain [dt]
       call readfloat(funit,lcount,mat%dn)           ! half width of gaussian in time domain [dt]
    else	
       call readfloat(funit,lcount,mat%a0)         ! gaussian start value as fraction of peak
       call readfloat(funit,lcount,mat%dlambdainv) ! half width of vac.wave in units of [2 pi c]
    end if
       
    ! read regions and output structures

    })

    M4_WRITE_DBG(". exit ReadMatTfsfObj")

  end subroutine ReadMatTfsfObj

!----------------------------------------------------------------------

  subroutine InitializeMatTfsf

    M4_MODLOOP_DECL({MATTFSF},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))

    M4_WRITE_DBG(". enter InitializeMatTfsf")
    M4_MODLOOP_EXPR({MATTFSF},mat,{

    ! check whether region is a face
    reg = regobj(mat%regidx)
    
    ! center frequency
    mat%omega0 = 2. * PI * mat%lambdainv0

    if ( mat%tmode ) then
       ! time mode: n0, dn given
       mat%n1 =  mat%n0 + mat%dn    
       mat%gamma = 1./DT * sqrt( log(2.) ) / mat%dn   
       mat%domega = log(2.) * mat%gamma
       mat%omega1 = mat%omega0 - mat%domega
       mat%a0 = exp(- mat%gamma**2 * ( mat%n0*DT )**2 )
    else
       ! wavelength mode: a0, dlambda given
       mat%omega1 = 2. * PI * (  mat%lambdainv0 + mat%dlambdainv  )
       mat%domega = mat%omega0 - mat%omega1
       mat%gamma =  mat%domega / log(2.)
       mat%dn = sqrt( log(2.) )/ ( mat%gamma * DT )    
       mat%n0 = 1./DT *  sqrt ( - log(mat%a0) / mat%gamma**2 )
       mat%n1 =  mat%n0 + mat%dn
    end if
      
    mat%nend = mat%n0 + mat%ncw + mat%n0 
    
    M4_IFELSE_DBG({call EchoMatTfsfObj(mat)},{call DisplayMatTfsfObj(mat)})
    
    })

    M4_WRITE_DBG(". exit InitializeMatTfsf")

  end subroutine InitializeMatTfsf

!----------------------------------------------------------------------

  subroutine FinalizeMatTfsf

    M4_WRITE_DBG(". enter FinalizeMatTfsf")
    
    M4_WRITE_DBG(". exit FinalizeMatTfsf")

  end subroutine FinalizeMatTfsf

!----------------------------------------------------------------------

  subroutine StepEMatTfsf(ncyc)

    integer :: ncyc
    
    M4_MODLOOP_DECL({MATTFSF},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))
    real(kind=8) :: ncyc0, nend0, es, wavefct

    M4_MODLOOP_EXPR({MATTFSF},mat,{

       M4_MODOBJ_GETREG(mat,reg)

       ncyc0 = 1.0*ncyc - mat%noffs
       
       if ( ncyc0 .ge. 0. .and. ncyc0 .le. mat%nend ) then 

         ! in between attack and decay there is a period of length ncw with cw operation. 
         mat%amp = 1.0
         ! attack phase
         if ( ncyc0 .le. mat%n0 ) then
            mat%amp =  exp ( - mat%gamma**2 * ( ( ncyc0 - mat%n0 ) * DT )**2 )
         end if
         ! decay phase
         if ( ncyc0 .ge. mat%n0 + mat%ncw ) then         	
            mat%amp =  exp ( - mat%gamma**2 * ( ( ncyc0 - mat%n0 - mat%ncw ) * DT )**2 )
         end if
         
         mat%wavefct = mat%amp * sin(mat%omega0*ncyc0*DT) 
            
         M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

M4_IFELSE_TM({
         Ex(i,j,k) = Ex(i,j,k) + DT * ( w(5)/M4_SZ(i,j,k) - w(6)/M4_SY(i,j,k) ) * epsinvx(i,j,k) * mat%wavefct
         Ey(i,j,k) = Ey(i,j,k) + DT * ( w(6)/M4_SX(i,j,k) - w(4)/M4_SZ(i,j,k) ) * epsinvy(i,j,k) * mat%wavefct
})
M4_IFELSE_TE({
         Ez(i,j,k) = Ez(i,j,k) + DT * ( w(4)/M4_SY(i,j,k) - w(5)/M4_SX(i,j,k) ) * epsinvz(i,j,k) * mat%wavefct
})
         })
            
       end if      
    })

  end subroutine StepEMatTfsf

!----------------------------------------------------------------------

  subroutine StepHMatTfsf(ncyc)

    integer :: ncyc

    M4_MODLOOP_DECL({MATTFSF},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))
    real(kind=8) :: ncyc0, nend0

    M4_MODLOOP_EXPR({MATTFSF},mat,{

    M4_MODOBJ_GETREG(mat,reg)

    ncyc0 = 1.0*ncyc - mat%noffs + 0.5 ! step H is half a timestep behind
       
    if ( ncyc0 .ge. 0. .and. ncyc0 .le. mat%nend ) then 

       ! in between attack and decay there is a period of length ncw with cw operation. 
       mat%amp = 1.0
       ! attack phase
       if ( ncyc0 .le. mat%n0 ) then
          mat%amp =  exp ( - mat%gamma**2 * ( ( ncyc0 - mat%n0 ) * DT )**2 )
       end if
       ! decay phase
       if ( ncyc0 .ge. mat%n0 + mat%ncw ) then         	
          mat%amp =  exp ( - mat%gamma**2 * ( ( ncyc0 - mat%n0 - mat%ncw ) * DT )**2 )
       end if
       
       mat%wavefct = mat%amp * sin(mat%omega0*ncyc0*DT)
         
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
          
M4_IFELSE_TE({
       Hx(i,j,k) = Hx(i,j,k) + DT * ( w(3)/M4_SY(i,j,k) - w(2)/M4_SZ(i,j,k) ) * M4_MUINVX(i,j,k) * mat%wavefct
       Hy(i,j,k) = Hy(i,j,k) + DT * ( w(1)/M4_SZ(i,j,k) - w(3)/M4_SX(i,j,k) ) * M4_MUINVY(i,j,k) * mat%wavefct
})
M4_IFELSE_TM({
       Hz(i,j,k) = Hz(i,j,k) + DT * ( w(2)/M4_SX(i,j,k) - w(1)/M4_SY(i,j,k) ) * M4_MUINVZ(i,j,k) * mat%wavefct
})
   
       })
          
    end if
    })
    
    
  end subroutine StepHMatTfsf

!----------------------------------------------------------------------

  real(kind=8) function SumJEMatTfsf(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
    real(kind=8) :: Jx, Jy, Jz
   
    M4_MODLOOP_DECL({MATTFSF},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))

    sum = 0

    M4_MODLOOP_EXPR({MATTFSF},mat,{

    ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

       if ( mask(i,j,k) ) then

M4_IFELSE_TM({
          Jx = - ( w(5)/M4_SZ(i,j,k) - w(6)/M4_SY(i,j,k) ) * mat%wavefct
          Jy = - ( w(6)/M4_SX(i,j,k) - w(4)/M4_SZ(i,j,k) ) * mat%wavefct
})
M4_IFELSE_TE({
          Jz = - ( w(4)/M4_SY(i,j,k) - w(5)/M4_SX(i,j,k) ) * mat%wavefct
})

          sum = sum + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * real(Ex(i,j,k)) * Jx + },{0. +}) &
M4_IFELSE_TM({ M4_VOLEY(i,j,k) * real(Ey(i,j,k)) * Jy + },{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * real(Ez(i,j,k)) * Jz   },{0.  }) &
               )
             
       endif

       })      

    })
    
    SumJEMatTfsf = sum
    
  end function SumJEMatTfsf

!----------------------------------------------------------------------

  real(kind=8) function SumKHMatTfsf(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
    real(kind=8) :: Kx, Ky, Kz
   
    M4_MODLOOP_DECL({MATTFSF},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))

    sum = 0

    M4_MODLOOP_EXPR({MATTFSF},mat,{

    ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

       if ( mask(i,j,k) ) then

M4_IFELSE_TE({
          Kx = - ( w(3)/M4_SY(i,j,k) - w(2)/M4_SZ(i,j,k) ) * mat%wavefct
          Ky = - ( w(1)/M4_SZ(i,j,k) - w(3)/M4_SX(i,j,k) ) * mat%wavefct
})
M4_IFELSE_TM({
          Kz = - ( w(2)/M4_SX(i,j,k) - w(1)/M4_SY(i,j,k) ) * mat%wavefct
})
          sum = sum + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * real(Hx(i,j,k)) * Kx +},{0. +}) &
M4_IFELSE_TM({ M4_VOLEY(i,j,k) * real(Hy(i,j,k)) * Ky +},{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * real(Hz(i,j,k)) * Kz  },{0.  }) &
               )
             
       endif

       })      

    })
    
    SumKHMatTfsf = sum

  end function SumKHMatTfsf


 !----------------------------------------------------------------------

  subroutine DisplayMatTfsfObj(mat)

    type(T_MATTFSF) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" ON=",TRIM(i2str(int(mat%noffs))),&
    	" ATT/DEC=",TRIM(i2str(int(mat%n0))),&
    	" SUS=",TRIM(i2str(int(mat%ncw))),&
    	" OFF=",TRIM(i2str(int(mat%nend))) })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatTfsfObj
  
 !----------------------------------------------------------------------

   subroutine EchoMatTfsfObj(mat)

    type(T_MATTFSF) :: mat

    M4_WRITE_INFO({"--- mat # ",TRIM(i2str(mat%idx))})
    M4_WRITE_INFO({"lambdainv0 = ",mat%lambdainv0   })
    M4_WRITE_INFO({"tmode  = ",mat%tmode   })
    M4_WRITE_INFO({"noffs/ncw  = ",mat%noffs,mat%ncw   })
    M4_WRITE_INFO({"n0 = ",mat%n0 })
    M4_WRITE_INFO({"dn = ",mat%dn })
    M4_WRITE_INFO({"dlambdainv = ",mat%dlambdainv })
    M4_WRITE_INFO({"gamma = ",mat%gamma})
    M4_WRITE_INFO({"a0 = ",mat%a0})
    
    M4_WRITE_INFO({"omega0/omega1 = ", &
         mat%omega0, mat%omega1 })
    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    
  end subroutine EchoMatTfsfObj

!----------------------------------------------------------------------

end module mattfsf

!
! Authors:  J.Hamm
! Modified: 19/02/2007
!
!======================================================================


