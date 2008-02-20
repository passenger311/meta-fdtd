!-*- F90 -*------------------------------------------------------------
!
!  module: matsource / meta
!
!  feed electromagnetic field with source 
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatSource module allows to add a hard source current to the 
! electromagnetic field equations. The time dependence of the current 
! consists of three phases A,S,D
!
! A = exp(-gamma**2 * ( ncyc - noffs - n0 )**2 * sin(omega*(ncyc-noffs))
! S = 1.0 * sin(omega*(ncyc-noffs))
! D = exp(-gamma**2 * ( ncyc - noffs - n0 - ncw )**2 * sin(omega*(ncyc-noffs))
!
! which are the attack, sustain and decay phases. The attack sets in at 
! ncyc = noffs and lasts for n0 till the peak of the gaussian is reached.
! The sustain phase which lasts ncw follows, till the decay phase sets
! in at ncyc = noffs+n0+ncw. The decay lasts as long as the attack phase.
! After ncyc = noffs+n0+ncw+n0, the source is off.
!

module matsource

  use constant
  use parse
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save


  M4_MATHEAD_DECL({MATSOURCE},MAXMATOBJ,{

     real(kind=8) :: lambdainv0    ! inverse vacuum wavelength in units of [2 pi c]

     logical      :: tmode         ! time or wavelength mode
! read these is frequency mode     
     real(kind=8) :: dlambdainv    ! spectral width of vac.wave in units of [2 pi c]
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

  subroutine ReadMatSourceObj(funit,lcount)


    M4_MODREAD_DECL({MATSOURCE}, funit,lcount,mat,reg,out)
    integer :: v(2)

    M4_WRITE_DBG(". enter ReadMatSourceObj")

    M4_MODREAD_EXPR({MATSOURCE}, funit,lcount,mat,reg, 3, out,{ 

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
    

    })

    M4_WRITE_DBG(". exit ReadMatSourceObj")

  end subroutine ReadMatSourceObj

!----------------------------------------------------------------------

  subroutine InitializeMatSource

    M4_MODLOOP_DECL({MATSOURCE},mat)
    M4_WRITE_DBG(". enter InitializeMatSource")
    M4_MODLOOP_EXPR({MATSOURCE},mat,{
    
    ! center frequency
    mat%omega0 = 2. * PI * mat%lambdainv0

    if ( mat%tmode ) then
	    ! time mode: n0, dn given
        mat%n1 =  mat%n0 + mat%dn    
        mat%gamma = sqrt( log(2.) ) / ( mat%dn * DT )   
        mat%domega = log(2.) * mat%gamma
        mat%omega1 = mat%omega0 - mat%domega
        mat%a0 = exp(- mat%gamma**2 * ( mat%n0 * DT )**2 )
    else
        ! wavelength mode: a0, dlambda given
        mat%omega1 = 2. * PI * ( mat%lambdainv0 + mat%dlambdainv )
        mat%domega = mat%omega0 - mat%omega1
        mat%gamma =  mat%domega / log(2.)
        mat%dn =  1./DT * sqrt( log(2.) )/mat%gamma    
        mat%n0 =  1./DT * sqrt ( - log(mat%a0) / mat%gamma**2 )
        mat%n1 =  mat%n0 + mat%dn
      end if
      
      mat%nend = mat%n0 + mat%ncw + mat%n0 

      M4_IFELSE_DBG({call EchoMatSourceObj(mat)},{call DisplayMatSourceObj(mat)})
 
    })

    M4_WRITE_DBG(". exit InitializeMatSource")

  end subroutine InitializeMatSource

!----------------------------------------------------------------------

  subroutine FinalizeMatSource

    M4_MODLOOP_DECL({MATSOURCE},mat)
    M4_WRITE_DBG(". enter FinalizeMatSource")
    M4_MODLOOP_EXPR({MATSOURCE},mat,{

       ! finalize mat object here

    })

    M4_WRITE_DBG(". exit FinalizeMatSource")

  end subroutine FinalizeMatSource

!----------------------------------------------------------------------

  subroutine StepEMatSource(ncyc)

    integer :: ncyc
    
    M4_MODLOOP_DECL({MATSOURCE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8) :: ncyc0, nend0

    M4_MODLOOP_EXPR({MATSOURCE},mat,{

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

           M4_WRITE_DBG({"source @ ",i,j,k})

M4_IFELSE_TM({
           Ex(i,j,k) = Ex(i,j,k) + DT * w(1) * epsinvx(i,j,k) * mat%wavefct
           Ey(i,j,k) = Ey(i,j,k) + DT * w(2) * epsinvy(i,j,k) * mat%wavefct
})
M4_IFELSE_TE({
           Ez(i,j,k) = Ez(i,j,k) + DT * w(3) * epsinvz(i,j,k) * mat%wavefct
})
         })

       end if      
    })

  end subroutine StepEMatSource

!----------------------------------------------------------------------

  subroutine StepHMatSource(ncyc)

    integer :: ncyc
    
  end subroutine StepHMatSource


!----------------------------------------------------------------------

  real(kind=8) function SumJEMatSource(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
    real(kind=8) :: Jx, Jy, Jz
   
    M4_MODLOOP_DECL({MATSOURCE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))

    sum = 0

    M4_MODLOOP_EXPR({MATSOURCE},mat,{

    ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

       if ( mask(i,j,k) ) then

          Jx = - w(1) * mat%wavefct
          Jy = - w(2) * mat%wavefct
          Jz = - w(3) * mat%wavefct
          
          sum = sum + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * real(Ex(i,j,k)) * Jx + },{0. +}) &
M4_IFELSE_TM({ M4_VOLEY(i,j,k) * real(Ey(i,j,k)) * Jy + },{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * real(Ez(i,j,k)) * Jz   },{0.  }) &
               )
             
       endif

       })      

    })
    
    SumJEMatSource = sum    

  end function SumJEMatSource

!----------------------------------------------------------------------

  real(kind=8) function SumKHMatSource(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
   
    SumKHMatSource = 0.

  end function SumKHMatSource


 !----------------------------------------------------------------------

  subroutine DisplayMatSourceObj(mat)

    type(T_MATSOURCE) :: mat
 
    
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" ON=",TRIM(i2str(int(mat%noffs))),&
    	" ATT/DEC=",TRIM(i2str(int(mat%n0))),&
    	" SUS=",TRIM(i2str(int(mat%ncw))),&
    	" OFF=",TRIM(i2str(int(mat%nend))) })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatSourceObj
  
 !----------------------------------------------------------------------

   subroutine EchoMatSourceObj(mat)

    type(T_MATSOURCE) :: mat

    M4_WRITE_INFO({"--- mat # ",TRIM(i2str(mat%idx))})
    M4_WRITE_INFO({"lambdainv0 = ",mat%lambdainv0   })
    M4_WRITE_INFO({"tmode  = ",mat%tmode   })
    M4_WRITE_INFO({"noffs/ncw  = ",mat%noffs,mat%ncw   })
!    M4_WRITE_INFO({"vec(3) = ", mat%vec(1),mat%vec(2), mat%vec(3)})
    M4_WRITE_INFO({"n0 = ",mat%n0 })
    M4_WRITE_INFO({"dn = ",mat%dn })
    M4_WRITE_INFO({"dlambdainv = ",mat%dlambdainv })
    M4_WRITE_INFO({"gamma = ",mat%gamma})
    M4_WRITE_INFO({"a0 = ",mat%a0})
    
    M4_WRITE_INFO({"omega0/omega1 = ", &
         mat%omega0, mat%omega1 })
    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    
  end subroutine EchoMatSourceObj

!----------------------------------------------------------------------

end module matsource

!
! Authors:  J.Hamm
! Modified: 20/12/2007
!
!======================================================================


