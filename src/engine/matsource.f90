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

     logical      :: tmode         ! time or wavelength mode
! read these is frequency mode     
     real(kind=8) :: dlambda       ! spectral width of vac.wave in units of [dx]
     real(kind=8) :: a0            ! gaussian start value as fraction of peak [dt]
! read these is time mode     
     real(kind=8) :: dn            ! time width of gaussian [dt]
     real(kind=8) :: n0            ! time offset of peak [dt]

     real(kind=8) :: noffs,ncw     ! offset in time domain and duration of cw [dt]    
     real(kind=8) :: vec(3)        ! J vector

     real(kind=8) :: n1,nend       ! some values used internally ...
     real(kind=8) :: gamma
     real(kind=8) :: omega0, omega1, domega 
     real(kind=8) :: es

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatSourceObj(funit)


    M4_MODREAD_DECL({MATSOURCE}, funit,mat,reg,out)
   
    M4_WRITE_DBG(". enter ReadMatSourceObj")

    M4_MODREAD_EXPR({MATSOURCE}, funit,mat,reg, 1, out,{ 

    ! read parameters here, as defined in mat data structure

    read(funit,*) mat%lambda0        ! vacuum wavelength in units of [dx]
    read(funit,*) mat%noffs, mat%ncw ! time offset and cw activity [dt]
    read(funit,*) mat%tmode          ! time mode? (or wavelength mode)
    if ( mat%tmode ) then
    	read(funit,*) mat%n0      ! peak of gaussian in time domain [dt]
    	read(funit,*) mat%dn      ! half width of gaussian in time domain [dt]
	else	
    	read(funit,*) mat%a0        ! gaussian start value as fraction of peak
    	read(funit,*) mat%dlambda   ! half width of vac.wave in units of [dx]
	end if
    read(funit,*) mat%vec(1),mat%vec(2), mat%vec(3) ! vector components
   
    ! read regions and output structures

    })

    M4_WRITE_DBG(". exit ReadMatSourceObj")

  end subroutine ReadMatSourceObj

!----------------------------------------------------------------------

  subroutine InitializeMatSource

    M4_MODLOOP_DECL({MATSOURCE},mat)
    M4_WRITE_DBG(". enter InitializeMatSource")
    M4_MODLOOP_EXPR({MATSOURCE},mat,{
    
    ! center frequency
    mat%omega0 = 2. * PI * 1. / ( mat%lambda0 * DT )

    if ( mat%tmode ) then
	    ! time mode: n0, dn given
        mat%n1 =  mat%n0 + mat%dn    
        mat%gamma = sqrt( log(2.) ) / mat%dn    
        mat%domega = log(2.) * mat%gamma
		mat%omega1 = mat%omega0 - mat%domega
		mat%a0 = exp(- mat%gamma**2 * ( mat%n0 )**2 )
    else
        ! wavelength mode: a0, dlambda given
        mat%omega1 = 2. * PI * 1. / ( ( mat%lambda0 + mat%dlambda ) * DT )
		mat%domega = mat%omega0 - mat%omega1
        mat%gamma =  mat%domega / log(2.)
        mat%dn = sqrt( log(2.) )/mat%gamma    
        mat%n0 =  sqrt ( - log(mat%a0) / mat%gamma**2 )
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
    M4_REGLOOP_DECL(reg,p,i,j,k,w(1))
    real(kind=8) :: ncyc0, nend0, es, wavefct

    M4_MODLOOP_EXPR({MATSOURCE},mat,{

       M4_MODOBJ_GETREG(mat,reg)

       ncyc0 = 1.0*ncyc - mat%noffs
       
       if ( ncyc0 .ge. 0. .and. ncyc0 .le. mat%nend ) then 

         ! in between attack and decay there is a period of length ncw with cw operation. 
         es = 1.0
         ! attack phase
         if ( ncyc0 .le. mat%n0 ) then
            es =  exp ( - mat%gamma**2 * ( ncyc0 - mat%n0 )**2 )
		 end if
         ! decay phase
         if ( ncyc0 .ge. mat%n0 + mat%ncw ) then         	
            es =  exp ( - mat%gamma**2 * ( ncyc0 - mat%n0 - mat%ncw )**2 )
		 end if
		 
		 wavefct = es * sin(mat%omega0*ncyc0) * DT
         
         M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

         M4_WRITE_DBG({"source @ ",i,j,k})
           Ex(i,j,k) = Ex(i,j,k) + w(1) * M4_EPSINVX(i,j,k) * mat%vec(1) * wavefct
           Ey(i,j,k) = Ey(i,j,k) + w(1) * M4_EPSINVY(i,j,k) * mat%vec(2) * wavefct
           Ez(i,j,k) = Ez(i,j,k) + w(1) * M4_EPSINVZ(i,j,k) * mat%vec(3) * wavefct
         })

       end if      
    })

  end subroutine StepEMatSource

!----------------------------------------------------------------------

  subroutine StepHMatSource(ncyc)

    integer :: ncyc
    
  end subroutine StepHMatSource

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
    M4_WRITE_INFO({"lambda = ",mat%lambda0   })
    M4_WRITE_INFO({"tmode  = ",mat%tmode   })
    M4_WRITE_INFO({"noffs/ncw  = ",mat%noffs,mat%ncw   })
    M4_WRITE_INFO({"vec(3) = ", mat%vec(1),mat%vec(2), mat%vec(3)})
    M4_WRITE_INFO({"n0 = ",mat%n0 })
    M4_WRITE_INFO({"dn = ",mat%dn })
    M4_WRITE_INFO({"dlambda0 = ",mat%dlambda })
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

