!-*- F90 -*------------------------------------------------------------
!
!  module: mattfsf / meta
!
!  total-field-scattered-field source module. 
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------


! =====================================================================
!
! The mattfsf module is similar to to the matsource module but feeds
! an incident EH field instead of a driving current.
!
! A mattfsf must be defined over a box region with zero volume (the 
! plane of incidence). The incident field is fed in accoding to 
! [Tavlov p 207]. This means 4 field components must be passed in
! (two E and two H). Depending on the face type, the components that 
! must be given are:
!
! box is a k=k_0 face: Ex, Ey, Hx, Hy
! box is a j=j_0 face: Ex, Ez, Hx, Hz
! box is a i=i_0 face: Ey, Ez, Hy, Hz


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


  M4_MATHEAD_DECL({MATTFSF},100,{

     real(kind=8) :: lambda0       ! vacuum wavelength in units of [dx]

     logical      :: tmode         ! time or wavelength mode
! read these is frequency mode     
     real(kind=8) :: dlambda       ! spectral width of vac.wave in units of [dx]
     real(kind=8) :: a0            ! gaussian start value as fraction of peak [dt]
! read these is time mode     
     real(kind=8) :: dn            ! time width of gaussian [dt]
     real(kind=8) :: n0            ! time offset of peak [dt]

     real(kind=8) :: noffs,ncw     ! offset in time domain and duration of cw [dt]    

     real(kind=8) :: n1,nend       ! some values used internally ...
     real(kind=8) :: gamma
     real(kind=8) :: omega0, omega1, domega 
     real(kind=8) :: es
     integer :: face               ! 1 -> i=const, 2 => j=const, 3 => z=const

     real(kind=8), pointer, dimension(:,:,:,:) :: Finc

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatTfsfObj(funit)


    M4_MODREAD_DECL({MATTFSF}, funit,mat,reg,out)
   
    M4_WRITE_DBG(". enter ReadMatTfsfObj")

    M4_MODREAD_EXPR({MATTFSF}, funit,mat,reg, 4, out,{ 

    ! read parameters here, as defined in mat data structure

    read(funit,*) mat%lambda0        ! vacuum wavelength in units of [dx]
    read(funit,*) mat%noffs, mat%ncw ! time offset and cw activity [dt]
    read(funit,*) mat%tmode          ! time mode? (or wavelength mode)
    if ( mat%tmode ) then
       read(funit,*) mat%n0          ! peak of gaussian in time domain [dt]
       read(funit,*) mat%dn          ! half width of gaussian in time domain [dt]
    else	
       read(funit,*) mat%a0          ! gaussian start value as fraction of peak
       read(funit,*) mat%dlambda     ! half width of vac.wave in units of [dx]
    end if
   
    ! read regions and output structures

    })

    M4_WRITE_DBG(". exit ReadMatTfsfObj")

  end subroutine ReadMatTfsfObj

!----------------------------------------------------------------------

  subroutine InitializeMatTfsf

    M4_MODLOOP_DECL({MATTFSF},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(4))
    integer :: err

    M4_WRITE_DBG(". enter InitializeMatTfsf")
    M4_MODLOOP_EXPR({MATTFSF},mat,{

    ! check whether region is a face
    reg = regobj(mat%regidx)

    mat%face = -1

    if ( reg%isbox ) then
       if ( reg%is .eq. reg%ie ) mat%face = 1 
       M4_IFELSE_1D({},{if ( reg%js .eq. reg%je ) mat%face = 2})  
       M4_IFELSE_3D({if ( reg%ks .eq. reg%ke ) mat%face = 3})
       if ( mat%face .le. 0 ) then
          M4_FATAL_ERROR("MATTFSF REGION MUST DEFINE A VALID INTERFACE!")
       endif
    else
       M4_FATAL_ERROR("MATTFSF OBJECT MUST BE DEFINED OVER A BOX!")
    end if

    
    ! allocate and fill incident wave components
    
    allocate(mat%Finc(reg%is-1:reg%ie+1,reg%js-1:reg%je+1,reg%ks-1:reg%ke+1, 1:4), stat = err)
    mat%Finc = 0.0
    M4_ALLOC_ERROR(err,{"InitializeMatTfsf"})
    M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
    mat%Finc(i,j,k,1) = w(1)
    mat%Finc(i,j,k,2) = w(2)
    mat%Finc(i,j,k,3) = w(3)
    mat%Finc(i,j,k,4) = w(4)
    })

    ! center frequency
    mat%omega0 = 2. * PI  / mat%lambda0

    if ( mat%tmode ) then
       ! time mode: n0, dn given
       mat%n1 =  mat%n0 + mat%dn    
       mat%gamma = 1./DT * sqrt( log(2.) ) / mat%dn   
       mat%domega = log(2.) * mat%gamma
       mat%omega1 = mat%omega0 - mat%domega
       mat%a0 = exp(- mat%gamma**2 * ( mat%n0*DT )**2 )
    else
       ! wavelength mode: a0, dlambda given
       mat%omega1 = 2. * PI / ( mat%lambda0 + mat%dlambda  )
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

    M4_MODLOOP_DECL({MATTFSF},mat)
    M4_WRITE_DBG(". enter FinalizeMatTfsf")
    M4_MODLOOP_EXPR({mattfsf},mat,{

    ! finalize mat object here
    deallocate(mat%Finc)

    })

    M4_WRITE_DBG(". exit FinalizeMatTfsf")

  end subroutine FinalizeMatTfsf

!----------------------------------------------------------------------

  subroutine StepEMatTfsf(ncyc)

    integer :: ncyc
    
    M4_MODLOOP_DECL({MATTFSF},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(4))
    real(kind=8) :: ncyc0, nend0, es, wavefct

    M4_MODLOOP_EXPR({MATTFSF},mat,{

       M4_MODOBJ_GETREG(mat,reg)

       ncyc0 = 1.0*ncyc - mat%noffs
       
       if ( ncyc0 .ge. 0. .and. ncyc0 .le. mat%nend ) then 

         ! in between attack and decay there is a period of length ncw with cw operation. 
         es = 1.0
         ! attack phase
         if ( ncyc0 .le. mat%n0 ) then
            es =  exp ( - mat%gamma**2 * ( ( ncyc0 - mat%n0 ) * DT )**2 )
         end if
         ! decay phase
         if ( ncyc0 .ge. mat%n0 + mat%ncw ) then         	
            es =  exp ( - mat%gamma**2 * ( ( ncyc0 - mat%n0 - mat%ncw ) * DT )**2 )
         end if
         
         wavefct = es * sin(mat%omega0*ncyc0*DT) 

         ! i=const face: Finc(3) => Hy, Finc(4) => Hz 
         if ( mat%face .eq. 1 ) then
            
            M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

            Ey(i,j,k) = Ey(i,j,k) + mat%Finc(i,j,k,4) * DT/SX * epsinvy(i,j,k) * wavefct
            Ez(i,j,k) = Ez(i,j,k) + mat%Finc(i,j,k,3) * DT/SX * epsinvz(i,j,k) * wavefct
           
            })
            
         end if
         
         ! j=const face: Finc(3) => Hx, Finc(4) => Hz
         if ( mat%face .eq. 2 ) then
            
            M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
            
            Ex(i,j,k) = Ex(i,j,k) + mat%Finc(i,j,k,4) * DT/SY * epsinvx(i,j,k) * wavefct
            Ez(i,j,k) = Ez(i,j,k) + mat%Finc(i,j,k,3) * DT/SY * epsinvz(i,j,k) * wavefct
            
            })

         end if
         
         ! k=const face: Finc(3) => Hx, Finc(4) => Hy       
         if ( mat%face .eq. 3 ) then
            
            M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
            
            Ex(i,j,k) = Ex(i,j,k) + mat%Finc(i,j,k,4) * DT/SZ * epsinvx(i,j,k) * wavefct
            Ey(i,j,k) = Ey(i,j,k) + mat%Finc(i,j,k,3) * DT/SZ * epsinvy(i,j,k) * wavefct
            
            })

         end if

       end if      
    })

  end subroutine StepEMatTfsf

!----------------------------------------------------------------------

  subroutine StepHMatTfsf(ncyc)

    integer :: ncyc

    M4_MODLOOP_DECL({MATTFSF},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(4))
    real(kind=8) :: ncyc0, nend0, es, wavefct

    M4_MODLOOP_EXPR({MATTFSF},mat,{

    M4_MODOBJ_GETREG(mat,reg)

    ncyc0 = 1.0*ncyc - mat%noffs + 0.5 ! step H is half a timestep behind
       
    if ( ncyc0 .ge. 0. .and. ncyc0 .le. mat%nend ) then 

       ! in between attack and decay there is a period of length ncw with cw operation. 
       es = 1.0
       ! attack phase
       if ( ncyc0 .le. mat%n0 ) then
          es =  exp ( - mat%gamma**2 * ( ( ncyc0 - mat%n0 ) * DT )**2 )
       end if
       ! decay phase
       if ( ncyc0 .ge. mat%n0 + mat%ncw ) then         	
          es =  exp ( - mat%gamma**2 * ( ( ncyc0 - mat%n0 - mat%ncw ) * DT )**2 )
       end if
       
       wavefct = es * sin(mat%omega0*ncyc0*DT)
       
       ! i=const face: Finc(1) => Ey, Finc(2) => Ez 
       if ( mat%face .eq. 1 ) then
          
          M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
          
          Hy(i,j,k) = Hy(i,j,k) + mat%Finc(i,j,k,2) * DT/SX * M4_MUINVY(i,j,k) * wavefct
          Hz(i,j,k) = Hz(i,j,k) + mat%Finc(i,j,k,1) * DT/SX * M4_MUINVZ(i,j,k) * wavefct
          
          })
          
       end if
       
       ! j=const face: Finc(3) => Hx, Finc(4) => Hz
       if ( mat%face .eq. 2 ) then
          
          M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
          
          Hx(i,j,k) = Hx(i,j,k) + mat%Finc(i,j,k,2) * DT/SY * M4_MUINVX(i,j,k) * wavefct
          Hz(i,j,k) = Hz(i,j,k) + mat%Finc(i,j,k,1) * DT/SY * M4_MUINVZ(i,j,k) * wavefct
          
          })
          
       end if
       
       ! k=const face: Finc(3) => Hx, Finc(4) => Hy       
       if ( mat%face .eq. 3 ) then
          
          M4_REGLOOP_EXPR(reg,p,i,j,k,w,{  
          
          Hx(i,j,k) = Hx(i,j,k) + mat%Finc(i,j,k,2) * DT/SZ * M4_MUINVX(i,j,k) * wavefct
          Hy(i,j,k) = Hy(i,j,k) + mat%Finc(i,j,k,1) * DT/SZ * M4_MUINVY(i,j,k) * wavefct
          
          })
          
       end if
       
    end if
    })
    
    
  end subroutine StepHMatTfsf

!----------------------------------------------------------------------

  real(kind=8) function SumJEMatTfsf(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
   
    SumJEMatTfsf = 0.

  end function SumJEMatTfsf

!----------------------------------------------------------------------

  real(kind=8) function SumKHMatTfsf(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
   
    SumKHMatTfsf = 0.

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
    M4_WRITE_INFO({"lambda = ",mat%lambda0   })
    M4_WRITE_INFO({"tmode  = ",mat%tmode   })
    M4_WRITE_INFO({"noffs/ncw  = ",mat%noffs,mat%ncw   })
    M4_WRITE_INFO({"n0 = ",mat%n0 })
    M4_WRITE_INFO({"dn = ",mat%dn })
    M4_WRITE_INFO({"dlambda0 = ",mat%dlambda })
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
! Modified: 25/12/2007
!
!======================================================================


