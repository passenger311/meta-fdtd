!-*- F90 -*------------------------------------------------------------
!
!  module: matjsource / meta
!
!  feed electromagnetic field with hard dipole current source 
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatJSource module allows to add a hard source current to the 
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

module matjsource

  use constant
  use parse
  use reglist
  use outlist
  use grid
  use signal
  use fdtd
  use fdtd_calc

  implicit none
  private
  save


  M4_MATHEAD_DECL({MATJSOURCE},MAXMATOBJ,{

     real(kind=8) :: lambdainv0                ! inverse vacuum wavelength in units of [2 pi c]
     real(kind=8) :: dn                        ! time width of gaussian [dt]

     real(kind=8) :: gamma                      
     real(kind=8) :: omega0, domega 

     real(kind=8) :: noffs, natt, nsus, ndcy   ! generic signal parameters [dt]
     real(kind=8) :: nend                      ! end of signal [dt]

     real(kind=8) :: theta, phi, psi           ! angles of incident wavefront

     real(kind=8) :: wavefct                   ! wave functions
     real(kind=8) :: finc(3)                   ! field components of incident field 

     logical :: planewave                      ! plane wave mode?

     real(kind=8) :: kinc(3)                   ! normed k vector of plane wave
     real(kind=8) :: orig(3)                   ! orgin of the coordinate system
     real(kind=8) :: pvel                      ! phase velocity

     real(kind=8) :: nrefr

     ! time delay field from point of origin
     real(kind=8), pointer, dimension(:) :: delay  

     ! component delay correction
     real(kind=8) :: cdelay(3)

     ! time signal lookup table
     integer :: tres
     real(kind=8), pointer, dimension(:) :: signal 
     integer :: signalp

     ! position of h field 
     integer, pointer, dimension(:,:) :: orient

     ! maximum delay: origin to point furthest away
     integer :: maxdelay


  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatJSourceObj(funit,lcount)


    M4_MODREAD_DECL({MATJSOURCE}, funit,lcount,mat,reg,out)
    real(kind=8) :: v(4)
    logical :: eof, err
    real(kind=8) :: angles(3)
    character(len=LINELNG) :: line

    M4_WRITE_DBG(". enter ReadMatJSourceObj")

    M4_MODREAD_EXPR({MATJSOURCE}, funit,lcount,mat,reg, 3, out,{ 

    ! read parameters here, as defined in mat data structure


    ! read parameters here, as defined in mat data structure

    call readfloat(funit, lcount, mat%lambdainv0)   ! inv. vacuum wavelength in units of [2 pi c]
    call readfloat(funit,lcount,mat%dn)             ! half width of gaussian in time domain [dt]

    call readfloats(funit,lcount,v,4)               ! generic signal parameters [dt] 
    mat%noffs = v(1)
    mat%natt =  v(2)
    mat%nsus =  v(3)
    mat%ndcy =  v(4)
    
    ! optional: configure plane wave source? 

    err = .false.
    angles = 0.
    mat%phi = 0.
    mat%theta = 0.
    mat%psi = 0.

    call readline(funit,lcount,eof,line)
    call getlogical(line,mat%planewave,err)
    call getfloats(line,angles,3,err)
    M4_PARSE_ERROR({err},lcount,{PLANE WAVE ANGLES})
    if ( mat%planewave ) then
       mat%phi = angles(1)
       mat%theta = angles(2)
       mat%psi = angles(3)
    end if

    })

    M4_WRITE_DBG(". exit ReadMatJSourceObj")

  end subroutine ReadMatJSourceObj

!----------------------------------------------------------------------

  subroutine InitializeMatJSource

    M4_MODLOOP_DECL({MATJSOURCE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8) :: r, rk(3), mdelay
    integer :: err

    M4_WRITE_DBG(". enter InitializeMatJSource")
    M4_MODLOOP_EXPR({MATJSOURCE},mat,{
    
    M4_MODOBJ_GETREG(mat,reg)
 
    ! center frequency
    mat%omega0 = 2. * PI * mat%lambdainv0

    ! derived gaussian parameters
    mat%gamma = sqrt( log(2.) ) / ( mat%dn * DT )   
    mat%domega = log(2.) * mat%gamma

    mat%maxdelay = 0
    
    if ( mat%planewave ) then

       ! normalised k vector
       mat%kinc(1) = sin(DEG*mat%theta)*cos(DEG*mat%phi)
       mat%kinc(2) = sin(DEG*mat%theta)*sin(DEG*mat%phi)
       mat%kinc(3) = cos(DEG*mat%theta)

       ! incident field/current amplitudes
       mat%finc(1) = cos(DEG*mat%psi)*sin(DEG*mat%phi) - sin(DEG*mat%psi)*cos(DEG*mat%theta)*cos(DEG*mat%phi)
       mat%finc(2) = -cos(DEG*mat%psi)*cos(DEG*mat%phi) - sin(DEG*mat%psi)*cos(DEG*mat%theta)*sin(DEG*mat%phi)
       mat%finc(3) = sin(DEG*mat%psi)*sin(DEG*mat%theta)

       ! decide on origin according to quadrant into which we emmit
       mat%theta = mod(mat%theta, 180.)
       mat%phi = mod(mat%phi, 360.)
       mat%psi = mod(mat%psi, 360.)

       if ( mat%theta .ge. 0. .and. mat%theta .lt. 90. ) then
          mat%orig(3) = reg%ks
       else
          mat%orig(3) = reg%ke
       endif

       if ( mat%phi .ge. 0. .and. mat%phi .lt. 90. ) then
          mat%orig(1) = reg%is
          mat%orig(2) = reg%js
       end if

       if ( mat%phi .ge. 90. .and. mat%phi .lt. 180 ) then
          mat%orig(1) = reg%ie
          mat%orig(2) = reg%js
       end if

       if ( mat%phi .ge. 180. .and. mat%phi .lt. 270. ) then
          mat%orig(1) = reg%ie
          mat%orig(2) = reg%je
       end if

       if ( mat%phi .ge. 270. .and. mat%phi .lt. 360. ) then
          mat%orig(1) = reg%is
          mat%orig(2) = reg%je
       end if

       ! determine numerical phase velocity
       
       mat%pvel = 1.

!       call NumericalPhaseVelocity(mat%kinc, mat%omega0, SX, SY, SZ, mat%pvel)

       ! optical delay of point-to-origin field

       i = mat%orig(1)
       j = mat%orig(2)
       k = mat%orig(3)

       ! calculate phase velocity

       mat%pvel = 1.0 / mat%nrefr
       call NumericalPhaseVelocity(mat%kinc, mat%omega0, mat%nrefr, SX, SY, SZ, mat%pvel)
       M4_WRITE_INFO({"phase velocity = ",mat%pvel})


       allocate(mat%delay(reg%numnodes), stat = err )
       M4_ALLOC_ERROR(err,{"InitializeMatJSource"})

       mdelay = 0
       
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
             ! project (i,j,k) location vector on kinc to create a distance field

             rk(1) = M4_DISTX(mat%orig(1),i)
             rk(2) = M4_DISTY(mat%orig(2),j)
             rk(3) = M4_DISTZ(mat%orig(3),k)
       
             ! distance projected to kinc

             r = rk(1)*mat%kinc(1) + rk(2)*mat%kinc(2) +  rk(3)*mat%kinc(3)

             ! time for signal at origin to arrive at (i,j,k)

             mat%delay(p) = r / ( mat%pvel * DT )

             ! maximum delay

             if ( mat%delay(p) .gt. mdelay ) mdelay = mat%delay(p)

       })


       mat%cdelay(1) = .5*SX*mat%kinc(1) /  ( mat%pvel * DT )
       mat%cdelay(2) = .5*SY*mat%kinc(2) /  ( mat%pvel * DT )
       mat%cdelay(3) = .5*SZ*mat%kinc(3) /  ( mat%pvel * DT )

       mat%maxdelay = int(mdelay + 1) + sqrt( SX*SX + SY*SY + SZ*SZ )/( mat%pvel * DT )

       mat%tres = 100. ! increased time resolution of delay buffer by this factor 

       allocate(mat%signal(0:mat%maxdelay * mat%tres -1) , stat = err )
       M4_ALLOC_ERROR(err,{"InitializeMatTfSource"})

       mat%signal = 0.

       mat%signalp = 0 ! initialise signal index pointers

    end if

    mat%nend = mat%natt + mat%nsus + mat%ndcy + mat%maxdelay

    M4_IFELSE_DBG({call EchoMatJSourceObj(mat)},{call DisplayMatJSourceObj(mat)})
      
    })

    M4_WRITE_DBG(". exit InitializeMatJSource")

  end subroutine InitializeMatJSource

!----------------------------------------------------------------------

  subroutine FinalizeMatJSource

    M4_MODLOOP_DECL({MATJSOURCE},mat)
    M4_WRITE_DBG(". enter FinalizeMatJSource")
    M4_MODLOOP_EXPR({MATJSOURCE},mat,{

      ! finalize mat object here
  
      if ( mat%planewave ) then
       
         deallocate(mat%signal, mat%delay)

      end if

    })

    M4_WRITE_DBG(". exit FinalizeMatJSource")

  end subroutine FinalizeMatJSource

!----------------------------------------------------------------------

  subroutine StepEMatJSource(ncyc)

    integer :: ncyc
    
    M4_MODLOOP_DECL({MATJSOURCE},mat)

    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8) :: ncyc1
    real(kind=8) :: wavefct(3), d
    integer :: di, l, n

    M4_MODLOOP_EXPR({MATJSOURCE},mat,{

      if ( ncyc .gt. mat%nend ) cycle 

      M4_MODOBJ_GETREG(mat,reg)

      ! enter e-field modulation loop

      if ( .not. mat%planewave ) then

        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

M4_IFELSE_TM({
           Ex(i,j,k) = Ex(i,j,k) + DT * w(1) * epsinvx(i,j,k) * mat%wavefct
           Ey(i,j,k) = Ey(i,j,k) + DT * w(2) * epsinvy(i,j,k) * mat%wavefct
})
M4_IFELSE_TE({
           Ez(i,j,k) = Ez(i,j,k) + DT * w(3) * epsinvz(i,j,k) * mat%wavefct
})

        })

      else ! plane wave!

         n = mat%signalp+mat%maxdelay

         M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

         do l = 1,3 
            
            d = mat%delay(p) + mat%cdelay(l)
            di = int(d * real(mat%tres) + 0.5)
            wavefct(l) = mat%signal(mod(n-di, mat%maxdelay*mat%tres ))

            w(l) = w(l) * mat%finc(l)

         end do

M4_IFELSE_TM({
           Ex(i,j,k) = Ex(i,j,k) + DT * w(1) * epsinvx(i,j,k) * wavefct(1) 
           Ey(i,j,k) = Ey(i,j,k) + DT * w(2) * epsinvy(i,j,k) * wavefct(2)
})
M4_IFELSE_TE({
           Ez(i,j,k) = Ez(i,j,k) + DT * w(3) * epsinvz(i,j,k) * wavefct(3)
})

         })

      endif
   })

  end subroutine StepEMatJSource

!----------------------------------------------------------------------

  subroutine StepHMatJSource(ncyc)

    integer :: ncyc

    M4_MODLOOP_DECL({MATJSOURCE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8) :: ncyc1

    M4_MODLOOP_EXPR({MATJSOURCE},mat,{

       if ( ncyc - mat%noffs .gt. mat%nend + mat%maxdelay ) cycle 

       M4_MODOBJ_GETREG(mat,reg)

       ! pre-calculate time signal for e-field modulation

       ncyc1 = 1.0*ncyc
       
       mat%wavefct = GaussianWave(ncyc1, mat%noffs, mat%natt, mat%nsus, mat%ndcy, &
            mat%gamma, mat%omega0)

      ! store time signal for delayed h-field modulation

      if ( mat%planewave ) then

         mat%signalp = mod(mat%signalp+1,mat%maxdelay)
         mat%signal(mat%signalp) =  mat%wavefct ! store signal function

      endif

   })
    
  end subroutine StepHMatJSource


!----------------------------------------------------------------------

  real(kind=8) function SumJEMatJSource(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc
    real(kind=8) :: Jx, Jy, Jz
    real(kind=8) :: wavefct(3), d
    integer :: di, l, n
   
    M4_MODLOOP_DECL({MATJSOURCE},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))

    sum = 0

    M4_MODLOOP_EXPR({MATJSOURCE},mat,{

    ! this loops over all mat structures, setting mat

       M4_MODOBJ_GETREG(mat,reg)

       if ( .not. mat%planewave ) then

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

       else ! plane wave!

          n = mat%signalp+mat%maxdelay

          M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

          if ( mask(i,j,k) ) then

             do l = 1,3 
                
                d = mat%delay(p) + mat%cdelay(l)
                di = int(d * real(mat%tres) + 0.5)
                wavefct(l) = mat%signal(mod(n-di, mat%maxdelay*mat%tres ))

                w(l) = w(l) * mat%finc(l)
                
             end do

M4_IFELSE_TM({
             Jx = - w(1) * epsinvx(i,j,k) * wavefct(1) 
             Jy = - w(2) * epsinvy(i,j,k) * wavefct(2)
})
M4_IFELSE_TE({
             Jz = - w(3) * epsinvz(i,j,k) * wavefct(3)
})
             sum = sum + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * real(Ex(i,j,k)) * Jx + },{0. +}) &
M4_IFELSE_TM({ M4_VOLEY(i,j,k) * real(Ey(i,j,k)) * Jy + },{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * real(Ez(i,j,k)) * Jz   },{0.  }) &
                )

          end if

         })

      end if

    })
    
    SumJEMatJSource = sum    

  end function SumJEMatJSource

!----------------------------------------------------------------------

  real(kind=8) function SumKHMatJSource(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
   
    SumKHMatJSource = 0.

  end function SumKHMatJSource


 !----------------------------------------------------------------------

  subroutine DisplayMatJSourceObj(mat)

    type(T_MATJSOURCE) :: mat
    
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
        " OMEGA=",TRIM(f2str(mat%omega0,4)),&
    	" OFFS=",TRIM(i2str(int(mat%noffs))),&
    	" ATT=",TRIM(i2str(int(mat%natt))),&
    	" SUS=",TRIM(i2str(int(mat%nsus))),&
    	" DCY=",TRIM(i2str(int(mat%ndcy))) })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatJSourceObj
  
 !----------------------------------------------------------------------

   subroutine EchoMatJSourceObj(mat)

    type(T_MATJSOURCE) :: mat

    M4_WRITE_INFO({"--- mat # ",TRIM(i2str(mat%idx))})
    M4_WRITE_INFO({"lambdainv0 = ",mat%lambdainv0   })
    M4_WRITE_INFO({"omega0  = ",mat%omega0 })
    M4_WRITE_INFO({"dn = ",mat%dn })
    M4_WRITE_INFO({"gamma = ",mat%gamma})
    M4_WRITE_INFO({"noffs/natt/nsus/ndcy  = ",mat%noffs,mat%natt,mat%nsus,mat%ndcy   })
    
    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    
  end subroutine EchoMatJSourceObj

!----------------------------------------------------------------------

end module matjsource

!
! Authors:  J.Hamm
! Modified: 20/12/2007
!
!======================================================================


