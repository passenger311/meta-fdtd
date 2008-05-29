!-*- F90 -*------------------------------------------------------------
!
!  module: srctfsfinj / meta
!
!  total-field-scattered-field box.
!
!----------------------------------------------------------------------


! =====================================================================
!
! The srctfsfinj module allows to setup a total-field scattered-field
! box region where a plane wave is fed in as incident wave. The total
! field then exists inside the box while the outside of the box 
! is associated with the scattered field.
! [see tavlov p 212 for calculation].
!

module srctfsfinj

  use constant
  use mpiworld
  use reglist
  use outlist
  use grid
  use signal
  use fdtd
  use fdtd_calc

  implicit none
  private
  save


  M4_SRCHEAD_DECL({SRCTFSFINJ},MAXSRCOBJ,{

     real(kind=8) :: lambdainv0                ! inverse vacuum wavelength in units of [2 pi c]
     real(kind=8) :: amp                       ! amplitude

     character(len=20) :: sigshape             ! signal shape

     real(kind=8) :: nhwhm                     ! time width of gaussian [dt]

     real(kind=8) :: gamma                      
     real(kind=8) :: omega0, domega 

     integer :: tres

     real(kind=8) :: noffs, natt, nsus, ndcy   ! generic signal parameters [dt]
     real(kind=8) :: nend                      ! end of signal [dt]

     real(kind=8) :: theta, phi, psi           ! angles of incident wavefront

     integer :: plane                          ! active plane

     real(kind=8) :: finc(6)                   ! field components of incident field 

     real(kind=8) :: kinc(3)                   ! normed k vector of plane wave
     real(kind=8) :: orig(3)                   ! orgin of the coordinate system

     real(kind=8) :: nrefr                     ! refractive index

     real(kind=8) :: epsinv,muinv              ! inverse material constant

     real(kind=8) :: wavefct                   ! wave function

     ! time delay field from point of origin

     real(kind=8), pointer, dimension(:,:,:) :: delay  

     ! component delay correction and inverse material constant (epsinv, muinv)

     real(kind=8) :: cdelay(6), pvel

     ! time signal lookup table

     real(kind=8), pointer, dimension(:) :: signal 
     integer :: signalp

     ! maximum delay: origin to point furthest away

     integer :: maxdelay

  })


contains

!----------------------------------------------------------------------

  subroutine ReadSrcTfsfInjObj(funit,lcount)

    M4_MODREAD_DECL({SRCTFSFINJ}, funit,lcount,src,reg,out)
    real(kind=8) :: v(4)
    logical :: eof, err
    character(len=LINELNG) :: line

    M4_WRITE_DBG(". enter ReadSrcTfsfInjObj")

    M4_MODREAD_EXPR({SRCTFSFINJ}, funit,lcount,src,reg, 2, out,{ 

    ! read parameters here, as defined in src data structure

    call readfloat(funit, lcount, src%lambdainv0)   ! inv. vacuum wavelength in units of [2 pi c]
    call readfloat(funit, lcount, src%amp)          ! amplitude 
    call readstring(funit, lcount, src%sigshape)     ! signal shape
    call readfloat(funit,lcount,src%nhwhm)          ! half width half max in time domain [dt]
    call readfloats(funit,lcount,v,4)               ! generic signal parameters [dt] 
    src%noffs = v(1)
    src%natt =  v(2)
    src%nsus =  v(3)
    src%ndcy =  v(4)
    
    call readline(funit,lcount,eof,line)

    err = .false.
    call getfloats(line,v,4,err)
    M4_SYNTAX_ERROR({err},lcount,{"[FLOATS]: phi,theta,psi,nrefr"})

    src%phi = v(1)
    src%theta = v(2)
    src%psi = v(3)
    src%nrefr = v(4)


    })

    M4_WRITE_DBG(". exit ReadSrcTfsfInjObj")

  end subroutine ReadSrcTfsfInjObj

!----------------------------------------------------------------------

  subroutine InitializeSrcTfsfInj

    M4_MODLOOP_DECL({SRCTFSFINJ},src)
    type(T_REG) :: reg
    integer :: i,j,k
    real(kind=8) :: r, rk(3), ksq, mdelay
    integer :: err, l 

    M4_WRITE_DBG(". enter InitializeSrcTfsfInj")
    M4_MODLOOP_EXPR({SRCTFSFINJ},src,{
    
    M4_MODOBJ_GETREG(src,reg)
 
    if ( .not. reg%isbox ) then
       M4_FATAL_ERROR("Region must be a single box!")
    endif

    src%plane = 0

    if ( reg%is .eq. reg%ie ) src%plane = 1
    if ( reg%js .eq. reg%je .and. DIM .ge. 2 ) src%plane = 3
    if ( reg%ks .eq. reg%ke .and. DIM .ge. 3 ) src%plane = 5


    if ( src%plane .eq. 0 ) then
       M4_FATAL_ERROR("Box region must define a proper plane!")
    end if

    ! center frequency
    src%omega0 = 2. * PI * src%lambdainv0


    src%maxdelay = 0
    
    ! normalised k vector
    src%kinc(1) = sin(DEG*src%theta)*cos(DEG*src%phi)
    src%kinc(2) = sin(DEG*src%theta)*sin(DEG*src%phi)
    src%kinc(3) = cos(DEG*src%theta)

    ! incident field/current amplitudes
    src%finc(1) = cos(DEG*src%psi)*sin(DEG*src%phi) - sin(DEG*src%psi)*cos(DEG*src%theta)*cos(DEG*src%phi)
    src%finc(2) = -cos(DEG*src%psi)*cos(DEG*src%phi) - sin(DEG*src%psi)*cos(DEG*src%theta)*sin(DEG*src%phi)
    src%finc(3) = sin(DEG*src%psi)*sin(DEG*src%theta)
    src%finc(4) = sin(DEG*src%psi)*sin(DEG*src%phi) + cos(DEG*src%psi)*cos(DEG*src%theta)*cos(DEG*src%phi)
    src%finc(5) = -sin(DEG*src%psi)*cos(DEG*src%phi) + cos(DEG*src%psi)*cos(DEG*src%theta)*sin(DEG*src%phi)
    src%finc(6) = -cos(DEG*src%psi)*sin(DEG*src%theta)

    ! decide on origin according to quadrant into which we emmit

    src%theta = abs(src%theta)

    if ( src%theta .gt. 180 ) then
       M4_FATAL_ERROR("0 <= Theta <= 180!")
    endif

    src%phi = mod(src%phi+360., 360.)
    src%psi = mod(src%psi+360., 360.)

    if ( DIM .eq. 3 ) then
       if ( src%theta .ge. 0. .and. src%theta .lt. 90. ) then
          src%orig(3) = reg%ks-1
       else
          src%orig(3) = reg%ke+1
          if ( src%plane .eq. 5 ) src%plane = 6
       endif
    else
        src%orig(3) = reg%ks
    end if


    if ( src%phi .ge. 0. .and. src%phi .lt. 90. ) then
       src%orig(1) = reg%is-1
       src%orig(2) = reg%js-1
    end if

    if ( src%phi .ge. 90. .and. src%phi .lt. 180 ) then
       src%orig(1) = reg%ie+1
       src%orig(2) = reg%js-1
       if ( src%plane .eq. 1 ) src%plane = 2
    end if

    if ( src%phi .ge. 180. .and. src%phi .lt. 270. ) then
       src%orig(1) = reg%ie+1
       src%orig(2) = reg%je+1
       if ( src%plane .eq. 1 ) src%plane = 2
       if ( src%plane .eq. 3 ) src%plane = 4
    end if
    
    if ( src%phi .ge. 270. .and. src%phi .lt. 360. ) then
       src%orig(1) = reg%is-1
       src%orig(2) = reg%je+1
       if ( src%plane .eq. 3 ) src%plane = 4
    end if
    
    ! optical delay of point-to-origin field

    i = src%orig(1)
    j = src%orig(2) 
    k = src%orig(3)  


    ! plane

    M4_WRITE_INFO({"plane = ",src%plane})

    ! calculate phase velocity

    src%pvel = 1.0 / src%nrefr
    call NumericalPhaseVelocity(src%kinc, src%omega0, src%nrefr, SX, SY, SZ, src%pvel)
    M4_WRITE_INFO({"phase velocity = ",src%pvel})

    ! be lazy and wasteful and do the whole box rather than just the border ...

    allocate(src%delay(reg%is-1:reg%ie+1,reg%js-1:reg%je+1,reg%ks-1:reg%ke+1), stat = err )
    M4_ALLOC_ERROR(err,{"InitializeSrcTfsfInj"})

    mdelay = 0

    do k = reg%ks-1, reg%ke+1, reg%dk
       do j = reg%js-1, reg%je+1, reg%dj
          do i = reg%is-1, reg%ie+1, reg%di

             ! project (i,j,k) location vector on kinc to create a distance field

             rk(1) = M4_DISTX(src%orig(1),i)
             rk(2) = M4_DISTY(src%orig(2),j)
             rk(3) = M4_DISTZ(src%orig(3),k)
       
             ! distance projected to kinc

             r = rk(1)*src%kinc(1) + rk(2)*src%kinc(2) +  rk(3)*src%kinc(3)

             ! time for signal at origin to arrive at (i,j,k)

             src%delay(i,j,k) = r / ( src%pvel * DT )

             ! maximum delay

             if ( src%delay(i,j,k) .gt. mdelay ) mdelay = src%delay(i,j,k)

          end do
       end do
    end do

    ! additional component delay times due to location of field components on staggered grid

    src%cdelay(1) = .5*SX*src%kinc(1) /  ( src%pvel * DT )
    src%cdelay(2) = .5*SY*src%kinc(2) /  ( src%pvel * DT )
    src%cdelay(3) = .5*SZ*src%kinc(3) /  ( src%pvel * DT )
    
    src%cdelay(4) = src%cdelay(2) + src%cdelay(3) 
    src%cdelay(5) = src%cdelay(1) + src%cdelay(3) 
    src%cdelay(6) = src%cdelay(1) + src%cdelay(2) 

    src%maxdelay = int(mdelay + 1) + sqrt( SX*SX + SY*SY + SZ*SZ )/( src%pvel * DT )

    src%tres = 100. ! increased time resolution of delay buffer by this factor 

    allocate(src%signal(0:src%maxdelay * src%tres -1) , stat = err )
    M4_ALLOC_ERROR(err,{"InitializeSrcTfsfInj"})

    src%signal = 0.

    src%signalp = 0 ! initialise signal index pointers

    src%nend = src%noffs + src%natt + src%nsus + src%ndcy + src%maxdelay


    M4_IFELSE_DBG({call EchoSrcTfsfInjObj(src)},{call DisplaySrcTfsfInjObj(src)})
      
    })

    M4_WRITE_DBG(". exit InitializeSrcTfsfInj")

  end subroutine InitializeSrcTfsfInj

!----------------------------------------------------------------------

  subroutine FinalizeSrcTfsfInj
    
    M4_MODLOOP_DECL({SRCTFSFINJ},src)

    M4_WRITE_DBG(". enter FinalizeSrcTfsfInj")
    
    M4_MODLOOP_EXPR({SRCTFSFINJ},src,{

       deallocate(src%signal, src%delay)

    })

    M4_WRITE_DBG(". exit FinalizeSrcTfsfInj")

  end subroutine FinalizeSrcTfsfInj

!----------------------------------------------------------------------

  subroutine StepESrcTfsfInj(ncyc)

    integer :: ncyc
    
    M4_MODLOOP_DECL({SRCTFSFINJ},src)
    type (T_REG) :: reg
    real(kind=8) :: ncyc1, ddt
    real(kind=8) :: val, zero = 0.
    integer :: l, sp


    M4_MODLOOP_EXPR({SRCTFSFINJ},src,{

       if ( ncyc .lt. src%noffs .or. ncyc .gt. src%nend ) cycle 

       M4_MODOBJ_GETREG(src,reg)

       if ( src%plane .eq. 1 ) then 
          ! i = is
          val = + src%finc(6) ! + Hz_inc
          call CalcEComp(src, Ey, -1,0,0, 6, zero, zero, val) 
          val = - src%finc(5) ! - Hy_inc
          call CalcEComp(src, Ez, -1,0,0, 5, zero, val, zero) 
 
       endif

       if ( src%plane .eq. 2 ) then 
          ! i = ie 
          val = - src%finc(6) ! - Hz_inc
          call CalcEComp(src, Ey, 0,0,0, 6, zero, zero, val) 
          val = + src%finc(5) ! + Hy_inc
          call CalcEComp(src, Ez, 0,0,0, 5, zero, val, zero) 
       end if

       if ( src%plane .eq. 3 ) then 
          ! j = js
          val = - src%finc(6) ! - Hz_inc
          call CalcEComp(src, Ex, 0,-1,0, 6, zero, zero, val)
          val = + src%finc(4) ! + Hx_inc
          call CalcEComp(src, Ez, 0,-1,0, 4, val, zero, zero)
       end if

       if ( src%plane .eq. 4 ) then 
          ! j = je
          val = + src%finc(6) ! + Hz_inc
          call CalcEComp(src, Ex, 0,0,0, 6, zero, zero, val)
          val = - src%finc(4) ! - Hx_inc
          call CalcEComp(src, Ez, 0,0,0, 4, val, zero, zero)
       end if
       
       if ( src%plane .eq. 5 ) then 
          ! k = ks
          val = + src%finc(5) ! + Hy_inc
          call CalcEComp(src, Ex, 0,0,-1, 5, zero, val, zero) 
          val = - src%finc(4) ! - Hx_inc
          call CalcEComp(src, Ey, 0,0,-1, 4, val, zero, zero) 
       end if

       if ( src%plane .eq. 6 ) then 
          ! k = ke
          val = - src%finc(5) ! - Hy_inc
          call CalcEComp(src, Ex, 0,0,0, 5, zero, val, zero) 
          val = + src%finc(4) ! - Hx_inc
          call CalcEComp(src, Ey, 0,0,0, 4, val, zero, zero) 
       end if

    })


    contains

      subroutine CalcEComp(src, E, o1,o2,o3, l, fx, fy, fz)

        type(T_SRCTFSFINJ) :: src
        M4_FTYPE, dimension(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)) :: E
        integer :: l, o1,o2,o3
        real(kind=8) :: fx, fy, fz
        
        M4_REGLOOP_DECL(reg,p,i,j,k,w(2))
        real(kind=8) :: wavefct, d, dd
        integer :: n, di

        n = src%signalp +  src%maxdelay * src%tres 

        M4_MODOBJ_GETREG(src,reg)

        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

                 ! time delay of inj field component in (i,j,k) from origin

                 d = src%delay(i+o1,j+o2,k+o3) + src%cdelay(l)
                 di = int(d * real(src%tres) + 0.5)
                 wavefct = src%signal(mod(n-di, src%maxdelay*src%tres ))
                 
                 ! update field component F from either fx or fy or fz
                 
                 E(i,j,k) = E(i,j,k) +  DT  * wavefct * w(2) * ( &
                      epsinvx(i,j,k) * fx/M4_SX(i,j,k) + &
                      epsinvy(i,j,k) * fy/M4_SY(i,j,k) + &
                      epsinvz(i,j,k) * fz/M4_SZ(i,j,k)  )

         })


      end subroutine CalcEComp

  end subroutine StepESrcTfsfInj

!----------------------------------------------------------------------

  subroutine StepHSrcTfsfInj(ncyc)


    integer :: ncyc
    
    M4_MODLOOP_DECL({SRCTFSFINJ},src)
    type(T_REG) :: reg
    real(kind=8) :: ddt, ncyc1
    real(kind=8) :: val, zero = 0.
    integer :: l, sp

    M4_MODLOOP_EXPR({SRCTFSFINJ},src,{

       if ( ncyc .lt.src%noffs .or. ncyc .gt. src%nend ) cycle 

       M4_MODOBJ_GETREG(src,reg)

       ! pre-calculate time signal for e-field modulation

       ddt = 1.0/real(src%tres)

       do l = 1, src%tres

!          ncyc1 = 1.0*ncyc  - 1.0 + l * ddt
          ncyc1 = 1.0*ncyc  - 0.5 + l * ddt  ! signal: n-1/2 -> n+1/2
          
          src%wavefct = src%amp * GenericWave(src%sigshape, ncyc1, src%noffs, src%natt, src%nsus, src%ndcy, & 
               src%nhwhm, src%omega0)

          ! store time signal for delayed e-field modulation

          sp = mod(src%signalp+l ,src%maxdelay * src%tres)
          src%signal(sp) = src%wavefct                ! store signal function

       end do

       src%signalp = sp

       ! enter e-field modulation loop

       if ( src%plane .eq. 1 ) then 
          ! i = is
          val = + src%finc(2) ! + Ey_inc
          call CalcHComp(src, Hz, 1,0,0, 2, zero, val, zero) 
          val = - src%finc(3) ! - Ez_inc
          call CalcHComp(src, Hy, 1,0,0, 3, zero, zero, val) 
       end if

       if ( src%plane .eq. 2 ) then 
          ! i = ie
          val = - src%finc(2) ! - Ey_inc
          call CalcHComp(src, Hz, 0,0,0, 2, zero, val, zero) 
          val = + src%finc(3) ! + Ez_inc
          call CalcHComp(src, Hy, 0,0,0, 3, zero, zero, val) 
       end if

       if ( src%plane .eq. 3 ) then 
          ! j = js
          val = - src%finc(1) ! - Ex_inc
          call CalcHComp(src, Hz, 0,1,0, 1, val, zero, zero)
          val = + src%finc(3) ! + Ez_inc
          call CalcHComp(src, Hx, 0,1,0, 3, zero, zero, val)
       end if

       if ( src%plane .eq. 4 ) then 
          ! j = je
          val = + src%finc(1) ! + Ex_inc
          call CalcHComp(src, Hz, 0,0,0, 1, val, zero, zero)
          val = - src%finc(3) ! - Ez_inc
          call CalcHComp(src, Hx, 0,0,0, 3, zero, zero, val)
       end if

       if ( src%plane .eq. 5 ) then 
          ! k = ks
          val = - src%finc(2) ! - Ey_inc
          call CalcHComp(src, Hx, 0,0,1, 2, val, zero, zero) 
          val = + src%finc(1) ! + Ex_inc
          call CalcHComp(src, Hy, 0,0,1, 1, zero, val, zero) 
       end if
       
       if ( src%plane .eq. 6 ) then 
          ! k = ke
          val = + src%finc(2) ! - Ey_inc
          call CalcHComp(src, Hx, 0,0,0, 2, val, zero, zero) 
          val = - src%finc(1) ! + Ex_inc
          call CalcHComp(src, Hy, 0,0,0, 1, zero, val, zero) 
       end if

    })


    contains

      subroutine CalcHComp(src, H, o1,o2,o3, l, fx, fy, fz)

        type(T_SRCTFSFINJ) :: src
        M4_FTYPE, dimension(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)) :: H
        integer :: l, o1,o2,o3
        real(kind=8) :: fx, fy, fz
        
        M4_REGLOOP_DECL(reg,p,i,j,k,w(2))
        real(kind=8) :: wavefct, d, dd
        integer :: n, di

        n = src%signalp + src%maxdelay * src%tres - 0.5 * src%tres
         
        M4_MODOBJ_GETREG(src,reg)

        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

                 ! time delay of inj field component in (i,j,k) from origin

                 d = src%delay(i,j,k) + src%cdelay(l)
                 di = int(d * real(src%tres) + 0.5)
                 wavefct = src%signal(mod(n-di , src%maxdelay * src%tres ))

                 ! update field component F from either fx or fy or fz

                 H(i-o1,j-o2,k-o3) = H(i-o1,j-o2,k-o3) +  DT * wavefct * w(1) * ( &
                      M4_MUINVX(i-o1,j-o2,k-o3) * fx/M4_SX(i,j,k) + &
                      M4_MUINVY(i-o1,j-o2,k-o3) * fy/M4_SY(i,j,k) + &
                      M4_MUINVZ(i-o1,j-o2,k-o3) * fz/M4_SZ(i,j,k)    ) 

         })

      end subroutine CalcHComp

  end subroutine StepHSrcTfsfInj

!----------------------------------------------------------------------

  real(kind=8) function SumJESrcTfsfInj(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc

    SumJESrcTfsfInj = 0.
    
  end function SumJESrcTfsfInj

!----------------------------------------------------------------------

  real(kind=8) function SumKHSrcTfsfInj(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc

    SumKHSrcTfsfInj = 0.

  end function SumKHSrcTfsfInj


 !----------------------------------------------------------------------

  subroutine DisplaySrcTfsfInjObj(src)

    type(T_SRCTFSFINJ) :: src
 
    M4_WRITE_INFO({"#",TRIM(i2str(src%idx)),&
        " OMEGA=",TRIM(f2str(src%omega0,4)),&
    	" OFFS=",TRIM(i2str(int(src%noffs))),&
    	" ATT=",TRIM(i2str(int(src%natt))),&
    	" SUS=",TRIM(i2str(int(src%nsus))),&
    	" DCY=",TRIM(i2str(int(src%ndcy))) })
    call DisplayRegObj(regobj(src%regidx))
    	
  end subroutine DisplaySrcTfsfInjObj
  
 !----------------------------------------------------------------------

   subroutine EchoSrcTfsfInjObj(src)

    type(T_SRCTFSFINJ) :: src

    M4_WRITE_INFO({"--- src # ",TRIM(i2str(src%idx))})
    M4_WRITE_INFO({"lambdainv0 = ",src%lambdainv0   })
    M4_WRITE_INFO({"omega0  = ",src%omega0 })
    M4_WRITE_INFO({"nhwhm = ",src%nhwhm })
    M4_WRITE_INFO({"noffs/natt/nsus/ndcy  = ",src%noffs,src%natt,src%nsus,src%ndcy   })
    
    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(src%regidx))
    
  end subroutine EchoSrcTfsfInjObj

!----------------------------------------------------------------------

end module srctfsfinj

!
! Authors:  J.Hamm, E.Kirby
! Modified: 26/04/2008
!
!======================================================================


