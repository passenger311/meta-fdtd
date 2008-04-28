!-*- F90 -*------------------------------------------------------------
!
!  module: mattfsource / meta
!
!  Uses the total-field scattered-field method to create a region that
!  acts like a transparent source.
!
!----------------------------------------------------------------------


! =====================================================================
!
! 

module mattfsource

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


  M4_MATHEAD_DECL({MATTFSOURCE},MAXMATOBJ,{

     real(kind=8) :: lambdainv0                ! inverse vacuum wavelength in units of [2 pi c]
     real(kind=8) :: amp                       ! amplitude
     real(kind=8) :: dn                        ! time width of gaussian [dt]

     real(kind=8) :: gamma                      
     real(kind=8) :: omega0, domega 

     integer :: tres

     real(kind=8) :: noffs, natt, nsus, ndcy   ! generic signal parameters [dt]
     real(kind=8) :: nend                      ! end of signal [dt]

     real(kind=8) :: theta, phi, psi           ! angles of incident wavefront

     integer, dimension(6) :: planeactive      ! active box planes

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

  subroutine ReadMatTfSourceObj(funit,lcount)

    M4_MODREAD_DECL({MATTFSOURCE}, funit,lcount,mat,reg,out)
    real(kind=8) :: v(4)
    logical :: eof, err
    real(kind=8) :: angles(3)
    character(len=LINELNG) :: line

    M4_WRITE_DBG(". enter ReadMatTfSourceObj")

    M4_MODREAD_EXPR({MATTFSOURCE}, funit,lcount,mat,reg, 2, out,{ 

    ! read parameters here, as defined in mat data structure

    call readfloat(funit, lcount, mat%lambdainv0)   ! inv. vacuum wavelength in units of [2 pi c]
    call readfloat(funit, lcount, mat%amp)          ! amplitude 
    call readfloat(funit,lcount,mat%dn)             ! half width of gaussian in time domain [dt]

    call readfloats(funit,lcount,v,4)               ! generic signal parameters [dt] 
    mat%noffs = v(1)
    mat%natt =  v(2)
    mat%nsus =  v(3)
    mat%ndcy =  v(4)
    
    ! read angles: phi, theta, psi

    angles = 0.

    call readline(funit,lcount,eof,line)

    err = .false.
    call getfloats(line,angles,3,err)
    M4_PARSE_ERROR({err},lcount,{bad format for plane wave angles})

    mat%phi = angles(1)
    mat%theta = angles(2)
    mat%psi = angles(3)

    ! call readintvec(funit, lcount, mat%planeactive, 6)

    })

    M4_WRITE_DBG(". exit ReadMatTfSourceObj")

  end subroutine ReadMatTfSourceObj

!----------------------------------------------------------------------

  subroutine InitializeMatTfSource

    M4_MODLOOP_DECL({MATTFSOURCE},mat)
    type(T_REG) :: reg
    integer :: i,j,k
    real(kind=8) :: r, rk(3), ksq, mdelay
    integer :: err, l 

    M4_WRITE_DBG(". enter InitializeMatTfSource")
    M4_MODLOOP_EXPR({MATTFSOURCE},mat,{
    
    M4_MODOBJ_GETREG(mat,reg)
 

    if ( .not. reg%isbox ) then
       M4_FATAL_ERROR("Region must be a single box!")
    endif

    ! center frequency
    mat%omega0 = 2. * PI * mat%lambdainv0

    ! derived gaussian parameters
    mat%gamma = sqrt( log(2.) ) / ( mat%dn * DT )   
    mat%domega = log(2.) * mat%gamma

    mat%maxdelay = 0
    
    ! normalised k vector
    mat%kinc(1) = sin(DEG*mat%theta)*cos(DEG*mat%phi)
    mat%kinc(2) = sin(DEG*mat%theta)*sin(DEG*mat%phi)
    mat%kinc(3) = cos(DEG*mat%theta)

    ! incident field/current amplitudes
    mat%finc(1) = cos(DEG*mat%psi)*sin(DEG*mat%phi) - sin(DEG*mat%psi)*cos(DEG*mat%theta)*cos(DEG*mat%phi)
    mat%finc(2) = -cos(DEG*mat%psi)*cos(DEG*mat%phi) - sin(DEG*mat%psi)*cos(DEG*mat%theta)*sin(DEG*mat%phi)
    mat%finc(3) = sin(DEG*mat%psi)*sin(DEG*mat%theta)
    mat%finc(4) = sin(DEG*mat%psi)*sin(DEG*mat%phi) + cos(DEG*mat%psi)*cos(DEG*mat%theta)*cos(DEG*mat%phi)
    mat%finc(5) = -sin(DEG*mat%psi)*cos(DEG*mat%phi) + cos(DEG*mat%psi)*cos(DEG*mat%theta)*sin(DEG*mat%phi)
    mat%finc(6) = -cos(DEG*mat%psi)*sin(DEG*mat%theta)

    ! decide on origin according to quadrant into which we emmit
    ! also decide on which tfsf plane per yee cell to activate. 

    mat%planeactive = 0.

    mat%theta = abs(mat%theta)

    if ( mat%theta .gt. 180 ) then
       M4_FATAL_ERROR("0 <= Theta <= 180!")
    endif

    mat%phi = mod(mat%phi+360., 360.)
    mat%psi = mod(mat%psi+360., 360.)


    if ( DIM .eq. 3 ) then
       if ( mat%theta .ge. 0. .and. mat%theta .lt. 90. ) then
          mat%orig(3) = reg%ks-1
          mat%planeactive(5) = 1
          mat%planeactive(6) = 0
       else
          mat%orig(3) = reg%ke+1
          mat%planeactive(5) = 0
          mat%planeactive(6) = 1
       endif
    else
        mat%orig(3) = reg%ks
    end if


    if ( mat%phi .ge. 0. .and. mat%phi .lt. 90. ) then
       mat%orig(1) = reg%is-1
       mat%orig(2) = reg%js-1
       mat%planeactive(1) = 1
       mat%planeactive(2) = 0
       mat%planeactive(3) = 1
       mat%planeactive(4) = 0
    end if

    if ( mat%phi .ge. 90. .and. mat%phi .lt. 180 ) then
       mat%orig(1) = reg%ie+1
       mat%orig(2) = reg%js-1
       mat%planeactive(1) = 0
       mat%planeactive(2) = 1
       mat%planeactive(3) = 1
       mat%planeactive(4) = 0
    end if

    if ( mat%phi .ge. 180. .and. mat%phi .lt. 270. ) then
       mat%orig(1) = reg%ie+1
       mat%orig(2) = reg%je+1
       mat%planeactive(1) = 0
       mat%planeactive(2) = 1
       mat%planeactive(3) = 0
       mat%planeactive(4) = 1
    end if
    
    if ( mat%phi .ge. 270. .and. mat%phi .lt. 360. ) then
       mat%orig(1) = reg%is-1
       mat%orig(2) = reg%je+1
       mat%planeactive(1) = 1
       mat%planeactive(2) = 0
       mat%planeactive(3) = 0
       mat%planeactive(4) = 1
    end if
    
    ! from here on everythhing is same as in the mattfsf module.
   
    ! optical delay of point-to-origin field

    i = mat%orig(1)
    j = mat%orig(2) 
    k = mat%orig(3)  

    ! get inverse material indices and refractive index at origin.

    if ( epsinvx(i,j,k) .ne. epsinvy(i,j,k) .or. epsinvx(i,j,k) .ne. epsinvz(i,j,k)  .or. &
         M4_MUINVX(i,j,k) .ne. M4_MUINVY(i,j,k) .or. M4_MUINVX(i,j,k) .ne. M4_MUINVZ(i,j,k)  ) then
       M4_FATAL_ERROR("Box origin is not set within isotropic material!")
    end if


    mat%epsinv = epsinvx(i,j,k)
    mat%muinv = M4_MUINVX(i,j,k)
    mat%nrefr = 1./sqrt( mat%epsinv * mat%muinv )

    ! scale planewave components to couple into dielectric material

    do l = 1,3

       mat%finc(l)   =  mat%amp * mat%finc(l)   * 1./sqrt(mat%muinv)
       mat%finc(l+3) =  mat%amp * mat%finc(l+3) * 1./sqrt(mat%epsinv)

    end do

    ! calculate phase velocity

    mat%pvel = 1.0 / mat%nrefr
    call NumericalPhaseVelocity(mat%kinc, mat%omega0, mat%nrefr, SX, SY, SZ, mat%pvel)
    M4_WRITE_INFO({"phase velocity = ",mat%pvel})

    ! be lazy and wasteful and do the whole box rather than just the border ...

    allocate(mat%delay(reg%is-1:reg%ie+1,reg%js-1:reg%je+1,reg%ks-1:reg%ke+1), stat = err )
    M4_ALLOC_ERROR(err,{"InitializeMatTfSource"})

    mdelay = 0

    do k = reg%ks-1, reg%ke+1, reg%dk
       do j = reg%js-1, reg%je+1, reg%dj
          do i = reg%is-1, reg%ie+1, reg%di

             ! project (i,j,k) location vector on kinc to create a distance field

             rk(1) = M4_DISTX(mat%orig(1),i)
             rk(2) = M4_DISTY(mat%orig(2),j)
             rk(3) = M4_DISTZ(mat%orig(3),k)
       
             ! distance projected to kinc

             r = rk(1)*mat%kinc(1) + rk(2)*mat%kinc(2) +  rk(3)*mat%kinc(3)

             ! time for signal at origin to arrive at (i,j,k)

             mat%delay(i,j,k) = r / ( mat%pvel * DT )

             ! maximum delay

             if ( mat%delay(i,j,k) .gt. mdelay ) mdelay = mat%delay(i,j,k)

          end do
       end do
    end do

    ! additional component delay times due to location of field components on staggered grid

    mat%cdelay(1) = .5*SX*mat%kinc(1) /  ( mat%pvel * DT )
    mat%cdelay(2) = .5*SY*mat%kinc(2) /  ( mat%pvel * DT )
    mat%cdelay(3) = .5*SZ*mat%kinc(3) /  ( mat%pvel * DT )
    
    mat%cdelay(4) = mat%cdelay(2) + mat%cdelay(3) 
    mat%cdelay(5) = mat%cdelay(1) + mat%cdelay(3) 
    mat%cdelay(6) = mat%cdelay(1) + mat%cdelay(2) 

    mat%maxdelay = int(mdelay + 1) + 1

    mat%tres = 100. ! increased time resolution of delay buffer by this factor 

    allocate(mat%signal(0:mat%maxdelay * mat%tres -1) , stat = err )
    M4_ALLOC_ERROR(err,{"InitializeMatTfSource"})

    mat%signal = 0.

    mat%signalp = 0 ! initialise signal index pointers

    mat%nend = mat%noffs + mat%natt + mat%nsus + mat%ndcy + mat%maxdelay


    M4_IFELSE_DBG({call EchoMatTfSourceObj(mat)},{call DisplayMatTfSourceObj(mat)})
      
    })

    M4_WRITE_DBG(". exit InitializeMatTfSource")

  end subroutine InitializeMatTfSource

!----------------------------------------------------------------------

  subroutine FinalizeMatTfSource
    
    M4_MODLOOP_DECL({MATTFSOURCE},mat)

    M4_WRITE_DBG(". enter FinalizeMatTfSource")
    
    M4_MODLOOP_EXPR({MATTFSOURCE},mat,{

       deallocate(mat%signal, mat%delay)

    })

    M4_WRITE_DBG(". exit FinalizeMatTfSource")

  end subroutine FinalizeMatTfSource

!----------------------------------------------------------------------

  subroutine StepEMatTfSource(ncyc)

    integer :: ncyc
    
    M4_MODLOOP_DECL({MATTFSOURCE},mat)
    type (T_REG) :: reg
    real(kind=8) :: ncyc1, ddt
    real(kind=8) :: val, zero = 0.
    integer :: l, sp


    M4_MODLOOP_EXPR({MATTFSOURCE},mat,{

       if ( ncyc .gt. mat%nend ) cycle 

       M4_MODOBJ_GETREG(mat,reg)


       if ( mat%planeactive(1) .ne. 0 ) then 
          ! i = is
          val = + mat%finc(6) ! + Hz_inc
          call CalcEComp(mat, Ey, -1,0,0, 6, zero, zero, val, 1) 
          val = - mat%finc(5) ! - Hy_inc
          call CalcEComp(mat, Ez, -1,0,0, 5, zero, val, zero, 1) 
       endif

       if ( mat%planeactive(2) .ne. 0 ) then 
          ! i = ie 
          val = - mat%finc(6) ! - Hz_inc
          call CalcEComp(mat, Ey, 0,0,0, 6, zero, zero, val, 1) 
          val = + mat%finc(5) ! + Hy_inc
          call CalcEComp(mat, Ez, 0,0,0, 5, zero, val, zero, 1) 
       end if

       if ( mat%planeactive(3) .ne. 0 ) then 
          ! j = js
          val = - mat%finc(6) ! - Hz_inc
          call CalcEComp(mat, Ex, 0,-1,0, 6, zero, zero, val, 2)
          val = + mat%finc(4) ! + Hx_inc
          call CalcEComp(mat, Ez, 0,-1,0, 4, val, zero, zero, 2)
       end if

       if ( mat%planeactive(4) .ne. 0 ) then 
          ! j = je
          val = + mat%finc(6) ! + Hz_inc
          call CalcEComp(mat, Ex, 0,0,0, 6, zero, zero, val, 2)
          val = - mat%finc(4) ! - Hx_inc
          call CalcEComp(mat, Ez, 0,0,0, 4, val, zero, zero, 2)
       end if
       
       if ( mat%planeactive(5) .ne. 0 ) then 
          ! k = ks
          val = + mat%finc(5) ! + Hy_inc
          call CalcEComp(mat, Ex, 0,0,-1, 5, zero, val, zero, 3) 
          val = - mat%finc(4) ! - Hx_inc
          call CalcEComp(mat, Ey, 0,0,-1, 4, val, zero, zero, 3) 
       end if

       if ( mat%planeactive(6) .ne. 0 ) then 
          ! k = ke
          val = - mat%finc(5) ! - Hy_inc
          call CalcEComp(mat, Ex, 0,0,0, 5, zero, val, zero, 3) 
          val = - mat%finc(4) ! - Hx_inc
          call CalcEComp(mat, Ey, 0,0,0, 4, val, zero, zero, 3) 
       end if

    })


    contains

      subroutine CalcEComp(mat, F, o1,o2,o3, l, fx, fy, fz, c)

        type(T_MATTFSOURCE) :: mat
        integer :: is,ie,js,je,ks,ke
        M4_FTYPE, dimension(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)) :: F
        integer :: l, o1,o2,o3, c
        real(kind=8) :: fx, fy, fz
        
        M4_REGLOOP_DECL(reg,p,i,j,k,w(2))
        real(kind=8) :: wavefct, d, dd
        integer :: n, di

        M4_MODOBJ_GETREG(mat,reg)

        n = mat%signalp +  mat%maxdelay * mat%tres 

        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
        
          ! time delay of inj field component in (i,j,k) from origin

          d = mat%delay(i+o1,j+o2,k+o3) + mat%cdelay(l)
          di = int(d * real(mat%tres) + 0.5)
          wavefct = mat%signal(mod(n-di, mat%maxdelay*mat%tres ))
          
          ! update field component F from either fx or fy or fz
          
          F(i,j,k) = F(i,j,k) +  DT * w(2) * mat%kinc(c) * wavefct * mat%epsinv * ( &
               fx/M4_SX(i,j,k) + &
               fy/M4_SY(i,j,k) + &
               fz/M4_SZ(i,j,k)  )
          
        })


      end subroutine CalcEComp

  end subroutine StepEMatTfSource

!----------------------------------------------------------------------

  subroutine StepHMatTfSource(ncyc)


    integer :: ncyc
    
    M4_MODLOOP_DECL({MATTFSOURCE},mat)
    type(T_REG) :: reg
    real(kind=8) :: ddt, ncyc1
    real(kind=8) :: val, zero = 0.
    integer :: l, sp

    M4_MODLOOP_EXPR({MATTFSOURCE},mat,{

       if ( ncyc .gt. mat%nend ) cycle 

       M4_MODOBJ_GETREG(mat,reg)

       ! pre-calculate time signal for e-field modulation

       ddt = 1.0/real(mat%tres)

       do l = 1, mat%tres

!          ncyc1 = 1.0*ncyc  - 1.0 + l * ddt
          ncyc1 = 1.0*ncyc  - 0.5 + l * ddt  ! signal: n-1/2 -> n+1/2
          
          mat%wavefct = GaussianWave(ncyc1, mat%noffs, mat%natt, mat%nsus, mat%ndcy, & 
               mat%gamma, mat%omega0)

          ! store time signal for delayed e-field modulation

          sp = mod(mat%signalp+l ,mat%maxdelay * mat%tres)
          mat%signal(sp) = mat%wavefct                ! store signal function

       end do

       mat%signalp = sp

       ! enter e-field modulation loop

       if ( mat%planeactive(1) .ne. 0 ) then 
          ! i = is
          val = + mat%finc(2) ! + Ey_inc
          call CalcHComp(mat, Hz, -1,0,0, 2, zero, val, zero, 1) 
          val = - mat%finc(3) ! - Ez_inc
          call CalcHComp(mat, Hy, -1,0,0, 3, zero, zero, val, 1) 
       end if

       if ( mat%planeactive(2) .ne. 0 ) then 
          ! i = ie
          val = - mat%finc(2) ! - Ey_inc
          call CalcHComp(mat, Hz, 0,0,0, 2, zero, val, zero, 1) 
          val = + mat%finc(3) ! + Ez_inc
          call CalcHComp(mat, Hy, 0,0,0, 3, zero, zero, val, 1) 
       end if

       if ( mat%planeactive(3) .ne. 0 ) then 
          ! j = js
          val = - mat%finc(1) ! - Ex_inc
          call CalcHComp(mat, Hz, 0,-1,0, 1, val, zero, zero, 2)
          val = + mat%finc(3) ! + Ez_inc
          call CalcHComp(mat, Hx, 0,-1,0, 3, zero, zero, val, 2)
       end if

       if ( mat%planeactive(4) .ne. 0 ) then 
          ! j = je
          val = + mat%finc(1) ! + Ex_inc
          call CalcHComp(mat, Hz, 0,0,0, 1, val, zero, zero, 2)
          val = - mat%finc(3) ! - Ez_inc
          call CalcHComp(mat, Hx, 0,0,0, 3, zero, zero, val, 2)
       end if

       if ( mat%planeactive(5) .ne. 0 ) then 
          ! k = ks
          val = + mat%finc(1) ! + Ex_inc
          call CalcHComp(mat, Hx, 0,0,-1, 1, val, zero, zero, 3) 
          val = - mat%finc(2) ! - Ey_inc
          call CalcHComp(mat, Hy, 0,0,-1, 2, zero, val, zero, 3) 
       end if
       
       if ( mat%planeactive(6) .ne. 0 ) then 
          ! k = ke
          val = - mat%finc(1) ! + Ex_inc
          call CalcHComp(mat, Hx, 0,0,0, 1, val, zero, zero, 3) 
          val = + mat%finc(2) ! - Ey_inc
          call CalcHComp(mat, Hy, 0,0,0, 2, zero, val, zero, 3) 
       end if

    })


    contains

      subroutine CalcHComp(mat, F, o1,o2,o3, l, fx, fy, fz, c)

        type(T_MATTFSOURCE) :: mat
        integer :: is,ie,js,je,ks,ke
        M4_FTYPE, dimension(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)) :: F
        integer :: l, o1,o2,o3, c
        real(kind=8) :: fx, fy, fz
        
        M4_REGLOOP_DECL(reg,p,i,j,k,w(2))
        real(kind=8) :: wavefct, d, dd
        integer :: n, di

        M4_MODOBJ_GETREG(mat,reg)

        n = mat%signalp + mat%maxdelay * mat%tres - 0.5 * mat%tres
         
        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

          ! time delay of inj field component in (i,j,k) from origin

          d = mat%delay(i,j,k) + mat%cdelay(l)
          di = int(d * real(mat%tres) + 0.5)
          wavefct = mat%signal(mod(n-di , mat%maxdelay * mat%tres ))

          ! update field component F from either fx or fy or fz

          F(i+o1,j+o2,k+o3) = F(i+o1,j+o2,k+o3) +  DT * w(1) * mat%kinc(c) * wavefct * mat%muinv * ( &
               fx/M4_SX(i,j,k) + &
               fy/M4_SY(i,j,k) + &
               fz/M4_SZ(i,j,k)    ) 

        })


      end subroutine CalcHComp

  end subroutine StepHMatTfSource

!----------------------------------------------------------------------

  real(kind=8) function SumJEMatTfSource(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc

    SumJEMatTfSource = 0.
    
  end function SumJEMatTfSource

!----------------------------------------------------------------------

  real(kind=8) function SumKHMatTfSource(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc

    SumKHMatTfSource = 0.

  end function SumKHMatTfSource


 !----------------------------------------------------------------------

  subroutine DisplayMatTfSourceObj(mat)

    type(T_MATTFSOURCE) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
        " OMEGA=",TRIM(f2str(mat%omega0,4)),&
    	" OFFS=",TRIM(i2str(int(mat%noffs))),&
    	" ATT=",TRIM(i2str(int(mat%natt))),&
    	" SUS=",TRIM(i2str(int(mat%nsus))),&
    	" DCY=",TRIM(i2str(int(mat%ndcy))) })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatTfSourceObj
  
 !----------------------------------------------------------------------

   subroutine EchoMatTfSourceObj(mat)

    type(T_MATTFSOURCE) :: mat

    M4_WRITE_INFO({"--- mat # ",TRIM(i2str(mat%idx))})
    M4_WRITE_INFO({"lambdainv0 = ",mat%lambdainv0   })
    M4_WRITE_INFO({"omega0  = ",mat%omega0 })
    M4_WRITE_INFO({"dn = ",mat%dn })
    M4_WRITE_INFO({"gamma = ",mat%gamma})
    M4_WRITE_INFO({"noffs/natt/nsus/ndcy  = ",mat%noffs,mat%natt,mat%nsus,mat%ndcy   })
    
    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    
  end subroutine EchoMatTfSourceObj

!----------------------------------------------------------------------

end module mattfsource

!
! Authors:  J.Hamm
! Modified: 26/04/2008
!
!======================================================================

