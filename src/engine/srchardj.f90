!-*- F90 -*------------------------------------------------------------
!
!  module: srchardj / meta
!
!  feed electromagnetic field with hard dipole current source 
!
!----------------------------------------------------------------------


! =====================================================================
!
! The SrcHardJ module allows to add a hard source current to the 
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

module srchardj

  use constant
  use checkpoint
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


  M4_SRCHEAD_DECL({SRCHARDJ},MAXSRCOBJ,{

     real(kind=8) :: lambdainv0                ! inverse vacuum wavelength in units of [2 pi c]
     real(kind=8) :: amp                       ! amplitude

     character(len=20) :: sigshape             ! signal shape

     real(kind=8) :: nhwhm                     ! time width of gaussian [dt]

     real(kind=8) :: gamma                      
     real(kind=8) :: omega0, domega 

     real(kind=8) :: noffs, natt, nsus, ndcy   ! generic signal parameters [dt]
     real(kind=8) :: nend                      ! end of signal [dt]
     real(kind=8) :: alpha

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

  subroutine ReadSrcHardJObj(funit,lcount)


    M4_MODREAD_DECL({SRCHARDJ}, funit,lcount,src,reg,out)
    real(kind=8) :: v(4)
    logical :: eof, err
    character(len=LINELNG) :: line

    M4_WRITE_DBG(". enter ReadSrcHardJObj")

    M4_MODREAD_EXPR({SRCHARDJ}, funit,lcount,src,reg, 3, out,{ 

    ! read parameters here, as defined in src data structure


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
    
    ! optional: configure plane wave source? 

    err = .false.
    src%phi = 0.
    src%theta = 0.
    src%psi = 0.

    call readlogical(funit,lcount,src%planewave)

    call readline(funit,lcount,eof,line)

    err = .false.
    call getfloats(line,v,4,err)

    M4_SYNTAX_ERROR({err},lcount,{"[FLOATS] : phi,theta,psi,nrefr"})

    src%phi = v(1)
    src%theta = v(2)
    src%psi = v(3)
    src%nrefr = v(4)
    
    call readfloat(funit,lcount,src%alpha)
    
    })

    M4_WRITE_DBG(". exit ReadSrcHardJObj")

  end subroutine ReadSrcHardJObj

!----------------------------------------------------------------------

  subroutine InitializeSrcHardJ

    M4_MODLOOP_DECL({SRCHARDJ},src)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8) :: r, rk(3), mdelay
    integer :: err

    M4_WRITE_DBG(". enter InitializeSrcHardJ")
    M4_MODLOOP_EXPR({SRCHARDJ},src,{
    
    M4_MODOBJ_GETREG(src,reg)
 
    ! center frequency
    src%omega0 = 2. * PI * src%lambdainv0

    src%maxdelay = 0
    
    if ( src%planewave ) then

       ! normalised k vector
       src%kinc(1) = sin(DEG*src%theta)*cos(DEG*src%phi)
       src%kinc(2) = sin(DEG*src%theta)*sin(DEG*src%phi)
       src%kinc(3) = cos(DEG*src%theta)

       ! incident field/current amplitudes
       src%finc(1) = cos(DEG*src%psi)*sin(DEG*src%phi) - sin(DEG*src%psi)*cos(DEG*src%theta)*cos(DEG*src%phi)
       src%finc(2) = -cos(DEG*src%psi)*cos(DEG*src%phi) - sin(DEG*src%psi)*cos(DEG*src%theta)*sin(DEG*src%phi)
       src%finc(3) = sin(DEG*src%psi)*sin(DEG*src%theta)

       ! decide on origin according to quadrant into which we emmit

       src%theta = abs(src%theta)
       
       if ( src%theta .gt. 180 ) then
          M4_FATAL_ERROR("0 <= Theta <= 180!")
       endif
       
       src%phi = mod(real(src%phi+360.), real(360.))
       src%psi = mod(real(src%psi+360.), real(360.))

       if ( src%theta .ge. 0. .and. src%theta .lt. 90. ) then
          src%orig(3) = reg%ks
       else
          src%orig(3) = reg%ke
       endif

       if ( src%phi .ge. 0. .and. src%phi .lt. 90. ) then
          src%orig(1) = reg%is
          src%orig(2) = reg%js
       end if

       if ( src%phi .ge. 90. .and. src%phi .lt. 180 ) then
          src%orig(1) = reg%ie
          src%orig(2) = reg%js
       end if

       if ( src%phi .ge. 180. .and. src%phi .lt. 270. ) then
          src%orig(1) = reg%ie
          src%orig(2) = reg%je
       end if

       if ( src%phi .ge. 270. .and. src%phi .lt. 360. ) then
          src%orig(1) = reg%is
          src%orig(2) = reg%je
       end if

       ! calculate phase velocity

       src%pvel = 1.0 / src%nrefr
       call NumericalPhaseVelocity(src%kinc, src%omega0, src%nrefr, SX, SY, SZ, src%pvel)
       M4_WRITE_INFO({"phase velocity = ",src%pvel})


       allocate(src%delay(reg%numnodes), stat = err )
       M4_ALLOC_ERROR(err,{"InitializeSrcHardJ"})

       mdelay = 0
       
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
             ! project (i,j,k) location vector on kinc to create a distance field

             rk(1) = M4_DISTX(src%orig(1),i)
             rk(2) = M4_DISTY(src%orig(2),j)
             rk(3) = M4_DISTZ(src%orig(3),k)
       
             ! distance projected to kinc

             r = rk(1)*src%kinc(1) + rk(2)*src%kinc(2) +  rk(3)*src%kinc(3)

             ! time for signal at origin to arrive at (i,j,k)

             src%delay(p) = r / ( src%pvel * DT )

             ! maximum delay

             if ( src%delay(p) .gt. mdelay ) mdelay = src%delay(p)

       })


       src%cdelay(1) = .5*SX*src%kinc(1) /  ( src%pvel * DT )
       src%cdelay(2) = .5*SY*src%kinc(2) /  ( src%pvel * DT )
       src%cdelay(3) = .5*SZ*src%kinc(3) /  ( src%pvel * DT )

       src%maxdelay = int(mdelay + 1) + sqrt( SX*SX + SY*SY + SZ*SZ )/( src%pvel * DT )

       src%tres = 100. ! increased time resolution of delay buffer by this factor 

       allocate(src%signal(0:src%maxdelay * src%tres -1) , stat = err )
       M4_ALLOC_ERROR(err,{"InitializeSrcTfSource"})

       src%signal = 0.

       src%signalp = 0 ! initialise signal index pointers

       if ( load_state .and. ( detail_level .eq. 1 .or. detail_level .eq. 3 ) ) then

          read(UNITCHK) src%signal
          read(UNITCHK) src%signalp

       end if

    end if

    src%nend = src%noffs +src%natt + src%nsus + src%ndcy + src%maxdelay

    M4_IFELSE_DBG({call EchoSrcHardJObj(src)},{call DisplaySrcHardJObj(src)})
      
    })

    M4_WRITE_DBG(". exit InitializeSrcHardJ")

  end subroutine InitializeSrcHardJ

!----------------------------------------------------------------------

  subroutine FinalizeSrcHardJ

    M4_MODLOOP_DECL({SRCHARDJ},src)
    M4_WRITE_DBG(". enter FinalizeSrcHardJ")
    M4_MODLOOP_EXPR({SRCHARDJ},src,{

      ! finalize src object here
  
      if ( src%planewave ) then

         if ( save_state .and. ( detail_level .eq. 1 .or. detail_level .eq. 3 ) ) then

            write(UNITCHK) src%signal
            write(UNITCHK) src%signalp
 
         end if

         deallocate(src%signal, src%delay)

      end if

    })

    M4_WRITE_DBG(". exit FinalizeSrcHardJ")

  end subroutine FinalizeSrcHardJ

!----------------------------------------------------------------------

  subroutine StepESrcHardJ(ncyc)

    integer :: ncyc
    
    M4_MODLOOP_DECL({SRCHARDJ},src)

    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8) :: ncyc1
    real(kind=8) :: wavefct(3), d
    integer :: di, l, n

    M4_MODLOOP_EXPR({SRCHARDJ},src,{

      if ( ncyc .lt. src%noffs .or. ncyc .gt. src%nend ) cycle 

      M4_MODOBJ_GETREG(src,reg)

      ! enter e-field modulation loop

      if ( .not. src%planewave ) then

        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

M4_IFELSE_TM({
           Ex(i,j,k) = Ex(i,j,k) + DT * w(1) * epsinvx(i,j,k) * src%wavefct
           Ey(i,j,k) = Ey(i,j,k) + DT * w(2) * epsinvy(i,j,k) * src%wavefct
})
M4_IFELSE_TE({
           Ez(i,j,k) = Ez(i,j,k) + DT * w(3) * epsinvz(i,j,k) * src%wavefct
})

        })

      else ! plane wave!

         n = src%signalp+src%maxdelay*src%tres

         M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

         do l = 1,3 
            
            d = src%delay(p) + src%cdelay(l)
            di = int(d * real(src%tres) + 0.5)
            wavefct(l) = src%signal(mod(n-di, src%maxdelay*src%tres ))            

!            write(6,*) src%delay(p),src%cdelay(l),src%maxdelay 

!            w(l) = w(l) * src%finc(l)

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

  end subroutine StepESrcHardJ

!----------------------------------------------------------------------

  subroutine StepHSrcHardJ(ncyc)

    integer :: ncyc

    M4_MODLOOP_DECL({SRCHARDJ},src)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8) :: ncyc1, ddt
    integer :: sp, l

    M4_MODLOOP_EXPR({SRCHARDJ},src,{

       if ( ncyc - src%noffs .gt. src%nend + src%maxdelay ) cycle 

       M4_MODOBJ_GETREG(src,reg)

       ! pre-calculate time signal for e-field modulation

      

      ! store time signal for delayed h-field modulation

      if ( .not. src%planewave ) then

         ncyc1 = 1.0*ncyc
         
         src%wavefct = src%amp * GenericWave(src%sigshape, ncyc1, src%noffs, src%natt, src%nsus, src%ndcy, &
              src%nhwhm, src%omega0, src%alpha)
         
      else

         ddt = 1.0/real(src%tres)

         do l = 1, src%tres

            ncyc1 = 1.0*ncyc  + l * ddt 
          
            src%wavefct =  src%amp * GenericWave(src%sigshape, ncyc1, src%noffs, src%natt, src%nsus, src%ndcy, & 
                 src%nhwhm, src%omega0, src%alpha)

            ! store time signal for delayed e-field modulation

            sp = mod(src%signalp+l ,src%maxdelay * src%tres)
            src%signal(sp) = src%wavefct                ! store signal function

         end do

         src%signalp = sp

      endif

   })
    
  end subroutine StepHSrcHardJ


!----------------------------------------------------------------------

  real(kind=8) function SumJESrcHardJ(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc
    real(kind=8) :: Jx, Jy, Jz
    real(kind=8) :: wavefct(3), d
    integer :: di, l, n
   
    M4_MODLOOP_DECL({SRCHARDJ},src)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))

    sum = 0

    M4_MODLOOP_EXPR({SRCHARDJ},src,{

       if ( ncyc .lt. src%noffs .or. ncyc .gt. src%nend ) cycle 

    ! this loops over all src structures, setting src

       M4_MODOBJ_GETREG(src,reg)

       if ( .not. src%planewave ) then

          M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
          ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

          ! J(*,m) is P(n+1) and J(*,n) is P(n)      
          
          if ( mask(i,j,k) ) then

             Jx = - w(1) * src%wavefct
             Jy = - w(2) * src%wavefct
             Jz = - w(3) * src%wavefct
          
             sum = sum + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * dble(Ex(i,j,k)) * Jx + },{0. +}) &
M4_IFELSE_TM({ M4_VOLEY(i,j,k) * dble(Ey(i,j,k)) * Jy + },{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * dble(Ez(i,j,k)) * Jz   },{0.  }) &
                 )
             
          endif

          })      

       else ! plane wave!

          n = src%signalp+src%maxdelay*src%tres

          M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

          if ( mask(i,j,k) ) then

             do l = 1,3 
                
                d = src%delay(p) + src%cdelay(l)
                di = int(d * real(src%tres) + 0.5)
                wavefct(l) = src%signal(mod(n-di, src%maxdelay*src%tres ))

               ! w(l) = w(l) * src%finc(l)
                
             end do

M4_IFELSE_TM({
             Jx = - w(1) * epsinvx(i,j,k) * wavefct(1) 
             Jy = - w(2) * epsinvy(i,j,k) * wavefct(2)
})
M4_IFELSE_TE({
             Jz = - w(3) * epsinvz(i,j,k) * wavefct(3)
})
             sum = sum + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * dble(Ex(i,j,k)) * Jx + },{0. +}) &
M4_IFELSE_TM({ M4_VOLEY(i,j,k) * dble(Ey(i,j,k)) * Jy + },{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * dble(Ez(i,j,k)) * Jz   },{0.  }) &
                )

          end if

         })

      end if

    })
    
    SumJESrcHardJ = sum    

  end function SumJESrcHardJ

!----------------------------------------------------------------------

  real(kind=8) function SumKHSrcHardJ(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
   
    SumKHSrcHardJ = 0.

  end function SumKHSrcHardJ


 !----------------------------------------------------------------------

  subroutine DisplaySrcHardJObj(src)

    type(T_SRCHARDJ) :: src
    
    M4_WRITE_INFO({"#",TRIM(i2str(src%idx)),&
        " OMEGA=",TRIM(f2str(src%omega0,4)),&
    	" OFFS=",TRIM(i2str(int(src%noffs))),&
    	" ATT=",TRIM(i2str(int(src%natt))),&
    	" SUS=",TRIM(i2str(int(src%nsus))),&
    	" DCY=",TRIM(i2str(int(src%ndcy))) })
    call DisplayRegObj(regobj(src%regidx))
    	
  end subroutine DisplaySrcHardJObj
  
 !----------------------------------------------------------------------

   subroutine EchoSrcHardJObj(src)

    type(T_SRCHARDJ) :: src

    M4_WRITE_INFO({"--- src # ",TRIM(i2str(src%idx))})
    M4_WRITE_INFO({"lambdainv0 = ",src%lambdainv0   })
    M4_WRITE_INFO({"omega0  = ",src%omega0 })
    M4_WRITE_INFO({"hwhm = ",src%nhwhm })
    M4_WRITE_INFO({"noffs/natt/nsus/ndcy  = ",src%noffs,src%natt,src%nsus,src%ndcy   })
    
    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(src%regidx))
    
  end subroutine EchoSrcHardJObj

!----------------------------------------------------------------------

end module srchardj

!
! Authors:  J.Hamm
! Modified: 20/12/2007
!
!======================================================================


