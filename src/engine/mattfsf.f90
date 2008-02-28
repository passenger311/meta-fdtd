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
  use signal
  use fdtd
  use fdtd_calc

  implicit none
  private
  save


  M4_MATHEAD_DECL({MATTFSF},MAXMATOBJ,{

     real(kind=8) :: lambdainv0                ! inverse vacuum wavelength in units of [2 pi c]
     real(kind=8) :: dn                        ! time width of gaussian [dt]

     real(kind=8) :: gamma                      
     real(kind=8) :: omega0, domega 

     real(kind=8) :: noffs, natt, nsus, ndcy   ! generic signal parameters [dt]
     real(kind=8) :: nend                      ! end of signal [dt]

     real(kind=8) :: theta, phi, psi           ! angles of incident wavefront

     real(kind=8) :: wavefcte, wavefcth        ! wave functions
     real(kind=8) :: finc(6)                   ! field components of incident field 

     logical :: planewave                      ! plane wave mode?
     real(kind=8) :: kinc(3)                   ! normed k vector of plane wave
     real(kind=8) :: orig(3)                   ! orgin of the coordinate system
     real(kind=8) :: pvel                      ! phase velocity

     ! time delay field from point of origin
     real(kind=8), pointer, dimension(:) :: delay  

     ! component delay correction
     real(kind=8) :: cdelay(6)

     ! time signal lookup table
     real(kind=8), pointer, dimension(:) :: signale, signalh 
     integer :: signalep, signalhp

     ! position shof of h field 
     integer, pointer, dimension(:,:) :: orient

     ! maximum delay: origin to point furthest away
     integer :: maxdelay

  })


contains

!----------------------------------------------------------------------

  subroutine ReadMatTfsfObj(funit,lcount)


    M4_MODREAD_DECL({MATTFSF}, funit,lcount,mat,reg,out)
    real(kind=8) :: v(4)
    logical :: eof, err
    real(kind=8) :: angles(3)
    character(len=LINELNG) :: line

    M4_WRITE_DBG(". enter ReadMatTfsfObj")

    M4_MODREAD_EXPR({MATTFSF}, funit,lcount,mat,reg, 3, out,{ 

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

    M4_WRITE_DBG(". exit ReadMatTfsfObj")

  end subroutine ReadMatTfsfObj

!----------------------------------------------------------------------

  subroutine InitializeMatTfsf

   M4_MODLOOP_DECL({MATTFSF},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))
    real(kind=8) :: r, in(3), mdelay
    integer :: err

    M4_WRITE_DBG(". enter InitializeMatTfsf")
    M4_MODLOOP_EXPR({MATTFSF},mat,{
    
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
       mat%finc(4) = sin(DEG*mat%psi)*sin(DEG*mat%phi) + cos(DEG*mat%psi)*cos(DEG*mat%theta)*cos(DEG*mat%phi)
       mat%finc(5) = -sin(DEG*mat%psi)*cos(DEG*mat%phi) + cos(DEG*mat%psi)*cos(DEG*mat%theta)*sin(DEG*mat%phi)
       mat%finc(6) = -cos(DEG*mat%psi)*sin(DEG*mat%theta)

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

       in(1) = sqrt( epsinvx(i,j,k) * M4_MUINVX(i,j,k) )  ! assumed to be homogeneous over source loc.
       in(2) = sqrt( epsinvy(i,j,k) * M4_MUINVY(i,j,k) )
       in(3) = sqrt( epsinvz(i,j,k) * M4_MUINVZ(i,j,k) )

       allocate(mat%delay(reg%numnodes), mat%orient(reg%numnodes,3), stat = err )
       M4_ALLOC_ERROR(err,{"InitializeMatTfsf"})

       mdelay = 0

       mat%orient = 0 
       
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! project (i,j,k) location vector on kinc to create a distance field
       
       r = in(1)*M4_DISTX(mat%orig(1),i)*mat%kinc(1) + in(2)*M4_DISTY(mat%orig(2),j)*mat%kinc(2) +  in(3)*M4_DISTZ(mat%orig(3),k)*mat%kinc(3)

       mat%delay(p) = r / ( mat%pvel * DT )

       if ( mat%delay(p) .gt. mdelay ) mdelay = mat%delay(p)

       if ( w(2) .lt. 0 .and. w(3) .lt. 0 ) mat%orient(p,1) = 1 ! flip orientations
       if ( w(1) .lt. 0 .and. w(2) .lt. 0 ) mat%orient(p,3) = 1
       if ( w(1) .lt. 0 .and. w(3) .lt. 0 ) mat%orient(p,2) = 1

       })

       ! additional component delay times due to location of field components on staggered grid

       mat%cdelay(1) = in(1)*.5*SX*mat%kinc(1) / ( mat%pvel * DT )
       mat%cdelay(2) = in(2)*.5*SY*mat%kinc(2) / ( mat%pvel * DT )
       mat%cdelay(3) = in(3)*.5*SZ*mat%kinc(3) / ( mat%pvel * DT )

       mat%cdelay(4) = mat%cdelay(2) + mat%cdelay(3) 
       mat%cdelay(5) = mat%cdelay(1) + mat%cdelay(3) 
       mat%cdelay(6) = mat%cdelay(1) + mat%cdelay(2) 

       mdelay =  mdelay  + sqrt( SX*SX + SY*SY + SZ*SZ )

       mat%maxdelay = int(mdelay + 0.5) + 1

       allocate(mat%signale(0:mat%maxdelay-1),mat%signalh(0:mat%maxdelay-1) , stat = err )
       M4_ALLOC_ERROR(err,{"InitializeMatTfsf"})

       mat%signale = 0.
       mat%signalh = 0.

       mat%signalep = 0 ! initialise signal index pointers
       mat%signalhp = 0

    end if

    mat%nend = mat%natt + mat%nsus + mat%ndcy + mat%maxdelay


    M4_IFELSE_DBG({call EchoMatTfsfObj(mat)},{call DisplayMatTfsfObj(mat)})
      
    })

    M4_WRITE_DBG(". exit InitializeMatTfsf")

  end subroutine InitializeMatTfsf

!----------------------------------------------------------------------

  subroutine FinalizeMatTfsf
    
    M4_MODLOOP_DECL({MATTFSF},mat)

    M4_WRITE_DBG(". enter FinalizeMatTfsf")
    
    M4_MODLOOP_EXPR({MATTFSF},mat,{

    if ( mat%planewave ) then

       deallocate(mat%signale, mat%signalh, mat%delay, mat%orient)

    end if

    })

    M4_WRITE_DBG(". exit FinalizeMatTfsf")

  end subroutine FinalizeMatTfsf

!----------------------------------------------------------------------

  subroutine StepEMatTfsf(ncyc)

    integer :: ncyc
    
    M4_MODLOOP_DECL({MATTFSF},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))
    real(kind=8) :: ncyc1
    real(kind=8) :: wavefct(3), d(3), dd(3)
    integer :: di(3), l, o(3), n

    M4_MODLOOP_EXPR({MATTFSF},mat,{

       if ( ncyc .gt. mat%nend ) cycle 

       M4_MODOBJ_GETREG(mat,reg)

       ! pre-calculate time signal for h-field modulation

       ncyc1 = 1.0*ncyc + 0.5
       
       mat%wavefcth = GaussianWave(ncyc1, mat%noffs, mat%natt, mat%nsus, mat%ndcy, & 
            mat%gamma, mat%omega0)

       ! store time signal for delayed h-field modulation

       if ( mat%planewave ) then

         mat%signalhp = mod(mat%signalhp+1,mat%maxdelay)
         mat%signalh(mat%signalhp) =  mat%wavefcth                ! store signal function

      endif

      ! enter e-field modulation loop

      n = mat%signalep+mat%maxdelay

      if ( .not. mat%planewave ) then

        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

M4_IFELSE_TM({
         Ex(i,j,k) = Ex(i,j,k) + DT * ( w(5)/M4_SZ(i,j,k) - w(6)/M4_SY(i,j,k) ) * epsinvx(i,j,k) * mat%wavefcte
         Ey(i,j,k) = Ey(i,j,k) + DT * ( w(6)/M4_SX(i,j,k) - w(4)/M4_SZ(i,j,k) ) * epsinvy(i,j,k) * mat%wavefcte
})
M4_IFELSE_TE({
         Ez(i,j,k) = Ez(i,j,k) + DT * ( w(4)/M4_SY(i,j,k) - w(5)/M4_SX(i,j,k) ) * epsinvz(i,j,k) * mat%wavefcte
})
        })

      else ! plane wave!

         M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

         do l = 1,3 
            
            d(l) = mat%delay(p) + mat%cdelay(l) ! time delay from origin
            di(l) = int(d(l))
            dd(l) = d(l) - di(l)
            ! use linear interpolation to read time signal
            wavefct(l) = (1. - dd(l)) * mat%signale(mod(n-di(l),mat%maxdelay)) + &
                 dd(l) * mat%signale(mod(n-di(l)-1,mat%maxdelay))

            w(3+l) = w(3+l) * mat%finc(3+l)
            w(l) = abs(w(l))

         end do

M4_IFELSE_TM({
         Ex(i,j,k) = Ex(i,j,k) + w(1) * DT * ( w(5)/M4_SZ(i,j,k) - w(6)/M4_SY(i,j,k) ) * epsinvx(i,j,k) * wavefct(1)
         Ey(i,j,k) = Ey(i,j,k) + w(2) * DT * ( w(6)/M4_SX(i,j,k) - w(4)/M4_SZ(i,j,k) ) * epsinvy(i,j,k) * wavefct(2)
})
M4_IFELSE_TE({
         Ez(i,j,k) = Ez(i,j,k) + w(3) * DT * ( w(4)/M4_SY(i,j,k) - w(5)/M4_SX(i,j,k) ) * epsinvz(i,j,k) * wavefct(3)
})
         })

      endif        

    })

  end subroutine StepEMatTfsf

!----------------------------------------------------------------------

  subroutine StepHMatTfsf(ncyc)

    integer :: ncyc

    M4_MODLOOP_DECL({MATTFSF},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))
    real(kind=8) :: ncyc1
    real(kind=8) :: wavefct(3), d(3), dd(3)
    integer :: di(3), n, l, o(3)

    M4_MODLOOP_EXPR({MATTFSF},mat,{

      if ( ncyc - mat%noffs .gt. mat%nend + mat%maxdelay ) cycle 

      M4_MODOBJ_GETREG(mat,reg)

      ! pre-calculate time signal for e-field modulation

      ncyc1 = 1.0*ncyc
       
      mat%wavefcte = GaussianWave(ncyc1, mat%noffs, mat%natt, mat%nsus, mat%ndcy, &
           mat%gamma, mat%omega0)

      ! store time signal for delayed h-field modulation

      if ( mat%planewave ) then

         mat%signalep = mod(mat%signalep+1,mat%maxdelay)
         mat%signale(mat%signalep) =  mat%wavefcte ! store signal function

      endif

      ! enter h-field modulation loop

      n = mat%signalhp+mat%maxdelay

      if ( .not. mat%planewave ) then

         M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

M4_IFELSE_TE({
         Hx(i+o(1),j,k) = Hx(i+o(1),j,k) + DT * ( w(3)/M4_SY(i,j,k) - w(2)/M4_SZ(i,j,k) ) * M4_MUINVX(i+o(1),j,k) * mat%wavefcth
         Hy(i,j+o(2),k) = Hy(i,j+o(2),k) + DT * ( w(1)/M4_SZ(i,j,k) - w(3)/M4_SX(i,j,k) ) * M4_MUINVY(i,j+o(2),k) * mat%wavefcth
})
M4_IFELSE_TM({
         Hz(i,j,k+o(3)) = Hz(i,j,k+o(3)) + DT * ( w(2)/M4_SX(i,j,k) - w(1)/M4_SY(i,j,k) ) * M4_MUINVZ(i,j,k+o(3)) * mat%wavefcth
})
   
         })
  
      else ! plane wave!

         M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

         do l = 1,3 
          
            d(l) = mat%delay(p) + mat%cdelay(l+3) ! time delay from origin
            di(l) = int(d(l))
            dd(l) = d(l) - di(l)
            ! use linear interpolation to read time signal
            wavefct(l) = (1. - dd(l)) * mat%signalh(mod(n-di(l),mat%maxdelay)) + &
                 dd(l) * mat%signalh(mod(n-di(l)-1,mat%maxdelay))
            
            w(l) = w(l) * mat%finc(l)
            o(l) = mat%orient(p,l)
            w(3+l) = abs(w(3+l))

            if ( p .gt. reg%numnodes ) write(6,*) p
            if ( w(l) .gt. 1. .or. w(l) .lt. -1 ) write(6,*) p, w(l)
            
         end do

M4_IFELSE_TE({
         Hx(i+o(1),j,k) = Hx(i+o(1),j,k) + w(4) * DT * ( w(3)/M4_SY(i,j,k) - w(2)/M4_SZ(i,j,k) ) * M4_MUINVX(i+o(1),j,k) * wavefct(1)
         Hy(i,j+o(2),k) = Hy(i,j+o(2),k) + w(5) * DT * ( w(1)/M4_SZ(i,j,k) - w(3)/M4_SX(i,j,k) ) * M4_MUINVY(i,j+o(2),k) * wavefct(2)
})
M4_IFELSE_TM({
         Hz(i,j,k+o(3)) = Hz(i,j,k+o(3)) + w(6) * DT * ( w(2)/M4_SX(i,j,k) - w(1)/M4_SY(i,j,k) ) * M4_MUINVZ(i,j,k+o(3)) * wavefct(3)
})
   
         })

      end if
  
  })
    
    
  end subroutine StepHMatTfsf

!----------------------------------------------------------------------

  real(kind=8) function SumJEMatTfsf(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc
    real(kind=8) :: Jx, Jy, Jz
   
    M4_MODLOOP_DECL({MATTFSF},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))

    sum = 0

    M4_MODLOOP_EXPR({MATTFSF},mat,{

    ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       if ( .not. mat%planewave ) then

          M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
          ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

          ! J(*,m) is P(n+1) and J(*,n) is P(n)      

          if ( mask(i,j,k) ) then
             
M4_IFELSE_TM({
            Jx = - ( w(5)/M4_SZ(i,j,k) - w(6)/M4_SY(i,j,k) ) * mat%wavefcte
            Jy = - ( w(6)/M4_SX(i,j,k) - w(4)/M4_SZ(i,j,k) ) * mat%wavefcte
})
M4_IFELSE_TE({
            Jz = - ( w(4)/M4_SY(i,j,k) - w(5)/M4_SX(i,j,k) ) * mat%wavefcte
})

            sum = sum + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * real(Ex(i,j,k)) * Jx + },{0. +}) &
M4_IFELSE_TM({ M4_VOLEY(i,j,k) * real(Ey(i,j,k)) * Jy + },{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * real(Ez(i,j,k)) * Jz   },{0.  }) &
               )
             
          endif

         })      

      else ! planewave case!







      end if

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

     if ( .not. mat%planewave ) then

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

       if ( mask(i,j,k) ) then

M4_IFELSE_TE({
          Kx = - ( w(3)/M4_SY(i,j,k) - w(2)/M4_SZ(i,j,k) ) * mat%wavefcth
          Ky = - ( w(1)/M4_SZ(i,j,k) - w(3)/M4_SX(i,j,k) ) * mat%wavefcth
})
M4_IFELSE_TM({
          Kz = - ( w(2)/M4_SX(i,j,k) - w(1)/M4_SY(i,j,k) ) * mat%wavefcth
})
          sum = sum + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * real(Hx(i,j,k)) * Kx +},{0. +}) &
M4_IFELSE_TM({ M4_VOLEY(i,j,k) * real(Hy(i,j,k)) * Ky +},{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * real(Hz(i,j,k)) * Kz  },{0.  }) &
               )
             
       endif

       })      


    else ! plane wave!



    endif

    })
    
    SumKHMatTfsf = sum

  end function SumKHMatTfsf


 !----------------------------------------------------------------------

  subroutine DisplayMatTfsfObj(mat)

    type(T_MATTFSF) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
        " OMEGA=",TRIM(f2str(mat%omega0,4)),&
    	" OFFS=",TRIM(i2str(int(mat%noffs))),&
    	" ATT=",TRIM(i2str(int(mat%natt))),&
    	" SUS=",TRIM(i2str(int(mat%nsus))),&
    	" DCY=",TRIM(i2str(int(mat%ndcy))) })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatTfsfObj
  
 !----------------------------------------------------------------------

   subroutine EchoMatTfsfObj(mat)

    type(T_MATTFSF) :: mat

    M4_WRITE_INFO({"--- mat # ",TRIM(i2str(mat%idx))})
    M4_WRITE_INFO({"lambdainv0 = ",mat%lambdainv0   })
    M4_WRITE_INFO({"omega0  = ",mat%omega0 })
    M4_WRITE_INFO({"dn = ",mat%dn })
    M4_WRITE_INFO({"gamma = ",mat%gamma})
    M4_WRITE_INFO({"noffs/natt/nsus/ndcy  = ",mat%noffs,mat%natt,mat%nsus,mat%ndcy   })
    
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


