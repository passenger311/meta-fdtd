!-*- F90 -*------------------------------------------------------------
!
!  module: matfourlvl / meta
!
!  Effective Maxwell Bloch material module for isotropic 4 level system.
!
!  subs:
!
!    InitializeMatFourlvl
!    FinalizeMatFourlvl
!    ReadMatFourlvlObj
!    StepEMatFourlvl
!    StepHMatFourlvl
!    SumJEKHMatFourlvl
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatFourlvl module calculates the reponse of an isotropic effective 4 level bloch system
!
! d/dt d/dt Pa_i + 2 * gammala * d/dt Pa_i + omegala**2 Pa_i = 2 * omegara * conv1 * Ma^2 * E_i (N1 - N2)
! d/dt d/dt Pb_i + 2 * gammalb * d/dt Pb_i + omegalb**2 Pb_i = 2 * omegarb * conv1 * Mb^2 * E_i (N0 - N3)
! d/dt N3 = conv2 / omegarb * sum ( ( d/dt Pb_i + gammal * Pb_i ) * E_i ) - gammanr32 * N3
! d/dt N2 = conv2 / omegara * sum ( ( d/dt Pa_i + gammal * Pa_i ) * E_i ) - gammanr21 * N2 + gammanr32 * N3  
! d/dt N1 = - conv2 / omegara * sum( ( d/dt Pa_i + gammal * Pa_i ) * E_i ) - gammanr10 * N1 + gammanr21 * N2
! d/dt N0 = - conv2 / omegarb * sum( ( d/dt Pb_i + gammal * Pb_i ) * E_i ) + gammanr10 * N1
! d/dt E = d/dt E* - epsinv * N * ( Ma d/dt Pa + Mb d/dt Pb )  
!
! where E* is the electric field as calculated without the sources.  
!
! In natural HL units we choose hbar = c = 1, that is the elementary charge
! would be
!
! e = sqrt ( 4 PI alpha ) with alpha = 1/137.0360
! 
!
! StepHMatFourlvl: update eq. Pa_i(n+1) = c1a * Pa_i(n) + c2a * Pa_i(n-1) + c3a * E_i(n) * ( N1(n) - N2(n) )
!                  update eq. Pb_i(n+1) = c1b * Pb_i(n) + c2b * Pb_i(n-1) + c3b * E_i(n) * ( N0(n) - N3(n) )
!                and update N (see below)
! StepEMatFourlvl: update eq. E_i(n+1) = E_i(n+1)* - epsinv * (Pa_i(n+1) - Pa_i(n) + Pb_i(n+1) - Pb_i(n))
!
! c1a = ( 2. - omegala**2 * DT**2 ) / ( 1 + gammala*dt )
! c2a = ( gammala*dt - 1 ) / ( 1 + gammala*dt )
! c3a = 2 * omegara * Ma / hbar / ( 1/dt^2 + gammala/dt )
! 
! N update eq:
! N3(n+1) = c43 * N3(n) + c5b * ( Pb(n+1) - Pb(n) )*( E(n+1) + E(n) ) + c6b * ( Pb(n+1) + Pb(n) )*( E(n+1) + E(n) )
! N2(n+1) = c42 * N2(n) + c5a * ( Pa(n+1) - Pa(n) )*( E(n+1) + E(n) ) + c6a * ( Pa(n+1) + Pa(n) )*( E(n+1) + E(n) )
! N1(n+1) = c41 * N1(n) - c5a * ( Pa(n+1) - Pa(n) )*( E(n+1) + E(n) ) - c6a * ( Pa(n+1) + Pa(n) )*( E(n+1) + E(n) )
! N0(n+1) = c40 * N0(n) - c5b * ( Pb(n+1) - Pb(n) )*( E(n+1) + E(n) ) - c6b * ( Pb(n+1) + Pb(n) )*( E(n+1) + E(n) )
! is calculated in StepHMatFourlvl in two steps. The first bit depends on E(n+1) and must be 
! calculated *after* StepEMat updated E to the proper E(n+1). However,
! n+1 -> n, so that
!
! N3(n) = N3(n)_p1 + c5b * ( Pb(n) - Pb(n-1) )*( E(n) ) + c6b * ( Pb(n) + Pb(n-1) )*( E(n) )
! N2(n) = N2(n)_p1 + c5a * ( Pa(n) - Pa(n-1) )*( E(n) ) + c6a * ( Pa(n) + Pa(n-1) )*( E(n) )
! N1(n) = N1(n)_p1 - c5a * ( Pa(n) - Pa(n-1) )*( E(n) ) - c6a * ( Pa(n) + Pa(n-1) )*( E(n) )
! N0(n) = N0(n)_p1 - c5b * ( Pb(n) - Pb(n-1) )*( E(n) ) - c6b * ( Pb(n) + Pb(n-1) )*( E(n) )
!
! The second bit (after calculating P) is actually the first part of the N calculation:
!
! N3(n+1)_p1 = c43 * N3(n) + c5b * ( Pb(n+1) - Pb(n) )*( E(n) ) + c6b * ( Pb(n+1) + Pb(n) )*( E(n) ) 
! N2(n+1)_p1 = c42 * N2(n) + N3(n)*(1 - c43) + c5a * ( Pa(n+1) - Pa(n) )*( E(n) ) + c6a * ( Pa(n+1) + Pb(n) )*( E(n) )
! N1(n+1)_p1 = c41 * N1(n) + N2(n)*(1 - c42) - c5a * ( Pa(n+1) - Pa(n) )*( E(n) ) + c6a * ( Pa(n+1) + Pa(n) )*( E(n) )
! N0(n+1)_p1 = N0(n) + N1(n)*(1 - c41) - c5b * ( Pb(n+1) - Pb(n) )*( E(n) ) + c6b * ( Pb(n+1) + Pb(n) )*( E(n) )
!
! In all Nj update eqs. P * E denotes the scalar product of vectors P and E
!
! c43 = ( 2. - gamma32 * dt ) / ( 2. + gamma32 * DT )
! c5a = 1. / ( 2. + gamma21 * DT ) * 1. / ( omegara ) * conv2
! c6a = c5a  * DT * mat%gammala / 2.
!


module matfourlvl

  use constant
  use parse
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MATHEAD_DECL({MATFOURLVL},MAXMATOBJ,{

  ! input parameters
  real(kind=8) :: lambdarinva    ! inv. vac. plasma wavelength
  real(kind=8) :: gammala        ! resonance width
  real(kind=8) :: lambdarinvb
  real(kind=8) :: gammalb

  real(kind=8) :: Ma             ! dipole strength of transition 1 <-> 2 [1/dt]
  real(kind=8) :: Mb             ! dipole strength of transition 0 <-> 3 [1/dt]
  real(kind=8) :: Nstart(3)      ! start density
  real(kind=8) :: N              ! number of 4lvl systems per grid cell
  real(kind=8) :: gamma32, gamma21, gamma10 ! nonradiative transition rates
  integer :: napprox, cyc
  
  ! calculated

  real(kind=8) :: omegara        ! resonance frequency 1 <-> 2
  real(kind=8) :: omegala        ! lorentz frequency 1 <-> 2
  real(kind=8) :: omegarb        ! resonance frequency 0 <-> 3
  real(kind=8) :: omegalb        ! lorentz frequency 0 <-> 3

  ! coefficients
  real(kind=8) :: c1a, c2a, c3a, c5a, c6a, c1b, c2b, c3b, c5b, c6b
  real(kind=8) :: c43, c42, c41
  
  ! polarisation field 
  M4_FTYPE, dimension(:,:), pointer :: Pax, Pay, Paz
  M4_FTYPE, dimension(:,:), pointer :: Pbx, Pby, Pbz

  ! number density
  real(kind=8), dimension(:), pointer :: N0
  real(kind=8), dimension(:), pointer :: N1
  real(kind=8), dimension(:), pointer :: N2
  real(kind=8), dimension(:), pointer :: N3

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatFourlvlObj(funit,lcount)

    M4_MODREAD_DECL({MATFOURLVL}, funit,lcount,mat,reg,out)
    real(kind=8) :: v(2)
    real(kind=8) :: c(3)
    logical :: eof,err
    character(len=LINELNG) :: line

    M4_WRITE_DBG(". enter ReadMatFourlvlObj")

    M4_MODREAD_EXPR({MATFOURLVL},funit,lcount,mat,reg,3,out,{ 

    ! read mat parameters here, as defined in mat data structure
    call readfloats(funit,lcount, v, 2)
    mat%lambdarinva = v(1)
    mat%lambdarinvb = v(2)

    call readfloats(funit,lcount, v, 2)
    mat%gammala = v(1)
    mat%gammalb = v(2)

    call readfloat(funit,lcount, mat%Ma)
    call readfloat(funit,lcount, mat%Mb)

    call readfloat(funit,lcount,mat%N)

    call readfloats(funit,lcount,c,3)
    mat%Nstart(1) = c(1)
    mat%Nstart(2) = c(2)
    mat%Nstart(3) = c(3)

    call readfloats(funit,lcount,c,3) 
    mat%gamma32 = c(1)
    mat%gamma21 = c(2)
    mat%gamma10 = c(3)

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatFourlvlObj")

  end subroutine ReadMatFourlvlObj

!----------------------------------------------------------------------

  subroutine InitializeMatFourlvl

    type (T_REG) :: reg
    integer :: err
    real (kind=8) :: conv1, conv2
    M4_MODLOOP_DECL({MATFOURLVL},mat) 
    M4_WRITE_DBG(". enter InitializeMatFourlvl")
    M4_MODLOOP_EXPR({MATFOURLVL},mat,{
    
       ! initialize mat object here

       mat%omegara = 2. * PI * mat%lambdarinva
       mat%omegarb = 2. * PI * mat%lambdarinvb
       mat%omegala = sqrt( mat%omegara**2 + mat%gammala**2 )
       mat%omegalb = sqrt( mat%omegarb**2 + mat%gammalb**2 )

       reg = regobj(mat%regidx)

       allocate(mat%Pax(reg%numnodes,2), mat%Pay(reg%numnodes,2), mat%Paz(reg%numnodes,2), stat=err)
       allocate(mat%Pbx(reg%numnodes,2), mat%Pby(reg%numnodes,2), mat%Pbz(reg%numnodes,2), stat=err)
       allocate(mat%N0(reg%numnodes), mat%N1(reg%numnodes), mat%N2(reg%numnodes), mat%N3(reg%numnodes), stat=err)
       M4_ALLOC_ERROR(err,"InitializeMatBloch")

       mat%Pax = 0.
       mat%Pay = 0.
       mat%Paz = 0.
       mat%Pbx = 0.
       mat%Pby = 0.
       mat%Pbz = 0.
       mat%N0 = 1 - mat%Nstart(1) - mat%Nstart(2) - mat%Nstart(3)
       if (mat%N0(1) < 0) then
          write(6,*) "Error in configuration: N0 is smaller zero" 
          !exit(1)
       endif
       mat%N1 = mat%Nstart(1)
       mat%N2 = mat%Nstart(2)
       mat%N3 = mat%Nstart(3)

! From unit conversion of polarisation equation - SI to Computational units 
       conv1 = SI_4PIALPHA
! From unit conversion of population equation - SI to Computation units
       conv2 = ( ( REAL_DX * SI_C ) / ( SI_HBAR ) )

! for polarisation integration
       mat%c1a = ( 2. - mat%omegala**2 * DT**2 ) / ( 1 + mat%gammala * DT )
       mat%c2a = ( mat%gammala * DT - 1. ) / ( 1. + mat%gammala * DT )
       mat%c3a = 2. * mat%omegara * conv1 * mat%Ma**2 / ( 1/DT**2 + mat%gammala * DT )
       mat%c1b = ( 2. - mat%omegalb**2 * DT**2 ) / ( 1 + mat%gammalb * DT )
       mat%c2b = ( mat%gammalb * DT - 1 ) / ( 1 + mat%gammalb * DT )
       mat%c3b = 2. * mat%omegarb * conv1 * mat%Mb**2 / ( 1/DT**2 + mat%gammalb * DT )

! for density integration
       mat%c43 = ( 2. - mat%gamma32 * DT ) / ( 2. + mat%gamma32 * DT )
       mat%c42 = ( 2. - mat%gamma21 * DT ) / ( 2. + mat%gamma21 * DT )
       mat%c41 = ( 2. - mat%gamma10 * DT ) / ( 2. + mat%gamma10 * DT )
       mat%c5a = 1. / ( 2. + ( mat%gamma21 + mat%gamma10 ) * DT ) * 1. / ( mat%omegara ) * conv2 
       mat%c6a = mat%c5a  * DT * mat%gammala / 2.
       mat%c5b = 1. / ( 2. + ( mat%gamma32 + mat%gamma10 ) * DT ) * 1. / ( mat%omegarb ) * conv2
       mat%c6b = mat%c5b * DT * mat%gammalb / 2.

       mat%cyc = 1

       M4_IFELSE_DBG({call EchoMatFourlvlObj(mat)},{call DisplayMatFourlvlObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatFourlvl")

  end subroutine InitializeMatFourlvl

!----------------------------------------------------------------------

  subroutine FinalizeMatFourlvl

    M4_MODLOOP_DECL({MATFOURLVL},mat)
    M4_WRITE_DBG(". enter FinalizeMatFourlvl")
    M4_MODLOOP_EXPR({MATFOURLVL},mat,{

    ! finalize mat object here
    deallocate(mat%Pax, mat%Pay, mat%Paz, mat%Pbx, mat%Pby, mat%Pbz, mat%N0, mat%N1, mat%N2, mat%N3)

    })
    M4_WRITE_DBG(". exit FinalizeMatFourlvl")

  end subroutine FinalizeMatFourlvl

!----------------------------------------------------------------------

  subroutine StepHMatFourlvl(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATFOURLVL},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8) :: pema, pena, pemb, penb
    real(kind=8) :: lEx, lEy, lEz, n32, n21, n10, ninva, ninvb
    M4_MODLOOP_EXPR({MATFOURLVL},mat,{

    ! this loops over all mat structures, setting mat

      M4_MODOBJ_GETREG(mat,reg)

      n = mod(ncyc-1+2,2) + 1
      m = mod(ncyc+2,2) + 1

      mat%cyc = m
      
        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

        ! calculate local field enhancement due to dielectric
        lEx = Ex(i,j,k) * ( 2. + 1./epsinvx(i,j,k) ) / 3.
        lEy = Ey(i,j,k) * ( 2. + 1./epsinvy(i,j,k) ) / 3.
        lEz = Ez(i,j,k) * ( 2. + 1./epsinvz(i,j,k) ) / 3.

        ! calculate second part of the density response (after the new E field got calculated)
        
        pema = mat%Pax(p,m) * lEx + mat%Pay(p,m) * lEy + mat%Paz(p,m) * lEz
        pena = mat%Pax(p,n) * lEx + mat%Pay(p,n) * lEy + mat%Paz(p,n) * lEz
        pemb = mat%Pbx(p,m) * lEx + mat%Pby(p,m) * lEy + mat%Pbz(p,m) * lEz
        penb = mat%Pbx(p,m) * lEx + mat%Pby(p,n) * lEy + mat%Pbz(p,n) * lEz
      
        mat%N3(p) = mat%N3(p) + mat%c5b * ( penb - pemb ) + mat%c6b * ( penb + pemb )
        mat%N2(p) = mat%N2(p) + mat%c5a * ( pena - pema ) + mat%c6a * ( pena + pema )
        mat%N1(p) = mat%N1(p) - mat%c5a * ( pena - pema ) - mat%c6a * ( pena + pema )
        mat%N0(p) = mat%N0(p) - mat%c5b * ( penb - pemb ) - mat%c6b * ( penb + pemb )
        ! calculate P(n+1) from P(n),P(n-1),E(n) and N(n)
        
        ! before: P(*,m) is P(n-1), P(*,n) is P(n)
        ninva = mat%N1(p) - mat%N2(p)
        ninvb = mat%N0(p) - mat%N3(p)

        mat%Pax(p,m) = mat%c1a * mat%Pax(p,n) + mat%c2a * mat%Pax(p,m) + mat%c3a * lEx * ninva
        mat%Pay(p,m) = mat%c1a * mat%Pay(p,n) + mat%c2a * mat%Pay(p,m) + mat%c3a * lEy * ninva
        mat%Paz(p,m) = mat%c1a * mat%Paz(p,n) + mat%c2a * mat%Paz(p,m) + mat%c3a * lEz * ninva
        mat%Pbx(p,m) = mat%c1b * mat%Pbx(p,n) + mat%c2b * mat%Pbx(p,m) + mat%c3b * lEx * ninvb
        mat%Pby(p,m) = mat%c1b * mat%Pby(p,n) + mat%c2b * mat%Pby(p,m) + mat%c3b * lEy * ninvb
        mat%Pbz(p,m) = mat%c1b * mat%Pbz(p,n) + mat%c2b * mat%Pbz(p,m) + mat%c3b * lEz * ninvb

        ! calculate first part of the density response

        pema = mat%Pax(p,m) * lEx + mat%Pay(p,m) * lEy + mat%Paz(p,m) * lEz
        pemb = mat%Pbx(p,m) * lEx + mat%Pby(p,m) * lEy + mat%Pbz(p,m) * lEz

        n32 = (1 - mat%c43) * mat%N3(p)
        n21 = (1 - mat%c42) * mat%N2(p)
        n10 = (1 - mat%c41) * mat%N1(p)

        mat%N3(p) = mat%N3(p) - n32 + mat%c5b * ( pemb - penb ) + mat%c6b * ( pemb + penb )
        mat%N2(p) = mat%N2(p) - n21 + n32 + mat%c5a * ( pema - pena ) + mat%c6a * ( pema + pena )
        mat%N1(p) = mat%N1(p) - n10 + n21 - mat%c5a * ( pema - pena ) - mat%c6a * ( pema + pena )
        mat%N0(p) = mat%N0(p) + n10 - mat%c5b * ( pemb - penb ) - mat%c6b * ( pemb + penb )

        ! after: J(*,m) is now P(n+1)
        ! m and n will be flipped in the next timestep!

        })      


    })
  
  end subroutine StepHMatFourlvl


!----------------------------------------------------------------------

  subroutine StepEMatFourlvl(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATFOURLVL},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATFOURLVL},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

M4_IFELSE_TM({
       Ex(i,j,k) = Ex(i,j,k) - w(1) * epsinvx(i,j,k) * mat%N * (  ( mat%Pax(p,m) - mat%Pax(p,n) )  + &
           ( mat%Pbx(p,m) - mat%Pbx(p,n) ) )
       Ey(i,j,k) = Ey(i,j,k) - w(2) * epsinvy(i,j,k) * mat%N * (  ( mat%Pay(p,m) - mat%Pay(p,n) ) + &
           ( mat%Pby(p,m) - mat%Pby(p,n) ) )
})
M4_IFELSE_TE({
       Ez(i,j,k) = Ez(i,j,k) - w(3) * epsinvz(i,j,k) * mat%N * (  ( mat%Paz(p,m) - mat%Paz(p,n) )  + &
           ( mat%Pbz(p,m) - mat%Pbz(p,n) ) )
})
       })      

    })

  end subroutine StepEMatFourlvl

!----------------------------------------------------------------------

  real(kind=8) function SumJEMatFourlvl(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
   
    M4_MODLOOP_DECL({MATFOURLVL},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    sum = 0

    M4_MODLOOP_EXPR({MATFOURLVL},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

       if ( mask(i,j,k) ) then

          sum = sum + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * w(1) * Ex(i,j,k) * mat%N * (  ( mat%Pax(p,m) - mat%Pax(p,n) ) + &
 ( mat%Pbx(p,m) - mat%Pbx(p,n) ) ) / DT +}, {0. +}) &
M4_IFELSE_TM({ M4_VOLEY(i,j,k) * w(2) * Ey(i,j,k) * mat%N * (  ( mat%Pay(p,m) - mat%Pay(p,n) ) + &
 ( mat%Pby(p,m) - mat%Pby(p,n) ) ) / DT +}, {0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * w(3) * Ez(i,j,k) * mat%N * (  ( mat%Paz(p,m) - mat%Paz(p,n) ) + &
 ( mat%Pbz(p,m) - mat%Pbz(p,n) ) ) / DT  }, {0.  }) &
               )
       endif

       })      

    })
    
    SumJEMatFourlvl = sum
    
  end function SumJEMatFourlvl

!----------------------------------------------------------------------

  real(kind=8) function SumKHMatFourlvl(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc

    SumKHMatFourlvl = 0.

  end function SumKHMatFourlvl
 
!----------------------------------------------------------------------

  subroutine DisplayMatFourlvlObj(mat)

    type(T_MATFOURLVL) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" lambdarinva=",TRIM(f2str(mat%lambdarinva,5)),&
    	" gammala=",TRIM(f2str(mat%gammala,5)),&
        " lamdarinvb=",TRIM(f2str(mat%lambdarinvb,5)),&
        " gammalb=",TRIM(f2str(mat%gammalb,5))
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatFourlvlObj

!----------------------------------------------------------------------

   subroutine EchoMatFourlvlObj(mat)

    type(T_MATFOURLVL) :: mat

    M4_WRITE_INFO({"--- matfourlvl # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"lambdarinva = ",mat%lambdarinva })
    M4_WRITE_INFO({"omegara = ",mat%omegara })
    M4_WRITE_INFO({"gammala = ",mat%gammala })
    M4_WRITE_INFO({"lambdarinvb = ",mat%lambdarinvb })
    M4_WRITE_INFO({"omegarb = ",mat%omegarb })
    M4_WRITE_INFO({"gammalb = ",mat%gammalb })
    M4_WRITE_INFO({"c1a = ",mat%c1a })
    M4_WRITE_INFO({"c2a = ",mat%c2a })
    M4_WRITE_INFO({"c3a = ",mat%c3a })
    M4_WRITE_INFO({"c1b = ",mat%c1b })
    M4_WRITE_INFO({"c2b = ",mat%c2b })
    M4_WRITE_INFO({"c3b = ",mat%c3b })
    M4_WRITE_INFO({"c5a = ",mat%c5a })
    M4_WRITE_INFO({"c5b = ",mat%c5b })
    M4_WRITE_INFO({"c6a = ",mat%c6a })
    M4_WRITE_INFO({"c6b = ",mat%c6b })
    M4_WRITE_INFO({"gamma32 = ",mat%gamma32 })
    M4_WRITE_INFO({"gamma21 = ",mat%gamma21 })
    M4_WRITE_INFO({"gamma10 = ",mat%gamma10 })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatFourlvlObj

  
!----------------------------------------------------------------------

end module matfourlvl

! Authors:  A.Pusch, J.Hamm 
! Modified: 11/01/2010
!
! =====================================================================


