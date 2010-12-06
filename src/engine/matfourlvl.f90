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
! d/dt N3 = conv2 / omegarb * sum ( ( d/dt Pb_i + gammal * Pb_i ) * E_i ) - gammanr32 * N3 - gammanr30 * N3
! d/dt N2 = conv2 / omegara * sum ( ( d/dt Pa_i + gammal * Pa_i ) * E_i ) - gammanr21 * N2 + gammanr32 * N3  
! d/dt N1 = - conv2 / omegara * sum( ( d/dt Pa_i + gammal * Pa_i ) * E_i ) - gammanr10 * N1 + gammanr21 * N2
! d/dt N0 = - conv2 / omegarb * sum( ( d/dt Pb_i + gammal * Pb_i ) * E_i ) + gammanr10 * N1 + gammanr30 * N3
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
! N update eq: ( 4 x 4 matrix equation  A * N(n+1) = B * N(n) + C )
!
! A11 * N3(n+1) = B11 * N3(n) + x1
! A21 * N3(n+1) + A22 * N2(n+1) = B21 * N3(n) + B22 * N2(n) + x2
! A31 *  N3(n+1) + A32 * N2(n+1) + A33 * N1(n+1) =
!                     B31 * N3(n) + B32 * N2(n) + B33 * N1(n) - x2
! A41 * N3(n+1) + A42 * N2(n+1) + A43 * N1(n+1) + A44 * N0(n+1) = 
!                     B41 * N3(n) + B42 * N2(n) + B43 * N1(n) + B44 * N0(n) - x1
! 
! With Matrices A and B and polarisation couplings x1 and x2
! These equations have to be transformed to a diagonal form in N(n+1) by multiplying with the inverse of A.
! The resulting equations are calculated in StepHMatFourlvl in two steps. The first bit depends on E(n+1) and must be 
! calculated *after* StepEMat updated E to the proper E(n+1). However,
! n+1 -> n, so that
! 
! AC = inv(A)*C with C = {{x1,0,0,0},{0,x2,0,0},{0,0,-x2,0},{0,0,0,-x1}}
! x1 = ( 1. / 2. / omegarb * conv2 / DT * ( Pb(n) - Pb(n-1) ) + gammalb / 4 . / omegarb * conv2 * ( Pb(n) - Pb(n-1) ) ) * 
!      * E(n)
! x2 = ( 1. / 2. / omegara * conv2 / DT * ( Pa(n) - Pa(n-1) ) + gammala / 4 . / omegara * conv2 * ( Pa(n) - Pa(n-1) ) ) * 
!      * E(n)
! N3(n) = N3(n)_p1 + AC11 * x1
! N2(n) = N2(n)_p1 + AC21 * x1 + AC22 * x2
! N1(n) = N1(n)_p1 + AC31 * x1 + AC32 * x2 + AC33 * x2
! N0(n) = N0(n)_p1 + AC41 * x1 + AC42 * x2 + AC43 * x2 + AC44 * x1
!
! The second bit (after calculating P) is actually the first part of the N calculation:
! 
! AB = inv(A)*B
!
! x1 = E(n) * (1. / 2. / omegarb * conv2 / DT * ( Pb(n+1) - Pb(n) ) + gammalb / 4 . / omegalb * conv2 * ( Pb(n+1) - Pb(n) ))
! x2 = E(n) * (1. / 2. / omegara * conv2 / DT * ( Pa(n+1) - Pa(n) ) + gammala / 4 . / omegala * conv2 * ( Pa(n+1) - Pa(n) ))
! N3(n+1)_p1 = AB11 * N3(n) + AC11 * x1
! N2(n+1)_p1 = AB21 * N3(n) + AB22 * N2(n) + AC21 * x1 + AC22 * x2
! N1(n+1)_p1 = AB31 * N3(n) + AB32 * N2(n) + AB33 * N1(n) + AC31 * x1 + AC32 * x2 + AC33 * x2 
! N0(n+1)_p1 = AB41 * N3(n) + AB42 * N2(n) + AB43 * N1(n) + AB44 * N0(n) + AC41 * x1 + AC42 * x2 + AC43 * x2 + AC44 * x1
!
! In all Nj update eqs. P * E denotes the scalar product of vectors P and E
!



module matfourlvl

  use constant
  use checkpoint
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
  real(kind=8) :: gamma32, gamma21, gamma10, gamma30 ! nonradiative transition rates
  integer :: napprox, cyc
  
  ! calculated

  real(kind=8) :: omegara        ! resonance frequency 1 <-> 2
  real(kind=8) :: omegala        ! lorentz frequency 1 <-> 2
  real(kind=8) :: omegarb        ! resonance frequency 0 <-> 3
  real(kind=8) :: omegalb        ! lorentz frequency 0 <-> 3

  ! coefficients
  real(kind=8) :: c1a, c2a, c3a, c1b, c2b, c3b
  real(kind=8) :: ab11, ab21, ab22, ab31, ab32, ab33, ab41, ab42, ab43, ab44
  real(kind=8) :: ac11, ac21, ac22, ac31, ac32, ac33, ac41, ac42, ac43, ac44
  real(kind=8) :: x1fac1, x1fac2, x2fac1, x2fac2
  
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
    real(kind=8) :: d(4)
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

    call readfloats(funit,lcount,d,4)
    mat%gamma30 = d(1) 
    mat%gamma32 = d(2)
    mat%gamma21 = d(3)
    mat%gamma10 = d(4)

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatFourlvlObj")

  end subroutine ReadMatFourlvlObj

!----------------------------------------------------------------------

  subroutine InitializeMatFourlvl

    integer :: err
    real (kind=8) :: conv1, conv2, g32, g30, g21, g10

    M4_MODLOOP_DECL({MATFOURLVL},mat) 
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

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
       mat%x1fac1 = 1. / 2. / mat%omegarb / DT * conv2
       mat%x1fac2 = mat%x1fac1 * DT * mat%gammalb / 2.
       mat%x2fac1 = 1. / 2. / mat%omegara / DT * conv2
       mat%x2fac2 = mat%x2fac2 * DT * mat%gammala / 2.
       g10 = mat%gamma10 * DT
       g21 = mat%gamma21 * DT
       g32 = mat%gamma32 * DT
       g30 = mat%gamma30 * DT
! relaxation terms
       mat%ab11 = - ( - 2. + g32 + g30 )/( 2. + g32 + g30)
       mat%ab21 = 4. * g32 / ( 2. + g32 + g30 ) / ( 2. + g21 )
       mat%ab31 = 4. * g32 * g21 / ( 2. + g10 ) / ( 2. + g21 ) / ( 2. + g32 + g30 )
       mat%ab41 = 2. * ( 4. * g30 + 2. * g30 * g21 + 2. * g30 * g10 + g30 * g10 * g21 + g10 * g21 * g32 ) &
                  / ( 2. + g10 ) / ( 2. + g21 ) / ( 2. + g32 + g30 )
       mat%ab22 = - ( - 2. + g21 ) / ( 2 + g21 )
       mat%ab32 = 4. * g21 / ( 2. + g10 ) / ( 2. + g21 )
       mat%ab42 = 2. * g10 * g21 / ( 2. + g10 ) / ( 2. + g21 )
       mat%ab33 = - ( - 2. + g10 ) / ( 2. + g10 )
       mat%ab43 = 2. * g10 / ( 2. + g10 )
       mat%ab44 = 1.
! stimulated emission and absorption
       mat%ac11 = 2. * DT / ( 2. + g32 + g30 )
       mat%ac21 = 2. * DT * g32 / ( 2. + g32 + g30 ) / ( 2. + g21 )
       mat%ac31 = 2. * DT * g21 * g32 / ( 2. + g21 ) / ( 2. + g10 ) / ( 2. + g32 + g30 )
       mat%ac41 = DT * ( 4. * g30 + 2. * g30 * g21 + 2. * g30 * g10 + g30 * g10 * g21 + g10 * g21 * g32 ) &
                   / ( 2. + g10 ) / ( 2. + g21 ) / ( 2. + g32 + g30 )
       mat%ac22 = 2. * DT / ( 2. + g21 )
       mat%ac32 = 2. * DT * g21 / ( 2. + g21 ) / ( 2. + g10 )
       mat%ac42 = g10 * DT * g21 / ( 2. + g10 ) / ( 2. + g21 )
       mat%ac33 = - 2. * DT / ( 2. + g10 )
       mat%ac43 = - g10 * DT / ( 2. + g10 )
       mat%ac44 = - DT

       mat%cyc = 1

! load from checkpoint file

       if ( load_state ) then

          M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

          read(UNITCHK) mat%Pax(p,1), mat%Pax(p,2)
          read(UNITCHK) mat%Pay(p,1), mat%Pay(p,2)
          read(UNITCHK) mat%Paz(p,1), mat%Paz(p,2)
          read(UNITCHK) mat%Pbx(p,1), mat%Pbx(p,2)
          read(UNITCHK) mat%Pby(p,1), mat%Pby(p,2)
          read(UNITCHK) mat%Pbz(p,1), mat%Pbz(p,2)
          read(UNITCHK) mat%Pbz(p,1), mat%Pbz(p,2)
          read(UNITCHK) mat%N0(p), mat%N1(p), mat%N2(p), mat%N3(p)

          })

       end if

       M4_IFELSE_DBG({call EchoMatFourlvlObj(mat)},{call DisplayMatFourlvlObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatFourlvl")

  end subroutine InitializeMatFourlvl

!----------------------------------------------------------------------

  subroutine FinalizeMatFourlvl

    M4_MODLOOP_DECL({MATFOURLVL},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_WRITE_DBG(". enter FinalizeMatFourlvl")
    M4_MODLOOP_EXPR({MATFOURLVL},mat,{

       M4_MODOBJ_GETREG(mat,reg)
! save to checkpoint file

       if ( save_state ) then

          M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

          write(UNITCHK) mat%Pax(p,1), mat%Pax(p,2)
          write(UNITCHK) mat%Pay(p,1), mat%Pay(p,2)
          write(UNITCHK) mat%Paz(p,1), mat%Paz(p,2)
          write(UNITCHK) mat%Pbx(p,1), mat%Pbx(p,2)
          write(UNITCHK) mat%Pby(p,1), mat%Pby(p,2)
          write(UNITCHK) mat%Pbz(p,1), mat%Pbz(p,2)
          write(UNITCHK) mat%Pbz(p,1), mat%Pbz(p,2)
          write(UNITCHK) mat%N0(p), mat%N1(p), mat%N2(p), mat%N3(p)

          })

       end if

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
    real(kind=8) :: lEx, lEy, lEz, n3, n2, n1, ninva, ninvb, x1 ,x2
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
        x1 = mat%x1fac1 * ( penb - pemb ) + mat%x1fac2 * ( penb + pemb ) 
        x2 = mat%x2fac1 * ( pena - pema ) + mat%x2fac2 * ( pena + pema )
        mat%N3(p) = mat%N3(p) + mat%ac11 * x1
        mat%N2(p) = mat%N2(p) + mat%ac21 * x1 + mat%ac22 * x2
        mat%N1(p) = mat%N1(p) + mat%ac31 * x1 + mat%ac32 * x2 + mat%ac33 * x2
        mat%N0(p) = mat%N0(p) + mat%ac41 * x1 + mat%ac42 * x2 + mat%ac43 * x2 + mat%ac44 * x1
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
        x1 = mat%x1fac1 * ( pemb - penb ) + mat%x1fac2 * ( penb + pemb ) 
        x2 = mat%x2fac1 * ( pema - pena ) + mat%x2fac2 * ( pena + pema )
        n3 = mat%N3(p)
        n2 = mat%N2(p)
        n1 = mat%N1(p)
        
        mat%N3(p) = mat%ab11 * n3 + mat%ac11 * x1
        mat%N2(p) = mat%ab21 * n3 + mat%ab22 * n2 + mat%ac21 * x1 + mat%ac22 * x2
        mat%N1(p) = mat%ab31 * n3 + mat%ab32 * n2 + mat%ab33 * n1 + mat%ac31 * x1 &
                    + mat%ac32 * x2 + mat%ac33 * x2
        mat%N0(p) = mat%ab41 * n3 + mat%ab42 * n2 + mat%ab43 * n1 + mat%ab44 * mat%N0(p) &
                    + mat%ac41 * x1 + mat%ac42 * x2 + mat%ac43 * x2 + mat%ac44 * x1

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
    M4_WRITE_INFO({"gamma30 = ",mat%gamma30 })
    M4_WRITE_INFO({"gamma32 = ",mat%gamma32 })
    M4_WRITE_INFO({"gamma21 = ",mat%gamma21 })
    M4_WRITE_INFO({"gamma10 = ",mat%gamma10 })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatFourlvlObj

  
!----------------------------------------------------------------------

end module matfourlvl

! Authors:  A.Pusch, J.Hamm 
! Modified: 22/03/2010
!
! =====================================================================


