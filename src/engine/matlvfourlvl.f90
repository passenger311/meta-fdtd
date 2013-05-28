!-*- F90 -*------------------------------------------------------------
!
!  module: matlvfourlvl / meta
!
!  Effective Maxwell Bloch material module for isotropic 4 level system with langevin noise terms.
!
!  subs:
!
!    InitializeMatLVFourlvl
!    FinalizeMatLVFourlvl
!    ReadMatLVFourlvlObj
!    StepEMatLVFourlvl
!    StepHMatLVFourlvl
!    SumJEKHMatLVFourlvl
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatLVFourlvl module calculates the reponse of an isotropic effective 4 level bloch system with langevin noise terms
!
! d/dt d/dt Pa_i + 2 * gammala * d/dt Pa_i + omegala**2 Pa_i = 2 * omegara * conv1 * Ma^2 * E_i (N1 - N2) - pfac * (omegara zeta2a_i sqrt(gammala N2) + d/dt (zeta1a_i sqrt(gammala N2)) ) 
! d/dt d/dt Pb_i + 2 * gammalb * d/dt Pb_i + omegalb**2 Pb_i = 2 * omegarb * conv1 * Mb^2 * E_i (N0 - N3) - pfac * (omegarb zeta2b_i sqrt(gammalb N2) + d/dt (zeta1b_i sqrt(gammalb N2)) ) 
! d/dt N3 = conv2 / omegarb * sum ( ( d/dt Pb_i + gammal * Pb_i ) * E_i ) - gammanr32 * N3 - gammanr30 * N3 + zeta32 sqrt(gammanr32/N) + zeta30 sqrt(gammanr30/N)
! d/dt N2 = conv2 / omegara * sum ( ( d/dt Pa_i + gammal * Pa_i ) * E_i ) - gammanr21 * N2 + gammanr32 * N3 - zeta32 sqrt(gammanr32/N) + zeta21 sqrt(gammanr21/N)  
! d/dt N1 = - conv2 / omegara * sum( ( d/dt Pa_i + gammal * Pa_i ) * E_i ) - gammanr10 * N1 + gammanr21 * N2 - zeta21 sqrt(gammanr21/N) + zeta10 sqrt(gammanr10/N)
! d/dt N0 = - conv2 / omegarb * sum( ( d/dt Pb_i + gammal * Pb_i ) * E_i ) + gammanr10 * N1 + gammanr30 * N3 - zeta30 sqrt(gammanr30/N) - zeta10 sqrt(gammanr10/N)
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
! StepHMatLVFourlvl: update eq. Pa_i(n+1) = c1a * Pa_i(n) + c2a * Pa_i(n-1) + c3a * E_i(n) * ( N1(n) - N2(n) ) +Noise
!                  update eq. Pb_i(n+1) = c1b * Pb_i(n) + c2b * Pb_i(n-1) + c3b * E_i(n) * ( N0(n) - N3(n) ) + Noise
!                and update N (see below)
! StepEMatLVFourlvl: update eq. E_i(n+1) = E_i(n+1)* - epsinv * (Pa_i(n+1) - Pa_i(n) + Pb_i(n+1) - Pb_i(n))
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
! The resulting equations are calculated in StepHMatLVFourlvl in two steps. The first bit depends on E(n+1) and must be 
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
! N3(n+1)_p1 = AB11 * N3(n) + AC11 * x1 + Noise
! N2(n+1)_p1 = AB21 * N3(n) + AB22 * N2(n) + AC21 * x1 + AC22 * x2 + Noise
! N1(n+1)_p1 = AB31 * N3(n) + AB32 * N2(n) + AB33 * N1(n) + AC31 * x1 + AC32 * x2 + AC33 * x2 + Noise
! N0(n+1)_p1 = AB41 * N3(n) + AB42 * N2(n) + AB43 * N1(n) + AB44 * N0(n) + AC41 * x1 + AC42 * x2 + AC43 * x2 + AC44 * x1 + Noise
!
! In all Nj update eqs. P * E denotes the scalar product of vectors P and E
!



module matlvfourlvl

  use constant
  use checkpoint
  use parse
  use reglist
  use outlist
  use grid
  use fdtd
  use ziggurat

  implicit none
  private
  save

  M4_MATHEAD_DECL({MATLVFOURLVL},MAXMATOBJ,{

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

  ! Lorenz-Lorentz field due to epsilon included?
  real(kind=8) :: epsLFE

  integer :: seed 		! random number seed
  real(kind=8) :: volfac	! factor between simulation volume dx^3 and real volume (==1 for 3D simulation)
  real(kind=8) :: linefac	! factor between homogeneous and inhomogeneous linewidth
  real(kind=8) :: lvPa1, lvPa2, lvPb1, lvPb2, lvN30, lvN32, lvN21, lvN10 ! factors for langevin noise magnitude
  real(kind=8), dimension(:), pointer :: lvtermax, lvtermay, lvtermaz, lvtermbx, lvtermby, lvtermbz ! old random numbers for d/dt P 

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatLVFourlvlObj(funit,lcount)

    M4_MODREAD_DECL({MATLVFOURLVL}, funit,lcount,mat,reg,out)
    real(kind=8) :: v(2)
    real(kind=8) :: c(3)
    real(kind=8) :: d(4)
    logical :: eof,err
    character(len=LINELNG) :: line

    M4_WRITE_DBG(". enter ReadMatLVFourlvlObj")

    M4_MODREAD_EXPR({MATLVFOURLVL},funit,lcount,mat,reg,3,out,{ 

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

    call readfloat(funit, lcount, mat%epsLFE)
    call readfloat(funit,lcount,mat%volfac)
    call readfloat(funit,lcount,mat%linefac)
    call readint(funit,lcount,mat%seed)

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatLVFourlvlObj")

  end subroutine ReadMatLVFourlvlObj

!----------------------------------------------------------------------

  subroutine InitializeMatLVFourlvl

    integer :: err
    real (kind=8) :: conv1, conv2, g32, g30, g21, g10, ninit, pfac
    real(kind=8) :: x1, x2, re, im, tipangle !only for Superfluorescence test

    M4_MODLOOP_DECL({MATLVFOURLVL},mat) 
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_WRITE_DBG(". enter InitializeMatLVFourlvl")
    M4_MODLOOP_EXPR({MATLVFOURLVL},mat,{
    
       ! initialize mat object here

       mat%omegara = 2. * PI * mat%lambdarinva
       mat%omegarb = 2. * PI * mat%lambdarinvb
       mat%omegala = sqrt( mat%omegara**2 + mat%gammala**2 )
       mat%omegalb = sqrt( mat%omegarb**2 + mat%gammalb**2 )

       reg = regobj(mat%regidx)

       allocate(mat%Pax(reg%numnodes,2), mat%Pay(reg%numnodes,2), mat%Paz(reg%numnodes,2), stat=err)
       allocate(mat%Pbx(reg%numnodes,2), mat%Pby(reg%numnodes,2), mat%Pbz(reg%numnodes,2), stat=err)
       allocate(mat%N0(reg%numnodes), mat%N1(reg%numnodes), mat%N2(reg%numnodes), mat%N3(reg%numnodes), stat=err)
       allocate(mat%lvtermax(reg%numnodes),mat%lvtermay(reg%numnodes),mat%lvtermaz(reg%numnodes),stat=err)
       allocate(mat%lvtermbx(reg%numnodes),mat%lvtermby(reg%numnodes),mat%lvtermbz(reg%numnodes),stat=err)
       M4_ALLOC_ERROR(err,"InitializeMatLVFOURLVL")

       mat%Pax = 0.
       mat%Pay = 0.
       mat%Paz = 0.
       mat%Pbx = 0.
       mat%Pby = 0.
       mat%Pbz = 0.
       ninit = 1 - mat%Nstart(1) - mat%Nstart(2) - mat%Nstart(3)
       if (ninit < 0) then
          write(6,*) "Error in configuration: N0 is smaller zero" 
          !exit(1)
       endif
       mat%N0 = ninit
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
       !mat%c3a = 2. * mat%omegara * conv1 * mat%Ma**2 / ( 1/DT**2 + mat%gammala * DT ) bugfix 28/05/2013 Influence of bug is minor
       mat%c3a = 2. * mat%omegara * conv1 * mat%Ma**2 / ( 1/DT**2 + mat%gammala / DT )
       mat%c1b = ( 2. - mat%omegalb**2 * DT**2 ) / ( 1 + mat%gammalb * DT )
       mat%c2b = ( mat%gammalb * DT - 1 ) / ( 1 + mat%gammalb * DT )
       !mat%c3b = 2. * mat%omegarb * conv1 * mat%Mb**2 / ( 1/DT**2 + mat%gammalb * DT ) bugfix 28/05/2013 Influence of bug is minor
       mat%c3b = 2. * mat%omegarb * conv1 * mat%Mb**2 / ( 1/DT**2 + mat%gammalb / DT )
       

! for density integration
       mat%x1fac1 = 1. / 2. / mat%omegarb / DT * conv2
       mat%x1fac2 = mat%x1fac1 * DT * mat%gammalb / 2.
       mat%x2fac1 = 1. / 2. / mat%omegara / DT * conv2
       !mat%x2fac2 = mat%x2fac2 * DT * mat%gammala / 2. bugfix 16/05/2013 Influence of bug is minor
       mat%x2fac2 = mat%x2fac1 * DT * mat%gammala / 2.
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

! for langevin noise
       call zigset(mat%seed)

       mat%lvtermax = 0.
       mat%lvtermay = 0.
       mat%lvtermaz = 0.
       mat%lvtermbx = 0.
       mat%lvtermby = 0.
       mat%lvtermbz = 0.
       pfac = 2*dsqrt(1./SI_EPS0/REAL_DX)*SI_E/SI_C ! factor between P and Re(rho_{ij})
       mat%lvPa2 = mat%omegara*dsqrt(DT/mat%volfac/mat%N)*pfac*mat%Ma
       mat%lvPa1 = dsqrt(1/DT/mat%volfac/mat%N)*pfac*mat%Ma
       mat%lvPb2 = mat%omegarb*dsqrt(DT*dabs(mat%gammalb/mat%linefac-0.5*mat%gamma30-0.5*mat%gamma32)/mat%volfac/mat%N)*pfac*mat%Mb
       mat%lvPb1 = dsqrt(dabs(mat%gammalb/mat%linefac-0.5*mat%gamma30-0.5*mat%gamma32)/DT/mat%volfac/mat%N)*pfac*mat%Mb
       mat%lvN30 = dsqrt(mat%gamma30*DT/mat%volfac/mat%N)
       mat%lvN32 = dsqrt(mat%gamma32*DT/mat%volfac/mat%N)
       mat%lvN21 = dsqrt(mat%gamma21*DT/mat%volfac/mat%N)
       mat%lvN10 = dsqrt(mat%gamma10*DT/mat%volfac/mat%N)


! load from checkpoint file

       if ( load_state .and. detail_level .ge. 2 ) then

          read(UNITCHK) mat%Pax, mat%Pay, mat%Paz
          read(UNITCHK) mat%Pbx, mat%Pby, mat%Pbz
          read(UNITCHK) mat%N0, mat%N1, mat%N2, mat%N3

       end if

       M4_IFELSE_DBG({call EchoMatLVFourlvlObj(mat)},{call DisplayMatLVFourlvlObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatLVFourlvl")

  end subroutine InitializeMatLVFourlvl

!----------------------------------------------------------------------

  subroutine FinalizeMatLVFourlvl

    M4_MODLOOP_DECL({MATLVFOURLVL},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_WRITE_DBG(". enter FinalizeMatLVFourlvl")
    M4_MODLOOP_EXPR({MATLVFOURLVL},mat,{

       M4_MODOBJ_GETREG(mat,reg)

! save to checkpoint file

       if ( save_state .and. detail_level .ge. 2 ) then

          write(UNITCHK) mat%Pax, mat%Pay, mat%Paz
          write(UNITCHK) mat%Pbx, mat%Pby, mat%Pbz
          write(UNITCHK) mat%N0, mat%N1, mat%N2, mat%N3

       end if

       ! finalize mat object here
       deallocate(mat%Pax, mat%Pay, mat%Paz, mat%Pbx, mat%Pby, mat%Pbz, mat%N0, mat%N1, mat%N2, mat%N3)

    })
    M4_WRITE_DBG(". exit FinalizeMatLVFourlvl")

  end subroutine FinalizeMatLVFourlvl

!----------------------------------------------------------------------

  subroutine StepHMatLVFourlvl(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATLVFOURLVL},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8) :: pema, pena, pemb, penb
    real(kind=8) :: zetaax1, zetaay1, zetaaz1, zetaax2, zetaay2, zetaaz2 ! random numbers for Pa
    real(kind=8) :: zetabx1, zetaby1, zetabz1, zetabx2, zetaby2, zetabz2 ! random numbers for Pb
    real(kind=8) :: zetaN32, zetaN21, zetaN10, zetaN30 ! random numbers for nonradiative relaxation 
    real(kind=8) :: lEx, lEy, lEz, n3, n2, n1, ninva, ninvb, x1 ,x2
    M4_MODLOOP_EXPR({MATLVFOURLVL},mat,{

    ! this loops over all mat structures, setting mat

      M4_MODOBJ_GETREG(mat,reg)

      n = mod(ncyc-1+2,2) + 1
      m = mod(ncyc+2,2) + 1

      mat%cyc = m
      
        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

        ! calculate local field enhancement due to dielectric if mat%epsLFE=1
        lEx = Ex(i,j,k) * ( mat%epsLFE * ( 2. + 1./epsinvx(i,j,k) ) / 3. + ( 1. - mat%epsLFE ) )
        lEy = Ey(i,j,k) * ( mat%epsLFE * ( 2. + 1./epsinvy(i,j,k) ) / 3. + ( 1. - mat%epsLFE ) )
        lEz = Ez(i,j,k) * ( mat%epsLFE * ( 2. + 1./epsinvz(i,j,k) ) / 3. + ( 1. - mat%epsLFE ) )

        ! calculate second part of the density response (after the new E field got calculated)
        pema = mat%Pax(p,m) * lEx + mat%Pay(p,m) * lEy + mat%Paz(p,m) * lEz
        pena = mat%Pax(p,n) * lEx + mat%Pay(p,n) * lEy + mat%Paz(p,n) * lEz
        pemb = mat%Pbx(p,m) * lEx + mat%Pby(p,m) * lEy + mat%Pbz(p,m) * lEz
        penb = mat%Pbx(p,m) * lEx + mat%Pby(p,n) * lEy + mat%Pbz(p,n) * lEz
        x1 = mat%x1fac1 * ( penb - pemb ) + mat%x1fac2 * ( penb + pemb ) 
        x2 = mat%x2fac1 * ( pena - pema ) + mat%x2fac2 * ( pena + pema )
	! get random numbers
	!zetaN30 = rnor()*dsqrt(dabs(mat%N3(p)))*mat%lvN30
	zetaN32 = rnor()*dsqrt(dabs(mat%N3(p)))*mat%lvN32
	!zetaN21 = rnor()*dsqrt(dabs(mat%N2(p)))*mat%lvN21
	zetaN10 = rnor()*dsqrt(dabs(mat%N1(p)))*mat%lvN10
        mat%N3(p) = mat%N3(p) + mat%ac11 * x1 + zetaN32 !+ zetaN30
        mat%N2(p) = mat%N2(p) + mat%ac21 * x1 + mat%ac22 * x2 - zetaN32 !+ zetaN21
        mat%N1(p) = mat%N1(p) + mat%ac31 * x1 + mat%ac32 * x2 + mat%ac33 * x2 + zetaN10 !- zetaN21
        mat%N0(p) = mat%N0(p) + mat%ac41 * x1 + mat%ac42 * x2 + mat%ac43 * x2 + mat%ac44 * x1 - zetaN10 !- zetaN30

        ! calculate P(n+1) from P(n),P(n-1),E(n) and N(n)
        
        ! before: P(*,m) is P(n-1), P(*,n) is P(n)
        ninva = mat%N1(p) - mat%N2(p)
        ninvb = mat%N0(p) - mat%N3(p)
	
	! get random numbers for polarization (zeros only for test purposes)
	zetaax1 = rnor()*dsqrt(dabs(mat%N2(p)*(mat%gammala/mat%linefac-0.5*mat%gamma21)+0.5*mat%N3(p)*mat%gamma32))*mat%lvPa1
	zetaax2 = rnor()*dsqrt(dabs(mat%N2(p)*(mat%gammala/mat%linefac-0.5*mat%gamma21)+0.5*mat%N3(p)*mat%gamma32))*mat%lvPa2
	zetaay1 = rnor()*dsqrt(dabs(mat%N2(p)*(mat%gammala/mat%linefac-0.5*mat%gamma21)+0.5*mat%N3(p)*mat%gamma32))*mat%lvPa1
	zetaay2 = rnor()*dsqrt(dabs(mat%N2(p)*(mat%gammala/mat%linefac-0.5*mat%gamma21)+0.5*mat%N3(p)*mat%gamma32))*mat%lvPa2
	zetaaz1 = rnor()*dsqrt(dabs(mat%N2(p)*(mat%gammala/mat%linefac-0.5*mat%gamma21)+0.5*mat%N3(p)*mat%gamma32))*mat%lvPa1
	zetaaz2 = rnor()*dsqrt(dabs(mat%N2(p)*(mat%gammala/mat%linefac-0.5*mat%gamma21)+0.5*mat%N3(p)*mat%gamma32))*mat%lvPa2
	zetabx1 = rnor()*dsqrt(dabs(mat%N3(p)))*mat%lvPb1
	zetabx2 = rnor()*dsqrt(dabs(mat%N3(p)))*mat%lvPb2
	zetaby1 = rnor()*dsqrt(dabs(mat%N3(p)))*mat%lvPb1
	zetaby2 = rnor()*dsqrt(dabs(mat%N3(p)))*mat%lvPb2
	zetabz1 = rnor()*dsqrt(dabs(mat%N3(p)))*mat%lvPb1
	zetabz2 = rnor()*dsqrt(dabs(mat%N3(p)))*mat%lvPb2
        mat%Pax(p,m) = mat%c1a * mat%Pax(p,n) + mat%c2a * mat%Pax(p,m) + mat%c3a * lEx * ninva &
             + zetaax1 + zetaax2 - mat%lvtermax(p)
        mat%Pay(p,m) = mat%c1a * mat%Pay(p,n) + mat%c2a * mat%Pay(p,m) + mat%c3a * lEy * ninva &
             + zetaay1 + zetaay2 - mat%lvtermay(p)
        mat%Paz(p,m) = mat%c1a * mat%Paz(p,n) + mat%c2a * mat%Paz(p,m) + mat%c3a * lEz * ninva &
             + zetaaz1 + zetaaz2 - mat%lvtermaz(p)
        mat%Pbx(p,m) = mat%c1b * mat%Pbx(p,n) + mat%c2b * mat%Pbx(p,m) + mat%c3b * lEx * ninvb &
             + zetabx1 + zetabx2 - mat%lvtermbx(p)
        mat%Pby(p,m) = mat%c1b * mat%Pby(p,n) + mat%c2b * mat%Pby(p,m) + mat%c3b * lEy * ninvb &
             + zetaby1 + zetaby2 - mat%lvtermby(p)
        mat%Pbz(p,m) = mat%c1b * mat%Pbz(p,n) + mat%c2b * mat%Pbz(p,m) + mat%c3b * lEz * ninvb &
             + zetabz1 + zetabz2 - mat%lvtermbz(p)
	! save old random numbers for next step
	mat%lvtermax(p) = zetaax1
	mat%lvtermay(p) = zetaay1
	mat%lvtermaz(p) = zetaaz1
	mat%lvtermbx(p) = zetabx1
	mat%lvtermby(p) = zetaby1
	mat%lvtermbz(p) = zetabz1
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
  
  end subroutine StepHMatLVFourlvl


!----------------------------------------------------------------------

  subroutine StepEMatLVFourlvl(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATLVFOURLVL},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATLVFOURLVL},mat,{

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

  end subroutine StepEMatLVFourlvl

!----------------------------------------------------------------------

  subroutine SumJEMatLVFourlvl(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH),sum1a,sum1b,sum2a,sum2b,ninva,ninvb,d34a,d34b
    integer :: ncyc, m, n, idx
   
    M4_MODLOOP_DECL({MatLVFourlvl},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))


    M4_MODLOOP_EXPR({MatLVFourlvl},mat,{

    sum1a = 0.
    sum1b = 0.
    sum2a = 0.
    sum2b = 0.

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

!       if ( mode ) then
!          d34a = mat%d4a
!          d34b = mat%d4b
!       else
!          d34a = mat%d3a
!          d34b = mat%d3b
!       endif

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
 
       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

       if ( mask(i,j,k) ) then
!
!          ninva = mat%N1(p) - mat%N2(p)
!          ninvb = mat%N0(p) - mat%N3(p)
!
!          sum1a = sum1a + ( &
!4_IFELSE_TM({ !4_VOLEX(i,j,k) * w(1) * d34a * Ex(i,j,k) * mat%N *( mat%Pax(p,m) - mat%Pax(p,n) ) / DT + &
!               !4_VOLEY(i,j,k) * w(2) * d34a * Ey(i,j,k) * mat%N *( mat%Pay(p,m) - mat%Pay(p,n) ) / DT +},{0. +}) &
!4_IFELSE_TE({ !4_VOLEZ(i,j,k) * w(3) * d34a * Ez(i,j,k) * mat%N *( mat%Paz(p,m) - mat%Paz(p,n) ) / DT  },{0.  }) &
!               )

!          if ( ninva .NE. 0 ) then ! if ninv=0 then P should also equal 0 -> term would give NAN rather than 0

!             sum1a = sum1a + ( &
!4_IFELSE_TM({ !4_VOLEX(i,j,k) * w(1) * ( mat%d1a * mat%Pax(p,m) + mat%d2a * mat%Pax(p,n) ) / ninva / &
!               ( mat%epsLFE * ( 2. + 1./epsinvx(i,j,k) ) / 3. + ( 1. - mat%epsLFE ) ) * &
!               mat%N * ( mat%Pax(p,m) - mat%Pax(p,n) ) / DT + &
!               !4_VOLEY(i,j,k) * w(2) * ( mat%d1a * mat%Pay(p,m) + mat%d2a * mat%Pay(p,n) ) / ninva / &
!               ( mat%epsLFE * ( 2. + 1./epsinvy(i,j,k) ) / 3. + ( 1. - mat%epsLFE ) ) * &
!               mat%N * ( mat%Pay(p,m) - mat%Pay(p,n) ) / DT +},{0. +}) &
!4_IFELSE_TE({ !4_VOLEZ(i,j,k) * w(3) * ( mat%d1a * mat%Paz(p,m) + mat%d2a * mat%Paz(p,n) ) / ninva / &
!               ( mat%epsLFE * ( 2. + 1./epsinvz(i,j,k) ) / 3. + ( 1. - mat%epsLFE ) ) * &
!               mat%N * ( mat%Paz(p,m) - mat%Paz(p,n) ) / DT  },{0.  }) &
!                  )

!          endif

!          sum1b = sum1b + ( &
!4_IFELSE_TM({ !4_VOLEX(i,j,k) * w(1) * d34b * Ex(i,j,k) * mat%N *( mat%Pbx(p,m) - mat%Pbx(p,n) ) / DT + &
               !4_VOLEY(i,j,k) * w(2) * d34b * Ey(i,j,k) * mat%N *( mat%Pby(p,m) - mat%Pby(p,n) ) / DT +},{0. +}) &
!4_IFELSE_TE({ !4_VOLEZ(i,j,k) * w(3) * d34b * Ez(i,j,k) * mat%N *( mat%Pbz(p,m) - mat%Pbz(p,n) ) / DT  },{0.  }) &
!               )

!          if ( ninvb .NE. 0 ) then ! if ninv=0 then P should also equal 0 -> term would give NAN rather than 0

!             sum1b = sum1b + ( &
!4_IFELSE_TM({ !4_VOLEX(i,j,k) * w(1) * ( mat%d1b * mat%Pbx(p,m) + mat%d2b * mat%Pbx(p,n) ) / ninvb / &
!               ( mat%epsLFE * ( 2. + 1./epsinvx(i,j,k) ) / 3. + ( 1. - mat%epsLFE ) ) * &
!               mat%N * ( mat%Pbx(p,m) - mat%Pbx(p,n) ) / DT + &
!               !4_VOLEY(i,j,k) * w(2) * ( mat%d1b * mat%Pby(p,m) + mat%d2b * mat%Pby(p,n) ) / ninvb / &
!               ( mat%epsLFE * ( 2. + 1./epsinvy(i,j,k) ) / 3. + ( 1. - mat%epsLFE ) ) * &
!               mat%N * ( mat%Pby(p,m) - mat%Pby(p,n) ) / DT +},{0. +}) &
!4_IFELSE_TE({ !4_VOLEZ(i,j,k) * w(3) * ( mat%d1b * mat%Pbz(p,m) + mat%d2b * mat%Pbz(p,n) ) / ninvb / &
!               ( mat%epsLFE * ( 2. + 1./epsinvz(i,j,k) ) / 3. + ( 1. - mat%epsLFE ) ) * &
!               mat%N * ( mat%Pbz(p,m) - mat%Pbz(p,n) ) / DT  },{0.  }) &
!                  )

!          endif

!          sum2a = sum2a + ( &
!4_IFELSE_TM({ !4_VOLEX(i,j,k) * w(1) * Ex(i,j,k) * mat%N * ( mat%Pax(p,m) - mat%Pax(p,n) ) / DT + &
               !4_VOLEY(i,j,k) * w(2) * Ey(i,j,k) * mat%N * ( mat%Pay(p,m) - mat%Pay(p,n) ) / DT +},{0. +}) &
!4_IFELSE_TE({ !4_VOLEZ(i,j,k) * w(3) * Ez(i,j,k) * mat%N * ( mat%Paz(p,m) - mat%Paz(p,n) ) / DT  },{0.  }) &
!               )

!          sum2b = sum2b + ( &
!4_IFELSE_TM({ !4_VOLEX(i,j,k) * w(1) * Ex(i,j,k) * mat%N * ( mat%Pbx(p,m) - mat%Pbx(p,n) ) / DT + &
               !4_VOLEY(i,j,k) * w(2) * Ey(i,j,k) * mat%N * ( mat%Pby(p,m) - mat%Pby(p,n) ) / DT +},{0. +}) &
!4_IFELSE_TE({ !4_VOLEZ(i,j,k) * w(3) * Ez(i,j,k) * mat%N * ( mat%Pbz(p,m) - mat%Pbz(p,n) ) / DT  },{0.  }) &
!               )

       endif

       })      

!    sum(idx) = sum(idx)  + sum1a + sum1b
!    sum(idx+1) = sum(idx+1) + sum2a !- sum1a
!    sum(idx+2) = sum(idx+2) + sum2b !- sum1b
    idx = idx + NUMEBALCH

    })


  end subroutine SumJEMatLVFourlvl

!----------------------------------------------------------------------

  subroutine SumKHMatLVFourlvl(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum(MAXEBALCH)
    logical :: mode
    integer :: ncyc, idx

    M4_MODLOOP_DECL({MatLVFourlvl},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MatLVFourlvl},mat,{

    idx = idx + NUMEBALCH

    })

  end subroutine SumKHMatLVFourlvl

!----------------------------------------------------------------------

  subroutine DisplayMatLVFourlvlObj(mat)

    type(T_MATLVFOURLVL) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" lambdarinva=",TRIM(f2str(mat%lambdarinva,5)),&
    	" gammala=",TRIM(f2str(mat%gammala,5)),&
        " lamdarinvb=",TRIM(f2str(mat%lambdarinvb,5)),&
        " gammalb=",TRIM(f2str(mat%gammalb,5))
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatLVFourlvlObj

!----------------------------------------------------------------------

   subroutine EchoMatLVFourlvlObj(mat)

    type(T_MATLVFOURLVL) :: mat

    M4_WRITE_INFO({"--- matlvfourlvl # ",&
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
    

  end subroutine EchoMatLVFourlvlObj

  
!----------------------------------------------------------------------

end module matlvfourlvl

! Authors:  A.Pusch, J.Hamm 
! Modified: 08/08/2011
!
! =====================================================================


