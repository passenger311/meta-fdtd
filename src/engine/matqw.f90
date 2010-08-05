!-*- F90 -*------------------------------------------------------------
!
!  module: matqw / meta
!
!  Effective Maxwell Bloch material module for a quantum well.
!
!  subs:
!
!    InitializeMatQW
!    FinalizeMatQW
!    ReadMatQWObj
!    StepEMatQW
!    StepHMatQW
!    SumJEKHMatQW
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatQW module calculates the reponse of a quantum well in free carrier and effective mass approximation
! For each transition between states with the same momentum k we have the following equations of motion:
!
! d/dt d/dt P_k + 2 * gammal_k * d/dt P_k + omegal_k**2 P_k = 2 * omegar_k * conv1 * M^2 * E (n^e_k + n^h_k - 1)
! d/dt n_k = 2 M * E / (hbar * omegar_k) * (d/dt P_k + gammal_k * P_k)
! d/dt E = d/dt E* - epsinv * N * sum_k( M d/dt P_k )  
! 
! where E* is the electric field as calculated without the sources.  
! 
! To account for carrier-carrier interactions and recombinations, additional terms are introduced into the equations
! (see PhD thesis of Klaus Boehringer page 26 ff.)
! We also introduce electrical pumping terms and assume the tendency of the carriers to return to thermodynamic
! quasiequilibrium if disturbed.
!
! In natural HL units we choose hbar = c = 1, that is the elementary charge
! would be
!
! e = sqrt ( 4 PI alpha ) with alpha = 1/137.0360
!
! StepEMat updates the electric field
! StepHMat updates polarisations and densities:
! First the densities n^i_k are updated using the old and new electric fields. This gives the densities of step n-1/2
! Second the polarisations P_k are updated  using the new electric field and the new densities. 
! This gives the polarisations of step n, which are needed for the update of the electric fields.
! To get the full equations and their discretisations see PhD thesis of Klaus Boehringer.

module matqw

  use constant
  use parse
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MATHEAD_DECL({MATQW},MAXMATOBJ,{

  ! input parameters
  real(kind=8) :: kmax            ! maximal momentum [1/m]
  real(kind=8) :: dk              ! momentum step [1/m]
  real(kind=8) :: k0              ! minimal momentum [1/m]
  real(kind=8) :: Ebandgap        ! bandgap [eV]
  real(kind=8) :: me, mh          ! effective masses of electrons and holes [SI_ME]
  real(kind=8) :: Eefac, Ehfac    ! factors for conversion of momentum to energy

  real(kind=8) :: dcv             ! dipole strength of transitions [dx]
  real(kind=8) :: nr              ! number of active layers per dx
  real(kind=8) :: ninitial        ! start density [dx^(-3)]
  real(kind=8) :: T               ! temperature [K]
  integer      :: kCount         ! number of momentum states
  real(kind=8) :: pgamma          ! homogeneous broadening [1/dt]
  real(kind=8) :: egamma, hgamma  ! intraband relaxation rates [1/dt]
  real(kind=8) :: Lambda          ! pumping rate [1/dt]
  real(kind=8) :: gamma_nr        ! nonradiative recombination rate [1/dt]
  real(kind=8) :: gamma_sp        ! radiative recombination rate ( spontaneous emission) [1/dt^2]
  real(kind=8) :: gamma_aug       ! Auger recombination [1/dt]
  real(kind=8) :: hbar, delta_t, delta_s, epsilon0
  real(kind=8) :: conv1, conv2
  integer :: napprox, cyc

  ! calculated
  real(kind=8), dimension(:), pointer :: omegar        ! resonance frequencies
  real(kind=8), dimension(:), pointer :: omegal        ! lorentz frequencies
  real(kind=8), dimension(:), pointer :: pW            ! weight of state

  ! coefficients
  real(kind=8), dimension(:), pointer :: pPa
  real(kind=8), dimension(:), pointer :: pPb
  real(kind=8), dimension(:), pointer :: pPc
  
  ! polarisation field
  real(kind=8), dimension(:,:,:), pointer :: Px, Py, Pz
  real(kind=8), dimension(:), pointer :: Exo, Eyo, Ezo

  ! macroscopic values
  real(kind=8), dimension(:), pointer :: n
  real(kind=8), dimension(:), pointer :: phie
  real(kind=8), dimension(:), pointer :: phih
  real(kind=8), dimension(:), pointer :: G
  real(kind=8), dimension(:), pointer :: P_ma

  ! microscopic values
  real(kind=8), dimension(:,:), pointer :: gnx, gny, gnz
  real(kind=8), dimension(:,:), pointer :: fe, feo
  real(kind=8), dimension(:,:), pointer :: fh, fho
  real(kind=8), dimension(:,:), pointer :: fe_eq
  real(kind=8), dimension(:,:), pointer :: fh_eq

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatQWObj(funit,lcount)

    M4_MODREAD_DECL({MATQW}, funit,lcount,mat,reg,out)
    real(kind=8) :: v(2)
    real(kind=8) :: c(3)
    logical :: eof,err
    character(len=LINELNG) :: line

    M4_WRITE_DBG(". enter ReadMatQWObj")

    M4_MODREAD_EXPR({MATQW},funit,lcount,mat,reg,3,out,{ 

    ! read mat parameters here, as defined in mat data structure
    call readfloats(funit,lcount, v, 2)
    mat%k0 = v(1)
    mat%kmax = v(2)
    call readint(funit,lcount,mat%kCount)

    call readfloat(funit,lcount,mat%Ebandgap)

    call readfloats(funit,lcount, v, 2)
    mat%me = v(1)
    mat%mh = v(2)

    call readfloat(funit,lcount, mat%dcv)
    call readfloat(funit,lcount, mat%nr)
    call readfloat(funit,lcount, mat%ninitial)
    call readfloat(funit,lcount, mat%T)
    call readfloat(funit,lcount, mat%Lambda)
    call readfloats(funit,lcount, c, 3)
    mat%pgamma = c(1)
    mat%egamma = c(2)
    mat%hgamma = c(3)

    call readfloats(funit,lcount,c,3)
    mat%gamma_nr = c(1)
    mat%gamma_sp = c(2)
    mat%gamma_aug = c(3)


    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatQWObj")

  end subroutine ReadMatQWObj

!----------------------------------------------------------------------

  subroutine InitializeMatQW

    integer :: err, pk
    real (kind=8) :: conv1, conv2, kB, mu_help, mu_e, mu_h, Kval, Ee, Eh
    M4_MODLOOP_DECL({MATQW},mat) 
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    M4_WRITE_DBG(". enter InitializeMatQW")
    M4_MODLOOP_EXPR({MATQW},mat,{
    
       ! initialize mat object here
       !reg = regobj(mat%regidx)
       M4_MODOBJ_GETREG(mat,reg)
       allocate(mat%Px(reg%numnodes,mat%kCount,2), mat%Py(reg%numnodes,mat%kCount,2),&
                mat%Pz(reg%numnodes,mat%kCount,2), stat=err)
       allocate(mat%Exo(reg%numnodes), mat%Eyo(reg%numnodes), mat%Ezo(reg%numnodes), stat=err)
       allocate(mat%omegar(mat%kCount), mat%omegal(mat%kCount), mat%pW(mat%kCount), stat=err)
       allocate(mat%pPa(mat%kCount), mat%pPb(mat%kCount), mat%pPc(mat%kCount), stat=err)
       allocate(mat%n(reg%numnodes), mat%phie(reg%numnodes), mat%phih(reg%numnodes), stat=err)
       allocate(mat%G(reg%numnodes), mat%P_ma(reg%numnodes), stat=err)
       allocate(mat%gnx(reg%numnodes,mat%kCount), mat%gny(reg%numnodes,mat%kCount),&
            mat%gnz(reg%numnodes,mat%kCount), stat=err) 
       allocate(mat%fe(reg%numnodes,mat%kCount), mat%fh(reg%numnodes,mat%kCount), stat=err)
       allocate(mat%feo(reg%numnodes,mat%kCount), mat%fho(reg%numnodes,mat%kCount), stat=err)
       allocate(mat%fe_eq(reg%numnodes,mat%kCount), mat%fh_eq(reg%numnodes,mat%kCount), stat=err)
       M4_ALLOC_ERROR(err,"InitializeMatQW")
       
! From unit conversion of polarisation equation - SI to Computational units 
       mat%conv1 = SI_4PIALPHA
! From unit conversion of population equation - SI to Computation units
       mat%conv2 = ( REAL_DX * SI_C ) / SI_HBAR
       write(6,*) "conversion factor =", mat%conv1, mat%conv2
       ! convert important constants into computation units
       mat%delta_s = REAL_DX
       mat%delta_t = mat%delta_s / SI_C * dt
       kB = SI_KB * mat%delta_t * mat%delta_t / mat%delta_s / mat%delta_s
       mat%hbar = SI_HBAR * mat%delta_t / mat%delta_s / mat%delta_s
       mat%epsilon0 = SI_EPS0 * mat%delta_s * mat%delta_s * mat%delta_s / mat%delta_t / mat%delta_t &
                     / mat%delta_t / mat%delta_t
       
       mat%dk = (mat%kmax - mat%k0)/(mat%kCount-1)
       ! calculate factors for conversion of momentum to energy
       mat%Eefac = SI_HBAR*SI_HBAR/(2.*SI_ME*mat%me)
       mat%Eefac = mat%Eefac * mat%delta_t * mat%delta_t/mat%delta_s/mat%delta_s/mat%delta_s/mat%delta_s
       mat%Ehfac = SI_HBAR*SI_HBAR/(2.*SI_ME*mat%mh)
       mat%Ehfac = mat%Ehfac * mat%delta_t * mat%delta_t/mat%delta_s/mat%delta_s/mat%delta_s/mat%delta_s
       ! convert bandgap from eV to computation units
       mat%EBandGap = mat%EBandGap * mat%delta_t * mat%delta_t / mat%delta_s / mat%delta_s * SI_E
       ! convert auger decay rate
       mat%gamma_aug = mat%gamma_aug / mat%delta_s / mat%delta_s / mat%delta_s / mat%delta_s
       ! convert matrix element
       !mat%dcv = mat%dcv * SI_E
       write(6,*) "dcv =", mat%dcv
       ! convert temperature
       mat%T = kB * mat%T
       write(6,*) "T =",mat%T
       write(6,*) mat%delta_s, mat%delta_t, mat%hbar, mat%epsilon0, mat%Eefac, mat%EBandGap
       do pk = 1,mat%kCount
             Kval = (pk-1) * mat%dk + mat%k0
             Ee = mat%Eefac * Kval * Kval
             Eh = mat%Ehfac * Kval * Kval
             mat%omegar(pk) = (Ee + Eh + mat%EBandGap) / mat%hbar ! hbar is normalised to simulation units
             mat%omegal(pk) = sqrt(mat%omegar(pk)*mat%omegar(pk) + mat%pgamma * mat%pgamma)
             mat%pPa(pk) = + (2. - mat%omegal(pk) * mat%omegal(pk) * dt * dt)/(1. + mat%pgamma * dt)
             mat%pPb(pk) = - (1. - mat%pgamma * dt )/( 1. + mat%pgamma * dt )
             mat%pPc(pk) = - ( dt*dt )/(1. + mat%pgamma * dt )*mat%omegar(pk)*mat%dcv*mat%conv1
             mat%pW(pk) = mat%dk * Kval / PI
             if ( pk == 1 .OR. pk == mat%kCount ) then
                write(6,*) "k(", pk,") =", Kval, "dk =", mat%dk, "k0 =", mat%k0 
                write(6,*) "Ee =", Ee, "Eh =", Eh
                write(6,*) "omegar =", mat%omegar(pk)
                write(6,*) "omegal =", mat%omegal(pk)
                write(6,*) "Pa =", mat%pPa(pk), "Pb =", mat%pPb(pk), "Pc =", mat%pPc(pk)
                write(6,*) "pW =", mat%pW(pk)
             endif
       end do
       mat%pW(1) = mat%pW(1) * 0.5
       mat%pW(mat%kCount) = mat%pW(mat%kCount) * 0.5


       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
            mat%n(p) = mat%ninitial
            mat%G(p) = 0.

            mu_help = PI*mat%hbar*mat%hbar/mat%T/SI_ME
            mu_e = log(exp( mu_help * mat%n(p) / mat%me ) - 1.) 
            mu_h = log(exp( mu_help * mat%n(p) / mat%mh ) - 1.)
            write(6,*) "mu_help =", mu_help
            do pk = 1,mat%kCount
                Kval = (pk-1) * mat%dk + mat%k0
                Ee = mat%Eefac * Kval * Kval
                Eh = mat%Ehfac * Kval * Kval
                mat%Px(p,pk,1) = 0.
                mat%Py(p,pk,1) = 0.
                mat%Pz(p,pk,1) = 0.
                mat%Px(p,pk,2) = 0.
                mat%Py(p,pk,2) = 0.
                mat%Pz(p,pk,2) = 0.
                mat%gnx(p,pk) = 0.
                mat%gny(p,pk) = 0.
                mat%gnz(p,pk) = 0.
                mat%fe_eq(p,pk) = 1./(exp(Ee/mat%T - mu_e) + 1.)
                mat%fh_eq(p,pk) = 1./(exp(Eh/mat%T - mu_h) + 1.)
                if ( pk == 1 .OR. pk == mat%kCount ) then
                   write(6,*) "fe =", mat%fe_eq(p,pk) 
                   write(6,*) "fh =", mat%fh_eq(p,pk)
                endif
                mat%fe(p,pk) = mat%fe_eq(p,pk)
                mat%fh(p,pk) = mat%fh_eq(p,pk)
               
            end do
            mat%phie(p) = 0.
            mat%phih(p) = 0.

            do pk = 1, mat%kCount
               mat%phie(p) = mat%phie(p) + mat%fe_eq(p,pk) * (1. -  mat%fe(p,pk)) * mat%pW(pk)
               mat%phih(p) = mat%phih(p) + mat%fh_eq(p,pk) * (1. -  mat%fh(p,pk)) * mat%pW(pk)
            end do
    })
       mat%cyc = 1

       M4_IFELSE_DBG({call EchoMatQWObj(mat)},{call DisplayMatQWObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatQW")

  end subroutine InitializeMatQW

!----------------------------------------------------------------------

  subroutine FinalizeMatQW

    M4_MODLOOP_DECL({MATQW},mat)
    M4_WRITE_DBG(". enter FinalizeMatQW")
    M4_MODLOOP_EXPR({MATQW},mat,{

    ! finalize mat object here
    deallocate(mat%Px, mat%Py, mat%Pz, mat%omegar, mat%omegal, mat%pPa, mat%pPb, mat%pPc)
    deallocate(mat%n, mat%phie, mat%phih, mat%G, mat%P_ma)
    deallocate(mat%gnx, mat%gny, mat%gnz, mat%fe, mat%fh, mat%fe_eq, mat%fh_eq) 

    })
    M4_WRITE_DBG(". exit FinalizeMatQW")

  end subroutine FinalizeMatQW

!----------------------------------------------------------------------

  subroutine StepHMatQW(ncyc)

    integer :: ncyc, m, n, pk
    real :: Pxsum, Pysum, Pzsum, dP, lEx, lEy, lEz, wspk, Wsp, mu_help, mu_e, mu_h
    real :: factor_e, factor_h, pFa, pFb, pFc, pFd, pFe, pFf, factor, Kval, Ee, Eh
    real :: gain, c_g, c_cx, c_cy, c_cz, lambdae, lambdah
    M4_MODLOOP_DECL({MATQW},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATQW},mat,{
    

      M4_MODOBJ_GETREG(mat,reg)

        n = mod(ncyc-1+2,2) + 1
        m = mod(ncyc+2,2) + 1

        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
          lEx = Ex(i,j,k)
          lEy = Ey(i,j,k)
          lEz = Ez(i,j,k)
          !calculate new densities
          c_g = (mat%dcv/mat%epsilon0)/mat%hbar
          !double c_c1 = 2. * dt
          !double c_c2 = c_g * (pE[i]+pEo[i]) * c_c1;
          c_cx = 2. * dt * c_g * ( lEx + mat%Exo(p) )
          c_cy = 2. * dt * c_g * ( lEy + mat%Eyo(p) )
          c_cz = 2. * dt * c_g * ( lEz + mat%Ezo(p) )
          lambdae = mat%Lambda/mat%phie(p)
          lambdah = mat%Lambda/mat%phih(p)
      
          Wsp = 0.
          mat%G(p) = 0.
          do pk = 1,mat%kCount
             mat%G(p) = mat%G(p) + ( mat%gnx(p,pk) * c_cx/2/dt + mat%gny(p,pk) * c_cy/2/dt + &
                  mat%gnz(p,pk) * c_cy/2/dt ) * mat%pW(pk)
             ! Weisskopf-Wigner   
             wspk = mat%gamma_sp * mat%omegar(pk)*mat%omegar(pk)*mat%omegar(pk) * mat%fe(p,pk) * mat%fh(p,pk)
             factor_e = ( mat%gamma_nr + mat%egamma + lambdae*mat%fe_eq(p,pk) ) * dt
             pFa = 1./( 2. + factor_e )
             pFb = 2. * dt * ( lambdae + mat%egamma ) * mat%fe_eq(p,pk) * pFa
             pFc = ( 2. - factor_e ) * pFa
             factor_h = ( mat%gamma_nr + mat%hgamma + lambdah*mat%fh_eq(p,pk) ) * dt
             pFd = 1./( 2. + factor_h )
             pFe = 2. * dt * ( lambdah + mat%hgamma ) * mat%fh_eq(p,pk) * pFd;
             pFf = ( 2. - factor_h ) * pFd;
             ! integration steps
             gain = ( mat%gnx(p,pk) * c_cx + mat%gny(p,pk) * c_cy + mat%gnz(p,pk) * c_cz )
             mat%fe(p,pk) = gain * pFa + pFb + pFc * mat%fe(p,pk)&
                  - 2. * dt * pFa * wspk
             mat%fh(p,pk) = gain * pFd + pFe + pFf * mat%fh(p,pk)&
                  - 2. * dt * pFd * wspk
	
             Wsp = Wsp + wspk * mat%pW(pk);
          end do
    
          ! macroscopic level
          factor = mat%gamma_nr*dt
          mat%n(p) = (2.-factor)/(2.+factor)*mat%n(p)&
               + 2.*dt/(2.+factor)*( mat%G(p) + mat%Lambda - Wsp )&
               - 2.*dt/(2.+factor)*mat%gamma_aug*mat%n(p)*mat%n(p)*mat%n(p);
    
          mu_help = PI*mat%hbar*mat%hbar/mat%T/SI_ME
          mu_e = (exp( mu_help * mat%n(p) / mat%me ) - 1.) 
          mu_h = (exp( mu_help * mat%n(p) / mat%mh ) - 1.)
          mat%phie(p) = 0.
          mat%phih(p) = 0.
          do pk = 1, mat%kCount
             Kval = (pk-1) * mat%dk + mat%k0
             Ee = mat%Eefac * Kval * Kval
             Eh = mat%Ehfac * Kval * Kval
             mat%fe_eq(p,pk) = 1./(exp(Ee/mat%T)/mu_e + 1.)
             mat%fh_eq(p,pk) = 1./(exp(Eh/mat%T)/mu_h + 1.)
      
             mat%phie(p) = mat%phie(p) + mat%fe_eq(p,pk)*(1. -  mat%fe(p,pk)) * mat%pW(pk)
             mat%phih(p) = mat%phih(p) + mat%fh_eq(p,pk)*(1. -  mat%fh(p,pk)) * mat%pW(pk)
          end do

          mat%Exo(p) = lEx
          mat%Eyo(p) = lEy
          mat%Ezo(p) = lEz
          
          ! calculate new polarisations
          dP = 0.
          !mat%G(p) = 0.
    
          do pk = 1, mat%kCount
             mat%Px(p,pk,m) = mat%pPa(pk) * mat%Px(p,pk,n) + mat%pPb(pk) * mat%Px(p,pk,m) + mat%pPc(pk) * lEx * &
                  ( mat%fe(p,pk) + mat%fh(p,pk) - 1. )
             mat%Py(p,pk,m) = mat%pPa(pk) * mat%Py(p,pk,n) + mat%pPb(pk) * mat%Py(p,pk,m) + mat%pPc(pk) * lEy * &
                  ( mat%fe(p,pk) + mat%fh(p,pk) - 1. )
             mat%Pz(p,pk,m) = mat%pPa(pk) * mat%Pz(p,pk,n) + mat%pPb(pk) * mat%Pz(p,pk,m) + mat%pPc(pk) * lEz * &
                  ( mat%fe(p,pk) + mat%fh(p,pk) - 1. )
             dP = dP + ( mat%Px(p,pk,m) + mat%Py(p,pk,m) + mat%Pz(p,pk,m) ) * mat%pW(pk);
             mat%gnx(p,pk) = ( (mat%Px(p,pk,m)-mat%Px(p,pk,n))/dt + 0.5*mat%pgamma*(mat%Px(p,pk,m)+mat%Px(p,pk,n) ))&
                    /mat%omegar(pk) * mat%conv2
             mat%gny(p,pk) =  ( (mat%Py(p,pk,m)-mat%Py(p,pk,n))/dt + 0.5*mat%pgamma*(mat%Py(p,pk,m)+&
                    mat%Py(p,pk,n) ) ) /mat%omegar(pk) * mat%conv2 
             mat%gnz(p,pk) = ( (mat%Pz(p,pk,m)-mat%Pz(p,pk,n))/dt + 0.5*mat%pgamma*(mat%Pz(p,pk,m)+&
                    mat%Pz(p,pk,n) ) ) /mat%omegar(pk) * mat%conv2
             !mat%G(p) = mat%G(p) + ( mat%gnx(p,pk) + mat%gny(p,pk) + mat%gnz(p,pk) ) * mat%pW(pk)
          end do
    
          !mat%G(p) = mat%G(p) * mat%dcv / mat%epsilon0 / mat%hbar
          mat%P_ma(p) = ( 2. * mat%nr * mat%dcv ) * dP;
       
          })
    })
  end subroutine StepHMatQW


!----------------------------------------------------------------------

  subroutine StepEMatQW(ncyc)

    integer :: ncyc, m, n, pk
    real :: Pxsum, Pysum, Pzsum, dP, lEx, lEy, lEz, wspk, Wsp, mu_help, mu_e, mu_h
    real :: factor_e, factor_h, pFa, pFb, pFc, pFd, pFe, pFf, factor, Kval, Ee, Eh
    M4_MODLOOP_DECL({MATQW},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATQW},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
      
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)
       ! J(*,m) is P(n+1) and J(*,n) is P(n)
       Pxsum = 0
       Pysum = 0
       Pzsum = 0      
       do pk=1, mat%kCount
            Pxsum = Pxsum + mat%pW(pk) * ( mat%Px(p,pk,m) - mat%Px(p,pk,n) )
            Pysum = Pysum + mat%pW(pk) * ( mat%Py(p,pk,m) - mat%Py(p,pk,n) )
            Pzsum = Pzsum + mat%pW(pk) * ( mat%Pz(p,pk,m) - mat%Pz(p,pk,n) )
       end do   
M4_IFELSE_TM({
       Ex(i,j,k) = Ex(i,j,k) - w(1) * epsinvx(i,j,k) * Pxsum * mat%dcv * 2. * mat%nr
       !write(6,*) "Ey =", Ey(i,j,k), "Py =", Pysum * mat%dcv * 2.
       Ey(i,j,k) = Ey(i,j,k) - w(2) * epsinvy(i,j,k) * Pysum * mat%dcv * 2. * mat%nr
})
M4_IFELSE_TE({
       Ez(i,j,k) = Ez(i,j,k) - w(3) * epsinvz(i,j,k) * Pzsum * mat%dcv * 2. * mat%nr
})
       })      

 
    })

  end subroutine StepEMatQW

!----------------------------------------------------------------------

  real(kind=8) function SumJEMatQW(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n, pk
   
    M4_MODLOOP_DECL({MATQW},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    sum = 0

    M4_MODLOOP_EXPR({MATQW},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

       if ( mask(i,j,k) ) then
          do pk=1,mat%kCount
          sum = sum + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * w(1) * Ex(i,j,k) * mat%nr * mat%pW(pk) * ( mat%Px(p,pk,m) - mat%Px(p,pk,n) )& 
/ DT +}, {0. +}) &
M4_IFELSE_TM({ M4_VOLEY(i,j,k) * w(2) * Ey(i,j,k) * mat%nr * mat%pW(pk) * ( mat%Py(p,pk,m) - mat%Py(p,pk,n) )&
 / DT +}, {0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * w(3) * Ez(i,j,k) * mat%nr * mat%pW(pk) * ( mat%Pz(p,pk,m) - mat%Pz(p,pk,n) )&
 / DT  }, {0.  }) &
               )
          end do
       endif

       })      

    })
    
    SumJEMatQW = sum
    
  end function SumJEMatQW

!----------------------------------------------------------------------

  real(kind=8) function SumKHMatQW(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc

    SumKHMatQW = 0.

  end function SumKHMatQW
 
!----------------------------------------------------------------------

  subroutine DisplayMatQWObj(mat)

    type(T_MATQW) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" omegar0 =",TRIM(f2str(mat%omegar(1),5)),&
        " omegarmax =",TRIM(f2str(mat%omegar(mat%kCount),5)),&
    	" pgamma =",TRIM(f2str(mat%pgamma,5)),&
        " me =",TRIM(f2str(mat%me,5)),&
        " mh =",TRIM(f2str(mat%mh,5))
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatQWObj

!----------------------------------------------------------------------

   subroutine EchoMatQWObj(mat)

    type(T_MATQW) :: mat

    M4_WRITE_INFO({"--- matqw # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"omegar0 = ",mat%omegar(1) })
    M4_WRITE_INFO({"omegarmax = ",mat%omegar(mat%kCount) })
    M4_WRITE_INFO({"pgamma = ",mat%pgamma })
    M4_WRITE_INFO({"gamma_nr = ",mat%gamma_nr })
    M4_WRITE_INFO({"pump = ",mat%Lambda })
    M4_WRITE_INFO({"gammae = ",mat%egamma })
    M4_WRITE_INFO({"gammah = ",mat%hgamma })
    M4_WRITE_INFO({"me = ",mat%me })
    M4_WRITE_INFO({"mh = ",mat%mh })
    M4_WRITE_INFO({"T = ",mat%T })
    M4_WRITE_INFO({"nr = ",mat%nr })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatQWObj

  
!----------------------------------------------------------------------

end module matqw

! Authors:  A.Pusch, J.Hamm 
! Modified: 17/04/2010
!
! =====================================================================


