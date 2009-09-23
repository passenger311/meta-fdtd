!-*- F90 -*------------------------------------------------------------
!
!  module: matthreelvl / meta
!
!  three level maxwell-bloch equations
!
!  subs:
!
!    InitializeMatThreeLvl
!    FinalizeMatThreeLvl
!    ReadMatThreeLvlObj
!    StepEMatThreeLvl
!    StepHMatThreeLvl
!    SumJEKHMatThreeLvl
!
!----------------------------------------------------------------------


! =====================================================================
!
! The Matthreelvl module calculates the reponse of an effective  3 lvl 
! bloch system with omega1 <= omega2 <= omega3.
! All possible 3lvl schemes (lambda,cascade,V,etc.) can be simulated.
!
! The equations for the density matrix are given by
! d/dt rho = i/hbar * [H,rho]
! with H/hbar = [{omega1, -Omega12 , -Omega13},{-Omega12*, omega2, -Omega23}
!  ,{-Omega13*,-Omega23*,omega3}]
! and Omega_{ij} = conj(M_{ij}) * E / hbar
! 
! E = D - eps_{-1} * P
! P = (M12*rho12 + M13*rho13 + M23*rho23) * n
! n = number of 3lvl systems per grid cell.
! 
! The equations are solved by the 4th order Runge-Kutta method in 
! combination with an interpolation scheme for the electric field E, 
! which relies on the prior calculation of an effective electric 
! displacement D taking all other polarisations at the same grid cell 
! into account.
!
! The equations are formulated in SI-units where meters are replaced by
! dx and seconds by dx/clight. Therefore the electric field has to be
! converted from NHL-units into these grid units first and the
! polarisation has to be converted as well, so that E=D-eps^{-1}P.
! 
! All material calculations and the calculation of the new electric
! field E are being done in StepEMatThreeLvl.
! In StepHMatThreeLvl the electric field is replaced by the effective
! displacement D in order to calculate the new D out of d/dt D = rot(H)
!
! =====================================================================
module matthreelvl

  use constant
  use parse
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save


  M4_MATHEAD_DECL({MATTHREELVL},MAXMATOBJ,{

  ! input parameters
  real(kind=8) :: lambdainv12, lambdainv13, lambdainv23    
  real(kind=8) :: gamma12, gamma13, gamma23       ! polarisation dephasing 
  real(kind=8) :: sigma12, sigma13, sigma23       ! nonradiative recombination rates

  complex, dimension(3) :: M12(3), M13(3), M23(3) ! dipole lengths
  
  ! calculated

  real(kind=8) :: omega12, omega13, omega23 ! resonance frequencies

  ! coefficients
  
  ! polarizations
  complex(kind=8), dimension(:), pointer :: rho12, rho13, rho23 

  ! densities
  real(kind=8), dimension(:), pointer :: rho11, rho22, rho33
  ! starting values (have to add up to 1)
  real(kind=8) :: rho11_0, rho22_0, rho33_0

  ! number of 3-level atoms per unit-cell
  real(kind=8) :: n

  ! local field effect is included unless LFE==0
  real(kind=8) :: LFE
  
  ! E = D - lfeval*P (lfeval = 2/3 for LFE and lfeval=1 for no LFE)
  real(kind=8) :: lfeval

  ! field conversion factor
  real(kind=8) :: conv

  ! old electric field
  M4_FTYPE, dimension(:,:), pointer :: Dold
  
  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatthreelvlObj(funit,lcount)

    M4_MODREAD_DECL({MATTHREELVL}, funit,lcount,mat,reg,out)
    real(kind=8) :: v(3)
    complex(kind=8) :: c(3)
    logical :: eof,err
    character(len=LINELNG) :: line
    integer :: m


    M4_WRITE_DBG(". enter ReadMatthreelvlObj")
    
    M4_MODREAD_EXPR({MATTHREELVL},funit,lcount,mat,reg,3,out,{ 


    ! read mat parameters here, as defined in mat data structure
    call readfloats(funit, lcount, v, 3)
    mat%lambdainv12 = v(1)
    mat%lambdainv13 = v(2)
    mat%lambdainv23 = v(3)

    call readfloats(funit, lcount, v, 3)
    mat%gamma12 = v(1)
    mat%gamma13 = v(2)
    mat%gamma23 = v(3)

    call readfloats(funit, lcount, v, 3)
    mat%sigma12 = v(1)
    mat%sigma13 = v(2)
    mat%sigma23 = v(3)
    
    do m=1,3

       call readcomplexs(funit,lcount, c, 3)
M4_IFELSE_CF({
       mat%M12(m) = c(1)
       mat%M13(m) = c(2)
       mat%M23(m) = c(3)
},{
       mat%M12(m) = complex( real(c(1)), 0 )
       mat%M13(m) = complex( real(c(2)), 0 )
       mat%M23(m) = complex( real(c(3)), 0 )
})
    end do

    call readfloats(funit, lcount, v, 3)
    mat%rho11_0 = v(1)
    mat%rho22_0 = v(2)
    mat%rho33_0 = v(3)
    
    call readfloat(funit, lcount, mat%n)

    call readfloat(funit, lcount, mat%LFE)

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatthreelvlObj")

  end subroutine ReadMatthreelvlObj

!----------------------------------------------------------------------

  subroutine InitializeMatthreelvl

    type (T_REG) :: reg
    integer :: err
    M4_MODLOOP_DECL({MATTHREELVL},mat) 
    M4_WRITE_DBG(". enter InitializeMatthreelvl")
    M4_MODLOOP_EXPR({MATTHREELVL},mat,{
    
       ! initialize mat object here

       mat%omega12 = 2. * PI * mat%lambdainv12
       mat%omega13 = 2. * PI * mat%lambdainv13
       mat%omega23 = 2. * PI * mat%lambdainv23

       ! field conversion factor 

       mat%conv = REAL_DX**(3.5)/sqrt(SI_EPS0) * SI_E / SI_HBAR
       
       ! Is local field effect included
       if (mat%LFE==0) THEN
          mat%lfeval = 1.
       else
          mat%lfeval = 2./3.
       end if
          

       reg = regobj(mat%regidx)

       allocate(mat%rho12(reg%numnodes), mat%rho13(reg%numnodes), &
            mat%rho23(reg%numnodes), stat = err)

       M4_ALLOC_ERROR(err,"InitializeMatthreelvl")

       ! set initial polarisations to zero
       mat%rho12(:) = 0
       mat%rho13(:) = 0
       mat%rho23(:) = 0

       allocate(mat%rho11(reg%numnodes), mat%rho22(reg%numnodes), &
            mat%rho33(reg%numnodes), stat = err)

       M4_ALLOC_ERROR(err,"InitializeMatthreelvl")
       
       ! set initial densities
       mat%rho11(:) = mat%rho11_0
       mat%rho22(:) = mat%rho22_0
       mat%rho33(:) = mat%rho33_0

       allocate(mat%Dold(reg%numnodes,3), stat = err)

       M4_ALLOC_ERROR(err,"InitializeMatthreelvl")

       mat%Dold(:,:) = 0
       
       M4_IFELSE_DBG({call EchoMatthreelvlObj(mat)},{call DisplayMatthreelvlObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatthreelvl")

  end subroutine InitializeMatthreelvl

!----------------------------------------------------------------------

  subroutine FinalizeMatthreelvl

    M4_MODLOOP_DECL({MATTHREELVL},mat)
    M4_WRITE_DBG(". enter FinalizeMatthreelvl")
    M4_MODLOOP_EXPR({MATTHREELVL},mat,{

    ! finalize mat object here
    deallocate(mat%rho12,mat%rho13,mat%rho23)
    deallocate(mat%rho11,mat%rho22,mat%rho33)
    deallocate(mat%Dold)

    })
    M4_WRITE_DBG(". exit FinalizeMatthreelvl")

  end subroutine FinalizeMatthreelvl

!----------------------------------------------------------------------

  subroutine StepHMatthreelvl(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATTHREELVL},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATTHREELVL},mat,{

    ! this loops over all mat structures, setting mat

      M4_MODOBJ_GETREG(mat,reg)

      M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
      ! replace electric field with effective electric displacement
        M4_IFELSE_TM({
           Ex(i,j,k) = mat%Dold(p,1) / mat%conv
           Ey(i,j,k) = mat%Dold(p,2) / mat%conv
        })
        M4_IFELSE_TE({
           Ez(i,j,k) = mat%Dold(p,3) / mat%conv
        })

      })

    })
  
  end subroutine StepHMatthreelvl


!----------------------------------------------------------------------

  subroutine StepEMatthreelvl(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({MATTHREELVL},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    M4_FTYPE, dimension(3) :: Dold, Dnew, D, Etmp
    complex(kind=8) :: ra12, ra13, ra23
    real(kind=8) :: ei(3)
    integer :: m, m1, m2
    complex(kind=8) :: k12(4), k13(4), k23(4)
    real(kind=8) :: k11(4), k22(4), k33(4)
    complex(kind=8) :: rho12, rho13, rho23
    complex(kind=8) :: rho12o, rho13o, rho23o
    real(kind=8) :: rho11, rho22, rho33
    real(kind=8) :: rho11o, rho22o, rho33o
    real(kind=8) :: om12, om13, om23, eps_lx, eps_ly, eps_lz

    D = 0
    ! determine which spatial field components have to be calculated
    m1 = 1
    m2 = 3

M4_IFELSE_TM({}, {
       m1 = 3
})
M4_IFELSE_TE({}, {
       m2 = 2
})



    M4_MODLOOP_EXPR({MATTHREELVL},mat,{

       ! this loops over all mat structures, setting mat

    om12 = mat%omega12
    om13 = mat%omega13
    om23 = mat%omega23


    M4_MODOBJ_GETREG(mat,reg)

    M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! If the electric permittivity is not equal to 1 the electric field acting on the particle is
       ! strengthened by a factor l = (eps + 2)/3. In effect Etmp(diel.) = Etmp * l.
       eps_lx = (1./epsinvx(i,j,k) + 2.)/3.
       eps_ly = (1./epsinvy(i,j,k) + 2.)/3.
       eps_lz = (1./epsinvz(i,j,k) + 2.)/3.

       ! set factor ei for E = D - ei*sum(real(Mij*rhoij))
       ! ei includes unit conversion and number of systems per unit cell
       ei(1) = 2 * w(1) * epsinvx(i,j,k) * mat%n * SI_4PIALPHA
       ei(2) = 2 * w(2) * epsinvy(i,j,k) * mat%n * SI_4PIALPHA
       ei(3) = 2 * w(3) * epsinvz(i,j,k) * mat%n * SI_4PIAlPHA

       ! save old values of rho
       rho12o = mat%rho12(p)
       rho13o = mat%rho13(p)
       rho23o = mat%rho23(p)
       rho11o = mat%rho11(p)
       rho22o = mat%rho22(p)
       rho33o = mat%rho33(p)

       ! get new values of displacement and save old ones
       M4_IFELSE_TM({
          Dold(1) = mat%Dold(p,1)
          Dold(2) = mat%Dold(p,2)
          Dnew(1) = mat%conv * Ex(i,j,k)
          Dnew(2) = mat%conv * Ey(i,j,k)
          mat%Dold(p,1) = Dnew(1) 
          mat%Dold(p,2) = Dnew(2)
       })
       M4_IFELSE_TE({
          Dold(3) = mat%Dold(p,3)
          Dnew(3) = mat%conv * Ez(i,j,k)
          mat%Dold(p,3) = Dnew(3)
       })


       ! rk4: first step
       ! D(t) needed to calculate E(t)
       D(:) = Dold(:)
       
       rho12 = rho12o
       rho13 = rho13o
       rho23 = rho23o
       rho11 = rho11o
       rho22 = rho22o
       rho33 = rho33o
       
       do m= m1, m2

       ! Local field effect E(loc) = E(macroscopic) + P/3 = D - 2/3*P 
          Etmp(m) = D(m) - mat%lfeval * ei(m) * real( mat%M12(m)*rho12 + &
              mat%M13(m)*rho13 + mat%M23(m)*rho23 )

       end do
       
       ! calculate rabi frequencies at time t
       ra12 = conjg(mat%M12(1))*Etmp(1)*eps_lx + conjg(mat%M12(2))*Etmp(2)*eps_ly + &
            conjg(mat%M12(3))*Etmp(3)*eps_lz
       ra13 = conjg(mat%M13(1))*Etmp(1)*eps_lx + conjg(mat%M13(2))*Etmp(2)*eps_ly + & 
            conjg(mat%M13(3))*Etmp(3)*eps_lz
       ra23 = conjg(mat%M23(1))*Etmp(1)*eps_lx + conjg(mat%M23(2))*Etmp(2)*eps_ly + & 
            conjg(mat%M23(3))*Etmp(3)
       

       ! calculate first ks
       k12(1) = IMAG*( -rho12*om12 + ra12*(rho11-rho22) + conjg(ra23)*rho13 - ra13*conjg(rho23) ) &
            - mat%gamma12 * rho12
       k13(1) = IMAG*( -rho13*om13 + ra13*(rho11-rho33) + ra23*rho12 - ra12*rho23 ) &
            - mat%gamma13 * rho13
       k23(1) = IMAG*( -rho23*om23 + ra23*(rho22-rho33) - conjg(ra12)*rho13 + ra13*conjg(rho12) ) &
            - mat%gamma23 * rho23
       k11(1) = - 2*aimag(rho12*conjg(ra12)) - 2*aimag(rho13*conjg(ra13)) &
            + mat%sigma12 * rho22 + mat%sigma13 * rho33
       k22(1) = 2*aimag(rho12*conjg(ra12))  - 2*aimag(rho23*conjg(ra23)) &
            - mat%sigma12 * rho22 + mat%sigma23 * rho33
       k33(1) = 2*aimag(rho13*conjg(ra13)) + 2*aimag(rho23*conjg(ra23)) &
            - (mat%sigma13 + mat%sigma23) * rho33


       ! rk4: second step
       ! D interpolated to t+dt/2 in order to calculate E(t+dt/2)
       D(:) = (Dold(:) + Dnew(:))/2
       
       ! rho interpolated to t+dt/2 according to runge-kutta scheme
       rho12 = rho12o + k12(1) * DT / 2
       rho13 = rho13o + k13(1) * DT / 2
       rho23 = rho23o + k23(1) * DT / 2
       rho11 = rho11o + k11(1) * DT / 2
       rho22 = rho22o + k22(1) * DT / 2
       rho33 = rho33o + k33(1) * DT / 2


       do m= m1, m2

       ! Local field effect E(loc) = E(macroscopic) + P/3 = D - 2/3*P
          Etmp(m) = D(m) - mat%lfeval * ei(m) * real( mat%M12(m)*rho12 + &
              mat%M13(m)*rho13 + mat%M23(m)*rho23 )

       end do

       ! calculate rabi frequencies at t+dt/2
       ra12 = conjg(mat%M12(1))*Etmp(1)*eps_lx + conjg(mat%M12(2))*Etmp(2)*eps_ly + &
            conjg(mat%M12(3))*Etmp(3)*eps_lz
       ra13 = conjg(mat%M13(1))*Etmp(1)*eps_lx + conjg(mat%M13(2))*Etmp(2)*eps_ly + &
            conjg(mat%M13(3))*Etmp(3)*eps_lz
       ra23 = conjg(mat%M23(1))*Etmp(1)*eps_lx + conjg(mat%M23(2))*Etmp(2)*eps_ly + &
            conjg(mat%M23(3))*Etmp(3)*eps_lz       


       k12(2) = IMAG*( -rho12*om12 + ra12*(rho11-rho22) + conjg(ra23)*rho13 - ra13*conjg(rho23) ) &
            - mat%gamma12 * rho12
       k13(2) = IMAG*( -rho13*om13 + ra13*(rho11-rho33) + ra23*rho12 - ra12*rho23 ) &
            - mat%gamma13 * rho13
       k23(2) = IMAG*( -rho23*om23 + ra23*(rho22-rho33) - conjg(ra12)*rho13 + ra13*conjg(rho12) ) &
            - mat%gamma23 * rho23
       k11(2) = - 2*aimag(rho12*conjg(ra12)) - 2*aimag(rho13*conjg(ra13)) &
            + mat%sigma12 * rho22 + mat%sigma13 * rho33
       k22(2) = 2*aimag(rho12*conjg(ra12))  - 2*aimag(rho23*conjg(ra23)) &
            - mat%sigma12 * rho22 + mat%sigma23 * rho33
       k33(2) = 2*aimag(rho13*conjg(ra13)) + 2*aimag(rho23*conjg(ra23)) &
            - (mat%sigma13 + mat%sigma23) * rho33


       ! rk4: third step
       ! rho interpolated to t+dt/2 for a second time. D stays the same
       
       rho12 = rho12o + k12(2) * DT / 2
       rho13 = rho13o + k13(2) * DT / 2
       rho23 = rho23o + k23(2) * DT / 2
       rho11 = rho11o + k11(2) * DT / 2
       rho22 = rho22o + k22(2) * DT / 2
       rho33 = rho33o + k33(2) * DT / 2

       do m= m1, m2

       ! Local field effect E(loc) = E(macroscopic) + P/3 = D - 2/3*P
          Etmp(m) = D(m) - mat%lfeval * ei(m) * real( mat%M12(m)*rho12 + &
              mat%M13(m)*rho13 + mat%M23(m)*rho23 )

       end do
       
       ! calculate rabi frequencies at t+dt/2
       ra12 = conjg(mat%M12(1))*Etmp(1)*eps_lx + conjg(mat%M12(2))*Etmp(2)*eps_ly + & 
            conjg(mat%M12(3))*Etmp(3)*eps_lz
       ra13 = conjg(mat%M13(1))*Etmp(1)*eps_lx + conjg(mat%M13(2))*Etmp(2)*eps_ly + & 
            conjg(mat%M13(3))*Etmp(3)*eps_lz
       ra23 = conjg(mat%M23(1))*Etmp(1)*eps_lx + conjg(mat%M23(2))*Etmp(2)*eps_ly + &
            conjg(mat%M23(3))*Etmp(3)*eps_lz       


       k12(3) = IMAG*( -rho12*om12 + ra12*(rho11-rho22) + conjg(ra23)*rho13 - ra13*conjg(rho23) ) &
            - mat%gamma12 * rho12
       k13(3) = IMAG*( -rho13*om13 + ra13*(rho11-rho33) + ra23*rho12 - ra12*rho23 ) &
            - mat%gamma13 * rho13
       k23(3) = IMAG*( -rho23*om23 + ra23*(rho22-rho33) - conjg(ra12)*rho13 + ra13*conjg(rho12) ) &
            - mat%gamma23 * rho23
       k11(3) = - 2*aimag(rho12*conjg(ra12)) - 2*aimag(rho13*conjg(ra13)) &
            + mat%sigma12 * rho22 + mat%sigma13 * rho33
       k22(3) = 2*aimag(rho12*conjg(ra12))  - 2*aimag(rho23*conjg(ra23)) &
            - mat%sigma12 * rho22 + mat%sigma23 * rho33
       k33(3) = 2*aimag(rho13*conjg(ra13)) + 2*aimag(rho23*conjg(ra23)) &
            - (mat%sigma13 + mat%sigma23) * rho33


       ! rk4: fourth step
       ! D(t+dt) for E(t+dt)
       D(:) = Dnew(:)

       ! rho(t+dt)
       rho12 = rho12o + k12(3) * DT
       rho13 = rho13o + k13(3) * DT
       rho23 = rho23o + k23(3) * DT
       rho11 = rho11o + k11(3) * DT
       rho22 = rho22o + k22(3) * DT
       rho33 = rho33o + k33(3) * DT   

       do m= m1, m2

       ! Local field effect E(loc) = E(macroscopic) + P/3 = D - 2/3*P
       Etmp(m) = D(m) - mat%lfeval * ei(m) * real( mat%M12(m)*rho12 + &
              mat%M13(m)*rho13 + mat%M23(m)*rho23 )

       end do

       ! calculate rabi frequencies at t+dt
       ra12 = conjg(mat%M12(1))*Etmp(1)*eps_lx + conjg(mat%M12(2))*Etmp(2)*eps_ly + &
            conjg(mat%M12(3))*Etmp(3)*eps_lz
       ra13 = conjg(mat%M13(1))*Etmp(1)*eps_lx + conjg(mat%M13(2))*Etmp(2)*eps_ly + &
            conjg(mat%M13(3))*Etmp(3)*eps_lz
       ra23 = conjg(mat%M23(1))*Etmp(1)*eps_lx + conjg(mat%M23(2))*Etmp(2)*eps_ly + &
            conjg(mat%M23(3))*Etmp(3)*eps_lz
 

       k12(4) = IMAG*( -rho12*om12 + ra12*(rho11-rho22) + conjg(ra23)*rho13 - ra13*conjg(rho23) ) &
            - mat%gamma12 * rho12
       k13(4) = IMAG*( -rho13*om13 + ra13*(rho11-rho33) + ra23*rho12 - ra12*rho23 ) &
            - mat%gamma13 * rho13
       k23(4) = IMAG*( -rho23*om23 + ra23*(rho22-rho33) - conjg(ra12)*rho13 + ra13*conjg(rho12) ) &
            - mat%gamma23 * rho23
       k11(4) = - 2*aimag(rho12*conjg(ra12)) - 2*aimag(rho13*conjg(ra13)) &
            + mat%sigma12 * rho22 + mat%sigma13 * rho33
       k22(4) = 2*aimag(rho12*conjg(ra12))  - 2*aimag(rho23*conjg(ra23)) &
            - mat%sigma12 * rho22 + mat%sigma23 * rho33
       k33(4) = 2*aimag(rho13*conjg(ra13)) + 2*aimag(rho23*conjg(ra23)) &
            - (mat%sigma13 + mat%sigma23) * rho33

       
       ! value at end point

       rho12 = rho12o + DT/6.*( k12(1) + 2.*k12(2) + 2.*k12(3) + k12(4) ) 
       rho13 = rho13o + DT/6.*( k13(1) + 2.*k13(2) + 2.*k13(3) + k13(4) ) 
       rho23 = rho23o + DT/6.*( k23(1) + 2.*k23(2) + 2.*k23(3) + k23(4) ) 
       rho11 = rho11o + DT/6.*( k11(1) + 2.*k11(2) + 2.*k11(3) + k11(4) ) 
       rho22 = rho22o + DT/6.*( k22(1) + 2.*k22(2) + 2.*k22(3) + k22(4) ) 
       rho33 = rho33o + DT/6.*( k33(1) + 2.*k33(2) + 2.*k33(3) + k33(4) ) 


       ! update fields

       do m= m1, m2
       
          Etmp(m) = D(m) - ei(m) * &
            real( mat%M12(m)*rho12 + mat%M13(m)*rho13 + mat%M23(m)*rho23 )

       end do
       
       ! update density matrix
       mat%rho12(p) = rho12
       mat%rho13(p) = rho13
       mat%rho23(p) = rho23

       mat%rho11(p) = rho11
       mat%rho22(p) = rho22
       mat%rho33(p) = rho33

       ! convert electric field back into NHL-units and write it back
M4_IFELSE_TM({
       Ex(i,j,k) = Etmp(1) / mat%conv
       Ey(i,j,k) = Etmp(2) / mat%conv
})
M4_IFELSE_TE({
       Ez(i,j,k) = Etmp(3) / mat%conv
})



       }) 

    })

  end subroutine StepEMatthreelvl

!----------------------------------------------------------------------

  real(kind=8) function SumJEMatthreelvl(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
   
    M4_MODLOOP_DECL({MATTHREELVL},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    sum = 0
    
    SumJEMatthreelvl = sum
    
  end function SumJEMatthreelvl

!----------------------------------------------------------------------

  real(kind=8) function SumKHMatthreelvl(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc

    SumKHMatthreelvl = 0.

  end function SumKHMatthreelvl
 
!----------------------------------------------------------------------

  subroutine DisplayMatthreelvlObj(mat)

    type(T_MATTHREELVL) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" lambdainv12=",TRIM(f2str(mat%lambdainv12,5)),&
    	" lambdainv13=",TRIM(f2str(mat%lambdainv13,5)),&
    	" lambdainv23=",TRIM(f2str(mat%lambdainv23,5))
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatthreelvlObj

!----------------------------------------------------------------------

   subroutine EchoMatthreelvlObj(mat)

    type(T_MATTHREELVL) :: mat

    M4_WRITE_INFO({"--- matthreelvl # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"lambdainv12 = ",mat%lambdainv12 })
    M4_WRITE_INFO({"lambdainv13 = ",mat%lambdainv13 })
    M4_WRITE_INFO({"lambdainv23 = ",mat%lambdainv23 })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatthreelvlObj

  
!----------------------------------------------------------------------

end module matthreelvl

! Authors: J.Hamm, A.Pusch 
! Modified: 9/9/2009
!
! =====================================================================


