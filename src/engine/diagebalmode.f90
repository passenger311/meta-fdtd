!-*- F90 -*------------------------------------------------------------
!
!  module: diagebalmode / meta
!
!  Energy balance diagnostics module.
!
!----------------------------------------------------------------------


! =====================================================================
!
! Edit ME!!!

module diagebalmode

  use constant
  use checkpoint
  use mpiworld
  use reglist
  use outlist
  use grid
  use parse
  use fdtd
  use mat
  use matlorentz

  implicit none
  private
  save

! modulename, max_nr_modules
  M4_MODHEAD_DECL({DIAGEBALMODE},100,{

  integer :: ns, ne, dn  ! time stepping 
  ! Array storage order is the other way around in FORTRAN!!!
  ! indices filt_history,i,j,k,vec_component
  complex(kind=8), pointer, dimension(:,:,:,:,:) :: buf_E
  complex(kind=8), pointer, dimension(:,:,:,:,:) :: buf_H
  complex(kind=8), pointer, dimension(:,:,:,:,:,:) :: buf_P
  ! filter parameters
  real(kind=8) :: omega
  real(kind=8) :: beta
  complex(kind=8), pointer, dimension(:) :: alpha
  integer :: p
  real(kind=8) :: freq_sep
  integer :: h_pos_H, h_pos_E, h_pos_P

  ! copied from diagebal.f90
  ! spatially integrated energy contributions
  real(kind=8) :: dudt, ds, je, kh, res
  real(kind=8) :: dsx, dsy, dsz
  ! spatially and time integrated energy contributions
  real(kind=8) :: sumdudt, sumds, sumje, sumkh, sumres
  real(kind=8) :: sumdsx, sumdsy, sumdsz
  
  logical, pointer, dimension(:,:,:) :: mask
  
  ! partial contributions to energy terms
  real(kind=8), dimension(3) :: en, skh, sje

  real(kind=8), dimension(3) :: dsx1, dsx2
  real(kind=8), dimension(3) :: dsy1, dsy2
  real(kind=8), dimension(3) :: dsz1, dsz2

  })

contains

!----------------------------------------------------------------------

  subroutine ReadDiagEBalModeObj(funit,lcount)

    M4_MODREAD_DECL({DIAGEBALMODE}, funit,lcount,diag,reg,out)
    integer :: v(5)

    M4_WRITE_DBG(". enter ReadMatEBalModeObj")
    
    M4_MODREAD_EXPR({DIAGEBALMODE}, funit,lcount,diag,reg,0,out, {

    call readints(funit,lcount,v,3) 
    diag%ns = v(1)
    diag%ne = v(2)
    diag%dn = v(3)
    call readfloat(funit, lcount, diag%omega)
    call readint(funit, lcount, diag%p)
    call readfloat(funit, lcount, diag%freq_sep)

    if ( diag%ns .ge. diag%ne .or. diag%ns .lt. 0 .or. diag%dn .lt. 1 ) then
       M4_FATAL_ERROR({"BAD TIME WINDOW!"})
    end if

    if ( diag%p .le. 0 ) then
       M4_FATAL_ERROR({"BAD BUFFER SIZE!"})
    end if
    if ( diag%omega .lt. 0 .or. diag%omega .gt. PI ) then
       M4_FATAL_ERROR({"BAD OMEGA VALUE"})
    end if
    if ( diag%freq_sep .lt. 0 .or. diag%freq_sep .gt. PI ) then
       M4_FATAL_ERROR({"BAD FREQUENCY SEPERATION VALUE"})
    end if})

    M4_WRITE_DBG(". exit ReadMatEBalModeObj")

  end subroutine ReadDiagEBalModeObj

!----------------------------------------------------------------------

  subroutine InitializeDiagEBalMode

    integer :: i,j,k,q,c,r
    integer :: err
    real(kind=8) :: variance
    real(kind=8) :: gamma
    real(kind=8) :: prefactor
    real(kind=8) :: ndb
    real(kind=8) :: A
    integer :: shift
    M4_MODLOOP_DECL({DIAGEBALMODE},diag)
    type (T_REG) :: reg, reg_outset
    
    ! this is used later when calculating gamma from freq-separation
    ! it's basically the value of the normal distribution at the 
    ! 2*sigma point
    ndb = ( 1. / sqrt( 2. * PI ) ) * exp( -.5 * ( 2. ** 2. ) )

    M4_WRITE_DBG(". enter InitializeMatEBalMode")
    M4_MODLOOP_EXPR({DIAGEBALMODE},diag,{
    diag%dudt = 0.
    diag%ds = 0.
    diag%dsx = 0.
    diag%dsy = 0.
    diag%dsz = 0.
    diag%je = 0.
    diag%res = 0.

    diag%sumdudt = 0.
    diag%sumds = 0.
    diag%sumdsx = 0.
    diag%sumdsy = 0.
    diag%sumdsz = 0.
    diag%sumje = 0.
    diag%sumres = 0.
! this will load the region specified in the input file
    reg = regobj(diag%regidx)
    reg_outset = reg
    call OutsetRegion( reg_outset )

    M4_WRITE_DBG("allocate mask")
    allocate(diag%mask(reg_outset%is:reg_outset%ie,reg_outset%js:reg_outset%je, &
         reg_outset%ks:reg_outset%ke), stat = err)
    M4_WRITE_DBG("initialize mask")
    call SetMaskRegObj(reg,diag%mask,reg_outset%is,reg_outset%ie,reg_outset%js,reg_outset%je, &
         reg_outset%ks,reg_outset%ke)
    
    A = ndb ** ( 2. / dble(diag%p))
    gamma = ( -1. + ( A * cos( 2. * (diag%freq_sep/(2. * PI)) ) ) + &
         sqrt( (-4.* ((A-1.)**2.) ) + &
         (2. - ( 2. * A * cos( 2. * (diag%freq_sep/(2. * PI)) ) ))**2. ) ) &
         /(A-1.)
    shift = NINT( ( dble(diag%p) * gamma ) / ( 1. - gamma ) )
    

    M4_WRITE_INFO({"--- diagebalmode: p = ", diag%p})
    M4_WRITE_INFO({"--- diagebalmode: gamma = ", gamma})
    M4_WRITE_INFO({"--- diagebalmode: omega = ", diag%omega})
    M4_WRITE_INFO({"--- diagebalmode: approx. history len = ", 4 * shift})
    M4_WRITE_INFO({"--- diagebalmode: output must be shifted by = ", shift})
    
    allocate(diag%buf_E(0:diag%p-1,reg_outset%is:reg_outset%ie,reg_outset%js:reg_outset%je, &
         reg_outset%ks:reg_outset%ke,0:2), stat = err)
    M4_ALLOC_ERROR( err, "Unable to allocate filter buffer for E" )
    allocate(diag%buf_H(0:diag%p-1,reg_outset%is:reg_outset%ie,reg_outset%js:reg_outset%je, &
         reg_outset%ks:reg_outset%ke,0:2), stat = err)
    M4_ALLOC_ERROR( err, "Unable to allocate filter buffer for H" )
    M4_WRITE_INFO({"--- diagebalmode: numMATLORENTZobj = ", numMATLORENTZobj})
    allocate(diag%buf_P(0:diag%p-1,reg_outset%is:reg_outset%ie,reg_outset%js:reg_outset%je, &
         reg_outset%ks:reg_outset%ke,0:2,0:(numMATLORENTZobj-1)), stat = err)
    M4_ALLOC_ERROR( err, "Unable to allocate filter buffer for P" )
    do c=0,2
       do k=reg%ks,reg%ke
          do j=reg%js,reg%je
             do i=reg%is,reg%ie
                do q=0,diag%p-1
                   diag%buf_E(q,i,j,k,c) = (0,0)
                   diag%buf_H(q,i,j,k,c) = (0,0)
                   do r=0,numMATLORENTZobj-1
                      diag%buf_P(q,i,j,k,c,r) = (0,0)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    diag%h_pos_H = 0
    diag%h_pos_P = 0
    diag%h_pos_E = diag%p - 1
    allocate(diag%alpha(0:diag%p-1), stat = err )
    M4_ALLOC_ERROR( err, "Unable to allocate alpha array" )
    diag%beta = -2.d0 * ((1.d0 - gamma)**dble(diag%p))
    M4_WRITE_INFO({"--- diagebalmode: beta = ", diag%beta})
    do i=0,diag%p-1,1
       prefactor = (-1.d0) * binomial( diag%p, i+1 ) &
            * ( (-1.d0*gamma)**dble(i+1) )
       diag%alpha(diag%p-i-1) = &
            complex(prefactor*cos(-1.d0 * diag%omega * dble(i+1) ), &
            prefactor*sin(-1.d0 * diag%omega * dble(i+1)))
       M4_WRITE_INFO({"--- diagebalmode: alpha = ", diag%alpha(diag%p-i-1)})
    enddo

    ! load from checkpoint file

    if ( load_state .and. ( detail_level .eq. 1 .or. detail_level .eq. 3 ) ) then

!       read(UNITCHK) diag%i

    end if

    M4_IFELSE_DBG({call EchoDiagEBalModeObj(diag)})

    })
    M4_WRITE_DBG(". exit InitializeMatEBalMode")

    contains
      
      function binomial( n_binom, k_binom ) result(c_binom)
        
        integer,intent(in) :: n_binom
        integer,intent(in) :: k_binom
        integer :: kp_binom
        integer :: incr_binom
        real(kind=8) :: c_binom
        
        kp_binom = k_binom
        if ( k_binom > ( n_binom - k_binom ) ) then
           kp_binom = n_binom - k_binom
        endif
        c_binom = 1.d0
        do incr_binom=0,kp_binom-1,1
           c_binom = c_binom * dble(n_binom-incr_binom)
           c_binom = c_binom / dble(incr_binom+1)
        enddo !incr
      end function binomial

  end subroutine InitializeDiagEBalMode

!----------------------------------------------------------------------

  subroutine FinalizeDiagEBalMode

    integer :: i
    M4_MODLOOP_DECL({DIAGEBALMODE},diag)
    M4_WRITE_DBG(". enter FinalizeMatEBalMode")

    M4_MODLOOP_EXPR({DIAGEBALMODE},diag,{


! save to checkpoint file
    deallocate(diag%alpha)
    deallocate(diag%buf_E)
    deallocate(diag%buf_H)
    if ( save_state .and. ( detail_level .eq. 1 .or. detail_level .eq. 3 ) ) then

!       write(UNITCHK) diag%i

    end if

    deallocate(diag%mask)
    

    })
    M4_WRITE_DBG(". exit FinalizeMatEBalMode")

  end subroutine FinalizeDiagEBalMode

!----------------------------------------------------------------------
  subroutine OutsetRegion(rgn)
    type (T_REG),intent(inout) :: rgn

    rgn%is = rgn%is - 1
    rgn%ie = rgn%ie + 1
    M4_IFELSE_1D({},{
    rgn%js = rgn%js - 1
    rgn%je = rgn%je + 1
    })
    M4_IFELSE_2D({},{
    rgn%ks = rgn%ks - 1
    rgn%ke = rgn%ke + 1
    })
  end subroutine OutsetRegion


  subroutine StepHDiagEBalMode(ncyc)

    integer :: ncyc, m
    integer :: mod_pos
    integer :: q
    complex(kind=8) :: out_value
    M4_MODLOOP_DECL({DIAGEBALMODE},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))
    type (T_REG) :: reg_outset
    
    m    = mod(ncyc*2+2,3) + 1

    M4_MODLOOP_EXPR({DIAGEBALMODE},diag,{

    if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle
    M4_MODOBJ_GETREG(diag,reg)
    reg_outset = reg
    call OutsetRegion( reg_outset )
    M4_CALC_FILTER_FIELD(H)

    diag%dsx1(m) = 0.
    diag%dsx2(m) = 0.
    diag%dsy1(m) = 0.
    diag%dsy2(m) = 0.
    diag%dsz1(m) = 0.
    diag%dsz2(m) = 0.

    M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
    ! poynting vector

    diag%dsx1(m) = diag%dsx1(m) &
         + 0.5 * ( realpart( M4_VOLEY(i,j,k)/M4_SX(i,j,k) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k,1) * & 
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,2) ) - &
           dconjg( diag%buf_H(diag%h_pos_H,i-1,j,k,2) ) ) )&
         - M4_VOLEZ(i,j,k)/M4_SX(i,j,k) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k,2) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,1) ) - &
         dconjg( diag%buf_H(diag%h_pos_H,i-1,j,k,1) ) ) ) ) )
    
    diag%dsx2(m) = diag%dsx2(m) &
         + 0.5 * ( realpart( M4_VOLHZ(i,j,k)/M4_HSX(i,j,k) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,2) ) *&
         ( diag%buf_E(diag%h_pos_E,i+1,j,k,1) - diag%buf_E(diag%h_pos_E,i,j,k,1) ) ) &
         - M4_VOLHY(i,j,k)/M4_HSX(i,j,k) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,1) ) * &
         ( diag%buf_E(diag%h_pos_E,i+1,j,k,2) - diag%buf_E(diag%h_pos_E,i,j,k,2) ) ) ) )
    
    M4_IFELSE_1D({},{
    
    
    diag%dsy1(m) = diag%dsy1(m) &
         + 0.5 * ( realpart( M4_VOLEZ(i,j,k)/M4_SY(i,j,k) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k,2) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,0) ) - &
         dconjg( diag%buf_H(diag%h_pos_H,i,j-1,k,0) ) ) ) &
         - M4_VOLEX(i,j,k)/M4_SY(i,j,k) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k,0) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,2) ) &
         - dconjg( diag%buf_H(diag%h_pos_H,i,j-1,k,2) ) ) ) ) )
    
    diag%dsy2(m) = diag%dsy2(m) &
         + 0.5 * ( realpart( M4_VOLHX(i,j,k)/M4_HSY(i,j,k) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,0) ) * &
         ( diag%buf_E(diag%h_pos_E,i,j+1,k,2) - diag%buf_E(diag%h_pos_E,i,j,k,2) ) ) &
         - M4_VOLHZ(i,j,k)/M4_HSY(i,j,k) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,2) ) * &
         ( diag%buf_E(diag%h_pos_E,i,j+1,k,0) - diag%buf_E(diag%h_pos_E,i,j,k,0) ) ) ) )
    
    })
    
    M4_IFELSE_3D({
    diag%dsz1(m) = diag%dsz1(m) &
         + 0.5 * ( realpart( M4_VOLEX(i,j,k)/M4_SZ(i,j,k) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k,0) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,1) ) - &
           dconjg( diag%buf_H(diag%h_pos_H,i,j,k-1,1) ) ) ) &
         - M4_VOLEY(i,j,k)/M4_SZ(i,j,k) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k,1) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,0) ) - &
           dconjg( diag%buf_H(diag%h_pos_H,i,j,k-1,0) ) ) ) ) )

    diag%dsz2(m) = diag%dsz2(m) &
         + 0.5 * ( realpart ( M4_VOLHY(i,j,k)/M4_HSZ(i,j,k) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,1) ) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k+1,0) - diag%buf_E(diag%h_pos_E,i,j,k,0) ) ) &
         - M4_VOLHX(i,j,k)/M4_HSZ(i,j,k) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,0) ) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k+1,1) - diag%buf_E(diag%h_pos_E,i,j,k,1) ) ) ) )
    
    })
    
    })



    })
  
  end subroutine StepHDiagEBalMode


!----------------------------------------------------------------------
  subroutine StepEDiagEBalMode(ncyc)

    integer :: ncyc, m, mo, moo, m_P, n_P
    integer :: mod_pos;
    integer :: q
    complex(kind=8) :: out_value
    M4_MODLOOP_DECL({DIAGEBALMODE},diag)
    type(T_MATLORENTZ) :: mat   
    type(T_REG) :: reg_outset
    integer :: m4_m
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    m    = mod(ncyc*2+3,3) + 1
    mo   = mod(ncyc*2+2,3) + 1 ! half timestep back
    moo  = mod(ncyc*2+1,3) + 1 ! full timestep back


    M4_MODLOOP_EXPR({DIAGEBALMODE},diag,{

    if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle
    
    M4_MODOBJ_GETREG(diag,reg)
    reg_outset = reg
    call OutsetRegion( reg_outset )
    
    diag%h_pos_E = MOD(diag%h_pos_E + 1, diag%p)
    
    M4_CALC_FILTER_FIELD(E)

    diag%en(m) = 0.
    diag%dsx1(m) = 0.
    diag%dsx2(m) = 0.
    diag%dsy1(m) = 0.
    diag%dsy2(m) = 0. 
    diag%dsz1(m) = 0.
    diag%dsz2(m) = 0.

    M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
    
       ! energy density
       diag%en(m) =  diag%en(m) + ( &
            M4_VOLEX(i,j,k) / epsinvx(i,j,k) * &
            (REALPART(diag%buf_E(diag%h_pos_E,i,j,k,0))**2+&
            IMAGPART(diag%buf_E(diag%h_pos_E,i,j,k,0))**2) + &
            M4_VOLEY(i,j,k) / epsinvy(i,j,k) * &
            (REALPART(diag%buf_E(diag%h_pos_E,i,j,k,1))**2+&
            IMAGPART(diag%buf_E(diag%h_pos_E,i,j,k,1))**2) + &
            M4_VOLEZ(i,j,k) / epsinvz(i,j,k) * &
            (REALPART(diag%buf_E(diag%h_pos_E,i,j,k,2))**2+&
            IMAGPART(diag%buf_E(diag%h_pos_E,i,j,k,2))**2) + &
            M4_VOLHX(i,j,k) / M4_MUINVX(i,j,k) * &
            (REALPART(diag%buf_H(diag%h_pos_H,i,j,k,0))**2+&
            IMAGPART(diag%buf_H(diag%h_pos_H,i,j,k,0))**2) + &
            M4_VOLHY(i,j,k) / M4_MUINVY(i,j,k) * &
            (REALPART(diag%buf_H(diag%h_pos_H,i,j,k,1))**2+&
            IMAGPART(diag%buf_H(diag%h_pos_H,i,j,k,1))**2) + &
            M4_VOLHZ(i,j,k) / M4_MUINVZ(i,j,k) * &
            (REALPART(diag%buf_H(diag%h_pos_H,i,j,k,2))**2+&
            IMAGPART(diag%buf_H(diag%h_pos_H,i,j,k,2))**2) )
    
    ! poynting vector
    diag%dsx1(m) = diag%dsx1(m) &
         + 0.5 * ( realpart( M4_VOLEY(i,j,k)/M4_SX(i,j,k) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k,1) * & 
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,2) ) - &
         dconjg( diag%buf_H(diag%h_pos_H,i-1,j,k,2) ) ) ) &
         - M4_VOLEZ(i,j,k)/M4_SX(i,j,k) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k,2) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,1) ) - &
         dconjg( diag%buf_H(diag%h_pos_H,i-1,j,k,1) ) ) ) ) )
    
    diag%dsx2(m) = diag%dsx2(m) &
         + 0.5 * ( realpart( M4_VOLHZ(i,j,k)/M4_HSX(i,j,k) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,2) ) *&
         ( diag%buf_E(diag%h_pos_E,i+1,j,k,1) - diag%buf_E(diag%h_pos_E,i,j,k,1) ) ) &
         - M4_VOLHY(i,j,k)/M4_HSX(i,j,k) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,1) ) * &
         ( diag%buf_E(diag%h_pos_E,i+1,j,k,2) - diag%buf_E(diag%h_pos_E,i,j,k,2) ) ) ) )
    
    M4_IFELSE_1D({},{
    
    
    diag%dsy1(m) = diag%dsy1(m) &
         + 0.5 * ( realpart( M4_VOLEZ(i,j,k)/M4_SY(i,j,k) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k,2) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,0) ) - &
         dconjg( diag%buf_H(diag%h_pos_H,i,j-1,k,0) ) ) ) &
         - M4_VOLEX(i,j,k)/M4_SY(i,j,k) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k,0) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,2) ) - &
         dconjg( diag%buf_H(diag%h_pos_H,i,j-1,k,2) ) ) ) ) )
    
    diag%dsy2(m) = diag%dsy2(m) &
         + 0.5 * ( realpart( M4_VOLHX(i,j,k)/M4_HSY(i,j,k) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,0) ) * &
         ( diag%buf_E(diag%h_pos_E,i,j+1,k,2) - diag%buf_E(diag%h_pos_E,i,j,k,2) ) ) &
         - M4_VOLHZ(i,j,k)/M4_HSY(i,j,k) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,2) ) * &
         ( diag%buf_E(diag%h_pos_E,i,j+1,k,0) - diag%buf_E(diag%h_pos_E,i,j,k,0) ) ) ) )
    
    })
    
    M4_IFELSE_3D({
    diag%dsz1(m) = diag%dsz1(m) &
         + 0.5 * ( realpart( M4_VOLEX(i,j,k)/M4_SZ(i,j,k) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k,0) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,1) ) -&
         dconjg( diag%buf_H(diag%h_pos_H,i,j,k-1,1) ) ) ) &
         - M4_VOLEY(i,j,k)/M4_SZ(i,j,k) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k,1) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,0) ) -&
         dconjg( diag%buf_H(diag%h_pos_H,i,j,k-1,0) ) ) ) ) )

    diag%dsz2(m) = diag%dsz2(m) &
         + 0.5 * ( realpart( M4_VOLHY(i,j,k)/M4_HSZ(i,j,k) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,1) ) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k+1,0) - diag%buf_E(diag%h_pos_E,i,j,k,0) ) ) &
         - M4_VOLHX(i,j,k)/M4_HSZ(i,j,k) * &
         ( dconjg( diag%buf_H(diag%h_pos_H,i,j,k,0) ) * &
         ( diag%buf_E(diag%h_pos_E,i,j,k+1,1) - diag%buf_E(diag%h_pos_E,i,j,k,1) ) ) ) )
    
    })
    
    })

    ! --- do the JE loop
    do m4_m = 1, numMATLORENTZobj
       mat = MATLORENTZobj(m4_m)
       M4_MODOBJ_GETREG(mat,reg)
       ! --- m is the last calculated polarization field
       M4_CALC_FILTER_FIELD_OF_MAT(m4_m-1)    
       MATLORENTZobj(m4_m) = mat
    enddo
    m_P = diag%h_pos_P
    if (m_P .eq. 0) then
       n_P = diag%p - 1
    else
       n_P = m_P - 1
    endif
    diag%h_pos_H = MOD(diag%h_pos_H + 1, diag%p)
    diag%h_pos_P = MOD(diag%h_pos_P + 1, diag%p)
    diag%sje(m) = 0
    
    do m4_m = 1, numMATLORENTZobj
       mat = MATLORENTZobj(m4_m)
       M4_MODOBJ_GETREG(mat,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       if ( diag%mask(i,j,k) ) then
          diag%sje(m) = diag%sje(m) + &
               ( &
               M4_IFELSE_TM({ M4_VOLEX(i,j,k) * w(1) * abs(diag%buf_E(diag%h_pos_E,i,j,k,0) *&
               dconjg( diag%buf_P(m_P,i,j,k,0,m4_m-1) - diag%buf_P(n_P,i,j,k,0,m4_m-1) ) ) / DT +},{0. +}) &
               M4_IFELSE_TM({ M4_VOLEY(i,j,k) * w(2) * abs(diag%buf_E(diag%h_pos_E,i,j,k,1) *&
               dconjg( diag%buf_P(m_P,i,j,k,1,m4_m-1) - diag%buf_P(n_P,i,j,k,1,m4_m-1) ) ) / DT +},{0. +}) &
               M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * w(3) * abs(diag%buf_E(diag%h_pos_E,i,j,k,2) *&
               dconjg( diag%buf_P(m_P,i,j,k,2,m4_m-1) - diag%buf_P(n_P,i,j,k,2,m4_m-1) ) ) / DT},{0. }) &
               )
       endif
       })
       MATLORENTZobj(m4_m) = mat
    enddo

    ! WE need an extra 0.5 compared to diagebal because the average over a
    ! sin^2 wave is 0.5*amplitude
    diag%dudt = .25/DT * ( diag%en(m) - diag%en(moo) )
    
    diag%dsx = .5 * ( diag%dsx1(m) + diag%dsx1(mo) + diag%dsx2(mo) + diag%dsx2(moo) )
    diag%dsy = .5 * ( diag%dsy1(m) + diag%dsy1(mo) + diag%dsy2(mo) + diag%dsy2(moo) )
    diag%dsz = .5 * ( diag%dsz1(m) + diag%dsz1(mo) + diag%dsz2(mo) + diag%dsz2(moo) )

    diag%je = 0.5 * ( diag%sje(m) + diag%sje(mo) )
    diag%res = - ( diag%dudt + diag%ds + diag%je )

    diag%ds   = diag%dsx + diag%dsy + diag%dsz 

    diag%sumdudt = diag%sumdudt + DT * diag%dudt
    diag%sumds = diag%sumds + DT * diag%ds
    diag%sumdsx = diag%sumdsx + DT * diag%dsx
    diag%sumdsy = diag%sumdsy + DT * diag%dsy
    diag%sumdsz = diag%sumdsz + DT * diag%dsz
    diag%sumje = diag%sumje + DT * diag%je

    diag%sumres = - ( diag%sumdudt + diag%sumds + diag%sumje )
    
    })

  end subroutine StepEDiagEBalMode

!----------------------------------------------------------------------

   subroutine EchoDiagEBalModeObj(diag)

    type(T_DIAGEBALMODE) :: diag
 
    M4_WRITE_INFO({"--- diagebalmode # ",&
         TRIM(i2str(diag%idx))," ", TRIM(diag%type)})

    ! -- write parameters to console 

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(diag%regidx))

  end subroutine EchoDiagEBalModeObj
  
!----------------------------------------------------------------------

end module diagebalmode

!
! Authors:  J.Hamm, E.Kirby, F.Renn
! Modified: 07/07/2011
!
! =====================================================================


