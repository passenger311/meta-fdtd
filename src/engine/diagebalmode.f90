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

  implicit none
  private
  save

! modulename, max_nr_modules
  M4_MODHEAD_DECL({DIAGEBALMODE},100,{

  integer :: ns, ne, dn  ! time stepping 
  logical, pointer, dimension(:,:,:) :: mask
  ! Array storage order is the other way around in FORTRAN!!!
  ! indices filt_history,i,j,k,vec_component
  complex(kind=8), pointer, dimension(:,:,:,:,:) :: buf_E
  complex(kind=8), pointer, dimension(:,:,:,:,:) :: buf_H
  ! filter parameters
  real(kind=8) :: omega
  real(kind=8) :: beta
  real(kind=8) :: energy_density, poynting_vector
  complex(kind=8), pointer, dimension(:) :: alpha
  integer :: p
  real(kind=8) :: freq_sep
  real(kind=8) :: sum
  integer :: h_pos
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

    integer :: i,j,k,q,c
    integer :: err
    real(kind=8) :: variance
    real(kind=8) :: gamma
    real(kind=8) :: prefactor
    real(kind=8) :: ndb
    real(kind=8) :: A
    integer :: shift
    M4_MODLOOP_DECL({DIAGEBALMODE},diag)
    type (T_REG) :: reg
    

    ! this is used later when calculating gamma from freq-separation
    ! it's basically the value of the normal distribution at the 
    ! 2*sigma point
    ndb = ( 1. / sqrt( 2. * PI ) ) * exp( -.5 * ( 2. ** 2. ) )

    M4_WRITE_DBG(". enter InitializeMatEBalMode")
    M4_MODLOOP_EXPR({DIAGEBALMODE},diag,{

! this will load the region specified in the input file
    reg = regobj(diag%regidx)

    M4_WRITE_DBG("allocate mask")
    allocate(diag%mask(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX), stat = err)
    M4_WRITE_DBG("initialize mask")
    call SetMaskRegObj(reg,diag%mask,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX)
    
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
    
    allocate(diag%buf_E(0:diag%p-1,IMIN:IMAX,JMIN:JMAX,KMIN:KMAX &
         ,0:2), stat = err)
    M4_ALLOC_ERROR( err, "Unable to allocate filter buffer" )
    allocate(diag%buf_H(0:diag%p-1,IMIN:IMAX,JMIN:JMAX,KMIN:KMAX &
         ,0:2), stat = err)
    M4_ALLOC_ERROR( err, "Unable to allocate filter buffer" )
    do c=0,2
       do k=KMIN,KMAX
          do j=JMIN,JMAX
             do i=IMIN,IMAX
                do q=0,diag%p-1
                   diag%buf_E(q,i,j,k,c) = (0,0)
                   diag%buf_H(q,i,j,k,c) = (0,0)
                enddo
             enddo
          enddo
       enddo
    enddo
    diag%h_pos = 0
    allocate(diag%alpha(0:diag%p-1), stat = err )
    M4_ALLOC_ERROR( err, "Unable to allocate alpha array" )
    diag%beta = -1.d0 * ((1.d0 - gamma)**dble(diag%p))
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

  subroutine StepHDiagEBalMode(ncyc)

    integer :: ncyc

!    M4_MODLOOP_DECL({DIAGEBALMODE},diag)
!    M4_MODLOOP_EXPR({DIAGEBALMODE},diag,{
 !   diag%i = diag%i + 1
 !   })
  
  end subroutine StepHDiagEBalMode


!----------------------------------------------------------------------
  subroutine StepEDiagEBalMode(ncyc)

    integer :: ncyc
    integer :: mod_pos;
    integer :: q
    complex(kind=8) :: out_value
    M4_MODLOOP_DECL({DIAGEBALMODE},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    M4_MODLOOP_EXPR({DIAGEBALMODE},diag,{

    diag%energy_density = 0
    diag%poynting_vector = 0
    if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle
    
    M4_MODOBJ_GETREG(diag,reg)
    
    M4_CALC_FILTER_FIELD(E)
    M4_CALC_FILTER_FIELD(H)

    diag%h_pos = MOD(diag%h_pos + 1, diag%p)
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


