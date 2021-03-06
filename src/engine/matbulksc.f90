!-*- F90 -*------------------------------------------------------------
!
!  module: matbulksc / meta
!
!  Effective Maxwell Bloch material module for isotropic 4 level system.
!
!  subs:
!
!    InitializeMatBulksc
!    FinalizeMatBulksc
!    ReadMatBulkscObj
!    StepEMatBulksc
!    StepHMatBulksc
!    SumJEMatBulksc
!    SumKHMatBulksc
!
!----------------------------------------------------------------------


! =====================================================================
!
! Bulk semiconductor model
!
! NOTE: calculations are performed in SI units. This means that the 
! electric field is converted before use into SI units (conv1) and the 
! integrated polarisation back into computational units (conv2).



module matbulksc

  use omp_lib
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

  M4_MATHEAD_DECL({MATBULKSC},MAXMATOBJ,{

  ! input parameters
  real(kind=8) :: gammap         ! dephasing rate [s]
  real(kind=8) :: M              ! dipole matrix length [m]
  real(kind=8) :: egap           ! bandgap [eV]
  real(kind=8) :: me             ! effective mass [electron mass]
  real(kind=8) :: mh             ! effective mass [electron mass]
  real(kind=8) :: N0             ! initial density [m^-3]
  real(kind=8) :: gammanr        ! nonradiative transition rates [s]
  real(kind=8) :: temp           ! temperature [K]
  real(kind=8) :: kmax           ! k cutoff value [1/m]
  real(kind=8) :: pump           ! pump rate [1/s m^-3]
  integer :: numk                ! number of k points for band discretisation

  integer :: cyc
  
  ! calculated

  real(kind=8) :: domega
  real(kind=8) :: mr0, me0, mh0     ! mass in [kg]
  real(kind=8) :: omegagap

  ! coefficients

  real(kind=8) :: cc, dd, c0, c1, c2, c3, K1, K2, K3
  real(kind=8) :: conv1, conv2
  real(kind=8) :: sigma, d2

  ! polarisation field 
  M4_FTYPE, dimension(:,:,:,:), pointer :: Pk

  ! sum of polarisations

  M4_FTYPE, dimension(:,:,:), pointer :: Qsum
  M4_FTYPE, dimension(:,:,:), pointer :: Psum

  ! index fields

  integer, dimension(:), pointer :: ipos,jpos,kpos


  ! density
  real(kind=8), dimension(:), pointer :: N

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatBulkSCObj(funit,lcount)

    M4_MODREAD_DECL({MATBULKSC}, funit,lcount,mat,reg,out)
    logical :: eof,err

    M4_WRITE_DBG(". enter ReadMatBulkSCObj")

    M4_MODREAD_EXPR({MATBULKSC},funit,lcount,mat,reg,3,out,{ 

    call readfloat(funit, lcount, mat%gammap)
    call readfloat(funit, lcount, mat%M)
    call readfloat(funit, lcount, mat%egap)
    call readfloat(funit, lcount, mat%me)
    call readfloat(funit, lcount, mat%mh)
    call readfloat(funit, lcount, mat%N0)
    call readfloat(funit, lcount, mat%gammanr)
    call readfloat(funit, lcount, mat%pump)
    call readfloat(funit, lcount, mat%temp) 
    call readfloat(funit, lcount, mat%kmax)
    call readint(funit, lcount, mat%numk)

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatBulkSCObj")

  end subroutine ReadMatBulkSCObj

!----------------------------------------------------------------------

  subroutine InitializeMatBulkSC

    integer :: err

    M4_MODLOOP_DECL({MATBULKSC},mat) 
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_WRITE_DBG(". enter InitializeMatBulkSC")
    M4_MODLOOP_EXPR({MATBULKSC},mat,{
    
       ! initialize mat object here

       reg = regobj(mat%regidx)

       allocate(M4_PK(3,2,0:mat%numk,reg%numnodes), stat=err)
       allocate(mat%N(reg%numnodes), M4_QSUM(3, 2, reg%numnodes), M4_PSUM(3, 2, reg%numnodes), stat=err)
       allocate(mat%ipos(reg%numnodes), mat%jpos(reg%numnodes), mat%kpos(reg%numnodes), stat=err)

       M4_ALLOC_ERROR(err,"InitializeMatBulkSC")

       mat%Pk = 0.
       mat%Psum = 0.
       mat%Qsum = 0.
       mat%N = mat%N0
  
! Initialize omegak(k)

       mat%me0 = mat%me * SI_ME
       mat%mh0 = mat%mh * SI_ME
       mat%omegagap = mat%egap * SI_E / SI_HBAR
       mat%mr0 = mat%me0 * mat%mh0 / ( mat%me0 + mat%mh0 )
       mat%domega = SI_HBAR*mat%kmax*mat%kmax / ( 2*mat%numk*mat%mr0) 

       mat%cyc = 1

! Initialize coefficients

        mat%cc = 2. + DT * REAL_DX / SI_C * mat%gammanr
        mat%c0 = 2 * DT * REAL_DX / mat%cc * mat%pump
        mat%c1 = ( 2. - DT * REAL_DX / SI_C * mat%gammanr ) / mat%cc
        mat%c2 = ( 2. + DT * REAL_DX / SI_C * mat%gammap ) / mat%cc
        mat%c3 = ( 2. - DT * REAL_DX / SI_C * mat%gammap ) / mat%cc

        mat%dd = 1. + DT * REAL_DX / SI_C * mat%gammap
        mat%sigma =  (DT * REAL_DX / SI_C )**2 / mat%dd  * 2. * SI_E**2 / SI_HBAR * mat%M**2
        mat%d2 = ( DT * REAL_DX / SI_C * mat%gammap - 1. ) / mat%dd

        mat%K1 = 4.8966851
        mat%K2 = 0.04496457
        mat%K3 = 0.1333760

        mat%conv1 = SI_C / ( REAL_DX**(3./2.) * sqrt(SI_EPS0) )
        mat%conv2 = REAL_DX**(3./2.) / ( SI_C * sqrt(SI_EPS0) )

! THIS IS A EXPERIMENTAL HOT-FIX! 
! reg.f90 will need to be rewritten to again enable list-mode rather than mask-mode.

        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
        
        mat%ipos(p) = i
        mat%jpos(p) = j
        mat%kpos(p) = k

        })

! load from checkpoint file

       if ( load_state .and. detail_level .ge. 2 ) then

          read(UNITCHK) mat%Pk 
          read(UNITCHK) mat%Psum, mat%Qsum, mat%N

       end if

       M4_IFELSE_DBG({call EchoMatBulkSCObj(mat)},{call DisplayMatBulkSCObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatBulkSC")

  end subroutine InitializeMatBulkSC

!----------------------------------------------------------------------

  subroutine FinalizeMatBulkSC

    M4_MODLOOP_DECL({MATBULKSC},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_WRITE_DBG(". enter FinalizeMatBulkSC")
    M4_MODLOOP_EXPR({MATBULKSC},mat,{

       M4_MODOBJ_GETREG(mat,reg)

! save to checkpoint file

       if ( save_state .and. detail_level .ge. 2 ) then

          write(UNITCHK) mat%Pk
          write(UNITCHK) mat%Psum, mat%Qsum, mat%N

       end if

       ! finalize mat object here
       deallocate(mat%Pk, mat%Psum, mat%Qsum, mat%N)
       deallocate(mat%ipos, mat%jpos, mat%kpos)


    })
    M4_WRITE_DBG(". exit FinalizeMatBulkSC")

  end subroutine FinalizeMatBulkSC

!----------------------------------------------------------------------

  subroutine StepHMatBulkSC(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATBULKSC},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8), dimension(3) :: lE, psum, qsum, arg, pk
    real(kind=8) :: qem, qen
    real(kind=8) :: d1, d2, d3, beta
    real(kind=8) :: ne0, nh0, nue, nuh
    real(kind=8) :: mue, muh
    real(kind=8) :: omega, omegabarsq, enew, enhw
    real(kind=8) :: fermiew, fermihw
    integer :: ki, t

    M4_MODLOOP_EXPR({MATBULKSC},mat,{

    ! this loops over all mat structures, setting mat

      M4_MODOBJ_GETREG(mat,reg)

      n = mod(ncyc+1,2) + 1
      m = mod(ncyc+2,2) + 1

      mat%cyc = m
  
!$OMP PARALLEL DO &
!$OMP& PRIVATE(i,j,k,lE,qem,qen,ne0,nh0,mue,muh,ki,omega,enew,enhw,fermiew,fermihw,psum,qsum,omegabarsq,d1,d3,pk,arg) &
!$OMP& SHARED(Ex,Ey,Ez,reg,mat)
      do p = 1, reg%numnodes 
 
        ! E(n) at current position
  
        i = mat%ipos(p)
        j = mat%jpos(p)
        k = mat%kpos(p)

        lE(1) = Ex(i,j,k) * mat%conv1
        lE(2) = Ey(i,j,k) * mat%conv1
        lE(3) = Ez(i,j,k) * mat%conv1

        ! calculate second part of the density response

        qem = lE(1) * M4_QSUM(1,m,p) + lE(2) * M4_QSUM(2,m,p) + lE(3) * M4_QSUM(3,m,p) 
        qen = lE(1) * M4_QSUM(1,n,p) + lE(2) * M4_QSUM(2,n,p) + lE(3) * M4_QSUM(3,n,p) 
        mat%N(p) = mat%N(p) + 0.5 * ( mat%c2 * qen - mat%c3 * qem )

        ! calculate P(n+1) from P(n),P(n-1),E(n) and N(n)
        
        ! before: P(*,m) is P(n-1), P(*,n) is P(n)

        beta = 1./( SI_KB * mat%temp )  

        ! calculate chemical potentials using the Pade approximations

        ne0 = 0.25 * ( 2. * mat%me0 / ( SI_HBAR**2 * pi * beta ) )**( 3./2. )
        nh0 = 0.25 * ( 2. * mat%mh0 / ( SI_HBAR**2 * pi * beta ) )**( 3./2. ) 

        nue = mat%N(p) / ne0 
        nuh = mat%N(p) / nh0 

        mue = log( nue ) + mat%K1 * log( mat%K2 * nue + 1 ) + mat%K3 * nue
        muh = log( nuh ) + mat%K1 * log( mat%K2 * nuh + 1 ) + mat%K3 * nuh 

        ! update the microscopic polarisation variables

        psum = 0.
        qsum = 0.

        do ki = 0, mat%numk
           
           omega = mat%omegagap + ki * mat%domega

           omegabarsq =  omega**2 + mat%gammap**2

           enew = beta * ( mat%mr0 / mat%me0 ) * ( SI_HBAR * ( omega - mat%omegagap ) )
           enhw = beta * ( mat%mr0 / mat%mh0 ) * ( SI_HBAR * ( omega - mat%omegagap ) )

           fermiew = 1./(dexp( enew - mue ) + 1.)  
           fermihw = 1./(dexp( enhw - muh ) + 1.)

           d1 = ( 2. - omegabarsq * ( DT * REAL_DX / SI_C )**2 ) / mat%dd
           d3 = mat%sigma * omega * ( 1. - fermiew - fermihw )

           pk(:) = d1 * M4_PK(:,n,ki,p) + mat%d2 * M4_PK(:,m,ki,p) + d3 * lE(:)

           arg =  ( SI_HBAR * ( omega - mat%omegagap ) * 2. * mat%mr0 )**(0.5) * mat%domega * mat%mr0 * pk / (SI_HBAR**2) 

           M4_PK(:,m,ki,p) = pk(:) ! store new polarisation variables

           psum = psum + arg

           qsum = qsum + arg / ( SI_HBAR * omega )

        end do

        ! remove half of first/last element (trapezian rule)

        psum = psum - 0.5 * arg
        qsum = qsum - 0.5 * arg / ( SI_HBAR * omega ) 

        M4_PSUM(:,m,p) = psum(:) / ( pi ** 2 ) 
        M4_QSUM(:,m,p) = qsum(:) / ( pi ** 2 )

        ! calculate first part of the density response

        qem = lE(1) * M4_QSUM(1,m,p) + lE(2) * M4_QSUM(2,m,p) + lE(3) * M4_QSUM(3,m,p) 
        qen = lE(1) * M4_QSUM(1,n,p) + lE(2) * M4_QSUM(2,n,p) + lE(3) * M4_QSUM(3,n,p) 
        mat%N(p) = mat%c1 * mat%N(p) + 0.5 * ( mat%c2 * qem - mat%c3 * qen ) + mat%c0

        ! after: J(*,m) is now P(n+1)
        ! m and n will be flipped in the next timestep!

     end do
!$OMP END PARALLEL DO
    })

  end subroutine StepHMatBulkSC


!----------------------------------------------------------------------

  subroutine StepEMatBulkSC(ncyc)

    integer :: ncyc, m, n, ii, jj, kk
    M4_MODLOOP_DECL({MATBULKSC},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    n = mod(ncyc-1+2,2) + 1
    m = mod(ncyc+2,2) + 1

    M4_MODLOOP_EXPR({MATBULKSC},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg) 

!$OMP PARALLEL DO &
!$OMP& PRIVATE(i,j,k,w) &
!$OMP& SHARED(reg,mat,Ex,Ey,Ez)
       do p=1, reg%numnodes 

          i = mat%ipos(p)
          j = mat%jpos(p)
          k = mat%kpos(p)

          w(:) = reg%val(:,p)

       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

M4_IFELSE_TM({
       Ex(i,j,k) = Ex(i,j,k) - w(1) * epsinvx(i,j,k) * mat%conv2 * ( M4_PSUM(1,m,p) - M4_PSUM(1,n,p) )
       Ey(i,j,k) = Ey(i,j,k) - w(2) * epsinvy(i,j,k) * mat%conv2 * ( M4_PSUM(2,m,p) - M4_PSUM(2,n,p) )
})
M4_IFELSE_TE({
       Ez(i,j,k) = Ez(i,j,k) - w(3) * epsinvz(i,j,k) * mat%conv2 * ( M4_PSUM(3,m,p) - M4_PSUM(3,n,p) )
})

       end do
!$OMP END PARALLEL DO

})

  end subroutine StepEMatBulkSC

!----------------------------------------------------------------------

  subroutine SumJEMatBulkSC(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc, idx
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH)


    M4_MODLOOP_DECL({MatBulkSC},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MatBulkSC},mat,{

    idx = idx + NUMEBALCH

    })

  end subroutine SumJEMatBulkSC

!----------------------------------------------------------------------

  subroutine SumKHMatBulkSC(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc, idx
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH)


    M4_MODLOOP_DECL({MatBulkSC},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MatBulkSC},mat,{

    idx = idx + NUMEBALCH

    })

  end subroutine SumKHMatBulkSC
 
!----------------------------------------------------------------------

  subroutine DisplayMatBulkSCObj(mat)

    type(T_MATBULKSC) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)), &
         " gammap=",TRIM(f2str(mat%gammap,5)), &
         " M=",TRIM(f2str(mat%M,5)), &
         " egap=",TRIM(f2str(mat%egap,5)), &
         " me=",TRIM(f2str(mat%me,5)), &
         " mh=",TRIM(f2str(mat%mh,5)), &
         " N0=",TRIM(f2str(mat%N0,5)), &
         " gammanr=",TRIM(f2str(mat%gammanr,5)), &
         " temp=",TRIM(f2str(mat%temp,5)), &
         " kmax=",TRIM(f2str(mat%kmax,5)), &
         " numk=",TRIM(i2str(mat%numk))
     })
    call DisplayRegObj(regobj(mat%regidx))

  end subroutine DisplayMatBulkSCObj

!----------------------------------------------------------------------

   subroutine EchoMatBulkSCObj(mat)

    type(T_MATBULKSC) :: mat

    M4_WRITE_INFO({"--- matbulksc # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type) })

    ! -- write parameters to console 
    M4_WRITE_INFO({"gammap = ", mat%gammap })
    M4_WRITE_INFO({"M = ", mat%M })
    M4_WRITE_INFO({"egap = ", mat%egap })
    M4_WRITE_INFO({"me = ", mat%me })
    M4_WRITE_INFO({"mh = ", mat%mh })
    M4_WRITE_INFO({"N0 = ", mat%N0 })
    M4_WRITE_INFO({"gammanr = ", mat%gammanr })
    M4_WRITE_INFO({"temp = ", mat%temp })
    M4_WRITE_INFO({"kmax = ", mat%kmax })
    M4_WRITE_INFO({"numk = ", mat%numk })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatBulkSCObj

  
!----------------------------------------------------------------------

end module matbulksc

! Authors: J.Hamm J. Wood S. Wuestner
! Modified: 02/07/2013
!
! =====================================================================


