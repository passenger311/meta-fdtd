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
  integer :: numk                ! number of k points for band discretisation

  !

  integer :: cyc
  
  ! calculated

  real(kind=8) :: dk
  real(kind=8) :: mr0, me0, mh0     ! mass in [kg]
  real(kind=8) :: omegagap

  ! polarisation field 
  M4_FTYPE, dimension(:,:,:,:), pointer :: Pk

  ! sum of polarisations

  M4_FTYPE, dimension(:,:,:), pointer :: Qsum
  M4_FTYPE, dimension(:,:,:), pointer :: Psum


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

    call readfloat(funit, lcount, mat%M)
    call readfloat(funit, lcount, mat%egap)
    call readfloat(funit, lcount, mat%me)
    call readfloat(funit, lcount, mat%mh)
    call readfloat(funit, lcount, mat%gammap)
    call readfloat(funit, lcount, mat%gammanr)
    call readfloat(funit, lcount, mat%N0)
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

       allocate(mat%Pk(3,2,0:mat%numk,reg%numnodes), stat=err)
       allocate(mat%N(reg%numnodes), mat%Qsum(3, 2, reg%numnodes), mat%Psum(3, 2, reg%numnodes), stat=err)

       M4_ALLOC_ERROR(err,"InitializeMatBulkSC")

       mat%Pk = 0.
       mat%Psum = 0.
       mat%Qsum = 0.
       mat%N = mat%N0
  
! Initialize omegak(k)

       mat%dk = mat%kmax / mat%numk
       mat%me0 = mat%me * SI_ME
       mat%mh0 = mat%mh * SI_ME
       mat%omegagap = mat%egap * SI_E / SI_HBAR
       mat%mr0 = mat%me0 * mat%mh0 / ( mat%me0 + mat%mh0 )

       mat%cyc = 1

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

    })
    M4_WRITE_DBG(". exit FinalizeMatBulkSC")

  end subroutine FinalizeMatBulkSC

!----------------------------------------------------------------------

  subroutine StepHMatBulkSC(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATBULKSC},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8), dimension(3) :: lE, psum, qsum, arg, pk
    real(kind=8) :: c1, c2, c3, dd
    real(kind=8) :: qem, qen
    real(kind=8) :: sigma, beta, ne0, nh0, nue, nuh
    real(kind=8) :: K1, K2, K3, mue, muh
    real(kind=8) :: kk, omegak, omegabarksq, enek, enhk
    real(kind=8) :: fermiek, fermihk
    real(kind=8) :: conv1
    integer :: ki

    M4_MODLOOP_EXPR({MATBULKSC},mat,{

    ! this loops over all mat structures, setting mat

      M4_MODOBJ_GETREG(mat,reg)

      n = mod(ncyc-1+2,2) + 1
      m = mod(ncyc+2,2) + 1

      mat%cyc = m
      
        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

        ! E(n) at current position
 
        conv1 = SI_C / ( REAL_DX**(3/2) * sqrt(SI_EPS0) )

        lE(1) = Ex(i,j,k) * conv1
        lE(2) = Ey(i,j,k) * conv1
        lE(3) = Ez(i,j,k) * conv1

        ! calculate second part of the density response

        dd = 2. + DT * mat%gammanr
        c1 = ( 2. - DT * mat%gammanr ) / dd
        c2 = ( 2. + DT * mat%gammap ) / dd
        c3 = ( 2. - DT * mat%gammap ) / dd

        qem = lE(1) * mat%Qsum(1,m,p) + lE(2) * mat%Qsum(2,m,p) + lE(3) * mat%Qsum(3,m,p) 
        qen = lE(1) * mat%Qsum(1,n,p) + lE(2) * mat%Qsum(2,n,p) + lE(3) * mat%Qsum(3,n,p) 

        mat%N(p) = mat%N(p) + 0.5 * ( c2 * qem - c3 * qen )

        ! calculate P(n+1) from P(n),P(n-1),E(n) and N(n)
        
        ! before: P(*,m) is P(n-1), P(*,n) is P(n)

        dd = 1 + mat%gammap * DT
        c2 = ( mat%gammap * DT - 1 ) / dd
        sigma =  DT**2 / dd  *  2 * SI_E**2 / SI_HBAR * mat%M**2
        beta = 1./( SI_KB * mat%temp )  

        ! calculate chemical potentials using the Pade approximations

        ne0 = 0.25 * ( 2. * mat%me0 / ( SI_HBAR**2 * pi * beta ) )**( 3./2. ) 
        nh0 = 0.25 * ( 2. * mat%mh0 / (SI_HBAR**2 * pi * beta ) )**( 3./2. ) 

        nue = mat%N(p) / ne0 
        nuh = mat%N(p) / nh0 

        K1 = 4.8966851
        K2 = 0.04496457
        K3 = 0.1333760

        mue = log( nue ) + K1 * log( K2 * nue + 1 ) + K3 * nue 
        muh = log( nue ) + K1 * log( K2 * nue + 1 ) + K3 * nue 

        ! update the microscopic polarisation variables

        psum = 0.
        qsum = 0.

        do ki = 0, mat%numk 
           
           kk = ki * mat%dk 

           omegak = mat%omegagap + SI_HBAR * kk**2 / ( 2*mat%mr0 )
           omegabarksq =  omegak**2 + mat%gammap**2
                     
           enek = beta * SI_HBAR**2 * kk**2 / ( 2*mat%me0 ) 
           enhk = beta * SI_HBAR**2 * kk**2 / ( 2*mat%mh0 ) 

           fermiek = 1./(exp( enek - mue ) + 1.)  
           fermihk = 1./(exp( enhk - muh ) + 1.)
           
           c1 = ( 2. - omegabarksq * DT**2 ) / dd
           c3 = sigma * omegak * ( 1- fermiek - fermihk )

           pk(:) = c1 * mat%Pk(:,n,ki,p) + c2 * mat%Pk(:,m,ki,p) + c3 * lE(:)
           
           ! sum up polarisations
           
           arg = kk**2 * mat%dk * pk

           psum = psum + arg

           qsum = qsum + arg / ( SI_HBAR * omegak )

           mat%Pk(:,m,ki,p) = pk(:) ! store new polarisation variables
           
        end do
  
        ! remove half of first/last element (trapezian rule)

        psum = psum - 0.5 * kk**2 * mat%dk * ( mat%Pk(:,m,0,p) + mat%Pk(:,m,mat%numk,p) )
        qsum = qsum - 0.5 * kk**2 * mat%dk * 1./SI_HBAR * ( & 
             mat%Pk(:,m,0,p) / mat%omegagap &
             + mat%Pk(:,m,mat%numk,p) / ( mat%omegagap + SI_HBAR * kk**2 / ( 2*mat%mr0 ) ) &
             )

        mat%Psum(:,m,p) = psum(:) / ( pi ** 2 ) 
        mat%Qsum(:,m,p) = qsum(:) / ( pi ** 2 )

        ! calculate first part of the density response

        dd = 2. + DT * mat%gammanr
        c1 = ( 2. - DT * mat%gammanr ) / dd
        c2 = ( 2. + DT * mat%gammap ) / dd
        c3 = ( 2. - DT * mat%gammap ) / dd

        qem = lE(1) * mat%Qsum(1,m,p) + lE(2) * mat%Qsum(2,m,p) + lE(3) * mat%Qsum(3,m,p) 
        qen = lE(1) * mat%Qsum(1,n,p) + lE(2) * mat%Qsum(2,n,p) + lE(3) * mat%Qsum(3,n,p) 

        mat%N(p) = c1 * mat%N(p) + 0.5 * ( c2 * qem - c3 * qen )

        ! after: J(*,m) is now P(n+1)
        ! m and n will be flipped in the next timestep!

        })      


    })
  
  end subroutine StepHMatBulkSC


!----------------------------------------------------------------------

  subroutine StepEMatBulkSC(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATBULKSC},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    real(kind=8) :: conv2

    conv2 = REAL_DX**(3/2) / ( SI_C * sqrt(SI_EPS0) )

    M4_MODLOOP_EXPR({MATBULKSC},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

M4_IFELSE_TM({
       Ex(i,j,k) = Ex(i,j,k) - w(1) * epsinvx(i,j,k) * conv2 * ( mat%Psum(1,m,p) - mat%Psum(1,n,p) )
       Ey(i,j,k) = Ey(i,j,k) - w(2) * epsinvy(i,j,k) * conv2 * ( mat%Psum(2,m,p) - mat%Psum(2,n,p) )
})
M4_IFELSE_TE({
       Ez(i,j,k) = Ez(i,j,k) - w(3) * epsinvz(i,j,k) * conv2 *( mat%Psum(3,m,p) - mat%Psum(3,n,p) )
})
       })      

    })

  end subroutine StepEMatBulkSC

!----------------------------------------------------------------------

  subroutine SumJEMatBulkSC(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc, idx
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH)

  end subroutine SumJEMatBulkSC

!----------------------------------------------------------------------

  subroutine SumKHMatBulkSC(mask, ncyc, sum, idx, mode)

   logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc, idx
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH)
  end subroutine SumKHMatBulkSC
 
!----------------------------------------------------------------------

  subroutine DisplayMatBulkSCObj(mat)

    type(T_MATBULKSC) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)), &
         " M=",TRIM(f2str(mat%M,5))
     })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatBulkSCObj

!----------------------------------------------------------------------

   subroutine EchoMatBulkSCObj(mat)

    type(T_MATBULKSC) :: mat

    M4_WRITE_INFO({"--- matbulksc # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type) })

    ! -- write parameters to console 
    M4_WRITE_INFO({"M = ", mat%M })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatBulkSCObj

  
!----------------------------------------------------------------------

end module matbulksc

! Authors: J.Hamm J. Wood
! Modified: 16/05/2013
!
! =====================================================================


