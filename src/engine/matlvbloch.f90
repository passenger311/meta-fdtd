!-*- F90 -*------------------------------------------------------------
!
!  module: matlvbloch / meta
!
!  Effective Maxwell Bloch material module.
!
!  subs:
!
!    InitializeMatLVBloch
!    FinalizeMatLVBloch
!    ReadMatLVBlochObj
!    StepEMatLVBloch
!    StepHMatLVBloch
!    SumJEKHMatLVBloch
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatlvBloch module calculates the reponse of an effective  2 level 
! bloch system with langevin noise for dephasing and relaxation
!
! d/dt d/dt P + 2 * gammal * d/dt P + omegal**2 P =  - 4 * omegar * conv1 * M (M * E) f(N,Ntr)
! d/dt N = conv2 * E/ (omegar ) * ( d/dt P + gammal * P ) - gammanr * N 
! d/dt E = d/dt E* - epsinv * d/dt P 
!
! where E* is the electric field as calculated without the sources.  
!
! In natural HL units we choose hbar = c = 1, that is the elementary charge
! would be
!
! e = sqrt ( 4 PI alpha ) with alpha = 1/137.0360
!
! 
!
! StepHMatLVBloch: update eq. P(n+1) = c1 * P(n) + c2 * P(n-1) + c3 * E(n)
!                and update N (see below)
! StepEMatLVBloch: update eq. E(n+1) = E(n+1)* - epsinv * (P(n+1) - P(n))
!
! N update eq:
! N(n+1) = c4 * N(n) + c5 * ( P(n+1) - P(n) )*( E(n+1) + E(n) ) + c6 * ( P(n+1) + P(n) )*( E(n+1) + E(n) )
!
! is calculated in StepHMatLVBloch in two steps. The first bit depends on E(n+1) and must be 
! calculated *after* StepEMat updated E to the proper E(n+1). However,
! n+1 -> n, so that
!
! N(n) = N(n)_p1 + c5 * ( P(n) - P(n-1) )*( E(n) ) + c6 * ( P(n) + P(n-1) )*( E(n) )
!
! The second bit (after calculating P) is actually the first part of the N calculation:
!
! N(n+1)_p1 = c4 * N(n) + c5 * ( P(n+1) - P(n) )*( E(n) ) + c6 * ( P(n+1) + P(n) )*( E(n) )
!

module matlvbloch

  use constant
  use parse
  use reglist
  use outlist
  use grid
  use fdtd
  use ziggurat

  implicit none
  private
  save

  M4_MATHEAD_DECL({MATLVBLOCH},MAXMATOBJ,{

  ! input parameters
  real(kind=8) :: lambdarinv    ! inv. vac. plasma wavelength
  real(kind=8) :: gammal        ! resonance width

  real(kind=8) :: M(3)          ! dipole matrix vector [1/dt]
  real(kind=8) :: N0, Ntr       ! transparency density
  real(kind=8) :: gammanr, pump
  integer :: napprox, cyc
  
  ! calculated

  real(kind=8) :: omegar        ! resonance frequency
  real(kind=8) :: omegal        ! lorentz frequency

  ! coefficients
  real(kind=8) :: c1, c2, c3, c4, c5, c6, c7

  ! langevin coefficients
  real(kind=8) :: lvP1, lvP2, lvP3
  
  ! polarisation field 
  M4_FTYPE, dimension(:,:), pointer :: P

  ! number density 
  real(kind=8), dimension(:), pointer :: N, Nold
  
  ! old langevin term for P
  real(kind=8), dimension(:), pointer :: lvterm
  
  ! factor to account for actual 3D volume of 2D or 1D representations
  real(kind=8) :: volfac

  ! seed for random number generator
  integer :: seed

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatlvblochObj(funit,lcount)

    M4_MODREAD_DECL({MATLVBLOCH}, funit,lcount,mat,reg,out)
    real(kind=8) :: v(2)
    real(kind=8) :: c(3)
    logical :: eof,err
    character(len=LINELNG) :: line

    M4_WRITE_DBG(". enter ReadMatlvblochObj")
    
    M4_MODREAD_EXPR({MATLVBLOCH},funit,lcount,mat,reg,3,out,{ 


    ! read mat parameters here, as defined in mat data structure
    call readfloat(funit,lcount, mat%lambdarinv)
    call readfloat(funit,lcount, mat%gammal)
    call readfloats(funit,lcount, c, 3)

    mat%M(1) = c(1)
    mat%M(2) = c(2)
    mat%M(3) = c(3)


    call readline(funit, lcount, eof, line)
    M4_EOF_ERROR(eof, lcount)
    err = .false.
    call getfloat(line, mat%Ntr, err)
    M4_SYNTAX_ERROR({err},lcount,"Ntr [ N0 ]")
    err = .false.
    call getfloat(line, mat%N0, err)
    if ( err ) mat%N0 = mat%Ntr
    call readfloat(funit,lcount, mat%gammanr)
    call readfloat(funit,lcount, mat%pump)
    call readfloat(funit,lcount, mat%volfac)
    call readint(funit,lcount, mat%seed)

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatlvblochObj")

  end subroutine ReadMatlvblochObj

!----------------------------------------------------------------------

  subroutine InitializeMatlvbloch

    !type (T_REG) :: reg
    integer :: err
    real (kind=8) :: conv1, conv2, x1, x2, tipangle, re, im
    M4_MODLOOP_DECL({MATLVBLOCH},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3)) 
    M4_WRITE_DBG(". enter InitializeMatlvbloch")
    M4_MODLOOP_EXPR({MATLVBLOCH},mat,{
    
       ! initialize mat object here

       mat%omegar = 2. * PI * mat%lambdarinv
       
       mat%omegal = sqrt( mat%omegar**2 + mat%gammal**2 ) 

       reg = regobj(mat%regidx)

       allocate(mat%P(reg%numnodes,2), mat%N(reg%numnodes), mat%Nold(reg%numnodes), &
            mat%lvterm(reg%numnodes), stat = err)
       M4_ALLOC_ERROR(err,"InitializeMatlvbloch")

       mat%P = 0.

       mat%N = mat%N0
       mat%Nold = mat%N0

       mat%lvterm = 0.

! From unit conversion of polarisation equation -  SI to Computational units 
       conv1 = SI_4PIALPHA
! From unit conversion of population equation - SI to Computation units
       conv2 = ( ( REAL_DX * SI_C ) / ( SI_HBAR ) )

! for polarisation integration
       mat%c1 = ( 2. - mat%omegal**2 * DT**2 ) / ( 1. + DT * mat%gammal )
       mat%c2 = ( -1. + DT * mat%gammal ) / ( 1. + DT * mat%gammal )
       mat%c3 = ( - 4. ) * conv1 * DT**2 * mat%omegar  / ( 1. + DT * mat%gammal )

! for density integration
       mat%c4 = ( 2. - mat%gammanr * DT ) / ( 2. + mat%gammanr * DT )
       mat%c5 = 1. / ( 2. + mat%gammanr * DT ) * 1. / ( mat%omegar ) * conv2 
       mat%c6 = mat%c5  * DT * mat%gammal / 2.
       mat%c7 = ( ( 2. * DT ) / ( 2. + mat%gammanr * DT ) ) * mat%pump

       mat%cyc = 1 

! for langevin noise
       call zigset(mat%seed)
       ! set initial tipping angle
       if (mat%N0 == 2*mat%Ntr) then
         M4_MODOBJ_GETREG(mat,reg)
         M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
            tipangle = 2/sqrt(mat%Ntr*2*mat%volfac)*rnor()
            x1 = uni()
            x2 = uni()
            re = x1/sqrt(x1**2+x2**2)*sin(tipangle)*mat%Ntr/&
                 sqrt(REAL_DX)*8.983e-23
            im = x2/sqrt(x1**2+x2**2)*sin(tipangle)*mat%Ntr/&
                 sqrt(REAL_DX)*8.983e-23
            mat%P(p,2) = re
            mat%P(p,1) = re*(1-mat%gammal*DT) - im*mat%omegar*DT
            mat%N(p) = (cos(tipangle)+1)*mat%Ntr
            mat%Nold(p) = mat%N(p)
         })
       endif
       mat%lvP2 = mat%omegar*sqrt(DT*mat%gammal/mat%volfac)*8.982e-23/sqrt(REAL_DX)
       mat%lvP1 = sqrt(mat%gammal*DT/mat%volfac)*8.982e-23/sqrt(REAL_DX)
       mat%lvP3 = sqrt(mat%gammanr*DT/mat%volfac)
       
M4_IFELSE_1D({
       M4_WRITE_INFO({"1D -> forcing M(1)=0!"})
       mat%M(1) = 0.
})

! not TE -> no Ex/Ey coupling
M4_IFELSE_TM({},{
       M4_WRITE_INFO({"Not TM -> forcing M(3)=0!"})
       mat%M(3) = 0.  
})

! not TM -> no Ex/Ey coupling
M4_IFELSE_TM({},{
       M4_WRITE_INFO({"Not TE -> forcing M(1)=M(2)=0!"})
       mat%M(1) = 0.  
       mat%M(2) = 0.  
})


       M4_IFELSE_DBG({call EchoMatlvblochObj(mat)},{call DisplayMatlvblochObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatlvbloch")

  end subroutine InitializeMatlvbloch

!----------------------------------------------------------------------

  subroutine FinalizeMatlvbloch

    M4_MODLOOP_DECL({MATLVBLOCH},mat)
    M4_WRITE_DBG(". enter FinalizeMatlvbloch")
    M4_MODLOOP_EXPR({MATLVBLOCH},mat,{

    ! finalize mat object here
    deallocate(mat%P,mat%N)

    })
    M4_WRITE_DBG(". exit FinalizeMatlvbloch")

  end subroutine FinalizeMatlvbloch

!----------------------------------------------------------------------

  subroutine StepHMatlvbloch(ncyc)

    integer :: ncyc, m, n, zeta1, zeta2, zeta3
    M4_MODLOOP_DECL({MATLVBLOCH},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    M4_FTYPE :: me
    real(kind=8) :: pem, pen, ninv

    M4_MODLOOP_EXPR({MATLVBLOCH},mat,{

    ! this loops over all mat structures, setting mat

      M4_MODOBJ_GETREG(mat,reg)

      n = mod(ncyc-1+2,2) + 1
      m = mod(ncyc+2,2) + 1

      mat%cyc = m
      
        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{


        me = mat%M(1) * Ex(i,j,k) + mat%M(2) * Ey(i,j,k) + mat%M(3) * Ez(i,j,k)

        ! calculate second part of the density response (after the new E field got calculated)
        

        pem =  mat%P(p,m) * me
        pen =  mat%P(p,n) * me

      
        mat%N(p) = mat%N(p) + mat%c5 * ( pen - pem ) + mat%c6 * ( pen + pem )
        
        ! calculate P(n+1) from P(n),P(n-1),E(n) and N(n)
        
        ! before: P(*,m) is P(n-1), P(*,n) is P(n)
        
        ninv = ( mat%N(p) - mat%Ntr )
        ! get random numbers for polarization noise
        zeta1 = rnor()
        zeta2 = rnor()
        mat%P(p,m) = mat%c1 * mat%P(p,n) + mat%c2 * mat%P(p,m) + mat%c3 * me * ninv &
             - zeta2*sqrt(mat%N(p))*mat%lvP2 + ( zeta1*sqrt(mat%N(p)) - &
             mat%lvterm(p)*sqrt(mat%Nold(p)) )*mat%lvP1
        mat%lvterm(p) = zeta1
        mat%Nold(p) = mat%N(p)

        ! calculate first part of the density response
        

        pem =  mat%P(p,m) * me

        zeta3 = rnor()
        mat%N(p) = mat%c4 * mat%N(p) + mat%c5 * ( pem - pen ) + mat%c6 * ( pem + pen ) + mat%c7&
             + zeta3*sqrt(mat%N(p))*mat%lvP3

        ! after: J(*,m) is now P(n+1)
        ! m and n will be flipped in the next timestep!

        })      


    })
  
  end subroutine StepHMatlvbloch


!----------------------------------------------------------------------

  subroutine StepEMatlvbloch(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATLVBLOCH},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATLVBLOCH},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

M4_IFELSE_TM({
       Ex(i,j,k) = Ex(i,j,k) - w(1) * epsinvx(i,j,k) * mat%M(1) * ( mat%P(p,m) - mat%P(p,n) )
       Ey(i,j,k) = Ey(i,j,k) - w(2) * epsinvy(i,j,k) * mat%M(2) * ( mat%P(p,m) - mat%P(p,n) )
})
M4_IFELSE_TE({
       Ez(i,j,k) = Ez(i,j,k) - w(3) * epsinvz(i,j,k) * mat%M(3) * ( mat%P(p,m) - mat%P(p,n) )
})
       })      

    })

  end subroutine StepEMatlvbloch

!----------------------------------------------------------------------

  subroutine SumJEMatlvbloch(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH),sum1a,sum1b,sum2a,sum2b,ninva,ninvb,d34a,d34b
    integer :: ncyc, m, n, idx
   
    M4_MODLOOP_DECL({MATlvbloch},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))


    M4_MODLOOP_EXPR({MATlvbloch},mat,{

    sum1a = 0.
    sum1b = 0.
    sum2a = 0.
    sum2b = 0.

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

!       n = mod(ncyc-1+2,2) + 1
!       m = mod(ncyc+2,2) + 1

!       if ( mode ) then
!          d34a = mat%d4a
!          d34b = mat%d4b
!       else
!          d34a = mat%d3a
!          d34b = mat%d3b
!       endif

!       !4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

!       if ( mask(i,j,k) ) then
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
!               !4_VOLEY(i,j,k) * w(2) * Ey(i,j,k) * mat%N * ( mat%Pay(p,m) - mat%Pay(p,n) ) / DT +},{0. +}) &
!4_IFELSE_TE({ !4_VOLEZ(i,j,k) * w(3) * Ez(i,j,k) * mat%N * ( mat%Paz(p,m) - mat%Paz(p,n) ) / DT  },{0.  }) &
!               )

!          sum2b = sum2b + ( &
!4_IFELSE_TM({ !4_VOLEX(i,j,k) * w(1) * Ex(i,j,k) * mat%N * ( mat%Pbx(p,m) - mat%Pbx(p,n) ) / DT + &
!               !4_VOLEY(i,j,k) * w(2) * Ey(i,j,k) * mat%N * ( mat%Pby(p,m) - mat%Pby(p,n) ) / DT +},{0. +}) &
!4_IFELSE_TE({ !4_VOLEZ(i,j,k) * w(3) * Ez(i,j,k) * mat%N * ( mat%Pbz(p,m) - mat%Pbz(p,n) ) / DT  },{0.  }) &
!               )

!       endif

!       })      

!    sum1a = sum1a
!    sum1b = sum1b
!    sum2a = sum2a
!    sum2b = sum2b

!    sum(idx) = sum(idx)  + sum1a + sum1b
!    sum(idx+1) = sum(idx+1) + sum2a - sum1a
!    sum(idx+2) = sum(idx+2) + sum2b - sum1b
    idx = idx + NUMEBALCH

    })


  end subroutine SumJEMatlvbloch

!----------------------------------------------------------------------

  subroutine SumKHMatlvbloch(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum(MAXEBALCH)
    logical :: mode
    integer :: ncyc, idx

    M4_MODLOOP_DECL({MATlvbloch},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATlvbloch},mat,{

    idx = idx + NUMEBALCH

    })

  end subroutine SumKHMatlvbloch
 
!----------------------------------------------------------------------

  subroutine DisplayMatlvblochObj(mat)

    type(T_MATLVBLOCH) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" lambdarinv=",TRIM(f2str(mat%lambdarinv,5)),&
    	" gammal=",TRIM(f2str(mat%gammal,5))
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatlvblochObj

!----------------------------------------------------------------------

   subroutine EchoMatlvblochObj(mat)

    type(T_MATLVBLOCH) :: mat

    M4_WRITE_INFO({"--- matlvbloch # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"lambdarinv = ",mat%lambdarinv })
    M4_WRITE_INFO({"omegar = ",mat%omegar })
    M4_WRITE_INFO({"omegal = ",mat%omegal })
    M4_WRITE_INFO({"gammal = ",mat%gammal })
    M4_WRITE_INFO({"c1 = ",mat%c1 })
    M4_WRITE_INFO({"c2 = ",mat%c2 })
    M4_WRITE_INFO({"c3 = ",mat%c3 })
    M4_WRITE_INFO({"c4 = ",mat%c4 })
    M4_WRITE_INFO({"c5 = ",mat%c5 })
    M4_WRITE_INFO({"c6 = ",mat%c6 })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatlvblochObj

  
!----------------------------------------------------------------------

end module matlvbloch

! Authors:  A.Pusch, K.Boehringer, J.Hamm 
! Modified: 29/07/2011
!
! =====================================================================


