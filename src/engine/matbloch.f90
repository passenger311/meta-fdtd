!-*- F90 -*------------------------------------------------------------
!
!  module: matbloch / meta
!
!  Effective Maxwell Bloch material module.
!
!  subs:
!
!    InitializeMatBloch
!    FinalizeMatBloch
!    ReadMatBlochObj
!    StepEMatBloch
!    StepHMatBloch
!    SumJEKHMatBloch
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatBloch module calculates the reponse of an effective  2 level 
! bloch system
!
! d/dt d/dt P + 2 * gammal * d/dt P + omegal**2 P =  - 4 * omegar * 1/hbar * M (M * E) f(N,Ntr)
! d/dt N = E/(hbar * omegar ) * ( d/dt P + gammal * P ) - gammanr * N 
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
! StepHMatBloch: update eq. P(n+1) = c1 * P(n) + c2 * P(n-1) + c3 * E(n)
!                and update N (see below)
! StepEMatBloch: update eq. E(n+1) = E(n+1)* - epsinv * (P(n+1) - P(n))
!
! N update eq:
! N(n+1) = c4 * N(n) + c5 * ( P(n+1) - P(n) )*( E(n+1) + E(n) ) + c6 * ( P(n+1) + P(n) )*( E(n+1) + E(n) )
!
! is calculated in StepHMatBloch in two steps. The first bit depends on E(n+1) and must be 
! calculated *after* StepEMat updated E to the proper E(n+1). However,
! n+1 -> n, so that
!
! N(n) = N(n)_p1 + c5 * ( P(n) - P(n-1) )*( E(n) ) + c6 * ( P(n) + P(n-1) )*( E(n) )
!
! The second bit (after calculating P) is actually the first part of the N calculation:
!
! N(n+1)_p1 = c4 * N(n) + c5 * ( P(n+1) - P(n) )*( E(n) ) + c6 * ( P(n+1) + P(n) )*( E(n) )
!

module matbloch

  use constant
  use parse
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  real(kind=8), parameter :: hbar = 1. ! planck constant


  M4_MATHEAD_DECL({MATBLOCH},MAXMATOBJ,{

  ! input parameters
  real(kind=8) :: lambdarinv    ! inv. vac. plasma wavelength
  real(kind=8) :: gammal        ! resonance width

  M4_FTYPE :: M(3)          ! dipole matrix vector [1/dt]
  real(kind=8) :: N0, Ntr       ! transparency density
  real(kind=8) :: gammanr, pump
  integer :: napprox, cyc
  
  ! calculated

  real(kind=8) :: omegar        ! resonance frequency
  real(kind=8) :: omegal        ! lorentz frequency

  ! coefficients
  real(kind=8) :: c1, c2, c3, c4, c5, c6, c7
  
  ! polarisation field 
  M4_FTYPE, dimension(:,:), pointer :: P

  ! number density 
  real(kind=8), dimension(:), pointer :: N
  

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatBlochObj(funit,lcount)

    M4_MODREAD_DECL({MATBLOCH}, funit,lcount,mat,reg,out)
    real(kind=8) :: v(2)
    complex(kind=8) :: c(3)
    logical :: eof,err
    character(len=LINELNG) :: line

    M4_WRITE_DBG(". enter ReadMatBlochObj")
    
    M4_MODREAD_EXPR({MATBLOCH},funit,lcount,mat,reg,3,out,{ 


    ! read mat parameters here, as defined in mat data structure
    call readfloat(funit,lcount, mat%lambdarinv)
    call readfloat(funit,lcount, mat%gammal)
    call readcomplexs(funit,lcount, c, 3)
M4_IFELSE_CF({
    mat%M(1) = c(1)
    mat%M(2) = c(2)
    mat%M(3) = c(3)
},{
    mat%M(1) = real(c(1))
    mat%M(2) = real(c(2))
    mat%M(3) = real(c(3))
})

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

    call readint(funit,lcount, mat%napprox) 
    ! napprox = 0 -> N - Ntr
    ! napprox = 1 -> Ntr * ln(N/Ntr)
    

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatBlochObj")

  end subroutine ReadMatBlochObj

!----------------------------------------------------------------------

  subroutine InitializeMatBloch

    type (T_REG) :: reg
    integer :: err
    M4_MODLOOP_DECL({MATBLOCH},mat) 
    M4_WRITE_DBG(". enter InitializeMatBloch")
    M4_MODLOOP_EXPR({MATBLOCH},mat,{
    
       ! initialize mat object here

       mat%omegar = 2. * PI * mat%lambdarinv
       
       mat%omegal = sqrt( mat%omegar**2 + mat%gammal**2 ) 

       reg = regobj(mat%regidx)

       allocate(mat%P(reg%numnodes,2), mat%N(reg%numnodes), stat = err)

       M4_ALLOC_ERROR(err,"InitializeMatBloch")

       mat%P = 0.

       mat%N = mat%N0

! for polarisation integration
       mat%c1 = ( 2. - mat%omegal**2 * DT**2 ) / ( 1. + DT * mat%gammal )
       mat%c2 = ( -1. + DT * mat%gammal ) / ( 1. + DT * mat%gammal )
       mat%c3 = - 4./hbar * DT**2 * mat%omegar  / ( 1. + DT * mat%gammal )

! for density integration
       mat%c4 = ( 2. - mat%gammanr * DT ) / ( 2. + mat%gammanr * DT )
       mat%c5 = 1. / ( 2. + mat%gammanr * DT ) * 1./(hbar * mat%omegar ) 
       mat%c6 = mat%c5  * DT * mat%gammal / 2.
       mat%c7 = ( ( 2. * DT ) / ( 2. + mat%gammanr * DT ) ) * mat%pump

       mat%cyc = 1 


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


       M4_IFELSE_DBG({call EchoMatBlochObj(mat)},{call DisplayMatBlochObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatBloch")

  end subroutine InitializeMatBloch

!----------------------------------------------------------------------

  subroutine FinalizeMatBloch

    M4_MODLOOP_DECL({MATBLOCH},mat)
    M4_WRITE_DBG(". enter FinalizeMatBloch")
    M4_MODLOOP_EXPR({MATBLOCH},mat,{

    ! finalize mat object here
    deallocate(mat%P,mat%N)

    })
    M4_WRITE_DBG(". exit FinalizeMatBloch")

  end subroutine FinalizeMatBloch

!----------------------------------------------------------------------

  subroutine StepHMatBloch(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATBLOCH},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
    M4_FTYPE :: me
    real(kind=8) :: pem, pen, ninv

    M4_MODLOOP_EXPR({MATBLOCH},mat,{

    ! this loops over all mat structures, setting mat

      M4_MODOBJ_GETREG(mat,reg)

      n = mod(ncyc-1+2,2) + 1
      m = mod(ncyc+2,2) + 1

      mat%cyc = m
    
      select case ( mat%napprox )
      
      case ( 1 ) ! ---- LOGARITHMIC DENSITY DEPENDENCY

        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

M4_IFELSE_CF({
        me = conjg(mat%M(1)) * Ex(i,j,k) + conjg(mat%M(2)) * Ey(i,j,k) + conjg(mat%M(3)) * Ez(i,j,k)
},{
        me = mat%M(1) * Ex(i,j,k) + mat%M(2) * Ey(i,j,k) + mat%M(3) * Ez(i,j,k)
})

        ! calculate second part of the density response (after the new E field got calculated)
        
M4_IFELSE_CF({
        pem =  real( mat%P(p,m) * conjg(me) )
        pen =  real( mat%P(p,n) * conjg(me) )
},{
        pem =  mat%P(p,m) * me
        pen =  mat%P(p,n) * me
})
        
        mat%N(p) = mat%N(p) + mat%c5 * ( pen - pem ) + mat%c6 * ( pen + pem )
        
        ! calculate P(n+1) from P(n),P(n-1),E(n) and N(n)
        
        ! before: P(*,m) is P(n-1), P(*,n) is P(n)
        
        ninv = mat%Ntr * log( mat%N(p)/mat%Ntr )
       
        mat%P(p,m) = mat%c1 * mat%P(p,n) + mat%c2 * mat%P(p,m) + mat%c3 * me * ninv 

        ! calculate first part of the density response
        
M4_IFELSE_CF({
       pem =  real( mat%P(p,m) * conjg(me) )
},{
       pem =  mat%P(p,m) * me
})
        mat%N(p) = mat%c4 * mat%N(p) + mat%c5 * ( pem - pen ) + mat%c6 * ( pem + pen ) + mat%c7

        ! after: J(*,m) is now P(n+1)
        ! m and n will be flipped in the next timestep!

        })      
   
     case default ! ---- LINEAR DENSITY DEPENDENCY
        
        M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

M4_IFELSE_CF({
        me = conjg(mat%M(1)) * Ex(i,j,k) + conjg(mat%M(2)) * Ey(i,j,k) + conjg(mat%M(3)) * Ez(i,j,k)
},{
        me = mat%M(1) * Ex(i,j,k) + mat%M(2) * Ey(i,j,k) + mat%M(3) * Ez(i,j,k)
})
        ! calculate second part of the density response (after the new E field got calculated)
        
M4_IFELSE_CF({
        pem =  real( mat%P(p,m) * conjg(me) )
        pen =  real( mat%P(p,n) * conjg(me) )
},{
        pem =  mat%P(p,m) * me
        pen =  mat%P(p,n) * me
})
      
        mat%N(p) = mat%N(p) + mat%c5 * ( pen - pem ) + mat%c6 * ( pen + pem )
        
        ! calculate P(n+1) from P(n),P(n-1),E(n) and N(n)
        
        ! before: P(*,m) is P(n-1), P(*,n) is P(n)
        
        ninv = ( mat%N(p) - mat%Ntr )

        mat%P(p,m) = mat%c1 * mat%P(p,n) + mat%c2 * mat%P(p,m) + mat%c3 * me * ninv 

        ! calculate first part of the density response
        
M4_IFELSE_CF({
        pem =  real( mat%P(p,m) * conjg(me) )
},{
        pem =  mat%P(p,m) * me
})

        mat%N(p) = mat%c4 * mat%N(p) + mat%c5 * ( pem - pen ) + mat%c6 * ( pem + pen ) + mat%c7

        ! after: J(*,m) is now P(n+1)
        ! m and n will be flipped in the next timestep!

        })      

     end select


    })
  
  end subroutine StepHMatBloch


!----------------------------------------------------------------------

  subroutine StepEMatBloch(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATBLOCH},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATBLOCH},mat,{

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

  end subroutine StepEMatBloch

!----------------------------------------------------------------------

  real(kind=8) function SumJEMatBloch(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
   
    M4_MODLOOP_DECL({MATBLOCH},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    sum = 0

    M4_MODLOOP_EXPR({MATBLOCH},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       n = mod(ncyc-1+2,2) + 1
       m = mod(ncyc+2,2) + 1

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! correct E(n+1) using E(n+1)_fdtd and P(n+1),P(n)

       ! J(*,m) is P(n+1) and J(*,n) is P(n)      

       if ( mask(i,j,k) ) then

          sum = sum + ( &
M4_IFELSE_TM({ M4_VOLEX(i,j,k) * w(1) * real(Ex(i,j,k)) * real( mat%M(1) * ( mat%P(p,m) - mat%P(p,n) ) ) / DT +},{0. +}) &
M4_IFELSE_TM({ M4_VOLEY(i,j,k) * w(2) * real(Ey(i,j,k)) * real( mat%M(2) * ( mat%P(p,m) - mat%P(p,n) ) ) / DT +},{0. +}) &
M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * w(3) * real(Ez(i,j,k)) * real( mat%M(3) * ( mat%P(p,m) - mat%P(p,n) ) ) / DT  },{0.  }) &
               )

       endif

       })      

    })
    
    SumJEMatBloch = sum
    
  end function SumJEMatBloch

!----------------------------------------------------------------------

  real(kind=8) function SumKHMatBloch(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc

    SumKHMatBloch = 0.

  end function SumKHMatBloch
 
!----------------------------------------------------------------------

  subroutine DisplayMatBlochObj(mat)

    type(T_MATBLOCH) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	" lambdarinv=",TRIM(f2str(mat%lambdarinv,5)),&
    	" gammal=",TRIM(f2str(mat%gammal,5))
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatBlochObj

!----------------------------------------------------------------------

   subroutine EchoMatBlochObj(mat)

    type(T_MATBLOCH) :: mat

    M4_WRITE_INFO({"--- matbloch # ",&
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
    

  end subroutine EchoMatBlochObj

  
!----------------------------------------------------------------------

end module matbloch

! Authors:  K.Boehringer, J.Hamm 
! Modified: 14/1/2008
!
! =====================================================================


