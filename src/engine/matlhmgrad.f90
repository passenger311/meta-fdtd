!-*- F90 -*------------------------------------------------------------
!
!  module: matlhmgrad / meta
!
!  Left-handed materials module.
!
!  subs:
!
!    InitializeMatlhmgrad
!    FinalizeMatlhmgrad
!    ReadMatlhmgradObj
!    StepEMatlhmgrad
!    StepHMatlhmgrad
!
!----------------------------------------------------------------------


! =====================================================================
!
! The Matlhmgrad calculates the material response of a Drude pole. 
!
! d/dt J + gammapl * J = omegapl**2 * E
! E = E* + J
!
! where E* is the electric field as calculated without the sources.  
!
! A first order (calculating J) or second order approach (calculating
! P) can be chosen.
!
! 1. order method
!
! StepHMatlhmgrad: update eq. J(n+1/2) = c1 * J(n-1/2) + c2 * E(n)
! StepEMatlhmgrad: update eq. E(n+1)* = E(n+1) - epsinv * DT * J(n+1/2)
!
! 2. order method (see matlorentz)
!
! StepHMatPdrude: update eq. P(n+1) = c1 * P(n) + c2 * P(n-1) + c3 * E(n)
! StepEMatPdrude: update eq. E(n+1)* = E(n+1) - epsinv * (P(n+1) - P(n))
!

module matlhmgrad

  use constant
  use parse
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MATHEAD_DECL({MATLHMGRAD},MAXMATOBJ,{

  ! input parameters
  real(kind=8) :: lambdaplinv ! vac. plasma wavelength [dx]
  real(kind=8) :: gammapl     ! current damping [1/dt]
  integer :: order = 1        ! use 1. or 2. order solver?
  character(len=80) :: filename
  integer, dimension(3) :: gradsp, gradep, gradv
  integer :: gradlen
  real(kind=8),dimension(:,:),pointer :: gradf

  real(kind=8) :: omegapl

  ! coefficients
  real(kind=8) :: c1, c2, c3
  

  ! current field: J (or Polarisation P) 
  ! magnetization current: K
  M4_FTYPE, dimension(:,:), pointer :: Jx, Jy, Jz
  M4_FTYPE, dimension(:,:), pointer :: Kx, Ky, Kz

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatlhmgradObj(funit,lcount)

    character(len=LINELNG) :: line
    logical :: eof, err

    M4_MODREAD_DECL({MATLHMGRAD},funit,lcount,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatlhmgradObj")
    
    M4_MODREAD_EXPR({MATLHMGRAD},funit,lcount,mat,reg,6,out,{ 

    ! read mat parameters here, as defined in mat data structure

    ! call readfloat(funit,lcount,mat%lambdaplinv)
    ! call readfloat(funit,lcount,mat%gammapl)
    ! call readint(funit,lcount,mat%order)

    call readstring(funit,lcount,mat%filename)
    
    err = .false.
    call readline(funit,lcount,eof,line)
    M4_EOF_ERROR({eof},lcount)
    call getintvec(line, mat%gradsp, 3, ":", err) 
    M4_SYNTAX_ERROR({err},lcount,{"BAD POINT FORMAT"})
    call readline(funit,lcount,eof,line)
    M4_EOF_ERROR({eof},lcount)
    call getintvec(line, mat%gradep, 3, ":", err) 
    M4_SYNTAX_ERROR({err},lcount,{"BAD POINT FORMAT"})

    if ( DIM .le. 2 ) then
       mat%gradsp(3) = KBEG
       mat%gradep(3) = KBEG
    end if

    if ( DIM .le. 1 ) then
       mat%gradsp(2) = JBEG
       mat%gradep(2) = JBEG
    end if

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatlhmgradObj")

  end subroutine ReadMatlhmgradObj

!----------------------------------------------------------------------

  subroutine InitializeMatlhmgrad

    type (T_REG) :: reg
    integer :: err, i, ios, l, j
    real(kind=8) :: v,w,f, t
    real(kind=8),dimension(:,:),allocatable :: fdata
    M4_MODLOOP_DECL({MATLHMGRAD},mat) 
    M4_WRITE_DBG(". enter InitializeMatlhmgrad")
    M4_MODLOOP_EXPR({MATLHMGRAD},mat,{
    
       ! initialize mat object here

       mat%omegapl = 2. * PI * mat%lambdaplinv
!       mat%gammapl = 2. / ( mat%abslenpl * DT )
       
       ! -- vector betwen gradsp and gradep
       
       mat%gradv(1) = mat%gradep(1) - mat%gradsp(1) 
       mat%gradv(2) = mat%gradep(2) - mat%gradsp(2) 
       mat%gradv(3) = mat%gradep(3) - mat%gradsp(3) 

       mat%gradlen = sqrt(real(mat%gradv(1))**2 + real(mat%gradv(2))**2 + real(mat%gradv(3))**2) + 0.5

       open(UNITTMP,FILE=mat%filename, STATUS="unknown", IOSTAT=ios)
       M4_OPEN_ERROR(ios,mat%filename)
       ios = 0
       l = 0
       do while ( .true. )
          read(UNITTMP,*,iostat = ios) v, w
          if ( ios .ne. 0 ) exit 
          l = l + 1
       end do
       l = l - 1
       M4_WRITE_INFO({"read ", TRIM(i2str(l+1)), " lines."})
       close(UNITTMP)
       allocate(fdata(2,0:l+1))
       open(UNITTMP,FILE=mat%filename, STATUS="unknown", IOSTAT=ios)
       do i = 0, l
          read(UNITTMP,*,iostat = ios) fdata(1,i), fdata(2,i)
       end do
       fdata(1,l+1) = fdata(1,l)
       fdata(2,l+1) = fdata(2,l)
       close(UNITTMP)
       M4_WRITE_INFO({"first: ",fdata(1,0),fdata(2,0)})
       M4_WRITE_INFO({"last: ",fdata(1,l),fdata(2,l)})
       
       allocate(mat%gradf(2,0:mat%gradlen+1), stat = err)
       M4_ALLOC_ERROR(err,"InitializeMatlhmgrad")
       
       f = real(l) / real(mat%gradlen)
       do i = 0, mat%gradlen
          j = int(i*f)
          t = i*f - j 
          mat%omegapl = 2. * PI * ( (1.-t)*fdata(1,j) + t*fdata(1,j+1) )
          mat%gammapl = (1.-t)*fdata(2,j) + t*fdata(2,j+1) 
          mat%gradf(1,i) = ( 2. - DT * mat%gammapl ) / ( 2. + DT * mat%gammapl )
          mat%gradf(2,i) = ( 2. * DT ) / ( 2. + DT * mat%gammapl ) * mat%omegapl**2
       end do
       mat%gradf(1,mat%gradlen+1) = mat%gradf(1,mat%gradlen) 
       mat%gradf(2,mat%gradlen+1) = mat%gradf(2,mat%gradlen) 
            
       deallocate(fdata)

       ! --

       reg = regobj(mat%regidx)

       if ( mat%order .lt. 1 .or. mat%order .gt. 2 ) then
          M4_FATAL_ERROR({"ORDER PARAMETER MUST BE 1 OR 2"})
       endif

       allocate(mat%Jx(reg%numnodes,mat%order),mat%Jy(reg%numnodes,mat%order),mat%Jz(reg%numnodes,mat%order), &
            mat%Kx(reg%numnodes,mat%order),mat%Ky(reg%numnodes,mat%order),mat%Kz(reg%numnodes,mat%order), stat = err)
       M4_ALLOC_ERROR(err,"InitializeMatlhmgrad")

       mat%Jx = 0.
       mat%Jy = 0.
       mat%Jz = 0.
       mat%Kx = 0.
       mat%Ky = 0.
       mat%Kz = 0.


       M4_IFELSE_DBG({call EchoMatlhmgradObj(mat)},{call DisplayMatlhmgradObj(mat)})

    })
    M4_WRITE_DBG(". exit InitializeMatlhmgrad")

  end subroutine InitializeMatlhmgrad

!----------------------------------------------------------------------

  subroutine FinalizeMatlhmgrad

    M4_MODLOOP_DECL({MATLHMGRAD},mat)
    M4_WRITE_DBG(". enter FinalizeMatlhmgrad")
    M4_MODLOOP_EXPR({MATLHMGRAD},mat,{

    ! finalize mat object here
    deallocate(mat%Jx,mat%Jy,mat%Jz,mat%Kx,mat%Ky,mat%Kz,mat%gradf)

    })
    M4_WRITE_DBG(". exit FinalizeMatlhmgrad")

  end subroutine FinalizeMatlhmgrad

!----------------------------------------------------------------------

  subroutine StepHMatlhmgrad(ncyc)

    integer :: ncyc, m, n, a
    real(kind=8) :: c1,c2,d,t
    M4_MODLOOP_DECL({MATLHMGRAD},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))

    M4_MODLOOP_EXPR({MATLHMGRAD},mat,{

    ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
      
       d = sqrt(real(i-mat%gradsp(1))**2 + real(j-mat%gradsp(2))**2 + real(k-mat%gradsp(3))**2)
       d = max(0.,d)
       d = min(d,real(mat%gradlen))
       
       a = int(d)
       t = d - a
       c1 = (1.-t)*mat%gradf(1,a) + t*mat%gradf(1,a+1) 
       c2 = (1.-t)*mat%gradf(2,a) + t*mat%gradf(2,a+1) 

       ! correct H(n+1/2)

M4_IFELSE_TE({
       Hx(i,j,k) = Hx(i,j,k) -  w(4) * M4_MUINVX(i,j,k) * DT * mat%Kx(p,1)
       Hy(i,j,k) = Hy(i,j,k) -  w(5) * M4_MUINVY(i,j,k) * DT * mat%Ky(p,1)
})
M4_IFELSE_TM({
       Hz(i,j,k) = Hz(i,j,k) -  w(6) * M4_MUINVZ(i,j,k) * DT * mat%Kz(p,1)
})

       ! calculate J(n+1/2) from J(n-1/2) and E(n)

M4_IFELSE_TM({
       mat%Jx(p,1) = c1 * mat%Jx(p,1) + c2 * Ex(i,j,k)
       mat%Jy(p,1) = c1 * mat%Jy(p,1) + c2 * Ey(i,j,k)
})
M4_IFELSE_TE({
       mat%Jz(p,1) = c1 * mat%Jz(p,1) + c2 * Ez(i,j,k)
})
       })      

   
    })
  
  end subroutine StepHMatlhmgrad


!----------------------------------------------------------------------


  subroutine StepEMatlhmgrad(ncyc)

    integer :: ncyc, m, n, a
    real(kind=8) :: c1,c2,d,t
    M4_MODLOOP_DECL({MATLHMGRAD},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))

    M4_MODLOOP_EXPR({MATLHMGRAD},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)
          
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       d = sqrt(real(i-mat%gradsp(1))**2 + real(j-mat%gradsp(2))**2 + real(k-mat%gradsp(3))**2)
       d = max(0.,d)
       d = min(d,real(mat%gradlen))
       
       a = int(d)
       t = d - a
       c1 = (1.-t)*mat%gradf(1,a) + t*mat%gradf(1,a+1) 
       c2 = (1.-t)*mat%gradf(2,a) + t*mat%gradf(2,a+1) 

       ! correct E(n+1)

M4_IFELSE_TM({
       Ex(i,j,k) = Ex(i,j,k) -  w(1) * epsinvx(i,j,k) * DT * mat%Jx(p,1)
       Ey(i,j,k) = Ey(i,j,k) -  w(2) * epsinvy(i,j,k) * DT * mat%Jy(p,1)
})
M4_IFELSE_TE({
       Ez(i,j,k) = Ez(i,j,k) -  w(3) * epsinvz(i,j,k) * DT * mat%Jz(p,1)
})

       ! calculate K(n+1) from K(n) and H(n+1/2)

M4_IFELSE_TE({
       mat%Kx(p,1) = c1 * mat%Kx(p,1) + c2 * Hx(i,j,k)
       mat%Ky(p,1) = c1 * mat%Ky(p,1) + c2 * Hy(i,j,k)
})
M4_IFELSE_TM({
       mat%Kz(p,1) = c1 * mat%Kz(p,1) + c2 * Hz(i,j,k)
})

       })      


    })

  end subroutine StepEMatlhmgrad

!----------------------------------------------------------------------

  real(kind=8) function SumJEMatlhmgrad(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
   
    M4_MODLOOP_DECL({MATLHMGRAD},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))

    sum = 0
    
    SumJEMatlhmgrad = sum
    
  end function SumJEMatlhmgrad

!----------------------------------------------------------------------

  real(kind=8) function SumKHMatlhmgrad(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
   
    M4_MODLOOP_DECL({MATLHMGRAD},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(6))

    sum = 0
    
    SumKHMatlhmgrad = sum
    
  end function SumKHMatlhmgrad

!----------------------------------------------------------------------

  subroutine DisplayMatlhmgradObj(mat)

    type(T_MATLHMGRAD) :: mat
 
   
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatlhmgradObj
  
!----------------------------------------------------------------------

   subroutine EchoMatlhmgradObj(mat)

    type(T_MATLHMGRAD) :: mat

    M4_WRITE_INFO({"--- matlhmgrad # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"omegapl = ",mat%omegapl })
    M4_WRITE_INFO({"gammapl = ",mat%gammapl })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatlhmgradObj

  
!----------------------------------------------------------------------

end module matlhmgrad

! Authors:  J.Hamm 
! Modified: 14/1/2008
!
! =====================================================================


