


























 













































































!-*- F90 -*------------------------------------------------------------
!
!  module: mattwolvl / meta
!
!  two level maxwell-bloch equations with 4th order Runge-Kutta
!
!  subs:
!
!    InitializeMattwoLvl
!    FinalizeMattwoLvl
!    ReadMattwoLvlObj
!    StepEMattwoLvl
!    StepHMattwoLvl
!    SumJEKHMattwoLvl
!
!----------------------------------------------------------------------


! =====================================================================
!
! The Mattwolvl module calculates the reponse of a 2 lvl 
! bloch system.
!
! The equations for the density matrix are given by
! d/dt rho = i/hbar * [H,rho]
! with H/hbar = [omega1,-Omega,-Omega*, omega2]
! and Omega_12 = conj(M) * E / hbar
! 
! E = D - eps_-1 * P
! P = M * rho12 * n
! n = number of 2lvl systems per grid cell.
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
! polarisation has to be converted as well, so that E=D-eps^-1P.
! 
! All material calculations and the calculation of the new electric
! field E are being done in StepEMattwoLvl.
! In StepHMattwoLvl the electric field is replaced by the effective
! displacement D in order to calculate the new D out of d/dt D = rot(H)
! 
! =====================================================================
module mattwolvl

  use constant
  use parse
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save


  
public :: SumJEMATTWOLVL
public :: SumKHMATTWOLVL


  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'MATTWOLVL'

  ! --- Public Methods

  public :: InitializeMATTWOLVL
  public :: FinalizeMATTWOLVL
  public :: StepEMATTWOLVL
  public :: StepHMATTWOLVL
  public :: ReadMATTWOLVLObj
  public :: EchoMATTWOLVLObj

  ! --- Public Data

  public :: T_MATTWOLVL
  public :: MATTWOLVLobj
  public :: numMATTWOLVLobj

  ! --- Constants

  integer, parameter :: MAXMATTWOLVLOBJ = MAXMATOBJ

  ! --- Types

  type T_MATTWOLVL

     character(len=STRLNG) :: type = "MATTWOLVL"
     integer :: regidx, idx       		! regobj index
     

  ! input parameters
  real(kind=8) :: lambdainv    
  real(kind=8) :: gamma       ! polarisation dephasing 
  real(kind=8) :: sigma       ! nonradiative recombination rate

  complex, dimension(3) :: M12(3) ! dipole length
  
  ! calculated

  real(kind=8) :: omega ! resonance frequency

  ! coefficients
  
  ! polarization
  complex(kind=8), dimension(:), pointer :: rho12 

  ! densities
  real(kind=8), dimension(:), pointer :: inversion ! rho22 - rho11
  ! starting values (have to add up to 1)
  real(kind=8) :: inversion_0

  ! number of 3-level atoms per unit-cell
  real(kind=8) :: n

  ! local field effect is included unless LFE==0
  real(kind=8) :: LFE
  
  ! Lorentz relation included unless Lorentz==0
  real(kind=8) :: Lorentz
  
  ! E = D - lfeval*P (lfeval = 2/3 for LFE and lfeval=1 for no LFE)
  real(kind=8) :: lfeval

  ! field conversion factor
  real(kind=8) :: conv

  ! old electric field
  real(kind=8), dimension(:,:), pointer :: Dold
  
  	

  end type T_MATTWOLVL

  ! --- Fields

  type(T_MATTWOLVL) :: MATTWOLVLobj(MAXMATTWOLVLOBJ) 
  integer :: numMATTWOLVLobj



contains

!----------------------------------------------------------------------

  subroutine ReadMattwolvlObj(funit,lcount)

    
    integer:: funit,lcount ! funit, lcount	
    type(T_MATTWOLVL) :: mat 
    type(T_REG) :: reg ! reg
    type(T_OUT) :: out ! out
    character(len=LINELNG) :: m4_string
    logical :: m4_eof


    complex(kind=8) :: c(3)
    logical :: eof,err
    character(len=LINELNG) :: line


    
    
    

 numMATTWOLVLobj = numMATTWOLVLobj + 1
 mat = MATTWOLVLobj(numMATTWOLVLobj)
 mat%idx = numMATTWOLVLobj	

  


    ! read mat parameters here, as defined in mat data structure
    call readfloat(funit, lcount, mat%lambdainv)

    call readfloat(funit, lcount, mat%gamma)

    call readfloat(funit, lcount, mat%sigma)
	
	call readcomplexs(funit,lcount, c, 3)
	mat%M12(:) = c(:)

    call readfloat(funit, lcount, mat%inversion_0)
    
    call readfloat(funit, lcount, mat%n)

    call readfloat(funit, lcount, mat%LFE)
	
	call readfloat(funit, lcount, mat%Lorentz)

    

! read regobj information

    call readline(funit,lcount, m4_eof, m4_string)
    

if ( m4_eof ) then
write(STDOUT,'(A,A,A,A,A,A,A)') " !ERR (",TRIM(modname),") ", &
"PARSE ERROR (@",TRIM(i2str(lcount-1)),") ",&
"UNEXPECTED EOF" 
stop

endif


    
    

if ( m4_string .ne. "(REG" ) then
write(STDOUT,'(A,A,A,A,A,A,A,A)') " !ERR (",TRIM(modname),") ", &
"PARSE ERROR (@",TRIM(i2str(lcount-1)),") ",&
"EXPECTED ",TRIM("(REG") 
stop

endif


    
    call ReadRegObj(reg, fdtdreg, funit, lcount, 3, .false.)
    mat%regidx = reg%idx

! get terminator

    do	
      call readline(funit,lcount, m4_eof, m4_string)
      

if ( m4_eof ) then
write(STDOUT,'(A,A,A,A,A,A,A)') " !ERR (",TRIM(modname),") ", &
"PARSE ERROR (@",TRIM(i2str(lcount-1)),") ",&
"UNEXPECTED EOF" 
stop

endif


      select case (m4_string)
      case("(OUT") 
	
       call ReadOutObj(out, reg, funit, lcount, modname, .true.)
       outobj(out%idx)%objidx = numMATTWOLVLobj
       mat%regidx = reg%idx
      case(")MATTWOLVL")
	exit
      case default
	

if ( .true. ) then
write(STDOUT,'(A,A,A,A,A,A,A,A)') " !ERR (",TRIM(modname),") ", &
"PARSE ERROR (@",TRIM(i2str(lcount-1)),") ",&
"EXPECTED ",TRIM(")MATTWOLVL") 
stop

endif


      end select
    enddo


! write back
    MATTWOLVLobj(numMATTWOLVLobj) = mat		


    call CompressValRegObj(reg) ! compress filling factors

    

  end subroutine ReadMattwolvlObj

!----------------------------------------------------------------------

  subroutine InitializeMattwolvl

    type (T_REG) :: reg
    integer :: err
    
    integer :: m4_n
    type(T_MATTWOLVL) :: mat   
 
    
    
  do m4_n = 1, numMATTWOLVLobj

       mat = MATTWOLVLobj(m4_n)
       
       
    
       ! initialize mat object here

       mat%omega = 2. * PI * mat%lambdainv

       ! field conversion factor 

       mat%conv = REAL_DX**(0.5)/sqrt(SI_EPS0) * SI_E / SI_HBAR
       
       ! Is local field effect included
       if (mat%LFE==0) THEN
          mat%lfeval = 1.
       else
          mat%lfeval = 2./3.
       end if
          

       reg = regobj(mat%regidx)

       allocate(mat%rho12(reg%numnodes), stat = err)

       
if ( err .ne. 0 ) then
write(STDOUT,*) "!ERR (",TRIM(modname),") ", &
"OUT OF MEMORY ","InitializeMattwolvl"
stop

endif


       ! set initial polarisations to zero
       mat%rho12(:) = 0

       allocate(mat%inversion(reg%numnodes), stat = err)

       
if ( err .ne. 0 ) then
write(STDOUT,*) "!ERR (",TRIM(modname),") ", &
"OUT OF MEMORY ","InitializeMattwolvl"
stop

endif

       
       ! set initial densities
       mat%inversion(:) = mat%inversion_0

       allocate(mat%Dold(reg%numnodes,3), stat = err)

       
if ( err .ne. 0 ) then
write(STDOUT,*) "!ERR (",TRIM(modname),") ", &
"OUT OF MEMORY ","InitializeMattwolvl"
stop

endif


       mat%Dold(:,:) = 0
       
       call DisplayMattwolvlObj(mat)

    

       MATTWOLVLobj(m4_n) = mat

  end do

    

  end subroutine InitializeMattwolvl

!----------------------------------------------------------------------

  subroutine FinalizeMattwolvl

    
    integer :: m4_n
    type(T_MATTWOLVL) :: mat   

    
    
  do m4_n = 1, numMATTWOLVLobj

       mat = MATTWOLVLobj(m4_n)
       
       

    ! finalize mat object here
    deallocate(mat%rho12,mat%inversion,mat%Dold)

    

       MATTWOLVLobj(m4_n) = mat

  end do

    

  end subroutine FinalizeMattwolvl

!----------------------------------------------------------------------

  subroutine StepHMattwolvl(ncyc)

    integer :: ncyc
    
    integer :: m4_n
    type(T_MATTWOLVL) :: mat   

    
 type(T_REG) :: reg
 integer :: p,i,j,k
 real(kind=8) :: w(3)


    
  do m4_n = 1, numMATTWOLVLobj

       mat = MATTWOLVLobj(m4_n)
       
       

    ! this loops over all mat structures, setting mat

      
  reg = regobj(mat%regidx)


      
! --- START M4: REGLOOP_EXPR
! reg,p,i,j,k,w,expr,spac1, spac2

if ( reg%numnodes .gt. 0 ) then

if ( reg%isbox ) then
p = 0
do k = reg%ks, reg%ke, reg%dk
do j = reg%js, reg%je, reg%dj
do i = reg%is, reg%ie, reg%di
p = p + 1	     

if ( .not. reg%isref ) then
if ( reg%compressval ) then 
w(1:reg%numval) = reg%val(1:reg%numval,reg%valptr(p))
else
w(1:reg%numval) = reg%val(1:reg%numval,p)
endif
endif 


      ! replace electric field with effective electric displacement
        
           Ex(i,j,k) = mat%Dold(p,1) / mat%conv
           Ey(i,j,k) = mat%Dold(p,2) / mat%conv
        
        
           Ez(i,j,k) = mat%Dold(p,3) / mat%conv
        

      
	   
enddo !i

enddo !j

enddo !k

else !isbox


do k = reg%ks, reg%ke, reg%dk ! -> k

do j = reg%js, reg%je, reg%dj ! -> j
!$OMP PARALLEL DO PRIVATE(p,w)
do i = reg%is, reg%ie, reg%di ! -> i
p = reg%mask(i,j,k) ! -> p
if ( p .gt. 0 ) then

if ( .not. reg%isref ) then
if ( reg%compressval ) then 
w(1:reg%numval) = reg%val(1:reg%numval,reg%valptr(p))
else
w(1:reg%numval) = reg%val(1:reg%numval,p)
endif
endif


      ! replace electric field with effective electric displacement
        
           Ex(i,j,k) = mat%Dold(p,1) / mat%conv
           Ey(i,j,k) = mat%Dold(p,2) / mat%conv
        
        
           Ez(i,j,k) = mat%Dold(p,3) / mat%conv
        

      
	   
endif 
enddo ! i
enddo ! j
enddo ! k

endif !isbox

endif !numnodes > 0

! --- END M4: REGLOOP_EXPR


    

       MATTWOLVLobj(m4_n) = mat

  end do

  
  end subroutine StepHMattwolvl


!----------------------------------------------------------------------

  subroutine StepEMattwolvl(ncyc)

    integer :: ncyc
    
    integer :: m4_n
    type(T_MATTWOLVL) :: mat   

    
 type(T_REG) :: reg
 integer :: p,i,j,k
 real(kind=8) :: w(3)

    real(kind=8), dimension(3) :: Dold, Dnew, D, Etmp
    complex(kind=8) :: ra12
    real(kind=8) :: ei(3)
    integer :: m, m1, m2
    complex(kind=8) :: k12(4)
    real(kind=8) :: kinversion(4)
    complex(kind=8) :: rho12, rho12o
    real(kind=8) :: inversion, inversion_o
    real(kind=8) :: om12, eps_lx, eps_ly, eps_lz

    D = 0
    ! determine which spatial field components have to be calculated
    m1 = 1
    m2 = 3






    
  do m4_n = 1, numMATTWOLVLobj

       mat = MATTWOLVLobj(m4_n)
       
       

       ! this loops over all mat structures, setting mat

    om12 = mat%omega

    
  reg = regobj(mat%regidx)


    
! --- START M4: REGLOOP_EXPR
! reg,p,i,j,k,w,expr,spac1, spac2

if ( reg%numnodes .gt. 0 ) then

if ( reg%isbox ) then
p = 0
do k = reg%ks, reg%ke, reg%dk
do j = reg%js, reg%je, reg%dj
do i = reg%is, reg%ie, reg%di
p = p + 1	     

if ( .not. reg%isref ) then
if ( reg%compressval ) then 
w(1:reg%numval) = reg%val(1:reg%numval,reg%valptr(p))
else
w(1:reg%numval) = reg%val(1:reg%numval,p)
endif
endif 


       

       ! If the electric permittivity is not equal to 1 the electric field acting on the particle is
       ! strengthened by a factor l = (eps + 2)/3. In effect Etmp(diel.) = Etmp * l.
       if (mat%Lorentz == 1) then
        	eps_lx = (1./epsinvx(i,j,k) + 2.)/3.
		eps_ly = (1./epsinvy(i,j,k) + 2.)/3.
		eps_lz = (1./epsinvz(i,j,k) + 2.)/3.
	else
		eps_lx = 1.
		eps_ly = 1.
		eps_lz = 1.
	end if
       ! set factor ei for E = D - ei*sum(real(Mij*rhoij))
       ! ei includes unit conversion and number of systems per unit cell
       ei(1) = 2 * w(1) * epsinvx(i,j,k) * mat%n * SI_4PIALPHA
       ei(2) = 2 * w(2) * epsinvy(i,j,k) * mat%n * SI_4PIALPHA
       ei(3) = 2 * w(3) * epsinvz(i,j,k) * mat%n * SI_4PIAlPHA
       
       ! save old values of rho
       rho12o = mat%rho12(p)
       inversion_o = mat%inversion(p)

       ! get new values of displacement and save old ones
       
          Dold(1) = mat%Dold(p,1)
          Dold(2) = mat%Dold(p,2)
          Dnew(1) = mat%conv * Ex(i,j,k)
          Dnew(2) = mat%conv * Ey(i,j,k)
          mat%Dold(p,1) = Dnew(1) 
          mat%Dold(p,2) = Dnew(2)
       
       
          Dold(3) = mat%Dold(p,3)
          Dnew(3) = mat%conv * Ez(i,j,k)
          mat%Dold(p,3) = Dnew(3)
       


       ! rk4: first step
       ! D(t) needed to calculate E(t)
       D(:) = Dold(:)
       
       rho12 = rho12o
       inversion = inversion_o
       
       do m= m1, m2

       ! Local field effect E(loc) = E(macroscopic) + P/3 = D - 2/3*P 
          Etmp(m) = D(m) - mat%lfeval * ei(m) * real(mat%M12(m)*rho12)

       end do
       
       ! calculate rabi frequencies at time t
       ra12 = mat%M12(1)*Etmp(1)*eps_lx + mat%M12(2)*Etmp(2)*eps_ly + &
            mat%M12(3)*Etmp(3)*eps_lz
       
       ! calculate first ks
       k12(1) = IMAG*( -rho12*om12 - ra12*inversion ) - mat%gamma * rho12
       kinversion(1) = 4*aimag(rho12*conjg(ra12)) - mat%sigma * (1 - inversion)

       ! rk4: second step
       ! D interpolated to t+dt/2 in order to calculate E(t+dt/2)
       D(:) = (Dold(:) + Dnew(:))/2
       
       ! rho interpolated to t+dt/2 according to runge-kutta scheme
       rho12 = rho12o + k12(1) * DT / 2
       inversion = inversion_o + kinversion(1) * DT / 2

       do m= m1, m2

       ! Local field effect E(loc) = E(macroscopic) + P/3 = D - 2/3*P
          Etmp(m) = D(m) - mat%lfeval * ei(m) * real( conjg(mat%M12(m))*rho12 )

       end do

       ! calculate rabi frequencies at t+dt/2
       ra12 = mat%M12(1)*Etmp(1)*eps_lx + mat%M12(2)*Etmp(2)*eps_ly + &
            mat%M12(3)*Etmp(3)*eps_lz

       k12(2) = IMAG*( -rho12*om12 - ra12*inversion ) - mat%gamma * rho12
       kinversion(2) = 4*aimag(rho12*conjg(ra12)) - mat%sigma * ( 1 - inversion )

       ! rk4: third step
       ! rho interpolated to t+dt/2 for a second time. D stays the same
       
       rho12 = rho12o + k12(2) * DT / 2
       inversion = inversion_o + kinversion(2) * DT / 2
	   
       do m= m1, m2

       ! Local field effect E(loc) = E(macroscopic) + P/3 = D - 2/3*P
          Etmp(m) = D(m) - mat%lfeval * ei(m) * real( conjg(mat%M12(m))*rho12 )

       end do
       
       ! calculate rabi frequencies at t+dt/2
       ra12 = mat%M12(1)*Etmp(1)*eps_lx + mat%M12(2)*Etmp(2)*eps_ly + & 
            mat%M12(3)*Etmp(3)*eps_lz


       k12(3) = IMAG*( -rho12*om12 - ra12*inversion ) - mat%gamma * rho12
       kinversion(3) = 4*aimag(rho12*conjg(ra12)) - mat%sigma * ( 1 - inversion )

       ! rk4: fourth step
       ! D(t+dt) for E(t+dt)
       D(:) = Dnew(:)

       ! rho(t+dt)
       rho12 = rho12o + k12(3) * DT
       inversion = inversion_o + kinversion(3) * DT

       do m= m1, m2

       ! Local field effect E(loc) = E(macroscopic) + P/3 = D - 2/3*P
       Etmp(m) = D(m) - mat%lfeval * ei(m) * real( conjg(mat%M12(m))*rho12 )

       end do

       ! calculate rabi frequencies at t+dt
       ra12 = mat%M12(1)*Etmp(1)*eps_lx + mat%M12(2)*Etmp(2)*eps_ly + &
            mat%M12(3)*Etmp(3)*eps_lz
 

       k12(4) = IMAG*( -rho12*om12 - ra12*inversion ) - mat%gamma * rho12
       kinversion(4) = 4*aimag(rho12*conjg(ra12)) - mat%sigma * ( 1 - inversion )

       ! value at end point

       rho12 = rho12o + DT/6.*( k12(1) + 2.*k12(2) + 2.*k12(3) + k12(4) ) 
       inversion = inversion_o + DT/6.*( kinversion(1) + 2.*kinversion(2) + &
				2.*kinversion(3) + kinversion(4) ) 

       ! update fields

       do m= m1, m2

          Etmp(m) = D(m) - ei(m) * &
            real( conjg(mat%M12(m))*rho12 )

       end do
       
       ! update density matrix
       mat%rho12(p) = rho12
       mat%inversion(p) = inversion

       ! convert electric field back into NHL-units and write it back

       Ex(i,j,k) = Etmp(1) / mat%conv
       Ey(i,j,k) = Etmp(2) / mat%conv


       Ez(i,j,k) = Etmp(3) / mat%conv




       
	   
enddo !i

enddo !j

enddo !k

else !isbox


do k = reg%ks, reg%ke, reg%dk ! -> k

do j = reg%js, reg%je, reg%dj ! -> j
!$OMP PARALLEL DO PRIVATE(p,w)
do i = reg%is, reg%ie, reg%di ! -> i
p = reg%mask(i,j,k) ! -> p
if ( p .gt. 0 ) then

if ( .not. reg%isref ) then
if ( reg%compressval ) then 
w(1:reg%numval) = reg%val(1:reg%numval,reg%valptr(p))
else
w(1:reg%numval) = reg%val(1:reg%numval,p)
endif
endif


       

       ! If the electric permittivity is not equal to 1 the electric field acting on the particle is
       ! strengthened by a factor l = (eps + 2)/3. In effect Etmp(diel.) = Etmp * l.
       if (mat%Lorentz == 1) then
        	eps_lx = (1./epsinvx(i,j,k) + 2.)/3.
		eps_ly = (1./epsinvy(i,j,k) + 2.)/3.
		eps_lz = (1./epsinvz(i,j,k) + 2.)/3.
	else
		eps_lx = 1.
		eps_ly = 1.
		eps_lz = 1.
	end if
       ! set factor ei for E = D - ei*sum(real(Mij*rhoij))
       ! ei includes unit conversion and number of systems per unit cell
       ei(1) = 2 * w(1) * epsinvx(i,j,k) * mat%n * SI_4PIALPHA
       ei(2) = 2 * w(2) * epsinvy(i,j,k) * mat%n * SI_4PIALPHA
       ei(3) = 2 * w(3) * epsinvz(i,j,k) * mat%n * SI_4PIAlPHA
       
       ! save old values of rho
       rho12o = mat%rho12(p)
       inversion_o = mat%inversion(p)

       ! get new values of displacement and save old ones
       
          Dold(1) = mat%Dold(p,1)
          Dold(2) = mat%Dold(p,2)
          Dnew(1) = mat%conv * Ex(i,j,k)
          Dnew(2) = mat%conv * Ey(i,j,k)
          mat%Dold(p,1) = Dnew(1) 
          mat%Dold(p,2) = Dnew(2)
       
       
          Dold(3) = mat%Dold(p,3)
          Dnew(3) = mat%conv * Ez(i,j,k)
          mat%Dold(p,3) = Dnew(3)
       


       ! rk4: first step
       ! D(t) needed to calculate E(t)
       D(:) = Dold(:)
       
       rho12 = rho12o
       inversion = inversion_o
       
       do m= m1, m2

       ! Local field effect E(loc) = E(macroscopic) + P/3 = D - 2/3*P 
          Etmp(m) = D(m) - mat%lfeval * ei(m) * real(mat%M12(m)*rho12)

       end do
       
       ! calculate rabi frequencies at time t
       ra12 = mat%M12(1)*Etmp(1)*eps_lx + mat%M12(2)*Etmp(2)*eps_ly + &
            mat%M12(3)*Etmp(3)*eps_lz
       
       ! calculate first ks
       k12(1) = IMAG*( -rho12*om12 - ra12*inversion ) - mat%gamma * rho12
       kinversion(1) = 4*aimag(rho12*conjg(ra12)) - mat%sigma * (1 - inversion)

       ! rk4: second step
       ! D interpolated to t+dt/2 in order to calculate E(t+dt/2)
       D(:) = (Dold(:) + Dnew(:))/2
       
       ! rho interpolated to t+dt/2 according to runge-kutta scheme
       rho12 = rho12o + k12(1) * DT / 2
       inversion = inversion_o + kinversion(1) * DT / 2

       do m= m1, m2

       ! Local field effect E(loc) = E(macroscopic) + P/3 = D - 2/3*P
          Etmp(m) = D(m) - mat%lfeval * ei(m) * real( conjg(mat%M12(m))*rho12 )

       end do

       ! calculate rabi frequencies at t+dt/2
       ra12 = mat%M12(1)*Etmp(1)*eps_lx + mat%M12(2)*Etmp(2)*eps_ly + &
            mat%M12(3)*Etmp(3)*eps_lz

       k12(2) = IMAG*( -rho12*om12 - ra12*inversion ) - mat%gamma * rho12
       kinversion(2) = 4*aimag(rho12*conjg(ra12)) - mat%sigma * ( 1 - inversion )

       ! rk4: third step
       ! rho interpolated to t+dt/2 for a second time. D stays the same
       
       rho12 = rho12o + k12(2) * DT / 2
       inversion = inversion_o + kinversion(2) * DT / 2
	   
       do m= m1, m2

       ! Local field effect E(loc) = E(macroscopic) + P/3 = D - 2/3*P
          Etmp(m) = D(m) - mat%lfeval * ei(m) * real( conjg(mat%M12(m))*rho12 )

       end do
       
       ! calculate rabi frequencies at t+dt/2
       ra12 = mat%M12(1)*Etmp(1)*eps_lx + mat%M12(2)*Etmp(2)*eps_ly + & 
            mat%M12(3)*Etmp(3)*eps_lz


       k12(3) = IMAG*( -rho12*om12 - ra12*inversion ) - mat%gamma * rho12
       kinversion(3) = 4*aimag(rho12*conjg(ra12)) - mat%sigma * ( 1 - inversion )

       ! rk4: fourth step
       ! D(t+dt) for E(t+dt)
       D(:) = Dnew(:)

       ! rho(t+dt)
       rho12 = rho12o + k12(3) * DT
       inversion = inversion_o + kinversion(3) * DT

       do m= m1, m2

       ! Local field effect E(loc) = E(macroscopic) + P/3 = D - 2/3*P
       Etmp(m) = D(m) - mat%lfeval * ei(m) * real( conjg(mat%M12(m))*rho12 )

       end do

       ! calculate rabi frequencies at t+dt
       ra12 = mat%M12(1)*Etmp(1)*eps_lx + mat%M12(2)*Etmp(2)*eps_ly + &
            mat%M12(3)*Etmp(3)*eps_lz
 

       k12(4) = IMAG*( -rho12*om12 - ra12*inversion ) - mat%gamma * rho12
       kinversion(4) = 4*aimag(rho12*conjg(ra12)) - mat%sigma * ( 1 - inversion )

       ! value at end point

       rho12 = rho12o + DT/6.*( k12(1) + 2.*k12(2) + 2.*k12(3) + k12(4) ) 
       inversion = inversion_o + DT/6.*( kinversion(1) + 2.*kinversion(2) + &
				2.*kinversion(3) + kinversion(4) ) 

       ! update fields

       do m= m1, m2

          Etmp(m) = D(m) - ei(m) * &
            real( conjg(mat%M12(m))*rho12 )

       end do
       
       ! update density matrix
       mat%rho12(p) = rho12
       mat%inversion(p) = inversion

       ! convert electric field back into NHL-units and write it back

       Ex(i,j,k) = Etmp(1) / mat%conv
       Ey(i,j,k) = Etmp(2) / mat%conv


       Ez(i,j,k) = Etmp(3) / mat%conv




       
	   
endif 
enddo ! i
enddo ! j
enddo ! k

endif !isbox

endif !numnodes > 0

! --- END M4: REGLOOP_EXPR
 

    

       MATTWOLVLobj(m4_n) = mat

  end do


  end subroutine StepEMattwolvl

!----------------------------------------------------------------------

  real(kind=8) function SumJEMattwolvl(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    real(kind=8) :: sum
    integer :: ncyc, m, n
   
    
    integer :: m4_n
    type(T_MATTWOLVL) :: mat   

    
 type(T_REG) :: reg
 integer :: p,i,j,k
 real(kind=8) :: w(3)


    sum = 0
    
    SumJEMattwolvl = sum
    
  end function SumJEMattwolvl

!----------------------------------------------------------------------

  real(kind=8) function SumKHMattwolvl(mask, ncyc)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc

    SumKHMattwolvl = 0.

  end function SumKHMattwolvl
 
!----------------------------------------------------------------------

  subroutine DisplayMattwolvlObj(mat)

    type(T_MATTWOLVL) :: mat
 
    write(STDOUT,*) "!INF (",TRIM(modname),") ", &
"#",TRIM(i2str(mat%idx)),&
    	" lambdainv=",TRIM(f2str(mat%lambdainv,5))
     
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMattwolvlObj

!----------------------------------------------------------------------

   subroutine EchoMattwolvlObj(mat)

    type(T_MATTWOLVL) :: mat

    write(STDOUT,*) "!INF (",TRIM(modname),") ", &
"--- mattwolvl # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type) 

    ! -- write parameters to console 
    write(STDOUT,*) "!INF (",TRIM(modname),") ", &
"lambdainv = ",mat%lambdainv  

    write(STDOUT,*) "!INF (",TRIM(modname),") ", &
"defined over:" 
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMattwolvlObj

  
!----------------------------------------------------------------------

end module mattwolvl

! Authors: A.Pusch 
! Modified: 13/8/2010
!
! =====================================================================


