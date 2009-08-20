!-*- F90 -*------------------------------------------------------------
!
!  module: lumped / meta
!
!  lumped elements. 
!
!----------------------------------------------------------------------

! ==========================================================
!
! The Lumped module provides treatment of basic Lumped Elements into the FDTD algorithm.
!
! 3 + 1 "Elements" are accounted for ::
! Resistor :: "R"		E(n+1) = (E(n+1) - V*E(n))/(1+V)					V = (dT*dL)/(2*dA*eps*R')   	
! Capacitor :: "C"	 	E(n+1) = (E(n+1) + V*E(n))/(1+V)					V = (dL*C')/(dA*eps) 	
! Inductor :: "L"		E(n+1) = E(n+1) - V*(2*(E(n)+E(n-1)+...+E(0)) - E(n))	V = (dL*dT^2)/(2*dA*L'*eps)	
! Conductivity :: "S" 	E(n+1) = (E(n+1) - V*E(n))/(1+V)					V = (dT*S')/(2*eps)
!
! The equations are constructed as such to provide a modification after standard FDTD steps
! and it is assumed that there is only one element at any unique location (one at Ex component of {i,j,k},
! one at Ey{i,j,k}, etc. Equations will fail if elements are defined at the same location. The inductor follows
! the equations set forth by Sui et al. [IEEE Trans. Micro. Theory Tech. v40 #4 p724 (1992)]; the equation
! in Taflove causes an instability. 
!
! Finally, the R', C', L' and S' values in the "V" terms above are modified values of R, C, L and S, taking into
! account unitary convertions. 
!
! 	R' = R*eps0*c			S' = (S*dL)/(eps0*c)
! 	L' = (L*eps0*c^2)/dL	C' =  C/(dL*eps0)
!
! config input Data Structure ::
!
! (LUMPED
! C						- Element Type	
! (REG
! (POINT
!	x y z : Cx Cy Cz		- Element Location and Value list. 
! )POINT
! )REG
! )LUMPED
!

module lumped
	
	USE strings
	USe constant
	Use reglist
	USe grid
	USE parse					
	USE fdtd					
	USE outlist					
	
	implicit none
	private
	save
	
	 ! --- Module Identifier

	character(len=STRLNG), parameter :: modname = 'LUMPED'

	! --- Public Methods

	PUBLIC :: ReadConfigLumped
	public :: InitializeLumped
	public :: FinalizeLumped
	public :: StepELumped
	public :: StepHLumped
	public :: DisplayLumpedObj
	public :: EchoLumpedObj

	! --- Public Data

	public :: T_LUMPED
	public :: LUMPEDobj
	public :: numLUMPEDobj

	! --- Constants

	integer, parameter :: MAXLUMPEDOBJ = MAXLUMPOBJ

	! --- Types

	type T_LUMPED

		character(len=STRLNG) :: type = "LUMPED"
		integer :: regidx, idx       		! regobj index
     						
		! input parameters
		character(len=1) :: element		! element type (R, C, L, etc).
		
		! Backup or Summation Fields, Linear Indexed
		REAL (8), allocatable, dimension(:,:) :: B
		real (8), allocatable, dimension(:,:) :: S
	
	end type T_LUMPED

	! --- Fields

	type(T_LUMPED) :: LUMPEDobj(MAXLUMPEDOBJ) 
	integer :: numLUMPEDobj
	
contains

!----------------------------------------------------------------

	subroutine ReadConfigLumped(funit, lcount)
				
		M4_MODREAD_DECL({LUMPED}, funit, lcount, lumped, reg, out)			
		
		M4_WRITE_DBG(". enter ReadLumpedObj")									
		
		M4_MODREAD_EXPR({LUMPED}, funit, lcount, lumped, reg, 3, out, {			

			! Read in lumped parameters, as per lumped data structure
			call readstring(funit,lcount,lumped%element)
					
		})		
				
		call CompressValRegObj(reg)
		
		M4_WRITE_DBG(". exit ReadLumpedObj")									

	end subroutine ReadConfigLumped
				
	subroutine InitializeLumped
		
		type (T_REG) :: reg
		integer :: err, i, j
		
		M4_MODLOOP_DECL({LUMPED}, lumped)									
		
		M4_WRITE_DBG(". enter InitializeLumped")								
		
		M4_MODLOOP_EXPR({LUMPED}, lumped, {									
		
			reg = regobj(lumped%regidx)
		
			allocate(lumped%B(reg%numnodes,1:3), stat = err)
			ALLOCATE(lumped%S(reg%numnodes,1:3), STAT = ERR)
			M4_ALLOC_ERROR(err,"InitializeLumped")
			
			do i = 1, reg%numnodes
				do j = 1, 3
					lumped%B(i,j) = 0.0d0
					lumped%S(i,j) = 0.0d0
				end do
			end do
			
			M4_IFELSE_DBG({call EchoLumpedObj(lumped)},{call DisplayLumpedObj(lumped)})		
		
		})
		
		M4_WRITE_DBG(". exit InitializeLumped") 								
							
	end subroutine
				
	subroutine FinalizeLumped
									
		M4_MODLOOP_DECL({LUMPED}, lumped)	
		
		M4_WRITE_DBG(". enter FinalizeLumped")
		
		M4_MODLOOP_EXPR({LUMPED}, lumped, {			

			deallocate(lumped%B,lumped%S)

		})
		
		M4_WRITE_DBG(". exit FinalizeLumped")
		
	end subroutine FinalizeLumped
	
	subroutine StepHLumped(ncyc)
				
	integer :: ncyc			
		
	M4_MODLOOP_DECL({LUMPED}, lumped)
		
	M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
	
	M4_MODLOOP_EXPR({LUMPED}, lumped, { 		

		M4_MODOBJ_GETREG(lumped,reg)			
	
		! Assign Backup / Summation Variables
		select case (TRIM(lumped%element))
		case ("R")										! Lumped Element Resistance
		
			M4_REGLOOP_EXPR(reg,p,i,j,k,w,{				
		
				lumped%B(p,1) = Ex(i,j,k)
				lumped%B(p,2) = Ey(i,j,k)
				lumped%B(p,3) = Ez(i,j,k)

			})

		CASE ("C")										! Lumped Element Capacitance
				
			M4_REGLOOP_EXPR(reg,p,i,j,k,w,{				
		
				lumped%B(p,1) = Ex(i,j,k)
				lumped%B(p,2) = Ey(i,j,k)
				lumped%B(p,3) = Ez(i,j,k)
				
			})
			
		CASE ("L")										! Lumped Element Inductance
			
			M4_REGLOOP_EXPR(reg,p,i,j,k,w,{				
		
				lumped%B(p,1) = Ex(i,j,k)		
				lumped%B(p,2) = Ey(i,j,k)
				lumped%B(p,3) = Ez(i,j,k)
		
				lumped%S(p,1) = lumped%S(p,1) + Ex(i,j,k)
				lumped%S(p,2) = lumped%S(p,2) + Ey(i,j,k)
				lumped%S(p,3) = lumped%S(p,3) + Ez(i,j,k)
				
			})
			
		CASE ("S")										! Conductivity Pseudo-'Element'
			
			M4_REGLOOP_EXPR(reg,p,i,j,k,w,{				
		
				lumped%B(p,1) = Ex(i,j,k)
				lumped%B(p,2) = Ey(i,j,k)
				lumped%B(p,3) = Ez(i,j,k)
				
			})
			
		case default 
			M4_FATAL_ERROR({"StepHLumped() :: Element Type Unclassified"})
		end select
		
	})

	end subroutine StepHLumped

	subroutine StepELumped(ncyc)
	
	integer :: ncyc			
	real(kind=8) :: Vx, Vy, Vz
		
	M4_MODLOOP_DECL({LUMPED}, lumped)
		
	M4_REGLOOP_DECL(reg,p,i,j,k,w(3))
		
	M4_MODLOOP_EXPR({LUMPED}, lumped, { 		
		
		select case (lumped%element)
		case("R")										! Lumped Element Resistance
			
			M4_MODOBJ_GETREG(lumped,reg)
			
			M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
			
				if ( w(1) .ne. 0.0d0 ) then
					Vx = epsinvx(i,j,k) * (DT * SX)/(2 * SY * SZ * w(1))
					Ex(i,j,k) = (Ex(i,j,k) - Vx * lumped%B(p,1))/(1 + Vx)
				end if
				
				if ( w(2) .ne. 0.0d0 ) then
					Vy = epsinvy(i,j,k) * (DT * SY)/(2 * SZ * SX * w(2))
					Ey(i,j,k) = (Ey(i,j,k) - Vy * lumped%B(p,2))/(1 + Vy)
				end if
				
				if ( w(3) .ne. 0.0d0 ) then
					Vz = epsinvz(i,j,k) * (DT * SZ)/(2 * SX * SY * w(3))
					Ez(i,j,k) = (Ez(i,j,k) - Vz * lumped%B(p,3))/(1 + Vz)
				end if
			
			})
			
		case("C")										! Lumped Element Capacitance
			
			M4_MODOBJ_GETREG(lumped,reg)
			
			M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
			
				if ( w(1) .ne. 0.0d0 ) then
					Vx = epsinvx(i,j,k) * (SX * w(1))/(SY * SZ) 
					Ex(i,j,k) = (Ex(i,j,k) + Vx * lumped%B(p,1))/(1 + Vx)
				end if
					
				if ( w(2) .ne. 0.0d0 ) then
					Vy = epsinvy(i,j,k) * (SY * w(2))/(SZ * SX) 
					Ey(i,j,k) = (Ey(i,j,k) + Vy * lumped%B(p,2))/(1 + Vy)
				end if
				
				if ( w(3) .ne. 0.0d0 ) then
					Vz = epsinvz(i,j,k) * (SZ * w(3))/(SX * SY) 	
					Ez(i,j,k) = (Ez(i,j,k) + Vz * lumped%B(p,3))/(1 + Vz)
				end if
				
			})	
			
		case("L")										! Lumped Element Inductance
			
			M4_MODOBJ_GETREG(lumped,reg)
			
			M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
			
				if ( w(1) .ne. 0.0d0 ) then
					Vx = epsinvx(i,j,k) * (DT**2 * SX)/(2 * SY * SZ * w(1))
					Ex(i,j,k) = Ex(i,j,k) - Vx * (2 * lumped%S(p,1) - lumped%B(p,1))
				end if
				
				if ( w(2) .ne. 0.0d0 ) then
					Vy = epsinvy(i,j,k) * (DT**2 * SY)/(2 * SZ * SX * w(2))
					Ey(i,j,k) = Ey(i,j,k) - Vy * (2 * lumped%S(p,2) - lumped%B(p,2))
			
				end if
				
				if ( w(3) .ne. 0.0d0 ) then
					Vz = epsinvz(i,j,k) * (DT**2 * SZ)/(2 * SX * SY * w(3))
					Ez(i,j,k) = Ez(i,j,k) - Vz * (2 * lumped%S(p,3) - lumped%B(p,3))
				end if
	
			})

		case("S")										! Conductivity
			
			M4_MODOBJ_GETREG(lumped,reg)
			
			M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
			
				if ( w(1) .ne. 0.0d0 ) then
					Vx = epsinvx(i,j,k) * (DT * w(1))/2 
					Ex(i,j,k) = (Ex(i,j,k) - Vx * lumped%B(p,1))/(1 + Vx)
				end if
				
				if ( w(2) .ne. 0.0d0 ) then
					Vy = epsinvy(i,j,k) * (DT * w(2))/2
					Ey(i,j,k) = (Ey(i,j,k) - Vy * lumped%B(p,2))/(1 + Vy)
				end if
				
				if ( w(3) .ne. 0.0d0 ) then
					Vz = epsinvz(i,j,k) * (DT * w(3))/2
					Ez(i,j,k) = (Ez(i,j,k) - Vz * lumped%B(p,3))/(1 + Vz)
				end if
			
			})	
			
		case default
			M4_FATAL_ERROR({"StepELumped() :: Element Type Unclassified"})
		end select
			
	})

	end subroutine StepELumped
				
	subroutine DisplayLumpedObj(lumped)
	
		type(T_LUMPED) :: lumped
			
		M4_WRITE_INFO({"#",TRIM(i2str(lumped%idx))," Element Type  = ",TRIM(lumped%element)})
		
		call DisplayRegObj(regobj(lumped%regidx))
		
	END SUBROUTINE DisplayLumpedObj
	
	subroutine EchoLumpedObj(lumped)
	
		type(T_LUMPED) :: lumped
		
		M4_WRITE_INFO({"--- lumped # ", TRIM(i2str(lumped%idx))," ", TRIM(lumped%type)})
		
		! write parameters to console		
		M4_WRITE_INFO({"Element Type = ", lumped%element})
		M4_WRITE_INFO({"defined over:"})
		
		call EchoRegObj(regobj(lumped%regidx))
		
	end subroutine EchoLumpedObj

!----------------------------------------------------------------

end module lumped

! Author: R.Crowter, J.Hamm
! Modified: 18/08/2009
!
! ==========================================================
